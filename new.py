import socket
import threading
import galois
import numpy as np
import json
import time
import copy
from datetime import datetime
import hashlib
from utils import (generator1, generator2, curve_order, normalize, validate_point, GPoint, SRS, numbers_to_hash, generate_challenge, patch_galois, dump_proof, load_proof, dump_circuit, load_circuit)


#node simulation
class Node:
	def __init__(self, node_id, crs, port, other_ports):
		self.node_id = node_id
		self.crs = crs
		self.port = port
		self.other_ports = other_ports
		self.ledger = []
		threading.Thread(target=self.listen_for_transactions, daemon=True).start()
	def listen_for_transactions(self):
		server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		server_socket.bind(('localhost', self.port))
		server_socket.listen()
		print(f"Node {self.node_id} listening on port {self.port}")
		
		while True:
			client_socket, _ = server_socket.accept()
			threading.Thread(target=self.handle_transaction, args=(client_socket,), daemon=True).start()
	def handle_transaction(self, client_socket):
		transaction_data = client_socket.recv(4096).decode()
		transaction = json.loads(transaction_data)
		client_socket.close()
		global checker
		deserialized_transaction = json_deserialize(transaction)
		
		if self.verify_transaction(deserialized_transaction):
			print(f"Node {self.node_id} verified transaction {transaction['tx_id']} successfully.")
			self.ledger.append(transaction)
		else:
			print(f"Node {self.node_id} failed to verify transaction {transaction['tx_id']}.")
			checker = False
	def broadcast_transaction(self, transaction):
		transaction_data = json.dumps(transaction).encode()
		
		for port in self.other_ports:
			with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
				try:
					s.connect(('localhost', port))
					s.sendall(transaction_data)
				except ConnectionRefusedError:
					print(f"Node {self.node_id} could not connect to Node on port {port}")
	def create_transaction(self, tx_id, tx_type, details):
		organ_type = int(hashlib.md5(details["organ_type"].encode()).hexdigest(), 16)
		name = int(hashlib.md5(details["name"].encode()).hexdigest(), 16)
		patient_id = int(hashlib.md5(details["patient_id"].encode()).hexdigest(), 16)
		blood_group = int(hashlib.md5(details["blood_group"].encode()).hexdigest(), 16)
		age = int(hashlib.md5(details["age"].encode()).hexdigest(), 16)
		organ_life = int(hashlib.md5(details["organ_life"].encode()).hexdigest(), 16)
		n, omega, roots = setup1(organ_type, name, patient_id, blood_group, age, organ_life)
		a, b, c, pi, ql, qr, qm, qc, qo = witnessGatesMain(n, organ_type, name, patient_id, blood_group, age, organ_life)
		sigmaFor1 = sigmaFinder(n)
		c1_roots, c2_roots, c3_roots, sigma1, sigma2, sigma3, k1, k2 = permutations(n, roots, sigmaFor1, a, b, c)
		QL, QR, QM, QC, QO, PI = gatePolynomials(roots, ql, qr, qm, qc, qo, pi)
		Zh, S1, S2, S3, I1, I2, I3 = permutationPolynomial(roots, sigma1, sigma2, sigma3, k1, k2, c1_roots, c2_roots, c3_roots)
		
		# Generate the proof
		proof = prove(n, roots, a, b, c, Zh, QL, QR, QM, QC, QO, PI, S1, S2, S3, I1, I2, I3, k1, k2, omega)
		# Add transaction with proof to ledger
		transaction = {
			"tx_id": tx_id,
			"tx_type": tx_type,
			"details": details,
			"proof": proof,
			"timestamp": datetime.utcnow().isoformat()
		}
		global checker
		checker = True
		serialized_transaction = json_serialize(copy.deepcopy(transaction))
		print(f"Node {self.node_id} created transaction {tx_id} and broadcasting...")
		self.broadcast_transaction(serialized_transaction)
		time.sleep(5)
		if (checker):
			self.ledger.append(serialized_transaction)
	def verify_transaction(self, transaction):
		proof = transaction['proof']
		return verify(proof)

def create_nodes(num_nodes, crs, base_port=5000):
	nodes = []
	ports = [base_port + i for i in range(num_nodes)]
	for i in range(num_nodes):
		other_ports = ports[:i] + ports[i+1:]
		node = Node(node_id=i, crs=crs, port=ports[i], other_ports=other_ports)
		nodes.append(node)
	return nodes

def json_serialize(transaction):
    """Convert a transaction to a JSON-serializable format with Galois Field elements as integers."""
    serialized_transaction = {}

    for key, value in transaction.items():
        if key == "proof" and isinstance(value, dict):
            # Process the "proof" dictionary
            serialized_proof = {}
            for sub_key, sub_value in value.items():
                if isinstance(sub_value, list):
                    # Convert each Galois field element in the list to integer
                    serialized_proof[sub_key] = [
                        int(item) if isinstance(item, galois.FieldArray) else item
                        for item in sub_value
                    ]
                else:
                    # Keep other proof elements as they are
                    serialized_proof[sub_key] = sub_value
            serialized_transaction[key] = serialized_proof
        else:
            # Directly add other fields (non-proof) to the serialized transaction
            serialized_transaction[key] = value

    return serialized_transaction

def json_deserialize(serialized_transaction):
    """Convert serialized transaction integers back to Galois Field elements where applicable."""
    Fp = galois.GF(521)  # Define the Galois Field for order 521
    
    deserialized_transaction = {}
    
    for key, value in serialized_transaction.items():
        if key == "proof" and isinstance(value, dict):
            # Process the "proof" dictionary
            deserialized_proof = {}
            for sub_key, sub_value in value.items():
                if sub_key == "pi":
                    # Keep 'pi' as a list of integers
                    deserialized_proof[sub_key] = [int(item) for item in sub_value]
                elif sub_key in ["round1", "round2", "round3", "round4", "round5"]:
                    # Convert each integer in the list back to a Galois Field element for specified rounds
                    deserialized_proof[sub_key] = [Fp(item) for item in sub_value]
                else:
                    # Keep other proof elements as they are
                    deserialized_proof[sub_key] = sub_value
            deserialized_transaction[key] = deserialized_proof
        else:
            # Directly add other fields (non-proof) to the deserialized transaction
            deserialized_transaction[key] = value
            
    return deserialized_transaction

def aggregate_commitments(proofs, field):
    """Aggregate polynomial commitments from multiple proofs using weighted commitments."""
    # Initialize accumulators for each round's commitments in the specified field
    aggregated_round1 = [field(0), field(0), field(0)]  # For [A, B, C]
    aggregated_round2 = field(0)                       # For [Z]
    aggregated_round3 = [field(0), field(0), field(0)]  # For [Tl, Tm, Th]
    aggregated_round4 = [field(0)] * 6                 # For [a_zeta, b_zeta, c_zeta, s1_zeta, s2_zeta, z_omega_zeta]
    aggregated_round5 = [field(0), field(0)]           # For [Wzeta, Womega_zeta]
    aggregated_pi = [field(0)] * len(proofs[0]["pi"])  # Initialize aggregated `pi`

    # Loop through each proof and aggregate each round's values
    for proof in proofs:
        chi_i = generate_challenge(proof, field)  # Unique challenge for this proof
        
        # Aggregate round1 (A, B, C)
        for j in range(3):  # A, B, C components in round1
            aggregated_round1[j] += proof["round1"][j] * chi_i

        # Aggregate round2 (Z)
        aggregated_round2 += proof["round2"][0] * chi_i

        # Aggregate round3 (Tl, Tm, Th)
        for j in range(3):  # Tl, Tm, Th components in round3
            aggregated_round3[j] += proof["round3"][j] * chi_i

        # Aggregate round4 (a_zeta, b_zeta, c_zeta, s1_zeta, s2_zeta, z_omega_zeta)
        for j in range(6):  # Each element in round4
            aggregated_round4[j] += proof["round4"][j] * chi_i

        # Aggregate round5 (Wzeta, Womega_zeta)
        for j in range(2):  # Wzeta and Womega_zeta in round5
            aggregated_round5[j] += proof["round5"][j] * chi_i

        # Aggregate unique `pi` values
        for j, pi_value in enumerate(proof["pi"]):
            aggregated_pi[j] += pi_value * chi_i

    # Return the aggregated commitments, including `pi`
    return {
        "round1": aggregated_round1,
        "round2": [aggregated_round2],
        "round3": aggregated_round3,
        "round4": aggregated_round4,
        "round5": aggregated_round5,
        "pi": aggregated_pi
    }

#setup
G1 = generator1()
G2 = generator2()
patch_galois(galois.Poly)

encrypted = False

p = 521 if not encrypted else curve_order

Fp = galois.GF(p)

def setup1(organ_type, name, patient_id, blood_group, age, organ_life):
	
	# We have 7 gates, next power of 2 is 8
	n = 5
	n = 2 ** int(np.ceil(np.log2(n)))
	assert n & n - 1 == 0, "n must be a power of 2"
	print(n)
	# Find primitive root of unity
	omega = Fp.primitive_root_of_unity(n)
	assert omega ** (n) == 1, f"omega (ω) {omega} is not a root of unity"
	
	roots = Fp([omega**i for i in range(n)])
	print(f"roots = {roots}")
	return n, omega, roots

#witness and gates
def pad_array(a, n):
    return a + [0] * (n - len(a))

checker = True

def witnessGates1(organ_type, name, patient_id, blood_group, age, organ_life):
	#witness vectors
	c1 = organ_type + name
	c2 = patient_id + blood_group
	c3 = age + organ_life
	c4 = c1 + c2
	c5 = c3 + c4
	a = [organ_type, patient_id, age, c1, c3]
	b = [name, blood_group, organ_life, c2, c4]
	c = [c1, c2, c3, c4, c5]
	pi = [0, 0, 0, 0, -c5]

	# gate vectors
	ql = [1, 1, 1, 1, 1]
	qm = [0, 0, 0, 0, 0]
	qr = [1, 1, 1, 1, 1]
	qc = [0, 0, 0, 0, 0]
	qo = [-1, -1, -1, -1, 0]
	
	return a, b, c, pi, ql, qr, qm, qc, qo

def witnessGatesMain(n, organ_type, name, patient_id, blood_group, age, organ_life):
	a, b, c, pi, ql, qr, qm, qc, qo = witnessGates1(organ_type, name, patient_id, blood_group, age, organ_life)
	# pad vectors to length n
	a = pad_array(a, n)
	b = pad_array(b, n)
	c = pad_array(c, n)
	pi = pad_array(pi, n)
	ql = pad_array(ql, n)
	qr = pad_array(qr, n)
	qm = pad_array(qm, n)
	qc = pad_array(qc, n)
	qo = pad_array(qo, n)
	
	print(f"a = {a}")
	print(f"b = {b}")
	print(f"c = {c}")
	print(f"ql = {ql}")
	print(f"qr = {qr}")
	print(f"qm = {qm}")
	print(f"qc = {qc}")
	print(f"qo = {qo}")
	print(f"pi = {pi}")
	
	return a, b, c, pi, ql, qr, qm, qc, qo


#permutations
def print_sigma(sigma, a, b, c, r):
    group_size = len(sigma) // 3
    padding = 6

    print(f"{' w'} | {'value':{padding}} | {'i':{padding}} | {'sigma(i)':{padding}}")

    for i in range(0, group_size):
        print(f"a{i} | {a[i]:{padding}} | {r[i]:{padding}} | {r[sigma[i]]:{padding}}")

    print(f"-- | {'--':{padding}} | {'--':{padding}} | {'--':{padding}}")

    for i in range(group_size, 2 * group_size):
        print(
            f"b{i - group_size} | {b[i - group_size]:{padding}} | {r[i]:{padding}} | {r[sigma[i]]:{padding}}"
        )

    print(f"-- | {'--':{padding}} | {'--':{padding}} | {'--':{padding}}")

    for i in range(2 * group_size, 3 * group_size):
        print(
            f"c{i - 2 * group_size} | {c[i - 2 * group_size]:{padding}} | {r[i]:{padding}} | {r[sigma[i]]:{padding}}"
        )

def permutations(n, roots, sigma, a, b, c):
	k1 = 2
	k2 = 4
	c1_roots = roots
	c2_roots = roots * k1
	c3_roots = roots * k2
	
	c_roots = np.concatenate((c1_roots, c2_roots, c3_roots))
	
	check = set()
	for r in c_roots:
		assert not int(r) in check, f"Duplicate root {r} in {c_roots}"
		check.add(int(r))
	
	sigma1 = Fp([c_roots[sigma[i]] for i in range(0, n)])
	sigma2 = Fp([c_roots[sigma[i + n]] for i in range(0, n)])
	sigma3 = Fp([c_roots[sigma[i + 2 * n]] for i in range(0, n)])
	
	print_sigma(sigma, a, b, c, c_roots)
	
	print("\n\n--- Cosest ---")
	print(f"c0 = {c1_roots}")
	print(f"c1 = {c2_roots}")
	print(f"c2 = {c3_roots}")
	
	print("\n\n--- Sigma ---")
	print(f"sigma1 = {sigma1}")
	print(f"sigma2 = {sigma2}")
	print(f"sigma3 = {sigma3}")
	
	return c1_roots, c2_roots, c3_roots, sigma1, sigma2, sigma3, k1, k2

def sigmaFinder(n):
	ai = range(0, n)
	bi = range(n, 2 * n)
	ci = range(2 * n, 3 * n)
	sigma = {
		ai[0]: ai[0],
		ai[1]: ai[1],
		ai[2]: ai[2],
		ai[3]: ci[0],
		ai[4]: ci[2],
		ai[5]: ai[5],
		ai[6]: ai[6],
		ai[7]: ai[7],
		bi[0]: bi[0],
		bi[1]: bi[1],
		bi[2]: bi[2],
		bi[3]: ci[1],
		bi[4]: ci[3],
		bi[5]: bi[5],
		bi[6]: bi[6],
		bi[7]: bi[7],
		ci[0]: ai[3],
		ci[1]: bi[3],
		ci[2]: ai[4],
		ci[3]: bi[4],
		ci[4]: ci[4],
		ci[5]: ci[5],
		ci[6]: ci[6],
		ci[7]: ci[7],
	}
	
	return sigma



#gate polynomials
def to_galois_array(vector, field):
    # normalize to positive values
    a = [x % field.order for x in vector]
    return field(a)

def to_poly(x, v, field):
    assert len(x) == len(v)
    y = to_galois_array(v, field) if type(v) == list else v
    return galois.lagrange_poly(x, y)

def evaluate_poly(poly, roots):
    """Evaluate the polynomial at given points x."""
    return [int(poly(root)) for root in roots]

def gatePolynomials(roots, ql, qr, qm, qc, qo, pi):
	QL = to_poly(roots, ql, Fp)
	QR = to_poly(roots, qr, Fp)
	QM = to_poly(roots, qm, Fp)
	QC = to_poly(roots, qc, Fp)
	QO = to_poly(roots, qo, Fp)
	PI = to_poly(roots, pi, Fp)
	
	print("--- Gate Polynomials ---")
	print(f"QL = {QL}")
	print(f"QR = {QR}")
	print(f"QM = {QM}")
	print(f"QC = {QC}")
	print(f"QO = {QO}")
	print(f"PI = {PI}")
	
	return QL, QR, QM, QC, QO, PI


#permutation polynomial
def to_vanishing_poly(roots, field):
    # Z^n - 1 = (Z - 1)(Z - w)(Z - w^2)...(Z - w^(n-1))
    return galois.Poly.Degrees([len(roots), 0], coeffs=[1, -1], field=field)

def permutationPolynomial(roots, sigma1, sigma2, sigma3, k1, k2, c1_roots, c2_roots, c3_roots):
	S1 = to_poly(roots, sigma1, Fp)
	S2 = to_poly(roots, sigma2, Fp)
	S3 = to_poly(roots, sigma3, Fp)
	
	I1 = to_poly(roots, c1_roots, Fp)
	I2 = to_poly(roots, c2_roots, Fp)
	I3 = to_poly(roots, c3_roots, Fp)
	
	padding = 3
	for i in range(0, len(roots)):
	    s = f"i = {i:{padding}} --> {roots[i]:{padding}} "
	    s += f"  I1({roots[i]:{padding}}) = {I1(roots[i]):{padding}} "
	    s += f"  I2({roots[i]:{padding}}) = {I2(roots[i]):{padding}} "
	    s += f"  I3({roots[i]:{padding}}) = {I3(roots[i]):{padding}} "
	    s += f"  S1({roots[i]:{padding}}) = {S1(roots[i]):{padding}} "
	    s += f"  S2({roots[i]:{padding}}) = {S2(roots[i]):{padding}} "
	    s += f"  S3({roots[i]:{padding}}) = {S3(roots[i]):{padding}} "
	    print(s)
	    
	    assert I1(roots[i]) == roots[i], f"I1({roots[i]}) != {roots[i]}"
	    assert I2(roots[i]) == k1 * roots[i], f"I2({roots[i]}) != {k1 * roots[i]}"
	    assert I3(roots[i]) == k2 * roots[i], f"I3({roots[i]}) != {k2 * roots[i]}"
	    
	    assert S1(roots[i]) == sigma1[i], f"S1({roots[i]}) != {sigma1[i]}"
	    assert S2(roots[i]) == sigma2[i], f"S2({roots[i]}) != {sigma2[i]}"
	    assert S3(roots[i]) == sigma3[i], f"S3({roots[i]}) != {sigma3[i]}"
	Zh = to_vanishing_poly(roots, Fp)
	for x in roots:
	    assert Zh(x) == 0
	print("--- Vanishing Polynomial ---")
	print(f"Zh = {Zh}")
	
	return Zh, S1, S2, S3, I1, I2, I3


#CRS construction
def generate_tau(encrypted=False):
    return SRS(Fp.Random(), n) if encrypted else Fp.Random()

def crsConstruct():
	tau = generate_tau(encrypted=encrypted)
	print(f"--- Tau ---")
	print(tau)
	return tau

tau = crsConstruct()

#prover
def shift_poly(poly: galois.Poly, omega: Fp):
    coeffs = poly.coeffs[::-1]
    coeffs = [c * omega**i for i, c in enumerate(coeffs)]
    return galois.Poly(coeffs[::-1], field=poly.field)

def prove(n, roots, a, b, c, Zh, QL, QR, QM, QC, QO, PI, S1, S2, S3, I1, I2, I3, k1, k2, omega):
	pi = evaluate_poly(PI, roots)
	#round1
	random_b = [Fp.Random() for i in range(0, 9)]
	
	bA = galois.Poly(random_b[:2], field=Fp)
	bB = galois.Poly(random_b[2:4], field=Fp)
	bC = galois.Poly(random_b[4:6], field=Fp)
	
	_A = to_poly(roots, a, Fp)
	_B = to_poly(roots, b, Fp)
	_C = to_poly(roots, c, Fp)
	
	A = _A + bA * Zh
	B = _B + bB * Zh
	C = _C + bC * Zh
	
	# gate constraints polynomial
	# g(x) = a(x)*ql(x) + b(x)*qr(x) + a(x)*b(x)*qm(x) + c(x)*qo(x) + qpi(x) + qc(x)
	G = A * QL + B * QR + A * B * QM + C * QO + QC + PI
	
	print("--- Gate Constraints Polynomial ---")
	print(f"G = {G}")
	for i in range(0, len(roots)):
	    print(
	        f"gate #{i} G({roots[i]}) = {G(roots[i])} --> {'OK' if G(roots[i]) == 0 else 'FAIL'}"
	    )
	    assert G(roots[i]) == 0, f"G({roots[i]}) != 0"
	
	assert G % Zh == 0, f"G(x) % Zh(x) != 0"
	
	padding = 3
	for i in range(0, len(roots)):
	    s = f"i = {i:{padding}} --> {roots[i]:{padding}} "
	    s += f"   A({roots[i]:{padding}}) = {A(roots[i]):{padding}} "
	    s += f"   B({roots[i]:{padding}}) = {B(roots[i]):{padding}} "
	    s += f"   C({roots[i]:{padding}}) = {C(roots[i]):{padding}} "
	    print(s)
	
	round1 = [A(tau), B(tau), C(tau)]
	print("\n\n--- Round 1 ---")
	print(f"Round 1 = {round1}")
	
	
	#round2
	beta = numbers_to_hash(round1 + [0], Fp)
	gamma = numbers_to_hash(round1 + [1], Fp)
	
	_F = (A + I1 * beta + gamma) * (B + I2 * beta + gamma) * (C + I3 * beta + gamma)
	_G = (A + S1 * beta + gamma) * (B + S2 * beta + gamma) * (C + S3 * beta + gamma)
	
	acc_eval = [Fp(1)]
	for i in range(0, n):
		acc_eval.append(acc_eval[-1] * (_F(roots[i]) / _G(roots[i])))
	assert acc_eval.pop() == Fp(1)
	ACC = galois.lagrange_poly(roots, Fp(acc_eval))
	print("\n\n--- Accumulator Polynomial ---")
	print(f"ACC(x) = {ACC}")
	
	bZ = galois.Poly(random_b[6:9], field=Fp)
	print(f"bZ = {bZ}")
	Z = bZ * Zh + ACC
	
	print("\n\n--- Z Polynomial ---")
	print(f"Z(x) = {Z}")
	
	for r in roots:
	    print(f"Z({r}) = {Z(r)}")
	
	assert Z(roots[0]) == 1
	assert Z(roots[-1]) == 1
	
	round2 = [Z(tau)]
	print("\n\n--- Round 2 ---")
	print(f"Round 2 = {round2}")
	
	
	#round3
	alpha = numbers_to_hash(round1 + round2, Fp)
	
	Zomega = shift_poly(Z, omega=omega)
	for r in roots:
	    print(f"Z({r:3}) = {Z(r):3} -> Zω({r:3}) = {Zomega(r):3}")
	
	print("\n\n--- Zω Polynomial ---")
	print(f"Z(ωx) = {Zomega}")
	
	L1 = galois.lagrange_poly(roots, Fp([1] + [Fp(0)] * (n - 1)))
	
	print("\n\n--- L1 Polynomial ---")
	print(f"L1(x) = {L1}")
	for i, r in enumerate(roots):
	    print(f"L1({r}) = {L1(r)}")
	    assert L1(r) == (Fp(1) if i == 0 else Fp(0))
	
	T0 = G
	assert T0 % Zh == 0, f"T0(x) % Zh(x) != 0"
	
	T1 = (_F * Z - _G * Zomega) * alpha
	assert T1 % Zh == 0, f"T1(x) % Zh(x) != 0"
	
	T2 = (Z - galois.Poly([1], field=Fp)) * L1 * alpha**2
	assert T2 % Zh == 0, f"T2(x) % Zh(x) != 0"
	
	T = T0 + T1 + T2
	assert T % Zh == 0, f"T(x) % Zh(x) != 0"
	
	for r in roots:
	    assert T(r) == 0, f"T({r}) != 0"
	
	T = T // Zh
	
	print("\n\n--- T Polynomial ---")
	print(f"T(x) = {T}")
	
	t_coeffs = T.coeffs[::-1]
	
	Tl = galois.Poly(t_coeffs[:n][::-1], field=Fp)
	Tm = galois.Poly(t_coeffs[n : 2 * (n)][::-1], field=Fp)
	Th = galois.Poly(t_coeffs[2 * (n) :][::-1], field=Fp)
	
	X_n = galois.Poly.Degrees([n, 0], coeffs=[1, 0], field=Fp)
	X_2n = galois.Poly.Degrees([2 * (n), 0], coeffs=[1, 0], field=Fp)
	# make sure that T was split correctly
	# T = TL + X^n * TM + X^2n * TH
	assert T == (Tl + X_n * Tm + X_2n * Th)
	assert T.degree == 3 * n + 5
	
	b10 = Fp.Random()
	b11 = Fp.Random()
	
	Tl = Tl + b10 * X_n
	Tm = Tm - b10 + b11 * X_n
	Th = Th - b11
	assert T == (Tl + X_n * Tm + X_2n * Th)
	
	print("\n\n--- T' ---")
	print(f"Tl(x) = {Tl}")
	print(f"Tm(x) = {Tm}")
	print(f"Th(x) = {Th}")
	
	round3 = [Tl(tau), Tm(tau), Th(tau)]
	print("\n\n--- Round 3 ---")
	print(f"Round 3 = {round3}")
	
	
	#round4
	zeta = numbers_to_hash(round1 + round2 + round3, Fp)

	a_zeta = A(zeta)
	b_zeta = B(zeta)
	c_zeta = C(zeta)
	s1_zeta = S1(zeta)
	s2_zeta = S2(zeta)
	z_omega_zeta = Zomega(zeta)

	round4 = [a_zeta, b_zeta, c_zeta, s1_zeta, s2_zeta, z_omega_zeta]
	print("\n\n--- Round 4 ---")
	print(f"Round 4 = {round4}")
	
	
	#round5
	v = numbers_to_hash(round1 + round2 + round3 + round4, Fp)
	
	pi_zeta = PI(zeta)
	
	R = QM * a_zeta * b_zeta + QL * a_zeta + QR * b_zeta + QO * c_zeta + QC + pi_zeta
	R += (
	    Z
	    * (a_zeta + beta * zeta + gamma)
	    * (b_zeta + beta * zeta * k1 + gamma)
	    * (c_zeta + beta * zeta * k2 + gamma)
	    * alpha
	)
	R -= (
	    z_omega_zeta
	    * (a_zeta + beta * s1_zeta + gamma)
	    * (b_zeta + beta * s2_zeta + gamma)
	    * (c_zeta + beta * S3 + gamma)
	    * alpha
	)
	R += (Z - Fp(1)) * L1(zeta) * alpha**2
	R -= Zh(zeta) * (Tl + zeta**n * Tm + zeta ** (2 * n) * Th)

	print("\n\n--- R ---")
	print(f"R(x) = {R}")
	
	X_minus_zeta = galois.Poly([1, -zeta], field=Fp)
	print(f"X - zeta = {X_minus_zeta}")
	
	Wzeta = (
	    R
	    + (A - a_zeta) * v
	    + (B - b_zeta) * v**2
	    + (C - c_zeta) * v**3
	    + (S1 - s1_zeta) * v**4
	    + (S2 - s2_zeta) * v**5
	)
	
	assert Wzeta % X_minus_zeta == 0, f"Wzeta(x) % X - zeta != 0"
	Wzeta = Wzeta // X_minus_zeta
	
	X_minus_omega_zeta = galois.Poly([1, -(omega * zeta)], field=Fp)
	print(f"X - ω*zeta = {X_minus_omega_zeta}")
	
	Womega_zeta = Z - z_omega_zeta
	assert Womega_zeta % X_minus_omega_zeta == 0, f"Womega_zeta(x) % X - ω*zeta != 0"
	Womega_zeta = Womega_zeta // X_minus_omega_zeta
	
	round5 = [Wzeta(tau), Womega_zeta(tau)]
	print("\n\n--- Round 5 ---")
	print(f"Round 5 = {round5}")

	u = numbers_to_hash(round1 + round2 + round3 + round4 + round5, Fp)
	proof = {
		"Fp": p,
		"n": n,
		"pi": pi,
		"round1": round1,
		"round2": round2,
		"round3": round3,
		"round4": round4,
		"round5": round5
	}
	#proof["round1"][0] = Fp.Random()
	pr = dump_proof(proof, "proof.json")
	
	circuit = {
	"add_patient": {
		"QM": QM.coeffs,
		"QL": QL.coeffs,
		"QR": QR.coeffs,
		"QO": QO.coeffs,
		"QC": QC.coeffs,
		"PI": PI.coeffs,
		"Zh": Zh.coeffs,
		"L1": L1.coeffs,
		"S1": S1.coeffs,
		"S2": S2.coeffs,
		"S3": S3.coeffs,
		"k1": k1,
		"k2": k2,
		"tau": tau.tau if encrypted else tau,
		"Fp": p,
		"omega": omega,
		"n": n,
		"encrypted": encrypted,
		}
	}

	dump_circuit(circuit, "circuit.json")
	
	return proof




#verifier
def verify(proof):
	# These evaluations are calculated beforehand during the setup phase
	circuits = load_circuit("circuit.json")
	circuit_config = circuits["add_patient"]
	
	QM = galois.Poly(circuit_config["QM"], field=Fp)
	QL = galois.Poly(circuit_config["QL"], field=Fp)
	QR = galois.Poly(circuit_config["QR"], field=Fp)
	QO = galois.Poly(circuit_config["QO"], field=Fp)
	QC = galois.Poly(circuit_config["QC"], field=Fp)
	S1 = galois.Poly(circuit_config["S1"], field=Fp)
	S2 = galois.Poly(circuit_config["S2"], field=Fp)
	S3 = galois.Poly(circuit_config["S3"], field=Fp)
	
	qm_exp = QM(tau)
	ql_exp = QL(tau)
	qr_exp = QR(tau)
	qo_exp = QO(tau)
	qc_exp = QC(tau)
	s1_exp = S1(tau)
	s2_exp = S2(tau)
	s3_exp = S3(tau)
	
	k1 = circuit_config["k1"]
	k2 = circuit_config["k2"]
	
	n = circuit_config["n"]
	
	omega = circuit_config["omega"]
	
	roots = Fp([omega**i for i in range(n)])
	
	L1 = galois.Poly(circuit_config["L1"], field=Fp)
	Zh = galois.Poly(circuit_config["Zh"], field=Fp)
	
	round1 = proof.get("round1")
	round2 = proof.get("round2")
	round3 = proof.get("round3")
	round4 = proof.get("round4")
	round5 = proof.get("round5")
	
	pi = proof.get("pi")
	pi = pad_array(pi, n)
	PI = to_poly(roots, pi, Fp)
	
	# Values provided by the prover (round 1 to 5) is a proof.
	a_exp = round1[0]
	b_exp = round1[1]
	c_exp = round1[2]

	z_exp = round2[0]

	tl_exp = round3[0]
	tm_exp = round3[1]
	th_exp = round3[2]

	# Note: verifier has to verify that the following values are in the correct Fp field
	a_zeta, b_zeta, c_zeta, s1_zeta, s2_zeta, z_omega_zeta = round4

	w_zeta_exp = round5[0]
	w_omega_zeta_exp = round5[1]

	# Note: verifier has to verify that the following values are on the curve
	if encrypted:
		validate_point(qm_exp)
		validate_point(ql_exp)
		validate_point(qr_exp)
		validate_point(qo_exp)
		validate_point(qc_exp)
		validate_point(z_exp)
		validate_point(s1_exp)
		validate_point(s2_exp)
		validate_point(s3_exp)
		validate_point(tl_exp)
		validate_point(tm_exp)
		validate_point(th_exp)
		validate_point(a_exp)
		validate_point(b_exp)
		validate_point(c_exp)
		validate_point(w_zeta_exp)
		validate_point(w_omega_zeta_exp)

	beta = numbers_to_hash(round1 + [0], Fp)
	gamma = numbers_to_hash(round1 + [1], Fp)
	alpha = numbers_to_hash(round1 + round2, Fp)
	zeta = numbers_to_hash(round1 + round2 + round3, Fp)
	v = numbers_to_hash(round1 + round2 + round3 + round4, Fp)
	u = numbers_to_hash(round1 + round2 + round3 + round4 + round5, Fp)
	Zh_z = Zh(zeta)
	L1_z = L1(zeta)
	PI_z = PI(zeta)
	
	r0 = (
		PI_z
		- L1_z * alpha**2
		- (a_zeta + beta * s1_zeta + gamma)
		* (b_zeta + beta * s2_zeta + gamma)
		* (c_zeta + gamma)
		* z_omega_zeta
		* alpha
	)

	D_exp = (
		qm_exp * a_zeta * b_zeta
		+ ql_exp * a_zeta
		+ qr_exp * b_zeta
		+ qo_exp * c_zeta
		+ qc_exp
	)

	D_exp += z_exp * (
		(a_zeta + beta * zeta + gamma)
		* (b_zeta + beta * zeta * k1 + gamma)
		* (c_zeta + beta * zeta * k2 + gamma)
		* alpha
		+ L1_z * alpha**2
		+ u
	)

	D_exp -= (
		s3_exp
		* (a_zeta + beta * s1_zeta + gamma)
		* (b_zeta + beta * s2_zeta + gamma)
		* alpha
		* beta
		* z_omega_zeta
	)

	D_exp -= (tl_exp + tm_exp * zeta**n + th_exp * zeta ** (2 * n)) * Zh_z
	F_exp = (
		D_exp
		+ a_exp * v
		+ b_exp * v**2
		+ c_exp * v**3
		+ s1_exp * v**4
		+ s2_exp * v**5
	)

	E_exp = (
		-r0
		+ v * a_zeta
		+ v**2 * b_zeta
		+ v**3 * c_zeta
		+ v**4 * s1_zeta
		+ v**5 * s2_zeta
		+ u * z_omega_zeta
	)

	if encrypted:
		E_exp = G1 * E_exp

	e1 = w_zeta_exp + w_omega_zeta_exp * u
	e2 = (
		w_zeta_exp * zeta
		+ w_omega_zeta_exp * (u * zeta * omega)
		+ F_exp
		+ (E_exp * Fp(p - 1))
	)

	if encrypted:
		pairing1 = tau.tau2.pair(e1)
		pairing2 = G2.pair(e2)

		print(f"pairing1 = {pairing1}")
		print(f"pairing2 = {pairing2}")

		return (pairing1 == pairing2)
	else:
		print("\n\n--- e1, e2 ---")
		print(f"e1 = {e1 * tau} = {e1} * tau")
		print(f"e2 = {e2}")
		return (e1 * tau == e2)




#working
print("Welcome to Organ transplantation system simulation")
nodes = create_nodes(3, tau, 7000)
#nodes[0].create_transaction("TX1001", "add_patient", {"patient_id": 123})

time.sleep(3)

c = "1"
while (c!="0"):
	print("Enter details for organ entry")
	organ_type = input("Enter organ type: ")
	name = input("Enter patient name: ")
	patient_id = input("Enter patient id: ")
	blood_group = input("Enter blood group: ")
	age = input("Enter patient age: ")
	organ_life = input("Enter organ life: ")
	nodes[0].create_transaction("TX100"+str(len(nodes[0].ledger)+1), "add_patient", {"organ_type":organ_type, "name":name, "patient_id": patient_id, "blood_group":blood_group, "age":age, "organ_life":organ_life})
	time.sleep(2)
	print("Enter 0 to exit")
	c = input()

print("Prover ledger transactions: ")
for i in range(len(nodes[0].ledger)):
	print(f"[{i}]: ", nodes[0].ledger[i])
print("\nVerifier ledger: ")
for i in range(len(nodes[1].ledger)):
	print(f"[{i}]: ", nodes[1].ledger[i])
