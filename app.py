#setup
import galois
import numpy as np
from utils import (
    generator1,
    generator2,
    curve_order,
    normalize,
    validate_point,
    GPoint,
    SRS,
    numbers_to_hash,
    patch_galois,
    dump_proof,
    dump_circuit
)

G1 = generator1()
G2 = generator2()

patch_galois(galois.Poly)


# initial values for x and y
x = 2
y = 3

# to switch between encrypted (ECC) and unencrypted mode (prime field p=241)
encrypted = False

# Prime field p
p = 241 if not encrypted else curve_order
# p = 241
Fp = galois.GF(p)

# 2x^2 - x^2y^2 + 3
out = 2 * x**2 - x**2 * y**2 + 3
print(f"out = {out}")

# We have 7 gates, next power of 2 is 8
n = 7
n = 2 ** int(np.ceil(np.log2(n)))
assert n & n - 1 == 0, "n must be a power of 2"

# Find primitive root of unity
omega = Fp.primitive_root_of_unity(n)
assert omega ** (n) == 1, f"omega (ω) {omega} is not a root of unity"

roots = Fp([omega**i for i in range(n)])
print(f"roots = {roots}")



#Witness and gates
def pad_array(a, n):
    return a + [0] * (n - len(a))


# witness vectors
a = [2, 2, 3, 4, 4, 8, -28]
b = [2, 2, 3, 0, 9, 36, 3]
c = [4, 4, 9, 8, 36, -28, -25]
pi = [0, 0, 0, 0, 0, 0, 25]

# gate vectors
ql = [0, 0, 0, 2, 0, 1, 1]
qr = [0, 0, 0, 0, 0, -1, 0]
qm = [1, 1, 1, 1, 1, 0, 0]
qc = [0, 0, 0, 0, 0, 0, 3]
qo = [-1, -1, -1, -1, -1, -1, 0]

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


ai = range(0, n)
bi = range(n, 2 * n)
ci = range(2 * n, 3 * n)

sigma = {
    ai[0]: ai[0],
    ai[1]: ai[1],
    ai[2]: ai[2],
    ai[3]: ci[0],
    ai[4]: ci[1],
    ai[5]: ci[3],
    ai[6]: ci[5],
    ai[7]: ai[7],
    bi[0]: bi[0],
    bi[1]: bi[1],
    bi[2]: bi[2],
    bi[3]: bi[3],
    bi[4]: ci[2],
    bi[5]: ci[4],
    bi[6]: bi[6],
    bi[7]: bi[7],
    ci[0]: ai[3],
    ci[1]: ai[4],
    ci[2]: bi[4],
    ci[3]: ai[5],
    ci[4]: bi[5],
    ci[5]: ai[6],
    ci[6]: ci[6],
    ci[7]: ci[7],
}

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



#gate polynomials
def to_galois_array(vector, field):
    # normalize to positive values
    a = [x % field.order for x in vector]
    return field(a)


def to_poly(x, v, field):
    assert len(x) == len(v)
    y = to_galois_array(v, field) if type(v) == list else v
    return galois.lagrange_poly(x, y)


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



#permutations polynomials
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
   
def to_vanishing_poly(roots, field):
    # Z^n - 1 = (Z - 1)(Z - w)(Z - w^2)...(Z - w^(n-1))
    return galois.Poly.Degrees([len(roots), 0], coeffs=[1, -1], field=field)


Zh = to_vanishing_poly(roots, Fp)
for x in roots:
    assert Zh(x) == 0

print("--- Vanishing Polynomial ---")
print(f"Zh = {Zh}")



#CRS construction
def generate_tau(encrypted=False):
    return SRS(Fp.Random(), n) if encrypted else Fp.Random()


tau = generate_tau(encrypted=encrypted)
print(f"--- Tau ---")
print(tau)






#prover

#round 1
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



#round 2
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



#round 3
def shift_poly(poly: galois.Poly, omega: Fp):
    coeffs = poly.coeffs[::-1]
    coeffs = [c * omega**i for i, c in enumerate(coeffs)]
    return galois.Poly(coeffs[::-1], field=poly.field)


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



#round 4
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



#round 5
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
    "A": round1[0],
    "B": round1[1],
    "C": round1[2],
    "Z": round2[0],
    "Tl": round3[0],
    "Tm": round3[1],
    "Th": round3[2],
    "Wzeta": round5[0],
    "Womega_zeta": round5[1],
    "a_zeta": round4[0],
    "b_zeta": round4[1],
    "c_zeta": round4[2],
    "s1_zeta": round4[3],
    "s2_zeta": round4[4],
    "z_omega_zeta": round4[5],
    "Fp": p,
}

dump_proof(proof, "proof.json")

circuit = {
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

dump_circuit(circuit, "circuit.json")





#verifier
# These evaluations are calculated beforehand during the setup phase
qm_exp = QM(tau)
ql_exp = QL(tau)
qr_exp = QR(tau)
qo_exp = QO(tau)
qc_exp = QC(tau)
s1_exp = S1(tau)
s2_exp = S2(tau)
s3_exp = S3(tau)

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

print(f"F_exp = {F_exp}")

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

print(f"E_exp = {E_exp}")
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

    assert pairing1 == pairing2, f"pairing1 != pairing2"
else:
    print("\n\n--- e1, e2 ---")
    print(f"e1 = {e1 * tau} = {e1} * tau")
    print(f"e2 = {e2}")
    assert e1 * tau == e2

