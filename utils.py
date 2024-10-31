"""
Utilities for the Plonk protocol.
"""
import json
from Crypto.Hash import keccak
import galois
from py_ecc.optimized_bn128 import (
    add,
    multiply,
    G1,
    G2,
    neg,
    pairing,
    eq,
    normalize,
    FQ,
    FQ2,
    curve_order,
    is_on_curve,
)

__all__ = [
    "GPoint",
    "generator1",
    "generator2",
    "validate_point",
    "normalize",
    "curve_order",
    "SRS",
    "numbers_to_hash",
    "patch_galois",
    "dump_proof",
    "dump_circuit",
]

GF241 = galois.GF(241)

class GPoint(tuple):

    """
    A point on the BN128 curve.
    This class is a wrapper G1 and G2 points to provide a more intuitive
    interface. For example, instead of writing `multiply(G1, 5)` you can
    write `G1 * 5` or `5 * G1`. Similarly, instead of writing `add(G1, G2)`
    you can write `G1 + G2`.
    """

    def __new__(cls, x, y, z):
        return tuple.__new__(cls, (x, y, z))

    def __str__(self):
        return f"({self[0]}, {self[1]}, {self[2]})"

    def __add__(self, other):
        return GPoint(*add(self, other))

    def __iadd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return GPoint(*add(self, neg(other)))

    def __isub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        return GPoint(*multiply(self, int(other)))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __neg__(self):
        return GPoint(*neg(self))

    def __eq__(self, other):
        return eq(self, other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def pair(self, other):
        """Pairing function."""

        return pairing(self, other)


def generator1():
    """Generator for G1."""
    return GPoint(*G1)


def generator2():
    """Generator for G2."""
    return GPoint(*G2)


def validate_point(pt):
    """
    Check if a point is on the curve.
    Used in Plonk's verifier. Weak curve attack mitigation.
    """
    if isinstance(pt[0], FQ):
        assert is_on_curve(pt, FQ(3))
    elif isinstance(pt[0], FQ2):
        assert is_on_curve(pt, FQ2([3, 0]) / FQ2([9, 1]))
    else:
        raise Exception("Invalid point")


class SRS:
    """Trusted Setup Class aka Structured Reference String"""

    def __init__(self, tau, n=2):
        self.tau = tau
        g1 = generator1()
        g2 = generator2()
        self.tau1 = [g1 * int(tau) ** i for i in range(0, n + 7)]
        self.tau2 = g2 * int(tau)

    def __str__(self):
        s = f"tau: {self.tau}\n"
        s += "".join(
            [
                f"[tau^{i}]G1: {str(normalize(point))}\n"
                for i, point in enumerate(self.tau1)
            ]
        )
        s += f"[tau]G2: {str(normalize(self.tau2))}\n"
        return s


def numbers_to_hash(numbers, field) -> int:
    """Hash a number."""
    engine = keccak.new(digest_bits=256)
    for number in numbers:
        if isinstance(number, tuple):
            x, y, z = number
            engine.update(bytes(hex(int(x)), "utf-8"))
            engine.update(bytes(hex(int(y)), "utf-8"))
            engine.update(bytes(hex(int(z)), "utf-8"))
        else:
            engine.update(bytes(hex(int(number)), "utf-8"))
    return field(int(engine.hexdigest(), 16) % field.order)


def patch_galois(Poly):
    def new_call(self, at, **kwargs):
        if isinstance(at, SRS):
            coeffs = self.coeffs[::-1]
            result = at.tau1[0] * coeffs[0]
            for i in range(1, len(coeffs)):
                result += at.tau1[i] * coeffs[i]
            return result

        return Poly.original_call(self, at, **kwargs)

    Poly.original_call = Poly.__call__
    Poly.__call__ = new_call


def dump_proof(proof, path):
    """Dump proof to file."""
    for k, v in proof.items():
        if isinstance(v, GPoint):  # Assuming GPoint is defined somewhere in your code
            proof[k] = str(v)
        elif isinstance(v, list):
            # Convert each item in the list, keeping GF(241) elements as is
            proof[k] = [
                str(item) if isinstance(item, GPoint) else item for item in v
            ]
        else:
            # Check for GF(241) type
            if isinstance(v, GF241):  # Check against the specific Galois Field instance
                proof[k] = v  # Keep it as is
            else:
                proof[k] = int(v)  # Convert other types to int

    # Saving the proof dictionary to a JSON file
    with open(path, "w") as f:
        json.dump(proof, f, indent=2, default=str)  # Use default=str to handle GF elements

def load_proof(path):
    """Load proof from file."""
    with open(path, "r") as f:
        proof = json.load(f)
    
    # List of proof keys that should be converted to GF(241)
    gf_keys = ["round1", "round2", "round3", "round4", "round5"]
    
    for k, v in proof.items():
        if isinstance(v, list):
            # Convert each item in the list back to Galois Field or other types as needed
            proof[k] = [
                galois.GF(241)(int(item)) if k in gf_keys and isinstance(item, (int, str)) else item
                for item in v
            ]
        else:
            # Convert back if the value was an int or str and should be GF(241)
            if isinstance(v, (int, str)) and k in gf_keys:
                proof[k] = galois.GF(241)(int(v))  # Convert str to int, then to GF element
    
    return proof

def dump_circuit(circuit, path):
    """Dump nested circuit to file, converting Galois field arrays and polynomials to integer lists."""
    for function_name, circuit_data in circuit.items():
        for k, v in circuit_data.items():
            if k in ["tau", "k1", "k2", "Fp", "omega", "n"]:
                circuit_data[k] = int(v)  # convert these specific keys to integers
            elif isinstance(v, bool):
                continue  # leave booleans as they are
            elif isinstance(v, galois.FieldArray):
                circuit_data[k] = v.tolist()  # store FieldArray as a list of integers
            elif isinstance(v, galois.Poly):
                circuit_data[k] = v.coeffs.tolist()  # store Poly coefficients as a list of integers
            elif isinstance(v, list):
                circuit_data[k] = [int(x) for x in v]  # convert lists of coefficients to integers
    with open(path, "w") as f:
        json.dump(circuit, f, indent=2)

def load_circuit(path):
    """Load circuit configuration from file, with lists of integers for polynomial coefficients."""
    with open(path, "r") as f:
        circuit = json.load(f)

    for function_name, circuit_data in circuit.items():
        # Initialize the Galois field based on the stored 'Fp' field value
        Fp = galois.GF(circuit_data["Fp"]) if "Fp" in circuit_data else None

        for k, v in circuit_data.items():
            if k in ["tau", "k1", "k2", "Fp", "n"]:
                circuit_data[k] = int(v)  # Restore integers directly
            elif k == "omega":
            	circuit_data[k] = Fp(int(v))
            elif isinstance(v, list):
                circuit_data[k] = v  # Keep lists as they are
            else:
                circuit_data[k] = v  # Leave other types unchanged

    return circuit
