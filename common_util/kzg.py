import py_ecc.bn128 as b
from common_util.curve import ec_lincomb, G1Point, G2Point, Scalar
from dataclasses import dataclass
from common_util.poly import Polynomial, Basis

# Recover the trusted setup from a file in the format used in
# https://github.com/iden3/snarkjs#7-prepare-phase-2
SETUP_FILE_G1_STARTPOS = 80
SETUP_FILE_POWERS_POS = 60


@dataclass
class Setup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( H,    xH,  ...,  x^{d-1}H ), where H is a generator of G_2
    powers_of_X2: list[G2Point]
    N: int

    @classmethod
    def from_file(cls, filename):
        contents = open(filename, "rb").read()
        # Byte 60 gives you the base-2 log of how many powers there are
        powers = 2 ** contents[SETUP_FILE_POWERS_POS]
        # Extract G1 points, which start at byte 80
        values = [
            int.from_bytes(contents[i : i + 32], "little")
            for i in range(
                SETUP_FILE_G1_STARTPOS, SETUP_FILE_G1_STARTPOS + 32 * powers * 2, 32
            )
        ]
        assert max(values) < b.field_modulus
        # The points are encoded in a weird encoding, where all x and y points
        # are multiplied by a factor (for montgomery optimization?). We can
        # extract the factor because we know the first point is the generator.
        factor = b.FQ(values[0]) / b.G1[0]
        values = [b.FQ(x) / factor for x in values]
        powers_of_x = [(values[i * 2], values[i * 2 + 1]) for i in range(powers)]
        print("Extracted G1 side, X^1 point: {}".format(powers_of_x[1]))
        # Search for start of G2 points. We again know that the first point is
        # the generator.
        pos = SETUP_FILE_G1_STARTPOS + 32 * powers * 2
        target = (factor * b.G2[0].coeffs[0]).n
        while pos < len(contents):
            v = int.from_bytes(contents[pos : pos + 32], "little")
            if v == target:
                break
            pos += 1
        # print("Detected start of G2 side at byte {}".format(pos))
        X2_encodings = [
            e
            for i in range(pos, pos + 32 * powers * 4, 128)
            for e in contents[i : i + 32 * 4]
        ]
        X2_values = [
            b.FQ(int.from_bytes(X2_encodings[i : i + 32], "little")) / factor
            for i in range(0, 32 * powers * 4, 32)
        ]
        powers_of_x2 = [
            (b.FQ2(X2_values[i * 2 : i * 2 + 2]), b.FQ2(X2_values[i * 2 + 2 : i * 2 + 4])) for i in range(0, powers * 2, 2)
        ]
        N = len(powers_of_x)
        print("Extracted G2 side, X^1 point: {}".format(powers_of_x2[1]))
        assert b.pairing(b.G2, powers_of_x[1]) == b.pairing(powers_of_x2[1], b.G1)
        # assert b.pairing(b.G2, powers_of_x[2]) == b.pairing(X2_all[2], b.G1)
        # print("X^1 points checked consistent")
        return cls(powers_of_x, powers_of_x2, N)

    # Encodes the KZG commitment that evaluates to the given values in the group
    def commit(self, values: Polynomial) -> G1Point:
        assert values.basis == Basis.LAGRANGE

        # inverse FFT from Lagrange basis to monomial basis
        coeffs = values.ifft().values
        if len(coeffs) > len(self.powers_of_x):
            raise Exception("Not enough powers in setup")
        return ec_lincomb([(s, x) for s, x in zip(self.powers_of_x, coeffs)])
    
if __name__ == "__main__":
    setup = Setup.from_file("powersOfTau28_hez_final_11.ptau")
    dummy_values = Polynomial(
        list(map(Scalar, [1, 2, 3, 4, 5, 6, 7, 8])), Basis.LAGRANGE
    )
    commitment = setup.commit(dummy_values)
    assert commitment == G1Point(
        (
            16120260411117808045030798560855586501988622612038310041007562782458075125622,
            3125847109934958347271782137825877642397632921923926105820408033549219695465,
        )
    )
    print("Pass Setup Test")