from dataclasses import dataclass

from KZG10 import TrustedSetup, GF, curve
from src.common_util.kzg_optimized import Setup as KZGSetup


@dataclass
class Setup:
    F: GF
    srs: TrustedSetup

    def __init__(self):
        self.F = GF(curve.curve_order)
        self.srs = TrustedSetup.generate(self.F, 5, True)

if __name__ == "__main__":
    kzgSetup = KZGSetup.manual_setup()
