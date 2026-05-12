import unittest

from spiceypy import ckobj, ckcov
from spiceypy.utils.support_types import SPICEDOUBLE_CELL


class TestCkcov(unittest.TestCase):
    def test_something(self):
        ScanPlatformCK="data/spice/ck/vg1_jup_version1_type1_iss_sedr.bc"
        ckids = ckobj(ScanPlatformCK)
        print(ckids[0])
        cover = SPICEDOUBLE_CELL(200000)
        cover = ckcov(ScanPlatformCK, ckids[0], False, "INTERVAL", 0.0, "SCLK", cover)
        # If the code above doesn't raise an exception, the test passes

if __name__ == '__main__':
    unittest.main()
