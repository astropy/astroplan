# Use the IERS-A table bundled with astropy-iers-data for all tests, so they
# never download IERS-A (and don't fail when the IERS servers are down).
from astropy.utils import iers

iers.earth_orientation_table.set(iers.IERS_A.open(iers.IERS_A_FILE))
