import sys

__version__ = '0.0.1'
__version_info__ = (0, 0, 1)

# check python version
if sys.version_info[0] < 3 or sys.version_info[1] < 6:
    raise ValueError("only supported is Python 3.6+. Please update your environment.")

#from pymcc.wtscans import XyProfile, PDD
#import pymcc.readmcc

__all__ = ["wtscans", "readmcc"]
