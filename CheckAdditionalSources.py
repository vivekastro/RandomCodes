import numpy as np
from astropy.io import fits
import scipy as sp
import matplotlib.pyplot as plt
#from pyspherematch import *
from pydl.pydlutils import yanny
from pydl.pydlutils.spheregroup import *
import os
from astropy.time import Time


baltargets = yanny.read_table_yanny(filename='master-BAL-targets-yanny-format1.dat.txt',tablename='TARGET')
print len(baltargets)

spAllfile = 'spAll-v5_10_7.fits'

spAll = fits.open(spAllfile)[1].data
xx=np.where(spAll['EBOSS_TARGET2'] & 2**25)[0]
balspAll = spAll[xx]


