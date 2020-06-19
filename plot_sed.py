import numpy as np
import astropy.io.fits as fits
import scipy as sp
import matplotlib.pyplot as plt
import astroML 
from matplotlib.ticker import NullFormatter
import pandas as pd
import os,glob
import urllib2
from scipy import   stats
#from pydl.pydl.pydlutils.spheregroup import *
from astroML.plotting import hist
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord
#from dustmaps.sfd import SFDQuery
#from specutils import extinction 
from matplotlib.ticker import MaxNLocator
import seaborn as sns

