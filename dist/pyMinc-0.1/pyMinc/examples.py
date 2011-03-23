from pyMinc import *
from numpy import *
from pylab import *

fname = 't1.mnc'
f = mincFile(fname)
imshow(f.getNumpyArray()[30])