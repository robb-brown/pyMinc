from mincFile import *
from pylab import *

fname = '/Volumes/Data/Main/bmt_002-OTT-1_chafe_1234_baseline1_MTR.mnc.gz'

minc = mincFile(fname)
data = minc.getNumpyArray()
print max(ravel(data))
imshow(data[30]); colorbar()
