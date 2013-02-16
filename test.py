from mincFile import *
from pylab import *

fname = '/Volumes/Data2/Data/MODELS/DEFAULT/MNI152/DEFAULT/icbm152_sym_t1w.mnc.gz'

io = VIOVolume()
status = io.read(fname)
data = io.getVolume()
d = data['data']
#d.shape = data['shape']
imshow(d[54])


#io.invent()
