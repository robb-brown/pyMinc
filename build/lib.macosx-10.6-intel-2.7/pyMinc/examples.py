from pyMinc import *
from pyMinctracc import *
from pylab import *

#fname = '/Volumes/Data2/Data/MODELS/DEFAULT/MNI152/DEFAULT/icbm152_sym_t1w.mnc.gz'

#io = VIOVolume()
#status = io.read(fname)
#data = io.getVolume()
#d = data['data']
#imshow(d[54])
#io.invent()


m = Minctracc()
m.minctracc('minctracc')
