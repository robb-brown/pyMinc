from pyMinc import *
from pyMinctracc import *
from pylab import *

sourceF = '/temp/CIS-B_001-HSC-1_pt_00001_v01_t1g.mnc.gz'
targetF = '/temp/icbm152_sym_t1w.mnc.gz'

#io = VIOVolume()
#status = io.read(fname)
#data = io.getVolume()
#d = data['data']
#imshow(d[54])
#io.invent()

source = VIOVolume()
print "Reading source: %d" % source.read(sourceF)
target = VIOVolume()
print "Reading target: %d" % target.read(sourceF)


m = Minctracc()
m.minctracc(source,target,command='-lsq12')