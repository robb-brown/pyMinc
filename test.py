#!/usr/local/bin/python

from pyMinc import *
from pyMinctracc import *
from pylab import *
from time import time
from scipy.ndimage.interpolation import *

sourceF = '/temp/CIS-B_001-HSC-1_pt_00001_v02_t1g.mnc.gz'
targetF = '/temp/icbm152_sym_t1w.mnc.gz'

#io = VIOVolume()
#status = io.read(fname)
#data = io.getVolume()
#d = data['data']
#imshow(d[54])
#io.invent()

source = VIOVolume(); source.read(sourceF)
target = VIOVolume(); target.read(targetF)

m = Minctracc()

if 0:			# Time comparison
	print "Simplex"
	t1 = time()
	xfm = m.minctracc(source,target,initialXFM=None,transformType='lsq12',debug=False,verbose=0)
	t2 = time()

	print "\n\nBFGS"
	t3 = time()
	xfm = m.minctracc(source,target,initialXFM=None,transformType='lsq12',optimizer='BFGS',debug=False,verbose=0)
	t4 = time()
	print "\nSimplex took %0.2fs.  BFGS took %0.2fs." % (t2-t1,t4-t3)


#transform = m.minctracc(source,target,initialXFM=None,transformType='lsq9',debug=False)
#transform.write('/temp/output.xfm')
transform = VIOGeneralTransform()
transform.read('/temp/output.xfm')
xfm = mat(transform.transform)

transformed = VIOVolume(); transformed.read('/temp/resampled.mnc.gz')


import lingo.mincUtils as minc


vol = source
metadata = vol.metadata
im = transpose(vol.data,metadata['spatialAxes'])
spacing = array(metadata['spacing'])[metadata['spatialAxes']]
starts = array(metadata['starts'])[metadata['spatialAxes']]
extents = array(metadata['shape'])[metadata['spatialAxes']]
figure(1); imshow(im[s])
for d in range(0,3):
	if spacing[d] < 0:
		starts[d] += extents[d] * spacing[d]
		spacing[d] *= -1
		slc = repeat(slice(None),3).tolist()
		slc[d] = slice(None,None,-1)
		im = im[slc]

figure(2); imshow(im[s])

sourceToWorld = mat(minc.xfmFromParameters(translations=starts,scales=spacing))

#sourceToWorld = mat(source.voxelToWorldTransform.transform)
targetToWorld = mat(target.voxelToWorldTransform.transform)

total = (sourceToWorld*xfm*targetToWorld.I).I
e1 = (total*mat([0,0,0,1]).T).T
e2 = (total*mat(list(array(source.metadata['shape'])[source.metadata['spatialAxes']])+[1]).T).T
s1 = (total*mat([1,1,1,1]).T).T - e1

#outputShape = (array(maximum(e1,e2)+1)[0,0:3]).astype(int16).tolist()   
outputShape = array(target.metadata['shape'])[target.metadata['spatialAxes']]

#im = transpose(source.data,source.metadata['spatialAxes'])
res = affine_transform(im,total[0:3,0:3],offset=array(total)[0:3,3],output_shape=outputShape,order=1)
res = transpose(res,target.metadata['spatialAxes'])
s = 75
figure(1); imshow(source.data[:,s])
figure(2); imshow(res[:,s])
figure(3); imshow(transformed.data[:,s])
