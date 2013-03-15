#!/usr/local/bin/python

from pyMinc import *
from pyMinctracc import *
from pylab import *
from time import time
from scipy.ndimage.interpolation import *
from scipy.ndimage import *
from numpy import ma as ma
import os.path as p
import lingo.mincUtils as minc
from lingo.useDB2 import bcolors



def linearResample(source,transform,like):
	vol = source; target = like; xfm = transform
	metadata = vol.metadata
	im = transpose(vol.data,metadata['spatialAxes'])			# Transforms are stored as x,y,z so our image should be too
	spacing = array(metadata['spacing'])[metadata['spatialAxes']]
	starts = array(metadata['starts'])[metadata['spatialAxes']]
	extents = array(metadata['shape'])[metadata['spatialAxes']]

	# Get the (flipped) source to world transform
	sourceToWorld = mat(minc.xfmFromParameters(translations=starts,scales=spacing))
	
	# Get the target to world transform.  Note that this is assuming the atlas, which doesn't need to be flipped
	targetToWorld = mat(target.voxelToWorldTransform.data)
	
	# Total source voxel to target voxel transform, inverted
	# It's important to multiply backwards, and don't forget about order of operations!
	total = (targetToWorld.I*(xfm*sourceToWorld)).I
	
	# These could be useful for doing transform input sampling
	e1 = (total*mat([0,0,0,1]).T).T
	e2 = (total*mat(list(array(source.metadata['shape'])[source.metadata['spatialAxes']])+[1]).T).T
	s1 = (total*mat([1,1,1,1]).T).T - e1
	#outputShape = (array(maximum(e1,e2)+1)[0,0:3]).astype(int16).tolist()   
	
	# Use this if there's a like (using the target for this one)
	outputShape = array(target.metadata['shape'])[target.metadata['spatialAxes']]
	
	# Do the actual (linear) transform
	linearResampled = affine_transform(im,total[0:3,0:3],offset=array(total)[0:3,3],output_shape=outputShape,order=1)
	
	# Put the transformed axes back in the MINC order.  Note that we should check for negative spacing and flip here too
	linearResampled = transpose(linearResampled,target.metadata['spatialAxes'])
	
	return linearResampled





sourceF = '/temp/source.mnc.gz'
targetF = '/temp/target.mnc.gz'
sourceMaskF = '/temp/sourceMask.mnc.gz'
targetMaskF = '/temp/targetMask.mnc.gz'

source = VIOVolume(sourceF,type=float32)
target = VIOVolume(targetF,type=float32)
sourceMask = VIOVolume(sourceMaskF,type=float32)
targetMask = VIOVolume(targetMaskF,type=float32)
(sdata,tdata,linearResampled,nonlinearResampled,mincResampled) = (None,None,None,None,None)

#source = VIOVolume(source.data.astype(uint8))
#target = VIOVolume(target.data.astype(uint8))
#sourceMask = VIOVolume(sourceMask.data.astype(uint8))
#targetMask = VIOVolume(targetMask.data.astype(uint8))


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


	
	
if 1:			# old minctracc registration
	t0 = time()
	if 1:
		print bcolors.BLUE + '\n\nDoing Blurring\n' + bcolors.END
		cmd = ''
		mSource32 = p.join('/temp/temp/source32')
#		cmd = 'mincblur -fwhm 32 -rect -clobber %s %s' % (sourceF,mSource32)
		minc.execute(cmd); mSource32 += '_blur.mnc'
		mSource16 = p.join('/temp/temp/source16')
#		cmd = 'mincblur -fwhm 16 -rect -clobber %s %s' % (sourceF,mSource16)
		minc.execute(cmd); mSource16 += '_blur.mnc'
		mSource8 = p.join('/temp/temp/source8')
#		cmd = 'mincblur -fwhm 8 -rect -clobber %s %s' % (sourceF,mSource8)
		minc.execute(cmd); mSource8 += '_blur.mnc'
		mSource4 = p.join('/temp/temp/source4')
#		cmd = 'mincblur -fwhm 4 -rect -clobber %s %s' % (sourceF,mSource4)
		minc.execute(cmd); mSource4 += '_blur.mnc'
		mSource2 = p.join('/temp/temp/source2')
#		cmd = 'mincblur -fwhm 2 -rect -clobber %s %s' % (sourceF,mSource2)
		minc.execute(cmd); mSource2 += '_blur.mnc'

		mTarget32 = p.join('/temp/temp/target32')
#		cmd = 'mincblur -fwhm 32 -rect -clobber %s %s' % (targetF,mTarget32)
		minc.execute(cmd); mTarget32 += '_blur.mnc'
		mTarget16 = p.join('/temp/temp/target16')
#		cmd = 'mincblur -fwhm 16 -rect -clobber %s %s' % (targetF,mTarget16)
		minc.execute(cmd); mTarget16 += '_blur.mnc'
		mTarget8 = p.join('/temp/temp/target8')
#		cmd = 'mincblur -fwhm 8 -rect -clobber %s %s' % (targetF,mTarget8)
		minc.execute(cmd); mTarget8 += '_blur.mnc'
		mTarget4 = p.join('/temp/temp/target4')
#		cmd = 'mincblur -fwhm 4 -rect -clobber %s %s' % (targetF,mTarget4)
		minc.execute(cmd); mTarget4 += '_blur.mnc'
		mTarget2 = p.join('/temp/temp/target2')
#		cmd = 'mincblur -fwhm 2 -rect -clobber %s %s' % (targetF,mTarget2)
		minc.execute(cmd); mTarget2 += '_blur.mnc'

	iterations = 20
	t1 = time()
	print bcolors.BLUE + '\n\nDoing linear registration\n' + bcolors.END
	mLinear = p.join('/temp/temp/linear.xfm')
	cmd = 'minctracc -clobber -lsq12 -step 4 4 4 -iterations %d -source_mask %s -model_mask %s %s %s %s' % (iterations,targetMaskF,sourceMaskF,mTarget4,mSource4,mLinear)	
	minc.execute(cmd)
	t2 = time()
	print bcolors.BLUE + '\n\nDoing non linear 32 registration\n' + bcolors.END
	mNonlinear32 = p.join('/temp/temp/nonlinear32.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 32 32 32 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mLinear,targetMaskF,sourceMaskF,mTarget32,mSource32,mNonlinear32)
	minc.execute(cmd)
	t3 = time()
	print bcolors.BLUE + '\n\nDoing non linear 16 registration\n' + bcolors.END
	mNonlinear16 = p.join('/temp/temp/nonlinear16.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 16 16 16 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear32,targetMaskF,sourceMaskF,mTarget16,mSource16,mNonlinear16)
	minc.execute(cmd)
	t4 = time()
	print bcolors.BLUE + '\n\nDoing non linear 8 registration\n' + bcolors.END
	mNonlinear8 = p.join('/temp/temp/nonlinear8.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 8 8 8 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear16,targetMaskF,sourceMaskF,mTarget8,mSource8,mNonlinear8)
	minc.execute(cmd)
	t5 = time()	
	iterations = 10
	print bcolors.BLUE + '\n\nDoing non linear 4 registration\n' + bcolors.END
#	mNonlinear4 = p.join('/temp/temp/nonlinear4.xfm')
#	cmd = 'minctracc -clobber -nonlinear -step 4 4 4 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear8,targetMaskF,sourceMaskF,mTarget4,mSource4,mNonlinear4)
#	minc.execute(cmd)
	t6 = time()
	print bcolors.BLUE + '\n\nDoing non linear 2 registration\n' + bcolors.END
	#	mNonlinear2 = p.join('/temp/temp/nonlinear2.xfm')
	#	cmd = 'minctracc -clobber -nonlinear -step 2 2 2 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear4,targetMaskF,sourceMaskF,mTarget2,mSource2,mNonlinear2)
#	minc.execute(cmd)
	t7 = time()

#	nonlinear2 = nonlinear2.transforms[1];
#	nonlinear2.flipInverseFlag()
#	nonlinear2 = VIOGeneralTransform([linear,nonlinear2])

	print "Executable Registration took: "
	print "	Blurring: %0.2f s" % (t1-t0)
	print "	Linear: %0.2f s" % (t2-t1)
	print "	Nonlinear 32: %0.2f s" % (t3-t2)
	print "	Nonlinear 16: %0.2f s" % (t4-t3)
	print "	Nonlinear 8: %0.2f s" % (t5-t4)
	print "	Nonlinear 4: %0.2f s" % (t6-t5)
	print "	Nonlinear 2: %0.2f s" % (t7-t6)


	mtransform = VIOGeneralTransform(mNonlinear8)

	mtransform.flipInverseFlag()

	mtransform.write('/temp/temp/nonlinear.xfm')
	VIOGeneralTransform(mLinear).write('/temp/temp/linear.xfm')
else:
	mtransform = VIOGeneralTransform('/temp/temp/nonlinear.xfm')


if 1:			#minctracclib registration
	sourceMeta = source.metadata
	targetMeta = target.metadata
	sourceArgs = {'spacing':sourceMeta['spacing'],'starts':sourceMeta['starts']}
	targetArgs = {'spacing':targetMeta['spacing'],'starts':targetMeta['starts']}

	t0 = time()
	print bcolors.BLUE + '\n\nDoing Blurring\n' + bcolors.END
	source32 = VIOVolume(filters.uniform_filter(source.data,size=32,mode='constant'),**sourceArgs)
	source16 = VIOVolume(filters.uniform_filter(source.data,size=16,mode='constant'),**sourceArgs)
	source8 = VIOVolume(filters.uniform_filter(source.data,size=8,mode='constant'),**sourceArgs)
	source4 = VIOVolume(filters.uniform_filter(source.data,size=4,mode='constant'),**sourceArgs)
	source2 = VIOVolume(filters.uniform_filter(source.data,size=2,mode='constant'),**sourceArgs)
	target32 = VIOVolume(filters.uniform_filter(target.data,size=32,mode='constant'),**targetArgs)
	target16 = VIOVolume(filters.uniform_filter(target.data,size=16,mode='constant'),**targetArgs)
	target8 = VIOVolume(filters.uniform_filter(target.data,size=8,mode='constant'),**targetArgs)
	target4 = VIOVolume(filters.uniform_filter(target.data,size=4,mode='constant'),**targetArgs)
	target2 = VIOVolume(filters.uniform_filter(target.data,size=2,mode='constant'),**targetArgs)
	
	source32 = VIOVolume(mSource32)
	source16 = VIOVolume(mSource16)
	source8 = VIOVolume(mSource8)
	source4 = VIOVolume(mSource4)
	source2 = VIOVolume(mSource2)
	target32 = VIOVolume(mTarget32)
	target16 = VIOVolume(mTarget16)
	target8 = VIOVolume(mTarget8)
	target4 = VIOVolume(mTarget4)
	target2 = VIOVolume(mTarget2)

	iterations = 20
	t1 = time()
	print bcolors.BLUE + '\n\nDoing linear registration\n' + bcolors.END
	linear = m.minctracc(target4,source4,sourceMask=targetMask,targetMask=sourceMask,initialXFM=None,iterations=iterations,transformType='lsq12',debug=False)
	t2 = time()
	print bcolors.BLUE + '\n\nDoing non linear 32 registration\n' + bcolors.END
	nonlinear32 = m.minctracc(target32,source32,sourceMask=targetMask,targetMask=sourceMask,initialXFM=linear,iterations=iterations,transformType='nonlinear',step=[32.0,32.0,32.0],debug=False)
	t3 = time()
	print bcolors.BLUE + '\n\nDoing non linear 16 registration\n' + bcolors.END
	nonlinear16 = m.minctracc(target16,source16,sourceMask=targetMask,targetMask=sourceMask,initialXFM=linear,iterations=iterations,transformType='nonlinear',step=[16.0,16.0,16.0],debug=False)
	t4 = time()
	print bcolors.BLUE + '\n\nDoing non linear 8 registration\n' + bcolors.END
	nonlinear8 = m.minctracc(target8,source8,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear16,iterations=iterations,transformType='nonlinear',step=[8.0,8.0,8.0],debug=False)
	t5 = time()	
	iterations = 10
	print bcolors.BLUE + '\n\nDoing non linear 4 registration\n' + bcolors.END
	#	nonlinear4 = m.minctracc(target4,source4,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear8,iterations=iterations,transformType='nonlinear',step=[4.0,4.0,4.0],debug=False)
	t6 = time()
	print bcolors.BLUE + '\n\nDoing non linear 2 registration\n' + bcolors.END
	#	nonlinear2 = m.minctracc(target2,source2,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear4,iterations=iterations,transformType='nonlinear',step=[2.0,2.0,2.0],debug=False)
	t7 = time()

#	nonlinear2 = nonlinear2.transforms[1];
#	nonlinear2.flipInverseFlag()
#	nonlinear2 = VIOGeneralTransform([linear,nonlinear2])

	print "Library Registration took: "
	print "	Blurring: %0.2f s" % (t1-t0)
	print "	Linear: %0.2f s" % (t2-t1)
	print "	Nonlinear 32: %0.2f s" % (t3-t2)
	print "	Nonlinear 16: %0.2f s" % (t4-t3)
	print "	Nonlinear 8: %0.2f s" % (t5-t4)
	print "	Nonlinear 4: %0.2f s" % (t6-t5)
	print "	Nonlinear 2: %0.2f s" % (t7-t6)


	transform = nonlinear8

	transform.transforms[0].calculateInverseLinearTransform()
	transform.flipInverseFlag()

	transform.write('/temp/nonlinear.xfm')
	linear.write('/temp/linear.xfm')

else:
	transform = VIOGeneralTransform('/temp/nonlinear.xfm')


linear = [i for i in transform.transforms if i.transformType == 'linear'][0]
xfm = mat(linear.data)

param = minc.xfmParameters(xfm); param.pop('shears')
xfm = mat(minc.xfmFromParameters(**param))

# Make a mincresample image
if not p.exists('/temp/sourceResampled.mnc.gz'):
	cmd = 'mincresample %s %s -transformation %s -like %s' % (sourceF,'/temp/sourceResampled.mnc.gz','/temp/temp/nonlinear.xfm',targetF)
	t1 = time()
	minc.execute(cmd,quiet=False)
	t2 = time()
	print "MINC Resampling took %0.2f s" % (t2-t1)

mincResampled = VIOVolume('/temp/sourceResampled.mnc.gz',type=float32)

# Resample source
linearResampled = linearResample(source,transform=xfm.I,like=target)


if 0:
	# Nonlinear resampling
	print "Doing nonlinear resampling"
	invert = True

	if 0:		# Alternate (faster?) way of making the full size deformation field
		xfm = [i for i in transform.transforms if i.transformType == 'grid'][0]
		t6 = time()
		deformation = xfm.getDeformation(invert=invert,coordinateMap=False)
		metadata = deformation.metadata; spacing = metadata['spacing'][1:]; starts = metadata['starts'][1:]; names = metadata['names'][1:]
		planes = [VIOVolume(deformation.data[i],spacing=spacing,starts=starts,names=names) for i in range(0,3)]
		big = array([linearResample(i,transform=identity(4),like=target) for i in planes])
		slc = [slice(0,i) for i in shape(big)[1:]]
		cMap = big + mgrid[slc]
		t7 = time()
		t8 = time()
		nonlinearResampled = map_coordinates(linearResampled,cMap,order=1)
		t9 = time()	

		figure(6); clf(); imshow(deformation.data[0,:,34]); colorbar()
		figure(7); clf(); imshow(big[0,:,34*4]); colorbar()
		figure(8); clf(); imshow(cMap[0,:,34*4]); colorbar()

	if 1:			# Same way minc does it?
		t6 = time()
		cMap = transform.getDeformation(invert=invert,coordinateMap=True,source=source,target=target).data
		t7 = time()
		t8 = time()
		nonlinearResampled = map_coordinates(source.data,cMap,order=1)
		t9 = time()
		figure(4); clf()
		s = array(shape(nonlinearResampled)) / 2
		subplot(221); imshow(nonlinearResampled[s[0]]);title('Nonlinear Resampled')
		subplot(222); imshow(nonlinearResampled[:,s[1]]);
		subplot(223); imshow(nonlinearResampled[:,:,s[2]]);


		try:
			deformation = xfm.getDeformation(invert=invert,coordinateMap=False)
			figure(6); clf(); imshow(deformation.data[0,:,34]); colorbar()
		except:
			pass

		slc = [slice(0,i) for i in shape(cMap)[1:]]
		big = cMap - mgrid[slc]	
		figure(7); clf(); imshow(big[0,:,34*4]); colorbar(); title('Coordinate change map')
		figure(8); clf(); imshow(cMap[0,:,34*4]); colorbar(); title('Coordinate Map')

	print "Time to do nonlinear resample: %0.2f s to get cMap, %0.2f s to do resample, %0.2f s total" % (t7-t6,t9-t8,t7-t6+t9-t8)


if 1:
	titles = ['Source','Target','Resampled Source (Linear)','Resampled Source (nonlinear)','Minc resampled source (nonlinear)']
	sdata = ma.array(source.data); sdata.mask = where(sourceMask.data > 0,0,1)
	tdata = ma.array(target.data); tdata.mask = where(targetMask.data > 0,0,1)
	for i,data in enumerate([sdata,tdata,linearResampled,nonlinearResampled,mincResampled.data]):
		figure(i+1);
		try:
			s = array(shape(data)) / 2
			subplot(221); imshow(data[s[0]]);title(titles[i])
			subplot(222); imshow(data[:,s[1]]);
			subplot(223); imshow(data[:,:,s[2]]);
		except:
			pass

