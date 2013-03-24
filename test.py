#!/usr/local/bin/python

from pyMinc import *
from pyMinctracc import *
from pylab import *
from time import time
from scipy.ndimage.interpolation import *
from scipy.ndimage import *
from numpy import ma as ma
import os.path as p
import mincUtils as minc
from mutualInformation import mi
from utilities import *


class bcolors:
	BLUE = '\033[94m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	RED = '\033[91m'
	END = '\033[0m'



sourceF = '/temp/source.mnc.gz'
targetF = '/temp/target.mnc.gz'
sourceMaskF = '/temp/sourceMask.mnc.gz'
targetMaskF = '/temp/targetMask.mnc.gz'

source = VIOVolume(sourceF,type=float64)
target = VIOVolume(targetF,type=float64)
sourceMask = VIOVolume(sourceMaskF)
targetMask = VIOVolume(targetMaskF)
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


	
	
if 0:			# old minctracc registration
	t0 = time()
	if 1:
		print bcolors.BLUE + 'Doing Blurring' + bcolors.END
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
		mSource2 = sourceF
		
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
		mTarget2 = targetF

	iterations = 20
	t1 = time()
	print bcolors.BLUE + 'Doing linear registration' + bcolors.END
	mLinear = p.join('/temp/temp/linear.xfm')
	cmd = 'minctracc -clobber -lsq12 -step 4 4 4 -iterations %d -source_mask %s -model_mask %s %s %s %s' % (iterations,targetMaskF,sourceMaskF,mTarget4,mSource4,mLinear)	
	minc.execute(cmd)
	t2 = time()
	print bcolors.BLUE + 'Doing non linear 32 registration' + bcolors.END
	mNonlinear32 = p.join('/temp/temp/nonlinear32.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 32 32 32 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mLinear,targetMaskF,sourceMaskF,mTarget32,mSource32,mNonlinear32)
	minc.execute(cmd)
	t3 = time()
	print bcolors.BLUE + 'Doing non linear 16 registration' + bcolors.END
	mNonlinear16 = p.join('/temp/temp/nonlinear16.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 16 16 16 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear32,targetMaskF,sourceMaskF,mTarget16,mSource16,mNonlinear16)
	minc.execute(cmd)
	t4 = time()
	print bcolors.BLUE + 'Doing non linear 8 registration' + bcolors.END
	mNonlinear8 = p.join('/temp/temp/nonlinear8.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 8 8 8 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear16,targetMaskF,sourceMaskF,mTarget8,mSource8,mNonlinear8)
	minc.execute(cmd)
	t5 = time()	
	iterations = 10
	print bcolors.BLUE + 'Doing non linear 4 registration' + bcolors.END
	mNonlinear4 = p.join('/temp/temp/nonlinear4.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 4 4 4 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear8,targetMaskF,sourceMaskF,mTarget4,mSource4,mNonlinear4)
	minc.execute(cmd)
	t6 = time()
	print bcolors.BLUE + 'Doing non linear 2 registration' + bcolors.END
	mNonlinear2 = p.join('/temp/temp/nonlinear2.xfm')
	cmd = 'minctracc -clobber -nonlinear -step 2 2 2 -iterations %d -transformation %s -source_mask %s -model_mask %s %s %s %s' % (iterations,mNonlinear4,targetMaskF,sourceMaskF,mTarget2,mSource2,mNonlinear2)
	minc.execute(cmd)
	t7 = time()

	print "Executable Registration took: "
	print "	Blurring: %0.2f s" % (t1-t0)
	print "	Linear: %0.2f s" % (t2-t1)
	print "	Nonlinear 32: %0.2f s" % (t3-t2)
	print "	Nonlinear 16: %0.2f s" % (t4-t3)
	print "	Nonlinear 8: %0.2f s" % (t5-t4)
	print "	Nonlinear 4: %0.2f s" % (t6-t5)
	print "	Nonlinear 2: %0.2f s" % (t7-t6)


	mtransform = VIOGeneralTransform(mNonlinear8)

	mtransform = mtransform.inverse

	mtransform.write('/temp/temp/nonlinear.xfm')
	mlinear = VIOGeneralTransform(mLinear); mlinear.write('/temp/temp/linear.xfm')
else:
	mtransform = VIOGeneralTransform('/temp/temp/nonlinear.xfm')
	mlinear = VIOGeneralTransform('/temp/temp/linear.xfm')


if 0:			#minctracclib registration
	sourceMeta = source.metadata
	targetMeta = target.metadata
	sourceArgs = {'spacing':sourceMeta['spacing'],'starts':sourceMeta['starts']}
	targetArgs = {'spacing':targetMeta['spacing'],'starts':targetMeta['starts']}

	t0 = time()
	print bcolors.BLUE + 'Doing Blurring' + bcolors.END
	if 1:
		# This gaussian filter seems even better than the MINC one
		size = 16 / 2.35
		source32 = VIOVolume(filters.gaussian_filter(source.data,sigma=size,mode='constant'),**sourceArgs); size /= 2
		source16 = VIOVolume(filters.gaussian_filter(source.data,sigma=size,mode='constant'),**sourceArgs); size /= 2
		source8 = VIOVolume(filters.gaussian_filter(source.data,sigma=size,mode='constant'),**sourceArgs); size /= 2
		source4 = VIOVolume(filters.gaussian_filter(source.data,sigma=size,mode='constant'),**sourceArgs); size /= 2
		source2 = VIOVolume(filters.gaussian_filter(source.data,sigma=size,mode='constant'),**sourceArgs); size = 16 / 2.35
		target32 = VIOVolume(filters.gaussian_filter(target.data,sigma=size,mode='constant'),**targetArgs); size /= 2
		target16 = VIOVolume(filters.gaussian_filter(target.data,sigma=size,mode='constant'),**targetArgs); size /= 2
		target8 = VIOVolume(filters.gaussian_filter(target.data,sigma=size,mode='constant'),**targetArgs); size /= 2
		target4 = VIOVolume(filters.gaussian_filter(target.data,sigma=size,mode='constant'),**targetArgs); size /= 2
		target2 = VIOVolume(filters.gaussian_filter(target.data,sigma=size,mode='constant'),**targetArgs); size /= 2
	else:
		size = 32
		source32 = VIOVolume(filters.uniform_filter(source.data,size=size,mode='constant'),**sourceArgs); size /= 2
		source16 = VIOVolume(filters.uniform_filter(source.data,size=size,mode='constant'),**sourceArgs); size /= 2
		source8 = VIOVolume(filters.uniform_filter(source.data,size=size,mode='constant'),**sourceArgs); size /= 2
		source4 = VIOVolume(filters.uniform_filter(source.data,size=size,mode='constant'),**sourceArgs); size /= 2
		source2 = VIOVolume(filters.uniform_filter(source.data,size=size,mode='constant'),**sourceArgs); size = 32
		target32 = VIOVolume(filters.uniform_filter(target.data,size=size,mode='constant'),**targetArgs); size /= 2
		target16 = VIOVolume(filters.uniform_filter(target.data,size=size,mode='constant'),**targetArgs); size /= 2
		target8 = VIOVolume(filters.uniform_filter(target.data,size=size,mode='constant'),**targetArgs); size /= 2
		target4 = VIOVolume(filters.uniform_filter(target.data,size=size,mode='constant'),**targetArgs); size /= 2
		target2 = VIOVolume(filters.uniform_filter(target.data,size=size,mode='constant'),**targetArgs); size /= 2
		
#	source32 = VIOVolume(mSource32,type=float64)	
#	source16 = VIOVolume(mSource16,type=float64)	
#	source8 = VIOVolume(mSource8,type=float64)	
#	source4 = VIOVolume(mSource4,type=float64)	
#	source2 = VIOVolume(mSource2,type=float64)	
#	target32 = VIOVolume(mTarget32,type=float64)	
#	target16 = VIOVolume(mTarget16,type=float64)	
#	target8 = VIOVolume(mTarget8,type=float64)	
#	target4 = VIOVolume(mTarget4,type=float64)	
#	target2 = VIOVolume(mTarget2,type=float64)	
		
		
	iterations = 20
	t1 = time()
	print bcolors.BLUE + 'Doing linear registration' + bcolors.END
	linear = m.minctracc(target4,source4,sourceMask=targetMask,targetMask=sourceMask,initialXFM=None,iterations=iterations,transformType='lsq12',debug=False)
	t2 = time()
	print bcolors.BLUE + 'Doing non linear 32 registration' + bcolors.END
	nonlinear32 = m.minctracc(target32,source32,sourceMask=targetMask,targetMask=sourceMask,initialXFM=linear,iterations=iterations,transformType='nonlinear',step=[32.0,32.0,32.0],debug=False)
	t3 = time()
	print bcolors.BLUE + 'Doing non linear 16 registration' + bcolors.END
	nonlinear16 = m.minctracc(target16,source16,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear32,iterations=iterations,transformType='nonlinear',step=[16.0,16.0,16.0],debug=False)
	t4 = time()
	print bcolors.BLUE + 'Doing non linear 8 registration' + bcolors.END
	nonlinear8 = m.minctracc(target8,source8,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear16,iterations=iterations,transformType='nonlinear',step=[8.0,8.0,8.0],debug=False)
	t5 = time()	
	iterations = 10
	print bcolors.BLUE + 'Doing non linear 4 registration' + bcolors.END
	nonlinear4 = m.minctracc(target4,source4,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear8,iterations=iterations,transformType='nonlinear',step=[4.0,4.0,4.0],debug=False)
	t6 = time()
	print bcolors.BLUE + 'Doing non linear 2 registration' + bcolors.END
	nonlinear2 = m.minctracc(target2,source2,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear4,iterations=iterations,transformType='nonlinear',step=[2.0,2.0,2.0],debug=False)
	t7 = time()

	print "Library Registration took: "
	print "	Blurring: %0.2f s" % (t1-t0)
	print "	Linear: %0.2f s" % (t2-t1)
	print "	Nonlinear 32: %0.2f s" % (t3-t2)
	print "	Nonlinear 16: %0.2f s" % (t4-t3)
	print "	Nonlinear 8: %0.2f s" % (t5-t4)
	print "	Nonlinear 4: %0.2f s" % (t6-t5)
	print "	Nonlinear 2: %0.2f s" % (t7-t6)


	transform = nonlinear8

	transform = transform.inverse

	transform.write('/temp/nonlinear.xfm')
	linear.write('/temp/linear.xfm')

else:
	transform = VIOGeneralTransform('/temp/nonlinear.xfm')
	linear = VIOGeneralTransform('/temp/linear.xfm')



f = 10
for xfm in [transform.transforms[0].data.data,mtransform.transforms[0].data.data]:
	figure(f); clf()
	s = array(shape(xfm)) / 2
	subplot(221); imshow(xfm[0,s[0]]); colorbar()
	subplot(222); imshow(xfm[0,:,s[0]]); colorbar()
	subplot(223); imshow(xfm[0,:,:,s[0]]); colorbar()
	f += 1




linear = [i for i in transform.transforms if i.transformType == 'linear'][0]
xfm = mat(linear.data)

param = minc.xfmParameters(xfm); param.pop('shears')
xfm = VIOGeneralTransform(minc.xfmFromParameters(**param))

# Make a mincresample image
if not p.exists('/temp/sourceResampled.mnc.gz'):
	cmd = 'mincresample %s %s -transformation %s -like %s' % (sourceF,'/temp/sourceResampled.mnc.gz','/temp/temp/nonlinear.xfm',targetF)
	t1 = time()
	minc.execute(cmd,quiet=False)
	t2 = time()
	print "MINC Resampling took %0.2f s" % (t2-t1)

mincResampled = VIOVolume('/temp/sourceResampled.mnc.gz',type=float64)

# Resample source
print "Linear Resampling"
linearResampled = resample(source,transform=linear,like=target)

draw()

if 1:
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
#		cMap = transform.getDeformation(invert=invert,coordinateMap=True,source=source,target=target).data

		nonlinearResampled = resample(source,transform,target)

		t7 = time()
		t8 = time()
#		nonlinearResampled = map_coordinates(source.data,cMap,order=1)
		t9 = time()



		try:
			deformation = xfm.getDeformation(invert=invert,coordinateMap=False)
			figure(6); clf(); imshow(deformation.data[0,:,34]); colorbar()
		except:
			pass

#		slc = [slice(0,i) for i in shape(cMap)[1:]]
#		big = cMap - mgrid[slc]	
#		figure(7); clf(); imshow(big[0,:,34*4]); colorbar(); title('Coordinate change map')
#		figure(8); clf(); imshow(cMap[0,:,34*4]); colorbar(); title('Coordinate Map')

	figure(4); clf()
	s = array(shape(nonlinearResampled.data)) / 2
	subplot(221); imshow(nonlinearResampled.data[s[0]]);title('Nonlinear Resampled')
	subplot(222); imshow(nonlinearResampled.data[:,s[1]]);
	subplot(223); imshow(nonlinearResampled.data[:,:,s[2]]);

	print "Time to do nonlinear resample: %0.2f s to get cMap, %0.2f s to do resample, %0.2f s total" % (t7-t6,t9-t8,t7-t6+t9-t8)


if 1:
	titles = ['Source','Target','Resampled Source (Linear)','Resampled Source (nonlinear)','Minc resampled source (nonlinear)']
	sdata = ma.array(source.data); sdata.mask = where(sourceMask.data > 0,0,1)
	tdata = ma.array(target.data); tdata.mask = where(targetMask.data > 0,0,1)
	for i,data in enumerate([sdata,tdata,linearResampled.data,nonlinearResampled.data,mincResampled.data]):
		figure(i+1);
		try:
			s = array(shape(data)) / 2
			subplot(221); imshow(data[s[0]]);title(titles[i])
			subplot(222); imshow(data[:,s[1]]);
			subplot(223); imshow(data[:,:,s[2]]);
		except:
			pass

