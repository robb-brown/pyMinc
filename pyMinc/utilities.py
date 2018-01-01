from scipy.ndimage.interpolation import *
from scipy.ndimage import *
from numpy import *
import subprocess as sub

from pyMinc import VIOVolume, VIOGeneralTransform
from . import myTempfile as tempfile
from . import xfmParam


# -----------------------------------  Execution ------------------------

def executeWithStdout(command):
	p = sub.Popen(command,shell=True,executable='/bin/bash',stdout=sub.PIPE)
	(stdout,stderr) = p.communicate(None)
	return (p.returncode,stdout)


def executeWithStdoutStderr(command):
	p = sub.Popen(command,shell=True,executable='/bin/bash',stdout=sub.PIPE,stderr=sub.PIPE)
	(stdout,stderr) = p.communicate(None)
	return (p.returncode,stdout,stderr)


def executeWithStdBoth(command):
	p = sub.Popen(command,shell=True,executable='/bin/bash',stdout=sub.PIPE,stderr=sub.STDOUT)
	(stdout,stderr) = p.communicate(None)
	return (p.returncode,stdout)


def executeWithBothOnStdOut(command):
	p = sub.Popen(command,shell=True,executable='/bin/bash',stderr=sub.STDOUT)
	ret = p.wait()
	return ret


def executeWithBothOnStdErr(command):
	p = sub.Popen(command,shell=True,executable='/bin/bash',stdout=sub.STDERR)
	ret = p.wait()
	return ret


def executeWithRedirect(command,stdout=None,stderr=None):
	closestderr = False
	closestdout = False
	if (stderr == stdout):
		stderr = sub.STDOUT
	if (stdout):
		stdout = open(stdout,'w')
		closestdout = True
	if (stderr and not stderr == sub.STDOUT):
		stderr = open(stderr,'w')
		closestderr = True
	p = sub.Popen(command,shell=True,executable='/bin/bash',stdout=stdout,stderr=stderr)
	ret = p.wait()
	if (closestderr):
		stderr.close()
	if (closestdout):
		stdout.close()
	return ret
	

def executeWithNoCapture(command):
	p = sub.Popen(command,shell=True,executable='/bin/bash')
	ret = p.wait()
	return ret
	
	
def executeWithNoWait(command):
	p = sub.Popen(command,shell=True,executable='/bin/bash')
	return


def execute(command,quiet=True,wait=True):
#	print command
#	print ""
	if not wait:
		ret = executeWithNoCapture(command)
	elif (not quiet):
		ret = executeWithNoCapture(command)
	else:
		(ret,stdout) = executeWithStdBoth(command)
	return ret

# -----------------------------------  XFM handling (from mincUtils.py) --------------------

def concatenateXFMs(xfms,output=None,linear=True,old=False):
	if not linear:
		old = True
	if not old:
#		import IPython; IPython.embed()
		xfm = copy.deepcopy(xfms[0])
		if xfm.__class__ == ''.__class__: xfm = readXFM(xfm)
		for xfmFname in xfms[1:]:
			xfm2 = xfmFname
			if xfm2.__class__ == ''.__class__: xfm2 = readXFM(xfm2)
			xfm = numpy.dot(xfm2,xfm)
		if output:
			writeXFM(xfm,output)
			return output
		else:
			return copy.deepcopy(xfm)
	else:			# Using minc
		cmd = "xfmconcat "
		for i in xfms:
			cmd += "%s " % i
		cmd += "%s" % output
	#	print cmd
		ret = execute(cmd,quiet=True)		
		if (not ret):
			return output
		else:
			return None


def sqrtXFM(fname,newfname=None):
	if not newfname:
		newfname = fname
	f = open(fname,'r')
	s = f.read()
	f.close()
	s2 = s[s.find('Linear_Transform ='):].replace(';','')
	xfm = s2.split()[2:]
	for i in range(0,len(xfm)):
		xfm[i] = float(xfm[i])
	xfm.extend([0,0,0,1])
	xfm = reshape(array(xfm),(4,4))
	xfm = sqrtm(xfm).real
	sout = s[:s.find('Linear_Transform =')]
	sout += 'Linear_Transform ='
	for r in range(0,3):
		sout += '\n'
		for c in range(0,4):
			sout += '%f '  % xfm[r,c]
	sout = sout[:-1] + ';\n'
	f = open(newfname,'w')
	f.write(sout)
	f.close()
	
def xfmParameters(xfm,degrees=False):
	return xfmParam.xfmParameters(xfm,degrees)


def xfmFromParameters(**args):
	return xfmParam.buildXFM(**args)


def printXFMParameters(xfm,degrees=True):
	params = xfmParameters(xfm,degrees=degrees)
	for p in ['center','translations','rotations','scales','shears']:
		print('%s: %0.4f %0.4f %0.4f' % (p.title(),params[p][0],params[p][1],params[p][2]))
	

# -----------------------------------------------------------------------------------------

def blurImage(image,level):
	size = level / 2. / 2.35			# sigma in mm
	spacing = array(image.metadata['spacing'])
	sizes = size / spacing
	return VIOVolume(filters.gaussian_filter(image.data,sigma=sizes,mode='constant'),**image.metadata)
	

def linearResample(source,transform,like,invert=False,order=3,originalSpacing=False):
	""" June 14, 2016 - this seems to work, both for the source==like and source != like cases"""
	if transform.__class__ == VIOGeneralTransform:
		xfm = transform.data
	else:
		xfm = transform

	if xfm.__class__ in (list,tuple) and not len(shape(xfm)) == 2:
		xfm = mat(concatenateXFMs(xfm))
	else:
		xfm = mat(xfm)
				
	if alltrue(xfm == identity(4)) and like is None:
		return VIOVolume(source)
		
	if invert:
		xfm = xfm.I
	vol = source; target = like;
	metadata = vol.metadata
	im = transpose(vol.data,metadata['spatialAxes'])			# Transforms are stored as x,y,z so our image should be too
	spacing = array(metadata['spacing'])[metadata['spatialAxes']]
	starts = array(metadata['starts'])[metadata['spatialAxes']]
	extents = array(metadata['shape'])[metadata['spatialAxes']]

	# Get the source to world transform
	sourceToWorld = mat(source.voxelToWorldTransform.data)
	
	# Get the target to world transform.
	targetToWorld = mat(target.voxelToWorldTransform.data)
	
	# Idea here is to change the targetToWorld transform so the output spacing is equal to the source spacing.
	# Might have to change the starts too
	# if originalSpacing:
	# 	params = xfmParameters(targetToWorld)
	# 	params['scales'] = spacing
	# 	targetToWorld = xfmFromParameters(**params)
	
	# Total source voxel to target voxel transform, inverted
	# It's important to multiply backwards, and don't forget about order of operations!
	total = (targetToWorld.I*(xfm*sourceToWorld)).I
	

	# These could be useful for doing transform input sampling
	e1 = (total*mat([0,0,0,1]).T).T
	e2 = (total*mat(list(array(source.metadata['shape'])[source.metadata['spatialAxes']])+[1]).T).T
	s1 = (total*mat([1,1,1,1]).T).T - e1
	#outputShape = (array(maximum(e1,e2)+1)[0,0:3]).astype(int16).tolist()   
	
	# Get the shape for our output image
	outputShape = array(target.metadata['shape'])[target.metadata['spatialAxes']]
	
	# Do the actual (linear) transform
	linearResampled = affine_transform(im,total[0:3,0:3],offset=array(total)[0:3,3],output_shape=outputShape,order=order)
	
	# Put the transformed axes back in the MINC order.  It doesn't appear to be necessary to flip if spacing is negative.
	linearResampled = transpose(linearResampled,target.metadata['spatialAxes'])
	
	metadata = dict(like.metadata)
	if originalSpacing:
		metadata['spacing'] = spacing[metadata['spatialAxes']]
	
	return VIOVolume(linearResampled,**metadata)
	
	
	
def mincResample(source,xfm,like,invert=False,quiet=True):
	print("WARNING: Using MINC resampling")
	import os.path as p
		
	tempdir = tempfile.getTempDir();
	if isinstance(source,VIOVolume):
		source.write(p.join(tempdir,'source.mnc')); source = p.join(tempdir,'source.mnc')
	if isinstance(like,VIOVolume):
		like.write(p.join(tempdir,'like.mnc')); like = p.join(tempdir,'like.mnc')
	if isinstance(xfm,VIOGeneralTransform):
		xfm.write(p.join(tempdir,'xfm.xfm')); xfm = p.join(tempdir,'xfm.xfm')
	
	params = {	'input':source, 
				'output':p.join(tempdir,'resampled.mnc'), 
				'xfm':xfm, 
				'like':like,
				'invert' : '-invert_transformation' if invert else '',
			}
	cmd = "mincresample -2 -clobber {invert} -like {like} {input} {output} -transformation {xfm}".format(**params)
	ret = execute(cmd,quiet=quiet)	
	#minc.resample(p.join(tempdir,'source.mnc.gz'),p.join(tempdir,'resampled.mnc.gz'),xfm=p.join(tempdir,'xfm.xfm'),like=p.join(tempdir,'like.mnc.gz'),simpleResample=True)
	resampled = VIOVolume(p.join(tempdir,'resampled.mnc'))
	tempfile.releaseTempDir(tempdir)								# DEBUG
	return resampled

	
def nonlinearResample(source,transform,like,invert=False,order=2):
	if invert:
		xfm = transform.inverse
	else:
		xfm = transform

	# DEBUG - using MINC resampling for now
	return mincResample(source,xfm,like,invert)
		
		
	cMap = xfm.getDeformation(invert=True,coordinateMap=True,source=source,target=like).data
	nonlinearResampled = map_coordinates(source.data,cMap,order=order)
	return VIOVolume(nonlinearResampled,**like.metadata)
	


def resample(source,transform,like,invert=False,order=3):
	if order.__class__ == str and order.startswith('nearest'):
		order = 0
	allLinear = len([1 for i in transform.transforms if i.transformType == 'grid']) == 0
	if allLinear:
		return linearResample(source,transform,like,invert,order=order)
	else:
		return nonlinearResample(source,transform,like,invert,order=order)
	
	
	