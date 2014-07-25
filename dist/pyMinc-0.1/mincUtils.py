import subprocess as sub
import os.path
from scipy.linalg.matfuncs import sqrtm, expm, inv
import random
import myTempfile; from myTempfile import *
import os.path as p
from ORO.Extras.MINC import *
from time import time
import sys, traceback
from time import time
from numpy import *; import numpy
import copy
import xfmParam
from pyMinc import VIOGeneralTransform, VIOVolume


def getID():
	rnd = random.Random()
	rnd.seed()
	return int(round(rnd.random()*100000))
	

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


def maskVolume(fname):
	ret = executeWithStdout("mincstats -quiet -volume %s -floor 0.1" % fname)[1]
	return float(ret)
	
	
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
	

def getCosines(im):
	if (im.__class__ == ''.__class__):
		im = OROMINCReader.read(im,headerOnly=True)
	h = im.getMINCHeaders()
	cosines = []
	for dim in ['zspace','yspace','xspace']: #h['dimensions']['order']:
		cosines.append([float(i) for i in h[dim]['direction_cosines'].split(', ')])
	cosines = array(cosines)
	temp = numpy.copy(cosines[0])
	cosines[0] = numpy.copy(cosines[2])
	cosines[2] = temp
	return cosines


def cosinesToTransform(cosines):
#	xfm = transpose(cosines)
	xfm = numpy.copy(cosines)
	xfm = concatenate((xfm,array([[0,0,0]]).transpose()),axis=1)
	xfm = concatenate((xfm,array([[0,0,0,1]])),axis=0)
	return xfm
	

def stripCosinesFromImage(im):
	dims = im.getMINCHeaders()['dimensions']['order']
	for dim in dims:
		cosines = im.getMINCHeaders()[dim].get('direction_cosines',None)
		if cosines:
			cosines = cosines.split(',')
			csns = []
			for c in cosines:
				csns.append(float(c))

			amax = argmax(csns)
			csns = zeros(len(csns),float32)
			csns[amax] = 1.0
			im.getMINCHeaders()[dim]['direction_cosines'] = str(list(csns)).strip('[').strip(']')
	return im


def readXFM(fname,linearOnly=False):
	f = open(fname,'r')
	s = f.read()
	f.close()
	begin = s.find('Linear_Transform ='); end = s[begin:].find(';')+begin
	s2 = s[begin:end]
	xfm = s2.split()[2:]
	for i in range(0,len(xfm)):
		xfm[i] = float(xfm[i])
	xfm.extend([0,0,0,1])
	xfm = reshape(array(xfm),(4,4))
	
	if not linearOnly:
		dvol = getDisplacementVolume(fname)
	
		if dvol:
			dvol = OROMINCReader.read(dvol)
			return (xfm,dvol)
		else:
			return xfm
	else:
		return xfm
	
	
def xfmParameters(xfm,degrees=False):
	return xfmParam.xfmParameters(xfm,degrees)
	

def xfmFromParameters(**args):
	return xfmParam.buildXFM(**args)
	
	
def printXFMParameters(xfm,degrees=True):
	params = xfmParameters(xfm,degrees=degrees)
	for p in ['center','translations','rotations','scales','shears']:
		print '%s: %0.4f %0.4f %0.4f' % (p.title(),params[p][0],params[p][1],params[p][2])


def isTransformLinear(transform):
	if (not transform):
		return None
	linear = True
	f = open(transform)
	s = f.read()
	f.close()
	s = s.split('\n')
	dvol = None
	for l in s:
		if (l.startswith('Displacement_Volume')):
			dvol = l.split()[-1][0:-1]
			dvol = p.join(p.dirname(transform),dvol)
			linear = False
	return linear


def getDisplacementVolume(transform):
	if transform.__class__ == VIOGeneralTransform:
		return transform.nonlinear[0]
	linear = True
	f = open(transform)
	s = f.read()
	f.close()
	s = s.split('\n')
	dvol = None
	for l in s:
		if (l.startswith('Displacement_Volume')):
			dvol = l.split()[-1][0:-1]
			dvol = p.join(p.dirname(transform),dvol)
			linear = False
	return dvol


def writeXFM(xfm,fname,dvol = None):
	s = 'MNI Transform File\n\n'
	s += 'Transform_Type = Linear;\n'
	s += 'Linear_Transform =\n\n'
	s += '%0.16f %0.16f %0.16f %0.0f\n' % tuple(xfm[0])
	s += '%0.16f %0.16f %0.16f %0.0f\n' % tuple(xfm[1])
	s += '%0.16f %0.16f %0.16f %0.0f;\n' % tuple(xfm[2])
	if dvol:
		s += 'Transform_Type = Grid_Transform;\n'
		s += 'Displacement_Volume = %s;' % p.basename(dvol)
	s += '\n'
	f = open(fname,'w')
	f.write(s)
	f.close()
	if dvol:
		if not p.exists(p.join(p.dirname(fname),p.basename(dvol))):
			execute('cp %s %s' % (dvol,p.join(p.dirname(fname),p.basename(dvol))))
	return fname


def invertXFMArray(xfm):
	return linalg.inv(xfm)
	

def stripScaleFromXFM(xfm):
	# THIS DOESN'T WORK.  WHY NOT??
	xfm = copy.deepcopy(xfm)
	for dim in [0,1,2]:
		xfm[:,dim] = xfm[:,dim] / numpy.sqrt(numpy.dot(xfm[:,dim],xfm[:,dim]))
	return xfm
		
	
	
def invertXFM(fname,newfname=None):
	if not fname.__class__ == ''.__class__:
		return invertXFMArray(fname)
	tempdir = None
	try:
		if not (newfname):
			tempdir = myTempfile.mkdtemp(fast=True)
			newfname = p.join(tempdir,'temp.xfm')
		execute('xfminvert %s %s' % (fname,newfname))
	except:
		printException("Error inverting XFM")
		myTempfile.releaseTempDir(tempdir)		
	if (tempdir):
		execute('cp %s %s',newfname,fname)
		myTempfile.releaseTempDir(tempdir)
	return newfname
	
	if (0):				# Old way
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
		xfm = inv(xfm)
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


def getXFMLevel(xfm,asString=False):
	if xfm.__class__ == str:
		xfm = VIOGeneralTransform(xfm)
	if xfm.__class__ == VIOGeneralTransform:
		temp = xfm
		if not temp.linear:
			return 'nonlinear'
		xfm = temp.data
	params = xfmParameters(xfm)
	if not sum(abs(params.get('shears',array([0,0,0])))) == 0:
		ret = 12
	elif not sum(abs(params.get('scales',array([0,0,0])))) == 0:
		ret = 9
	elif not sum(abs(params.get('rotations',array([0,0,0])))) == 0:
		ret = 6
	elif not sum(abs(params.get('translations',array([0,0,0])))) == 0:
		ret = 3
	else:
		ret = 0
	if asString:
		ret = xfmLevelToString(ret)
	return ret
		
		
def xfmLevelToInteger(level):
	if level.__class__ == str:
		if level == 'nonlinear':
			return 1000
		else:
			return int(level.replace('lsq',''))
	else:
		return level


def xfmLevelToString(level):
	if not level.__class__ == str:
		if level > 100:
			return 'nonlinear'
		else:
			return 'lsq%d' % level
	else:
		return level


		
def reduceXFMLevel(xfm,level):
	if xfm.__class__ == VIOGeneralTransform:
		params = xfmParameters(xfm.linear.data)
	else:
		params = xfmParameters(xfm)
	if level.__class__ == str:
		lev = int(level.replace('lsq',''))
	else:
		lev = int(level)
	if lev < 12:
		params.pop('shears')
	if lev < 9:
		params.pop('scales')
	if lev < 6:
		params.pop('rotations')
	if xfm.__class__ == VIOGeneralTransform:
		xfm = VIOGeneralTransform(xfmFromParameters(**params))
	else:
		xfm = xfmFromParameters(**params)
	return xfm


def changeDimensionOrder(image,order=['zspace','yspace','xspace']):
	currentOrder = image.getMINCHeaders()['dimensions']['order']
	swapper = [currentOrder.index(d) for d in order]
	image.setNumpyArray(data=image.getNumpyArray().transpose(swapper),labels=order,spacing=list(array(image.getSpacing())[swapper]))
	

def N3(image):
	tempdir = myTempfile.mkdtemp(fast=True)
	if image.__class__ == VIOVolume:
		image.write(p.join(tempdir,'nu-in.mnc'))
	else:
		OROMINCWriter.write(image,p.join(tempdir,'nu-in.mnc'))
	execute('cd %s; nu_correct -clobber nu-in.mnc nu-out.mnc' % tempdir)
	if image.__class__ == VIOVolume:
		im = VIOVolume(p.join(tempdir,'nu-out.mnc'),type=float32)		
	else:
		im = OROMINCReader.read(p.join(tempdir,'nu-out.mnc'))
		changeDimensionOrder(im,image.getMINCHeaders()['dimensions']['order'])
	myTempfile.releaseTempDir(tempdir)
	return im


def BET(image,quiet=True):
	image = copy.deepcopy(image)
	originalOrder = image.getMINCHeaders()['dimensions']['order']
	changeDimensionOrder(image)
	
	# Make all the spacing positive, just in case
	h = image.getMINCHeaders()
	for dim in h['dimensions']['order']:
		h[dim]['step'] = str(abs(float(h[dim]['step'])))
	
	tempdir = myTempfile.mkdtemp(fast=True)
	OROMINCWriter.write(image,p.join(tempdir,'bet-in.mnc'))
	bm = None
	if not execute('cd %s; mnc2nii -nii bet-in.mnc bet-in-nii.nii' % tempdir,quiet=quiet):
		if not execute('cd %s; bet bet-in-nii.nii bm.nii -R -S' % tempdir,quiet=quiet): 
			if not execute('cd %s; gunzip bm_mask.nii.gz; nii2mnc bm_mask.nii' % tempdir,quiet=quiet):
				bm = OROMINCReader.read(p.join(tempdir,'bm_mask.mnc'))
				
	# Put back the correct spacing
	h = bm.getMINCHeaders()
	for dim in h['dimensions']['order']:
		h[dim]['step'] = image.getMINCHeaders()[dim]['step']

	# Put the dimension order back to the original
	changeDimensionOrder(bm,originalOrder)	

	myTempfile.releaseTempDir(tempdir)
	return bm


def FAST(image,brainMask,quiet=True):
	tempdir = myTempfile.mkdtemp(fast=True)
	csf = None; wm = None; gm = None
	brain = copy.deepcopy(image)
	brain.setNumpyArray(image.getNumpyArray()*where(brainMask.getNumpyArray()>0,1,0))
	
	# FAST likes to flip things, so make all the spacing positive
	h = brain.getMINCHeaders()
	for dim in h['dimensions']['order']:
		h[dim]['step'] = str(abs(float(h[dim]['step'])))
		h[dim]['start'] = 0
	
	originalOrder = brain.getMINCHeaders()['dimensions']['order']
	changeDimensionOrder(brain)
	
	OROMINCWriter.write(brain,p.join(tempdir,'fast-in.mnc'))
	if not execute('cd %s; mnc2nii -nii fast-in.mnc fast-in-nii.nii' % tempdir,quiet=quiet):
		if not execute('cd %s; fast -t 1 -o segmentation -p fast-in-nii.nii' % tempdir,quiet=quiet):
			cmd = "cd %s; gunzip segmentation_prob_0.nii.gz; nii2mnc segmentation_prob_0.nii CSFProb.mnc; \
					gunzip segmentation_prob_1.nii.gz; nii2mnc segmentation_prob_1.nii GMProb.mnc; \
					gunzip segmentation_prob_2.nii.gz; nii2mnc segmentation_prob_2.nii WMProb.mnc" % tempdir
			execute(cmd,quiet=quiet)
			csf = OROMINCReader.read(p.join(tempdir,'CSFProb.mnc'))
			wm = OROMINCReader.read(p.join(tempdir,'WMProb.mnc'))
			gm = OROMINCReader.read(p.join(tempdir,'GMProb.mnc'))
	
	for prob in [csf,wm,gm]:
		h = prob.getMINCHeaders()
		for dim in h['dimensions']['order']:
			h[dim]['step'] = image.getMINCHeaders()[dim]['step']
			h[dim]['start'] = image.getMINCHeaders()[dim]['start']

		# FAST seems to flip the last dimension no matter what you do.  So flip it back.  Flipping FAST.
		prob.setNumpyArray((prob.getNumpyArray()[:,:,::-1]))

		# Make sure the dimension order is the same as the original
		changeDimensionOrder(prob,originalOrder)	
	
	myTempfile.releaseTempDir(tempdir)
	return (wm,gm,csf)





















def createMask(input, output, threshold1, threshold2):
	ret = execute('mincmath -clobber -segment -const2 %f %f %s %s' % (threshold1, threshold2,input,output))
	if not ret:
		return output
	else:
		return None

	
def intensityNormalizeFile(contrast,infile,outfile):
	contrastMap = {'T2' : 'NeuroRx_t2w_cerebrum-hist_range1023.2009-01-20.txt',
													'PD' : 'NeuroRx_pdw_cerebrum-hist_range1023.2009-01-20.txt',
													'T1' : 'NeuroRx_t1g_cerebrum-hist_range1023.2008-12-01.txt',
													'T1c' : 'NeuroRx_t1p_cerebrum-hist_range1023.2009-01-20.txt',
													'FLAIR' : 'NeuroRx_t2flair_cerebrum-hist_range1023.2008-12-01.txt',
												}
	model = p.join('models',contrastMap[contrast])
	ret = execute('cd normalization; ./build_string_and_invoke_int_std.pl %(infile)s %(outfile)s -text %(model)s -low 0.01 -high 0.99 -quartiles' %
				{'infile' : infile, 'outfile' : outfile, 'model' : model},quiet=True)
	ret = 0
	if (not ret):
		return outfile
	else:
		return None


def intensityNormalize(contrast,image):
	tempdir = myTempfile.mkdtemp()
	try:
		infile = p.join(tempdir,'infile.mnc')
		outfile = p.join(tempdir,'outfile.mnc')
		OROMINCWriter.write(image,p.join(tempdir,'infile.mnc'))
		outfile = intensityNormalizeFile(contrast,infile,outfile)
		image = OROMINCReader.read(outfile)
	except:
		printException("Error doing intensity normalization")
		myTempfile.releaseTempDir(tempdir)
		return None
	myTempfile.releaseTempDir(tempdir)
	return image
	
	
def registerToTalairach(source,xfm,model="icbm_avg_152_t1_tal_lin"):
	params = { 'model' : model,
						'source' : source,
						'xfm' : xfm,
					}
#	print 'mritotal -model %(model)s %(source)s %(xfm)s' % params
#	ret = execute('mritotal -clobber -model %(model)s %(source)s %(xfm)s' % params,quiet=False)
	
	try:
		t1 = time()
		tempdir = myTempfile.mkdtemp()
#		tempdir = '/Users/rbrown/projects/lingo/temp'			# DEBUG
		target = p.join("models","%s.mnc.gz" % model)
#		target = p.join("/opt/local/mni/bin/../share/mni_autoreg","%s.mnc.gz" % model)
#		target = removeCosines(target,p.join(tempdir,"target.mnc"))
		source = resample(source,p.join(tempdir,"source-resampled.mnc"),None,like=target)

		cmd = "mritoself -clobber -mi -far -lsq7 %(source)s %(target)s %(xfm)s" % ({'source':source,'target':target,'xfm':xfm})
#		print cmd

		ret = execute(cmd,quiet=True)

#		cmd = "mritotal -clobber -model %(model)s %(source)s %(xfm)s" % ({'model':model,'source':source,'xfm':xfm})
#		print cmd
#		ret = execute(cmd,quiet=False)

		t2 = time()
		print "Linear registration took %s s" % (t2-t1)
				
	except:
		exception = sys.exc_info()
		print "Exception: ", exception[0],' ', exception[1]
		traceback.print_tb(exception[2])
		print ""
		del exception		
		return None
	execute("rm -R %s" % tempdir)

	if (not ret):
		return xfm
	else:
		return None
		
		
def registerToTalairachNonlinear(source, xfmLinear, xfmNonLinear, xfmCombined, model="icbm_avg_152_t1_tal_lin", sourceMask = None, targetMask = None):
	xfm = {}
	if (sourceMask and not targetMask):
		targetMask = model + '_mask'
	try:
		tempdir = myTempfile.mkdtemp()
#		tempdir = '/Users/rbrown/projects/lingo/temp'			# DEBUG
		target = p.join("models","%s.mnc.gz" % model)
#		target = p.join("/opt/local/mni/bin/../share/mni_autoreg","%s.mnc.gz" % model)
#		target = removeCosines(target,p.join(tempdir,"target.mnc"))
		source = resample(source,p.join(tempdir,"source-resampled.mnc"),xfmLinear,like=target,convertToUInt=True)
		if (sourceMask):
			sourceMask = resample(sourceMask,p.join(tempdir,"sourcemask-resampled.mnc"),xfmLinear,like=source,convertToUInt=True)
		if (targetMask):
			targetMask = resample(targetMask,p.join(tempdir,"targetmask-resampled.mnc"),None,like=target,convertToUInt=True)

		# nonlinear reg part
		t1 = time()
		cmd = 'nl_register -clobber ' #-init_xfm %s' % xfm['linear']
		if (sourceMask):
			cmd += ' -source_mask %s' % sourceMask
		if (targetMask):
			cmd += ' -target_mask %s' % targetMask
		cmd += ' %(source)s %(target)s %(nl_xfm)s' % ({'source':source,'target':target,'nl_xfm':xfmNonLinear})
		print cmd
		ret = execute(cmd,quiet=False)	

#		ret = execute("mritotal -clobber -nonlinear -model %(model)s %(source)s %(xfm)s" % 
#			({'model':model,'source':source,'xfm':xfmNonLinear}),quiet=True)

		t2 = time()
		print "Nonlinear registration took %s s" % (t2-t1)
		
		xfmCombined = concatenateXFMs([xfmLinear,xfmNonLinear],xfmCombined)
		
		execute("rm -R %s" % tempdir)
		if (not ret):
			xfm['nonlinear'] = xfmNonLinear
			xfm['combined'] = xfmCombined
		return xfm
	except:
		exception = sys.exc_info()
		print "Exception: ", exception[0],' ', exception[1]
		traceback.print_tb(exception[2])
		print ""
		del exception		
		return xfm
		
	
	
def makeNonBrainMask(source,xfm,brainMask,output):
	"""Here are the instructions for running makeNonBrainMask.... """
	tempdir = myTempfile.mkdtemp()
	params = { 'model' : "icbm_avg_152_t1_tal_lin.mnc.gz",
						'brainMask' : brainMask,
						'source' : source,
						'xfm' : xfm,
						'output' : output,
						'temp1': p.join(tempdir,'temp1TMP.mnc'),
						'temp2': p.join(tempdir,'temp2TMP.mnc'),
						'temp3': p.join(tempdir,'temp3TMP.mnc'),
						'temp4': p.join(tempdir,'temp4TMP.mnc'),
						'temp5': p.join(tempdir,'temp5TMP.mnc'),
						'temp6': p.join(tempdir,'temp6TMP.mnc'),
						
					}
	params['model'] = os.path.join(os.path.dirname(executeWithStdout('which mritotal')[1]),'../share/mni_autoreg/icbm_avg_152_t1_tal_lin.mnc.gz')
	
	ret = execute('mincmorph -clobber -successive DDD %(brainMask)s %(temp1)s' % params)
	if (ret):
		print "*** makeNonBrainMask failed: mincmorph"
		return None
	ret += execute('mincmath -clobber -segment -const2 0.9 1.1 %(temp1)s %(temp2)s' % params)
	if (ret):
		print "*** makeNonBrainMask failed: mincmath1"
		return None
	s = 'mincmath -clobber -mult %(source)s %(temp2)s %(temp3)s' % params
	ret += execute('mincmath -clobber -mult %(source)s %(temp2)s %(temp3)s' % params)		
	if (ret):
		print "*** makeNonBrainMask failed: mincmath2: %s"
		return None
	ret += execute("""mincresample -clobber -nearest_neighbour \
									-like %(model)s -transformation %(xfm)s \
									%(temp3)s %(temp4)s""" % params)
	if (ret):
		print "*** makeNonBrainMask failed: mincresample"
		return None
	ret += execute('mincreshape -clobber -dimrange zspace=48,103 %(temp4)s %(temp5)s' % params)
	if (ret):
		print "*** makeNonBrainMask failed: mincreshape"
		return None
	ret += execute("""mincresample -clobber -nearest_neighbour -like %(source)s \
										%(temp5)s %(temp6)s \
										-invert_transformation -transformation %(xfm)s""" % params)
	if (ret):
		print "*** makeNonBrainMask failed: mincresample"
		return None
	ret += execute('mincmath -clobber -mult %(source)s %(temp6)s %(output)s' % params)
	if (ret):
		print "*** makeNonBrainMask failed: mincmath3"
		return None
	execute("rm -R %s" % tempdir)
	if (not ret):
		return output
	else:
		return None
	

def NUCorrect(source,output):
	ret = execute("nu_correct -clobber %(source)s %(output)s" % {'source':source,'output':output});
	if (not ret):
		return output
	else:
		return None
	
	
def fixCosines(input,like,extraArgs=""):
#	print "	Fixing cosines"
	tempdir = myTempfile.mkdtemp()
	try:
		temp = p.join(tempdir,'temp.mnc')
	#	print "Resampling %s to %s like %s" % (input,temp,like)
		ret = execute("mincresample -clobber %(extraArgs)s %(input)s %(temp)s -like %(like)s" % {'extraArgs':extraArgs,'input':input, 'temp':temp, 'like':like})
		if (p.splitext(input)[-1] == '.gz'):
	#		print "Gzipping"
			ret += execute("gzip %s" % temp)
			temp += '.gz'
	#	print "Copying"
		ret += execute("cp %s %s" % (temp, input))
	except:
		printException()
		myTempfile.releaseTempDir(tempdir)
		return False
	myTempfile.releaseTempDir(tempdir)
	return ret==0
	
	
def removeCosines(input, output,convertToUInt=False):
#	print "	Removing cosines"

	quiet = True
	
	if (1):
		tempInput = p.join(p.dirname(output),input.replace(input.split('.')[0],'workingRemoveCosines'))
		execute('cp %s %s' % (input,tempInput),quiet=quiet)
		if tempInput.endswith('.mnc.gz') and output.endswith('.mnc'):
			if (execute('gunzip %s' % (tempInput),quiet=quiet)):
				execute('mv %s %s' % (tempInput,tempInput.replace('.gz','')),quiet=quiet)
			tempInput = tempInput.replace('.gz','')

		if (convertToUInt):
			cmd = 'mincresample -unsigned -short %s %s' % (tempInput,output)
			execute(cmd,quiet=quiet)
		else:
			execute('mv %s %s' % (tempInput,output),quiet=quiet)

		im = OROMINCReader.read(output,headerOnly=True)
		try:
			h = im.getMINCHeaders()
		except:
			print "Error Removing cosines from %s" % input
			return None
		cmd = 'minc_modify_header'
		for dim in h['dimensions']['order']:
			if h[dim].has_key('direction_cosines'):
				cosines = around([float(i) for i in h[dim]['direction_cosines'].split(', ')])
				if sum(cosines) < 0:
					return 'BadCosines'
				cmd += ' -dinsert %s:direction_cosines=%f,%f,%f' % (dim,cosines[0],cosines[1],cosines[2])

		# Fix the stupid irregular slice spacing problem while we're at it
		cmd += " -sinsert zspace:spacing=regular__ "
		cmd += ' %s' % output
		execute(cmd,quiet=quiet)

	if (0): 			# Old way
		inD = OROMINCReader.read(input)
		h = inD.getMINCHeaders()
		badCosines = False
		for d in range(len(h['dimensions']['order'])):
			dim = h['dimensions']['order'][d]
			try:
				cosines = h[dim]['direction_cosines'].split(', ')
			except:
				cosines = [1.0,1.0,1.0]

			for i in range(len(cosines)):
				cosines[i] = float(cosines[i])
			cosines = array(cosines)
			m = argmax(cosines)
			for i in range(len(cosines)):
				cosines[i] = 0.
			cosines[m] = 1.
			h[dim]['direction_cosines'] = '%f, %f, %f' % (cosines[0],cosines[1],cosines[2])
			inD.setMINCHeaders(h)
			if (sum(cosines) < 0):
				badCosines = True
		if (badCosines):
			return 'BadCosines'
		if (convertToUInt):
			data = inD.getNumpyArray()
			mn = min(ravel(data))
			mx = max(ravel(data))
			data = (data + mn) / (mx-mn) * 64000
			data = data.astype(uint16)
			inD.setNumpyArray(data)
		OROMINCWriter.write(inD,output,nocosines=True)
	#	OROMINCWriter.write(inD,output,nocosines=False)
	return output


def addCosines(input, output, like):
#	print "	Adding cosines"
	inD = OROMINCReader.read(input)
	like = OROMINCReader.read(like)
	h = inD.getMINCHeaders()
	hlike = like.getMINCHeaders()
	for d in hlike['dimensions']['order']:
		try:
			cosines = hlike[d]['direction_cosines']
			h[d]['direction_cosines'] = cosines
		except:
			pass
	OROMINCWriter.write(inD,output,nocosines=False)
	return output


		
def resample(input,output=None,xfm=None,like=None,extraArgs="",convertToUInt=False,simpleResample=False,quiet=True,spaceRepresentative=None):
	if not input or (input.__class__ == ''.__class__ and not p.exists(input)):
		if not quiet: print "Input doesn't exist: %s" % input
		return None
	inputImage = None
	if not output:
		outputImage = True
	else:
		outputImage = False
		
	tempdir = myTempfile.mkdtemp(fast=True)
	try:
		if not input.__class__ in [str,unicode]:
			inputImage = input
			input = p.join(tempdir,'tempInput.mnc')
			OROMINCWriter.write(inputImage,input)
		if not output:
			output = p.join(tempdir,'tempOutput.mnc')
		
		if not like == None and not like.__class__ in [str,unicode]:
			likeImage = like
			like = p.join(tempdir,'tempLike.mnc')
			OROMINCWriter.write(likeImage,like)

		xfmMatrix = None
		if not xfm == None and not xfm.__class__ in [str,unicode]:
			xfmMatrix = xfm
			xfm = p.join(tempdir,'tempXFM.xfm')
			writeXFM(xfmMatrix,xfm)
		
	#	print "Resampling %s like %s output: %s" % (input,like,output)
		if (not simpleResample):
			temp = removeCosines(input,p.join(tempdir,'input-resampled.mnc'),convertToUInt=convertToUInt)
		
			if (temp == 'BadCosines'):
				if not spaceRepresentative:
					print "Cosines are broken but resample has no space representative... do a resample with a space representative to fix"
					return None
				print "Fixing cosines"
				fixCosines(input,spaceRepresentative,extraArgs)
				input = removeCosines(input,p.join(tempdir,'input-resampled.mnc'),convertToUInt=convertToUInt)
			else:
				input = temp
			if (like):
				like = removeCosines(like,p.join(tempdir,'like-resampled.mnc'),convertToUInt=convertToUInt)
			
		if (not extraArgs):
			extraArgs = ""
		if (like == None):
			if (xfm):
				cmd = "mincresample -clobber -tfm_input_sampling %(extraArgs)s %(input)s %(output)s -transformation %(xfm)s" % {'extraArgs':extraArgs,'input':input, 'output':output, 'xfm':xfm}
				ret = execute(cmd,quiet=quiet)
			else:
				cmd = "mincresample -tfm_input_sampling -clobber %(extraArgs)s %(input)s %(output)s " % {'extraArgs':extraArgs,'input':input, 'output':output}
				if not quiet:
					print cmd
				ret = execute(cmd)
		else:
			if (like == 'Talairach'):
				like = os.path.join(os.path.dirname(executeWithStdout('which mincresample')[1]),'../share/mni_autoreg/icbm_avg_152_t1_tal_lin.mnc.gz')
			if (xfm):
				cmd = "mincresample -clobber -like %(like)s %(extraArgs)s %(input)s %(output)s -transformation %(xfm)s" % {'extraArgs':extraArgs,'input':input, 'output':output, 'xfm':xfm, 'like':like}
				if not quiet:
					print cmd
				ret = execute(cmd,quiet=quiet)
			else:
				cmd = "mincresample -clobber -like %(like)s %(extraArgs)s %(input)s %(output)s " % {'extraArgs':extraArgs,'input':input, 'output':output, 'like':like}
				if not quiet:
					print cmd
				ret = execute(cmd,quiet=quiet)

		if not quiet:
			print cmd
		if outputImage:
			output = OROMINCReader.read(output)
	except:
		printException()
		myTempfile.releaseTempDir(tempdir)
		return None

	myTempfile.releaseTempDir(tempdir)
		
	if (not ret):
		return output
	else:
		return None

	
		
def rigidRegistration(source, target, xfm, distance='far', removeDirectionCosines=True, doResample=True,extraArgs='',tempdir=None):
	try:
		cmd = ""
		if not tempdir:
			tempdir = myTempfile.mkdtemp()
		t1 = time()
		if (removeDirectionCosines):
			source = removeCosines(source,p.join(tempdir,"source.mnc"))
			target = removeCosines(target,p.join(tempdir,"target.mnc"))
		if (doResample):
			source = resample(source,p.join(tempdir,"source-resampled.mnc"),None,like=target)
			
		cmd = "mritoself -clobber -mi -%(distance)s %(extraArgs)s %(source)s %(target)s %(xfm)s" % ({'distance':distance,'source':source,'target':target,'xfm':xfm,'extraArgs':extraArgs})
		(ret,output) = executeWithStdBoth(cmd)
		t2 = time()
		print "Rigid registration took %s s" % (t2-t1)
		executeWithStdBoth("rm -R %s" % tempdir)
		
		lines = output.splitlines()
		values = []
		for line in lines:
			if line.find('objective function val') > 0:
				val = line[line.find('=')+1:]
				val = float(val)
				values.append(val)
		initial = values[0]
		final = values[-1]		
		
		if (not ret):
			return (xfm,initial,final)
		else:
			print
			print "\nError in rigidRegistration: %s to %s" % (source,target)
			print cmd
			print output
			print
			return None
	except:
		try:
			print "\nError in rigidRegistration: %s to %s" % (source,target)
			print cmd
			print output
		except:
			pass
		traceback.print_exc()
		print ""
		return (None,None,None)

	
		
def linearRegistration(source, target, xfm, sourceMask=None, targetMask=None, removeDirectionCosines=True, doResample=True,metric='mi',type='lsq9',extraArgs=''):
	try:
		cmd = ""
		tempdir = myTempfile.mkdtemp()
		t1 = time()
		if (removeDirectionCosines):
			source = removeCosines(source,p.join(tempdir,"source.mnc"))
			target = removeCosines(target,p.join(tempdir,"target.mnc"))
		if (doResample):
			source = resample(source,p.join(tempdir,"source-resampled.mnc"),None,like=target)
			if sourceMask:
				sourceMask = resample(sourceMask,p.join(tempdir,"sourcemask-resampled.mnc"),None,like=source,convertToUInt=True)
			if targetMask:
				targetMask = resample(targetMask,p.join(tempdir,"targetmask-resampled.mnc"),None,like=target,convertToUInt=True)

		cmd = "minctracc -clobber -%s -%s %s " % (type,metric,extraArgs)
		if (sourceMask):
			cmd += '-source_mask %s ' % sourceMask
		if (targetMask):
			cmd += '-model_mask %s ' % targetMask
		cmd += "%(source)s %(target)s %(xfm)s" % ({'source':source,'target':target,'xfm':xfm})
#		execute(cmd,quiet=False)
#		output = ''
		(ret,output) = executeWithStdBoth(cmd)
		t2 = time()
		print "Linear registration took %s s" % (t2-t1)
		executeWithStdBoth("rm -R %s" % tempdir)

		lines = output.splitlines()
		values = []
		for line in lines:
			if line.find('objective function val') > 0:
				val = line[line.find('=')+1:]
				val = float(val)
				values.append(val)
		initial = values[0]
		final = values[-1]		

		if (not ret):
			return (xfm,initial,final)
		else:
			return None
	except:
		try:
			print "\nError in linearRegistration: %s to %s" % (source,target)
			print cmd
		except:
			pass
		traceback.print_exc()
		print ""
		return (None,None,None)
	
		

		
def nonlinearRegistration(source, target, xfmLinear, xfmNonLinear, xfmCombined, sourceMask=None, targetMask=None,type='-lsq9'):
	xfm = {'linear':None, 'nonlinear':None, 'combined':None}
	try:
#		tempdir = '/Users/rbrown/projects/lingo/temp'			# DEBUG 
		tempdir = myTempfile.mkdtemp()
		target = removeCosines(target,p.join(tempdir,"target.mnc"),convertToUInt=True)
		source = resample(source,p.join(tempdir,"source-resampled.mnc"),None,like=target,convertToUInt=True)
		if (sourceMask):
			sourceMask = resample(sourceMask,p.join(tempdir,"sourcemask-resampled.mnc"),None,like=source,convertToUInt=True)
		if (targetMask):
			targetMask = resample(targetMask,p.join(tempdir,"targetmask-resampled.mnc"),None,like=target,convertToUInt=True)

		# Linear registration part
		linearCalculated = False
		if not p.exists(xfmLinear):
			(xfm['linear'],initial,final) = linearRegistration(source,target,xfmLinear,sourceMask=sourceMask,targetMask=targetMask,removeDirectionCosines=False,doResample=False,type=type)
			linearCalculated = True
		else:
			xfm['linear'] = xfmLinear

		source = resample(source,p.join(tempdir,"source-linearregistered.mnc"),xfm['linear'],like=target)
		if (sourceMask):
			sourceMask = resample(sourceMask,p.join(tempdir,"sourcemask-linearregistered.mnc"),xfm['linear'],like=source)	
	
		# nonlinear reg part
		t1 = time()
		cmd = 'nl_register -clobber ' #-init_xfm %s' % xfm['linear']
		if (sourceMask):
			cmd += ' -source_mask %s' % sourceMask
		if (targetMask):
			cmd += ' -target_mask %s' % targetMask
		cmd += ' %(source)s %(target)s %(nl_xfm)s' % ({'source':source,'target':target,'nl_xfm':xfmNonLinear})
#		print cmd
		ret = execute(cmd,quiet=True)	
#		print ret
		t2 = time()
		print "Nonlinear registration took %s s" % (t2-t1)
		
		xfmCombined = concatenateXFMs([xfm['linear'],xfmNonLinear],xfmCombined,linear=False)
		
		if not linearCalculated:
			xfm.pop('linear')
		
		execute("rm -R %s" % tempdir)
		if (not ret):
			xfm['nonlinear'] = xfmNonLinear
			xfm['combined'] = xfmCombined
		return xfm
	except:
		exception = sys.exc_info()
		print "Exception: ", exception[0],' ', exception[1]
		traceback.print_tb(exception[2])
		print ""
		del exception		
		return xfm



def registerToTalairach(source, target, xfmLinear, xfmNonLinear, xfmCombined, sourceMask=None, targetMask=None,linear=False):
	xfm = {'linear':None, 'nonlinear':None, 'combined':None}
	try:
		tempdir = myTempfile.mkdtemp()
#		target = removeCosines(target,p.join(tempdir,"target.mnc"),convertToUInt=True)
		source = resample(source,p.join(tempdir,"source-resampled.mnc"),None,like=target,convertToUInt=True)
		if (sourceMask):
			sourceMask = resample(sourceMask,p.join(tempdir,"sourcemask-resampled.mnc"),None,like=source,convertToUInt=True)
#		if (targetMask):
#			targetMask = resample(targetMask,p.join(tempdir,"targetmask-resampled.mnc"),None,like=target,convertToUInt=True)

		# Linear registration part
		cmd = 'minctracc -clobber -est_center -est_translations -est_scales -lsq9 -mi -groups 256 '
		if (targetMask):
			cmd += '-model_mask %s ' % targetMask
		if (sourceMask):
			cmd += '-source_mask %s ' % sourceMask
		cmd += '%s %s %s' % (source,target,xfmLinear)
		ret = execute(cmd,quiet=True)
		xfm['linear'] = xfmLinear

		# nonlinear reg part
		if not linear:
			source = resample(source,p.join(tempdir,"source-linearregistered.mnc"),xfm['linear'],like=target)
			if (sourceMask):
				sourceMask = resample(sourceMask,p.join(tempdir,"sourcemask-linearregistered.mnc"),xfm['linear'],like=source)	

			t1 = time()
			cmd = 'nl_register -clobber ' #-init_xfm %s' % xfm['linear']
			if (sourceMask):
				cmd += ' -source_mask %s' % sourceMask
			if (targetMask):
				cmd += ' -target_mask %s' % targetMask
			cmd += ' %(source)s %(target)s %(nl_xfm)s' % ({'source':source,'target':target,'nl_xfm':xfmNonLinear})
			ret = execute(cmd,quiet=True)	
			t2 = time()
			print "Nonlinear registration took %s s" % (t2-t1)

		xfmCombined = concatenateXFMs([xfm['linear'],xfmNonLinear],xfmCombined)
		xfm['nonlinear'] = xfmNonLinear
		xfm['combined'] = xfmCombined

		execute("rm -R %s" % tempdir)
		return xfm
	except:
		exception = sys.exc_info()
		print "Exception: ", exception[0],' ', exception[1]
		traceback.print_tb(exception[2])
		print ""
		del exception		
		return xfm



def rigidRegistrationOld(source, target, xfm,extraArgs="",distance='far'):
	tempdir = myTempfile.mkdtemp()
	ret = execute('uniformize_minc --clobber %(target)s %(utarget)s' % {'target':target,'utarget':p.join(tempdir,'targetuniformized.mnc')})
	ret += execute('uniformize_minc --clobber %(source)s %(usource)s' % {'source':source,'usource':p.join(tempdir,'sourceuniformized.mnc')})
	ret += execute("mritoself -clobber -%(distance)s -mi %(usource)s %(utarget)s %(xfm)s" % {'distance':distance, 'source':source, 'target':target, 'xfm':xfm, 'utarget':p.join(tempdir,'targetuniformized.mnc'), 'usource':p.join(tempdir,'sourceuniformized.mnc')})
	execute("rm -R %s" % tempdir)
	if (not ret):
		return xfm
	else:
		return None
		
		

def rigidRegistrationWithNUCorrect(sources, target, xfms, distance='far'):
	tempdir = myTempfile.mkdtemp()
	if (sources.__class__ == "".__class__):
		sources = [sources]
	if (xfms.__class__ == "".__class__):
		xfms = [xfms]
	params = {	'source'			:	sources,
							'target'			:	target,
							'sourceNU'		:	p.join(tempdir,'sourceNUTMP.mnc'),
							'targetNU'		:	p.join(tempdir,'targetNUTMP.mnc'),
						}

#	print "Doing Nonuniformity correction of target"
	if (NUCorrect(params['target'],params['targetNU'])):	
		for i in range(0,len(sources)):
			params['source'] = sources[i]
#			print "Doing Nonuniformity correction of source %d" % i					
			if (NUCorrect(params['source'],params['sourceNU'])):
#				print "Doing rigid registration"
				(xfms[i],initial,final) = rigidRegistration(params['sourceNU'], params['targetNU'], xfms[i],distance=distance)
	execute("rm -R %s" % tempdir)
	return xfms
	


def nonbrainConstrainedRigidRegistration(source, target, xfm, targetBrainMask, sourceBrainMask, targetNonbrainMask, sourceNonbrainMask):
	tempdir = myTempfile.mkdtemp()
	params = {	'source' 							: source,
							'target' 							: target,
							'xfm' 								: xfm,
							'targetBrainMask' 		: targetBrainMask,
							'sourceBrainMask'			:	sourceBrainMask,
							'targetNonbrainMask'	: targetNonbrainMask,
							'sourceNonbrainMask'	: sourceNonbrainMask,
							'utarget' 					: p.join(tempdir,'uniformtarget.mnc'),
							'usource' 					: p.join(tempdir,'uniformsource.mnc'),
							'headToHeadXFM'				:	p.join(tempdir,'headToHeadTMP%d.xfm'),
							'brainToBrainXFM'			:	p.join(tempdir,'brainToBrainTMP%d.xfm'),
							'skullToSkullXFM'			:	p.join(tempdir,'skullToSkullTMP%d.xfm'),
	}
	

	ret = execute('uniformize_minc --clobber %(target)s %(utarget)s' % params)
	ret += execute('uniformize_minc --clobber %(source)s %(usource)s' % params)
	
	# Head to head registration
	ret += execute("mritoself -clobber -far -mi %(usource)s %(utarget)s %(headToHeadXFM)s" % params)
	
	# Brain to brain registration
	if (not ret):
		ret += execute("minctracc -quiet -clobber -mi -step 2 2 2 -simplex 1 -lsq12 -transformation %(headToHeadXFM)s -model_mask %(targetBrainMask)s %(usource)s %(utarget)s %(brainToBrainXFM)s" % params)
	
	# Skull to skull registration
	if (not ret):
		ret += execute("minctracc -quiet -clobber -mi -step 2 2 2 -simplex 1 -lsq12 -transformation %(brainToBrainXFM)s %(sourceNonbrainMask)s %(targetNonbrainMask)s %(skullToSkullXFM)s" % params);	

	# Brain to brain (final) registration
	if (not ret):
		ret += execute("minctracc -quiet -clobber -mi -step 2 2 2 -transformation %(skullToSkullXFM)s -simplex 1 -lsq12 -w_scales 0 0 0 -w_shear 0 0 0 -model_mask %(targetBrainMask)s %(usource)s %(utarget)s %(xfm)s" % params);	

	execute("rm -R %s" % tempdir)
	
	if (not ret):
		return xfm
	else:
		return None
	
	
	
def nonbrainConstrainedRigidRegistrationWithNUCorrect(source, target, xfm, targetBrainMask, sourceBrainMask, targetNonbrainMask, sourceNonbrainMask):	
	tempdir = myTempfile.mkdtemp()
	params = {	'source'			:	source,
							'target'			:	target,
							'sourceNU'		:	p.join(tempdir,'sourceNUTMP.mnc'),
							'targetNU'		:	p.join(tempdir,'targetNUTMP.mnc'),
							'targetBrainMask' : targetBrainMask,
							'sourceBrainMask' : sourceBrainMask,
							'targetNonbrainMask' : targetNonbrainMask,
							'sourceNonbrainMask' : sourceNonbrainMask,
						}
	
	if (not NUCorrect(params['source'],params['sourceNU'])):
		return None
	if (not NUCorrect(params['target'],params['targetNU'])):
		return None

	xfm = nonbrainConstrainedRigidRegistration(params['sourceNU'], params['targetNU'], xfm, targetBrainMask, sourceBrainMask, targetNonbrainMask, sourceNonbrainMask)

	execute("rm -R %s" % tempdir)
	return xfm
	
	
def nonlinearNonbrainConstrainedRegistrationWithNUCorrect(source, target, xfmLinear, xfmNonLinear, targetBrainMask, sourceBrainMask, targetNonbrainMask, sourceNonbrainMask, targetLesions, sourceLesions):
	tempdir = myTempfile.mkdtemp()
	params = {	'source'			:	source,
							'target'			:	target,
							'xfmNonLinear'	:	xfmNonLinear,
							'xfmLinear'			: xfmLinear,
							'sourceNU'		:	p.join(tempdir,'sourceNUTMP.mnc'),
							'targetNU'		:	p.join(tempdir,'targetNUTMP.mnc'),
							'targetBrainMask' : targetBrainMask,
							'sourceBrainMask' : sourceBrainMask,
							'targetNonbrainMask' : targetNonbrainMask,
							'sourceNonbrainMask' : sourceNonbrainMask,
							'sourceLesions' : sourceLesions,
							'targetLesions' : targetLesions,
							'linearRegisteredSource' : p.join(tempdir,'linearRegisteredSourceTMP.mnc'),
							'linearRegisteredSourceNormalMask'	: p.join(tempdir,'linearRegisteredSourceNormalMaskTMP.mnc'),
							'uniformizedTarget' : p.join(tempdir,'uniformizedTargetTMP.mnc'),
							'sourceBlur8' : p.join(tempdir,'sourceBlur8TMP.mnc'),
							'sourceBlur4' : p.join(tempdir,'sourceBlur4TMP.mnc'),
							'targetBlur8' : p.join(tempdir,'targetBlur8TMP.mnc'),
							'targetBlur4' : p.join(tempdir,'targetBlur4TMP.mnc'),
							'targetResampledLesionMask' : p.join(tempdir,'targetResampledLesionMaskTMP.mnc'),
							'sourceResampledLesionMask' : p.join(tempdir,'sourceResampledLesionMaskTMP.mnc'),
							'targetNormalMaskPreSeg' : p.join(tempdir,'targetNormalMaskPreSegTMP.mnc'),
							'sourceNormalMaskPreSeg' : p.join(tempdir,'sourceNormalMaskPreSegTMP.mnc'),
							'targetNormalMask' : p.join(tempdir,'targetNormalMaskTMP.mnc'),
							'sourceNormalMask' : p.join(tempdir,'sourceNormalMaskTMP.mnc'),
							'xfm4'	: p.join(tempdir,'xfm4TMP.xfm'),
							'weight'			: 0.8,
							'stiffness'		: 0.98,
							'similarity'	: 0.88,
							'iterations4'	: 40,
							'iterations2'	:	28,
						}


	for k in params.keys():
		if (params[k] == None):
			print "Could not calculate nonlinear registration because %s was not provided." % k
			return {'linear':None, 'nonlinear':None}

	if (not NUCorrect(params['source'],params['sourceNU'])):
		return None
	if (not NUCorrect(params['target'],params['targetNU'])):
		return None

	xfm = {}
	ret = 0

	#	Linear Registration Part
#	print "\n\nDoing linear registration\n"
	xfm['linear'] = nonbrainConstrainedRigidRegistration(params['sourceNU'], params['targetNU'], xfmLinear, targetBrainMask, sourceBrainMask, targetNonbrainMask, sourceNonbrainMask)
	
	# create masks w/o lesions here
	ret += execute("mincresample -clobber -nearest_neighbour %(targetLesions)s -like %(targetBrainMask)s %(targetResampledLesionMask)s" % params)
	ret += execute("mincmath -clobber -sub %(targetBrainMask)s %(targetResampledLesionMask)s %(targetNormalMaskPreSeg)s" % params)
	ret += execute("mincmath -clobber -segment -const2 0.9 1.1 %(targetNormalMaskPreSeg)s %(targetNormalMask)s" % params)

	ret += execute("mincresample -clobber -nearest_neighbour %(sourceLesions)s -like %(sourceBrainMask)s %(sourceResampledLesionMask)s" % params)
	ret += execute("mincmath -clobber -sub %(sourceBrainMask)s %(sourceResampledLesionMask)s %(sourceNormalMaskPreSeg)s" % params)
	ret += execute("mincmath -clobber -segment -const2 0.9 1.1 %(sourceNormalMaskPreSeg)s %(sourceNormalMask)s" % params)
	if ret:
		print "Failed to make masks for nonlinear registration"
	
	ret += execute('uniformize_minc --clobber --transform %(xfmLinear)s %(sourceNU)s %(linearRegisteredSource)s' % params)
	ret += execute('uniformize_minc --clobber %(targetNU)s %(uniformizedTarget)s' % params)

	params['linearRegisteredSourceNormalMask'] = resample(params['sourceNormalMask'],params['linearRegisteredSourceNormalMask'],xfm['linear'],params['linearRegisteredSource'],'-nearest_neighbour')

	ret += execute('nl_register -clobber -source_mask %(linearRegisteredSourceNormalMask)s -target_mask %(targetNormalMask)s %(linearRegisteredSource)s %(uniformizedTarget)s %(xfmNonLinear)s' % params)

	if (not ret):
		xfm['nonlinear'] = params['xfmNonLinear']
	else:
		xfm['nonlinear'] = None

	execute("rm -R %s" % tempdir)
	return xfm

	
	# OLD WAY - NOT USED

	params['linearRegisteredSource'] = resample(params['sourceNU'],params['linearRegisteredSource'],xfm['linear'])
	
	# Nonlinear registration part
	# Create blurred images
	ret += execute('mincblur -fwhm 8 -clobber %(source)s %(sourceBlur8)s' % params)
	ret += execute('mincblur -fwhm 4 -clobber %(source)s %(sourceBlur4)s' % params)
	ret += execute('mincblur -fwhm 8 -clobber %(target)s %(targetBlur8)s' % params)
	ret += execute('mincblur -fwhm 4 -clobber %(target)s %(targetBlur4)s' % params)
	if ret:
		print "Failed to make blurred images for nonlinear registration"
	if (not ret):
		# Level 8 registration
#		print "\n\nDoing Level 8 nonlinear registration\n"
		ret = execute("""minctracc -clobber -nonlinear xcorr \
							-source_mask %(linearRegisteredSourceNormalMask)s \
							-model_mask %(targetNormalMask)s \
							-weight %(weight)s \
							-stiffness %(stiffness)s \
							-similarity %(similarity)s \
							-iterations %(iterations4)s \
							-max_def_magnitude 20 \
							-step 4 4 4 \
							-sub_lattice 6 \
							-lattice_diam 12 12 12 \
							-identity \
							%(sourceBlur8)s_blur.mnc %(targetBlur8)s_blur.mnc \
							%(xfm4)s
					""" % params)
	# Level 4 registration
	if (not ret):
#		print "\n\nDoing Level 4 nonlinear registration\n"
		ret = execute("""minctracc -clobber -nonlinear xcorr \
							-source_mask %(linearRegisteredSourceNormalMask)s \
							-model_mask %(targetNormalMask)s \
							-weight %(weight)s \
							-stiffness %(stiffness)s \
							-similarity %(similarity)s \
							-iterations %(iterations2)s \
							-max_def_magnitude 20 \
							-step 2 2 2 \
							-sub_lattice 6 \
							-lattice_diam 6 6 6 \
							-transformation %(xfm4)s \
							%(sourceBlur4)s_blur.mnc %(targetBlur4)s_blur.mnc \
							%(xfmNonLinear)s
					""" % params)
	if (not ret):
		xfm['nonlinear'] = params['xfmNonLinear']
	else:
		xfm['nonlinear'] = None

	execute("rm -R %s" % tempdir)
	return xfm
					
					
def symmetricNonlinearNonbrainConstrainedRegistrationWithNUCorrect(
				source, target, xfmLinearSource, xfmLinearTarget, xfmNonLinearSource, xfmNonLinearTarget,
				targetBrainMask, sourceBrainMask, targetNonbrainMask, sourceNonbrainMask, 
				targetLesions, sourceLesions):
	tempdir = myTempfile.mkdtemp()
	params = {	'source'			:	source,
							'target'			:	target,
							'xfmNonLinearSource'	:	xfmNonLinearSource,
							'xfmLinearSource'			: xfmLinearSource,
							'xfmNonLinearTarget'	:	xfmNonLinearTarget,
							'xfmLinearTarget'			: xfmLinearTarget,
							'sourceNU'		:	p.join(tempdir,'sourceNUTMP.mnc'),
							'targetNU'		:		p.join(tempdir,'targetNUTMP.mnc'),
							'targetBrainMask' : targetBrainMask,
							'sourceBrainMask' : sourceBrainMask,
							'targetNonbrainMask' : targetNonbrainMask,
							'sourceNonbrainMask' : sourceNonbrainMask,
							'sourceLesions' : sourceLesions,
							'targetLesions' : targetLesions,
							'linearRegisteredSource' : 	p.join(tempdir,'linearRegisteredSourceTMP.mnc'),
							'linearRegisteredTarget' : 	p.join(tempdir,'linearRegisteredTargetTMP.mnc'),
							'sourceBlur8' : 	p.join(tempdir,'sourceBlur8TMP.mnc'),
							'sourceBlur4' : 	p.join(tempdir,'sourceBlur4TMP.mnc'),
							'targetBlur8' : 	p.join(tempdir,'targetBlur8TMP.mnc'),
							'targetBlur4' : 	p.join(tempdir,'targetBlur4TMP.mnc'),
							'targetResampledLesionMask' : 	p.join(tempdir,'targetResampledLesionMaskTMP.mnc'),
							'sourceResampledLesionMask' : 	p.join(tempdir,'sourceResampledLesionMaskTMP.mnc'),
							'targetNormalMaskPreSeg' : 	p.join(tempdir,'targetNormalMaskPreSegTMP.mnc'),
							'sourceNormalMaskPreSeg' : 	p.join(tempdir,'sourceNormalMaskPreSegTMP.mnc'),
							'targetNormalMask' : 	p.join(tempdir,'targetNormalMaskTMP.mnc'),
							'linearRegisteredTargetNormalMask' : 	p.join(tempdir,'linearRegisteredTargetNormalMaskTMP.mnc'),
							'sourceNormalMask' : 	p.join(tempdir,'sourceNormalMaskTMP.mnc'),
							'linearRegisteredSourceNormalMask' : 	p.join(tempdir,'linearRegisteredSourceNormalMaskTMP.mnc'),
							'xfm4'	: 	p.join(tempdir,'xfm4TMP.xfm'),
							'weight'			: 0.8,
							'stiffness'		: 0.98,
							'similarity'	: 0.88,
							'iterations4'	: 40,
							'iterations2'	:	28,
						}


	ret = 0
	for k in params.keys():
		if (params[k] == None):
			print "Could not calculate nonlinear registration because %s was not provided." % k
			return {'linear':None, 'nonlinear':None}

	if (not NUCorrect(params['source'],params['sourceNU'])):
		return None
	if (not NUCorrect(params['target'],params['targetNU'])):
		return None

	xfmSource = {}
	xfmTarget = {}

	#	Linear Registration Part
#	print "\n\nDoing linear registration\n"
	xfmSource['linear'] = nonbrainConstrainedRigidRegistration(params['sourceNU'], params['targetNU'], xfmLinearSource, targetBrainMask, sourceBrainMask, targetNonbrainMask, sourceNonbrainMask)
	sqrtXFM(xfmSource['linear'])
	xfmTarget['linear'] = xfmLinearTarget
	invertXFM(xfmSource['linear'],xfmTarget['linear'])
	
	# create masks w/o lesions here
	ret += execute("mincresample -clobber -nearest_neighbour %(targetLesions)s -like %(targetBrainMask)s %(targetResampledLesionMask)s" % params)
	ret += execute("mincmath -clobber -sub %(targetBrainMask)s %(targetResampledLesionMask)s %(targetNormalMaskPreSeg)s" % params)
	ret += execute("mincmath -clobber -segment -const2 0.9 1.1 %(targetNormalMaskPreSeg)s %(targetNormalMask)s" % params)

	ret += execute("mincresample -clobber -nearest_neighbour %(sourceLesions)s -like %(sourceBrainMask)s %(sourceResampledLesionMask)s" % params)
	ret += execute("mincmath -clobber -sub %(sourceBrainMask)s %(sourceResampledLesionMask)s %(sourceNormalMaskPreSeg)s" % params)
	ret += execute("mincmath -clobber -segment -const2 0.9 1.1 %(sourceNormalMaskPreSeg)s %(sourceNormalMask)s" % params)
	
	ret += execute('uniformize_minc --clobber --transform %(xfmLinearSource)s %(sourceNU)s %(linearRegisteredSource)s' % params)
	ret += execute('uniformize_minc --clobber --transform %(xfmLinearTarget)s %(targetNU)s %(linearRegisteredTarget)s' % params)

	params['linearRegisteredTargetNormalMask'] = resample(params['targetNormalMask'],params['linearRegisteredTargetNormalMask'],xfmTarget['linear'],params['linearRegisteredTarget'],'-nearest_neighbour')
	params['linearRegisteredSourceNormalMask'] = resample(params['sourceNormalMask'],params['linearRegisteredSourceNormalMask'],xfmSource['linear'],params['linearRegisteredSource'],'-nearest_neighbour')

	ret += execute('nl_register -clobber -target_mask %(linearRegisteredTargetNormalMask)s -source_mask %(linearRegisteredSourceNormalMask)s %(linearRegisteredSource)s %(linearRegisteredTarget)s %(xfmNonLinearSource)s' % params)

	if (not ret):
		xfmSource['nonlinear'] = params['xfmNonLinearSource']
		xfmTarget['nonlinear'] = None
	else:
		xfmSource['nonlinear'] = None
		xfmTarget['nonlinear'] = None

	execute("rm -R %s" % tempdir)
	return (xfmSource, xfmTarget)

	
	
	# Old way, not used anymore 
	params['linearRegisteredTarget'] = resample(params['targetNU'],params['linearRegisteredTarget'],xfmTarget['linear'])
	params['linearRegisteredSource'] = resample(params['sourceNU'],params['linearRegisteredSource'],xfmSource['linear'])

	
	# Nonlinear registration part
	# Create blurred images
	ret += execute('mincblur -fwhm 8 -clobber %(source)s %(sourceBlur8)s' % params)
	ret += execute('mincblur -fwhm 4 -clobber %(source)s %(sourceBlur4)s' % params)
	ret += execute('mincblur -fwhm 8 -clobber %(target)s %(targetBlur8)s' % params)
	ret += execute('mincblur -fwhm 4 -clobber %(target)s %(targetBlur4)s' % params)
	if (not ret):
		# Level 8 registration
#		print "\n\nDoing Level 8 nonlinear registration\n"
		ret = execute("""minctracc -clobber -nonlinear xcorr \
							-source_mask %(linearRegisteredSourceNormalMask)s \
							-model_mask %(linearRegisteredTargetNormalMask)s \
							-weight %(weight)s \
							-stiffness %(stiffness)s \
							-similarity %(similarity)s \
							-iterations %(iterations4)s \
							-max_def_magnitude 20 \
							-step 4 4 4 \
							-sub_lattice 6 \
							-lattice_diam 12 12 12 \
							-identity \
							%(sourceBlur8)s_blur.mnc %(targetBlur8)s_blur.mnc \
							%(xfm4)s
					""" % params)
	# Level 4 registration
	if (not ret):
#		print "\n\nDoing Level 4 nonlinear registration\n"
		ret = execute("""minctracc -clobber -nonlinear xcorr \
							-source_mask %(linearRegisteredSourceNormalMask)s \
							-model_mask %(linearRegisteredTargetNormalMask)s \
							-weight %(weight)s \
							-stiffness %(stiffness)s \
							-similarity %(similarity)s \
							-iterations %(iterations2)s \
							-max_def_magnitude 20 \
							-step 2 2 2 \
							-sub_lattice 6 \
							-lattice_diam 6 6 6 \
							-transformation %(xfm4)s \
							%(sourceBlur4)s_blur.mnc %(targetBlur4)s_blur.mnc \
							%(xfmNonLinearSource)s
					""" % params)
	if (not ret):
		xfmSource['nonlinear'] = params['xfmNonLinearSource']
		xfmTarget['nonlinear'] = None
	else:
		xfmSource['nonlinear'] = None
		xfmTarget['nonlinear'] = None

	execute("rm %(targetNormalMask)s %(sourceNormalMask)s %(sourceNU)s %(targetNU)s %(sourceBlur4)s_blur.mnc %(sourceBlur8)s_blur.mnc %(targetBlur4)s_blur.mnc %(targetBlur8)s_blur.mnc %(xfm4)s %(linearRegisteredSource)s %(linearRegisteredSourceNormalMask)s %(linearRegisteredTarget)s %(linearRegisteredTargetNormalMask)s" % params)
	return (xfmSource, xfmTarget)

