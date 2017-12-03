import random
import os
import os.path as p
import tempfile
import shutil
import traceback
import cStringIO as sIO


# creat a 2 GB RAM disk with:
# diskutil erasevolume HFS+ "lingoRamDisk" `hdiutil attach -nomount ram://4194304`

class bcolors:
	BLUE = '\033[94m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	RED = '\033[91m'
	END = '\033[0m'

assigned = []

def mkdtemp(fast=False):
	if fast:
		if p.exists('/Volumes/lingoRamDisk'):
			try:
				path = p.join('/Volumes/lingoRamDisk','lingo-%s' % `int(random.random()*100000000000)`)
				tries = 0
				while (p.exists(path)) and tries < 10:
					tries += 1
					path = p.join('/Volumes/lingoRamDisk','lingo-%s' % `int(random.random()*100000000000)`)
				if p.exists(path):
					print "Tried ten times but could not make a unique tempdir! Returning a slow one."
					path = tempfile.mkdtemp()
				else:
					os.mkdir(path)
			except:
				print "Exception making fast temp dir (RAM disk full?).  Defaulting to a slow one."
				path = tempfile.mkdtemp()
		else:
			path = tempfile.mkdtemp() 
	else:
		path = tempfile.mkdtemp()
	assigned.append(path)
	return path


def getTempDir(fast=False):
	return mkdtemp(fast)
	
	
def getTempfileName(extension=''):
	if not extension == '' and not extension.startswith('.'):
		extension = '.' + extension
	return 'lingo-%s%s' % (`int(random.random()*10000000)`,extension)
	
	
def releaseTempDir(tmp):
	if tmp in assigned:
		try:
			shutil.rmtree(tmp)
		except:
			pass
		assigned.remove(tmp)
		
		
def releaseAllTempDirs():
	for tmp in assigned:
		try:
			shutil.rmtree(tmp)
		except:
			pass
		assigned.remove(tmp)

		
def printException(s1 = None):
	s = sIO.StringIO()
	s.write(bcolors.RED)
	if s1:
		s.write(s1+'\n')
	traceback.print_exc(file=s)
	s.write(bcolors.END+'\n')
	print s.getvalue()
	s.close()
	
def printError(s1):
	s = sIO.StringIO()
	s.write(bcolors.RED)
	s.write(s1+'\n')
	s.write(bcolors.END+'\n')
	print s.getvalue()
	s.close()
	
def printWarning(s1):
	s = sIO.StringIO()
	s.write(bcolors.YELLOW)
	s.write(s1+'\n')
	s.write(bcolors.END+'\n')
	print s.getvalue()
	s.close()

def printMessage(s1):
	s = sIO.StringIO()
#	s.write(bcolors.RED)
	s.write(s1+'\n')
#	s.write(bcolors.END+'\n')
	print s.getvalue()
	s.close()
