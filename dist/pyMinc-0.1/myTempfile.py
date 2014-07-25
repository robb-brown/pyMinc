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
			path = p.join('/Volumes/lingoRamDisk','lingo-%s' % `int(random.random()*10000000)`)
			os.mkdir(path)
		else:
			print "Asked for a fast temp file but no RAM disk exists."
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
	