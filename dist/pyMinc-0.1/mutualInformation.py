#from useDB2 import *
#import matplotlib.mlab as mlab 
import numpy

def entropy(counts):
	'''Compute entropy.'''
	ps = counts/float(numpy.sum(counts))  # coerce to float and normalize
	ps = ps[numpy.nonzero(ps)]            # toss out zeros
	H = -numpy.sum(ps * numpy.log2(ps))   # compute entropy
	return H


def mi(x,y,bins=100):
	H_xy = entropy(numpy.histogram2d(x.ravel(),y.ravel(),bins=bins,normed=True)[0])
	H_x = entropy(numpy.histogram(x.ravel(), bins=bins,normed=True)[0])
	H_y = entropy(numpy.histogram(y.ravel(), bins=bins,normed=True)[0])
	mi = H_x + H_y - H_xy
	return mi
