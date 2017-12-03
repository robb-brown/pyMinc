from numpy import array, matrix, cos, sin, arccos, arcsin, identity, dot, sqrt, diag, float64
from numpy import finfo, double
from numpy.linalg import norm
import math
import copy

epsilon = finfo(double).eps

#xfm = array([
#		[  1.16217590e+00,   1.62585841e-02,   2.77861518e-02, 9.39389729e-01],
#		[ -1.87208174e-02,   1.08681377e+00,   1.24523637e-01, -3.22229175e+01],
#		[ -2.59336571e-02,  -1.40190082e-01,   1.21928711e+00, -3.82661589e+01],
#		[  0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 1.00000000e+00]])

# should transform into:

#-center         0.00000    0.00000    0.00000
#-translation    0.93939  -32.22292  -38.26616
#-rotation      -6.55799    1.26915   -0.95244
#-scale          1.16262    1.09408    1.22759
#-shear         -0.00049    0.00102   -0.00003


def translationMatrix(translation):
	translation = matrix(translation)
	T = matrix(identity(4))
	T[0:3,-1] = translation.T
	return T


def rotationMatrix(rotations):
	rx = rotateXMatrix(rotations[0])
	ry = rotateYMatrix(rotations[1])
	rz = rotateZMatrix(rotations[2])
	return rz*ry*rx

	
def shearMatrix(shears):
	Sh = matrix(identity(4))
	Sh[1,0] = shears[0]
	Sh[2,0] = shears[1]
	Sh[2,1] = shears[2]
	return Sh
	

def scaleMatrix(scales):
	S = matrix(identity(4))
	S[0:3,0:3] = matrix(diag(scales))
	return S


def buildXFM(translations=[0,0,0],scales=[1,1,1],shears=[0,0,0],rotations=[0,0,0],center=[0,0,0],degrees=False):
	center = array(center,float64); translations = array(translations,float64); scales = array(scales,float64)
	shears = array(shears,float64); rotations = array(rotations,float64)
	if degrees:
		rotations = rotations* math.pi / 180.

	T = translationMatrix(translations+center)
	R = rotationMatrix(rotations)
	Sh = shearMatrix(shears)
	S = scaleMatrix(scales)
	C = translationMatrix(-center)
	return array(T*S*Sh*R*C)


def rotateZMatrix(a):
	M = matrix(identity(4))
	M[0,0] = cos(a); M[0,1] = -sin(a)
	M[1,0] = sin(a); M[1,1] = cos(a)
	return M


def rotateYMatrix(a):
	M = matrix(identity(4))
	M[0,0] = cos(a); M[0,2] = sin(a)
	M[2,0] = -sin(a); M[2,2] = cos(a)
	return M


def rotateXMatrix(a):
	M = matrix(identity(4))
	M[1,1] = cos(a); M[1,2] = -sin(a)
	M[2,1] = sin(a); M[2,2] = cos(a)
	return M


def anglesFromMatrix(R):
	R = copy.deepcopy(R)
	R = matrix(R)

	# Find the RZ rotation required to bring the local x into the world XZ plane

	t = matrix(R[:,0]); t[3] = 1
	i,j,k,_ = array(t)[:,0]			# SHOULD CHECK HERE THAT i is positive and not too close to zero

	len = sqrt(i**2+j**2)			# length of a vector on the XY plane... check if this is positive and > 0

	rz = abs(arcsin(j/len)) if abs(i)>abs(j) else abs(arccos(i/len))
	if j > 0: rz = -rz				# Find the counterclockwise angle and correct if necessary

	# Find the RY rotation required to align the local x on the world X axis

	# make a rotate Z matrix and rotate x into the XZ plane
	Rz = rotateZMatrix(rz)
	s = Rz*t
	i,j,k,_ = array(s)[:,0]

	len = sqrt(i**2 + k**2)
	ry = abs(arcsin(k/len)) if abs(i)>abs(k) else abs(arccos(i/len))
	if k<0: ry = -ry

	# Rotate around RX to align the local y with Y and z with Z
	t = array(R[:,2]); t[3] = 1
	Ry = rotateYMatrix(ry)
	s = Rz*t
	t = Ry*s
	i,j,k,_ = array(t)[:,0]
	len = sqrt(j**2+k**2)			# ROBB Problem here.  This is zero and it looks like it shouldn't be		
	rx = abs(arcsin(j/len)) if abs(k)>abs(j) else abs(arccos(k/len))
	if j<0: rx = -rx
	
	return array([-rx,-ry,-rz])
	
	

	
def xfmParameters(xfm,degrees=False):
	xmat = matrix(xfm)

	# find the center of rotation
	center = matrix([0,0,0])
	center_of_rotation = matrix([0,0,0,1]).T; center_of_rotation[0:3] = center.T;
	result = (xmat * center_of_rotation) - center_of_rotation
	translations = array(result[0:3].T)[0]

	# scaling values
	# construct the inverse translation matrix
	Tinv = translationMatrix(-translations)
	# contruct the origin to center translation matrix
	C = translationMatrix(center)
	# construct the center to origin translation matrix
	Cinv = translationMatrix(-center)
	# get scaling*rotation*shear matrix and its inverse
	SRS = Cinv*Tinv*xmat*C
	SRSinv = SRS.I
	# find each scale by mapping a unit vector backwards and finding the magnitude
	Sinv = matrix(identity(4))
	scales = matrix([1.0,1.0,1.0])
	for i in range(0,3):
		unit_vec = matrix([0,0,0,1]).T; unit_vec[i] = 1.0
		result = SRSinv*unit_vec
		magnitude = norm(result[0:3])
		scales[0,i] = 1./magnitude if not magnitude == 0 else 1.0
		Sinv[i,i] = magnitude if not magnitude == 0 else 1.0

	scales = array(scales)[0]

	# Get the skews
	# Use Sinv (scale inverse) to make the shear-rotation matrix from the shear-rotation-scale matrix
	SR = Sinv*SRS; SRinv = SR.I
	x = SRinv[:,0]; y = SRinv[:,1]; z = SRinv[:,2]
	x[-1] = 1.0; y[-1] = 1.0; z[-1] = 1.0
	# get normalized z direction
	nz = z * 1./norm(z[0:3]); nz[-1] = 1.0
	# get a direction perpendicular to z in the yz plane
	y_on_z = nz*float(dot(y[0:3].T,nz[0:3])); y_on_z[-1] = 1.0
	ortho_y = y - y_on_z; ortho_y = ortho_y / norm(ortho_y[0:3]); ortho_y[-1] = 1.0
	# Get C for the skew matrix
	ci = float(dot(y[0:3].T,nz[0:3])) / sqrt(norm(y[0:3])**2 - float(dot(y[0:3].T,nz[0:3]))**2)
	#Get B for the skew matrix
	bi = float(dot(x[0:3].T,nz[0:3])) / sqrt(norm(x[0:3])**2 - float(dot(x[0:3].T,nz[0:3]))**2 - float(dot(x[0:3].T,ortho_y[0:3]))**2)
	# Get A for the skew matrix
	ai = float(dot(x[0:3].T,ortho_y[0:3])) / sqrt(norm(x[0:3])**2 - float(dot(x[0:3].T,nz[0:3]))**2 - float(dot(x[0:3].T,ortho_y[0:3]))**2)
	n1 = sqrt(1+ai**2 + bi**2)
	n2 = sqrt(1 + ci**2)
	ai = ai / n1; bi = bi / n1; ci = ci / n2
	shears = array([ai,bi,ci])

	# Now get the rotation angles
	# inverse scale matrix
	tmp1 = matrix(identity(4))
	tmp1[0,0] = 1./scales[0]; tmp1[1,1] = 1./scales[1]; tmp1[2,2] = 1./scales[2]
	#inverse normalized shear matrix
	tmp2 = matrix(identity(4))
	tmp2[0,0] = sqrt(1-ai**2-bi**2)
	tmp2[1,1] = sqrt(1-ci**2)
	tmp2[1,0] = ai
	tmp2[2,0] = bi
	tmp2[2,1] = ci
	# extract rotation matrix
	T = tmp2 * tmp1
	R = T * SRS
	rotations = anglesFromMatrix(R)

	# Finally, adjust the scale and skew paramaters
	Tinv = T.I
	scales[0] = Tinv[0,0]
	scales[1] = Tinv[1,1]
	scales[2] = Tinv[2,2]
	shears[0] = Tinv[1,0]/scales[1]
	shears[1] = Tinv[2,0]/scales[2]
	shears[2] = Tinv[2,1]/scales[2]

	center = array(center)[0]
	
	if degrees: rotations = rotations*180 / math.pi

	#print 'Center:				%s' % center
	#print 'Translations:			%s' % translations
	#print 'Rotations:			%s' % (rotations*180 / math.pi)
	#print 'Scales:				%s' % scales
	#print 'Shears:				%s' % shears
	
	return {'center':center,'translations':translations,'rotations':rotations,'scales':scales,'shears':shears}





