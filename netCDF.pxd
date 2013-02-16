from builtins cimport *

cdef extern from "netcdf.h":

	ctypedef int nc_type
	int ncvarid(int ncid, const_char_ptr name)
	int ncvarinq(int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp, int *dimidsp, int *nattsp)
	int ncdiminq(int ncid, int dimid, char *name, long *lenp)	

	cdef enum:
		NC_NAT =	0	#/* NAT = 'Not A Type' (c.f. NaN) */
		NC_BYTE =	1	#/* signed 1 byte integer */
		NC_CHAR =	2	#/* ISO/ASCII character */
		NC_SHORT =	3	#/* signed 2 byte integer */
		NC_INT =	4	#/* signed 4 byte integer */
		NC_FLOAT =	5	#/* single precision floating point number */
		NC_DOUBLE =	6	#/* double precision floating point number */
	
#		"""
#		 * 	Default fill values, used unless _FillValue attribute is set.
#		 * These values are stuffed into newly allocated space as appropriate.
#		 * The hope is that one might use these to notice that a particular datum
#		 * has not been set.
#		 """
		NC_FILL_BYTE=	(-127)
		NC_FILL_CHAR	=(0)
		NC_FILL_SHORT=	(-32767)
		NC_FILL_INT	=(-2147483647L)
	
	
	cdef float NC_FILL_FLOAT	= 9.9692099683868690e+36
	cdef double NC_FILL_DOUBLE	= 9.9692099683868690e+36
	
#		"""
#		 * The above values are defaults.
#		 * If you wish a variable to use a different value than the above
#		 * defaults, create an attribute with the same type as the variable
#		 * and the following reserved name. The value you give the attribute
#		 * will be used as the fill value for that variable.
#		 """


	cdef char* _FillValue	="_FillValue"

	cdef enum:
		NC_FILL	=	0	
		NC_NOFILL=	0x100	

#		"""
#		 * 'mode' flags for ncopen
#		 """
		NC_NOWRITE=	0
		NC_WRITE   = 	0x1

#		"""
#		 * 'mode' flags for nccreate
#		 """
		NC_CLOBBER=	0
		NC_NOCLOBBER=	0x4
		NC_64BIT_OFFSET =0x0200 

#		"""
#		 * 'mode' flags for nccreate and ncopen
#		 """
		NC_SHARE=	0x0800	
		NC_STRICT_NC3  =(0x8)

#		""" The following flag currently is ignored, but use in
#		 * nc_open() or nc_create() may someday support use of advisory
#		 * locking to prevent multiple writers from clobbering a file 
#		 """
		NC_LOCK	=	0x0400

#		"""
#		 * Starting with version 3.6, there were two different format netCDF
#		 * files.  netCDF-4 introduces the third one.
#		 """
		NC_FORMAT_CLASSIC =(1)
		NC_FORMAT_64BIT   =(2)
		NC_FORMAT_NETCDF4 =(3)
		NC_FORMAT_NETCDF4_CLASSIC = (4) 

#		"""
#		 * Let nc__create() or nc__open() figure out
#		 * as suitable chunk size.
#		 """
		NC_SIZEHINT_DEFAULT= 0

#		"""
#		 * In nc__enddef(), align to the chunk size.
#		 """
		NC_ALIGN_CHUNK= ((-1))

#		"""
#		 * 'size' argument to ncdimdef for an unlimited dimension
#		 """
		NC_UNLIMITED= 0L

#		"""
#		 * attribute id to put/get a global attribute
#		 """
		NC_GLOBAL= -1

#		"""
#		 * These maximums are enforced by the interface, to facilitate writing
#		 * applications and utilities.  However, nothing is statically allocated to
#		 * these sizes internally.
#		 """
		NC_MAX_DIMS	=1024
		NC_MAX_ATTRS	=8192
		NC_MAX_VARS	=8192	 
		NC_MAX_NAME	=256	 
		NC_MAX_VAR_DIMS	=NC_MAX_DIMS


		FILL_BYTE=	NC_FILL_BYTE
		FILL_CHAR	=NC_FILL_CHAR
		FILL_SHORT=	NC_FILL_SHORT
		FILL_LONG=	NC_FILL_INT

	
	cdef float FILL_FLOAT
	cdef double FILL_DOUBLE
	
	cdef enum:

		MAX_NC_DIMS	=NC_MAX_DIMS
		MAX_NC_ATTRS=	NC_MAX_ATTRS
		MAX_NC_VARS	=NC_MAX_VARS
		MAX_NC_NAME	=NC_MAX_NAME
		MAX_VAR_DIMS=	NC_MAX_VAR_DIMS


		NC_NOERR=	0	

		NC2_ERR     =    (-1) 
		NC_EBADID=	(-33)
		NC_ENFILE=	(-34)
		NC_EEXIST=	(-35)
		NC_EINVAL=	(-36)
		NC_EPERM=	(-37)
		NC_ENOTINDEFINE=	(-38)	
		NC_EINDEFINE	=(-39)
		NC_EINVALCOORDS=	(-40)
		NC_EMAXDIMS=	(-41)
		NC_ENAMEINUSE	=(-42)
		NC_ENOTATT=	(-43)	
		NC_EMAXATTS	=(-44)	
		NC_EBADTYPE=	(-45)	
		NC_EBADDIM	=(-46)	
		NC_EUNLIMPOS	=(-47)
		NC_EMAXVARS	=(-48)	
		NC_ENOTVAR=	(-49)	
		NC_EGLOBAL	=(-50)	
		NC_ENOTNC=	(-51)	
		NC_ESTS     =   	(-52)
		NC_EMAXNAME  =  	(-53)
		NC_EUNLIMIT   = 	(-54)
		NC_ENORECVARS  =	(-55)
		NC_ECHAR=	(-56)	
		NC_EEDGE=	(-57)	
		NC_ESTRIDE=	(-58)	
		NC_EBADNAME=	(-59)	

#		""" N.B. following must match value in ncx.h """
		NC_ERANGE	=(-60)
		NC_ENOMEM	=(-61)

		NC_EVARSIZE  =   (-62) 

		NC_EDIMSIZE   =  (-63)
		NC_ETRUNC     =  (-64)
		NC_EAXISTYPE  =  (-65)

#		""" Following errors are added for DAP """
		NC_EDAP      =   (-66)
		NC_ECURL     =   (-67)
		NC_EIO       =   (-68)
		NC_ENODATA   =   (-69)
		NC_EDAPSVC   =   (-70)
		NC_EDAS	=	(-71)  
		NC_EDDS	=	(-72)  
		NC_EDATADDS	=(-73) 
		NC_EDAPURL	=(-74) 
		NC_EDAPCONSTRAINT =(-75)
