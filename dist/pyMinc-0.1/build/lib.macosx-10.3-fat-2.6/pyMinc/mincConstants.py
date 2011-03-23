# Constant declarations

# NetCDF

NC_NAT =	0	#/* NAT = 'Not A Type' (c.f. NaN) */
NC_BYTE =	1	#/* signed 1 byte integer */
NC_CHAR =	2	#/* ISO/ASCII character */
NC_SHORT =	3	#/* signed 2 byte integer */
NC_INT =	4	#/* signed 4 byte integer */
NC_FLOAT =	5	#/* single precision floating point number */
NC_DOUBLE =	6	#/* double precision floating point number */



"""
 * 	Default fill values, used unless _FillValue attribute is set.
 * These values are stuffed into newly allocated space as appropriate.
 * The hope is that one might use these to notice that a particular datum
 * has not been set.
 """
NC_FILL_BYTE=	(-127)
NC_FILL_CHAR	=(0)
NC_FILL_SHORT=	(-32767)
NC_FILL_INT	=(-2147483647L)
NC_FILL_FLOAT	= 9.9692099683868690e+36
NC_FILL_DOUBLE	= 9.9692099683868690e+36


"""
 * The above values are defaults.
 * If you wish a variable to use a different value than the above
 * defaults, create an attribute with the same type as the variable
 * and the following reserved name. The value you give the attribute
 * will be used as the fill value for that variable.
 """
_FillValue	="_FillValue"
NC_FILL	=	0	
NC_NOFILL=	0x100	

"""
 * 'mode' flags for ncopen
 """
NC_NOWRITE=	0
NC_WRITE   = 	0x1

"""
 * 'mode' flags for nccreate
 """
NC_CLOBBER=	0
NC_NOCLOBBER=	0x4
NC_64BIT_OFFSET =0x0200 

"""
 * 'mode' flags for nccreate and ncopen
 """
NC_SHARE=	0x0800	
NC_STRICT_NC3  =(0x8)

""" The following flag currently is ignored, but use in
 * nc_open() or nc_create() may someday support use of advisory
 * locking to prevent multiple writers from clobbering a file 
 """
NC_LOCK	=	0x0400

"""
 * Starting with version 3.6, there were two different format netCDF
 * files.  netCDF-4 introduces the third one.
 """
NC_FORMAT_CLASSIC =(1)
NC_FORMAT_64BIT   =(2)
NC_FORMAT_NETCDF4 =(3)
NC_FORMAT_NETCDF4_CLASSIC = (4) 

"""
 * Let nc__create() or nc__open() figure out
 * as suitable chunk size.
 """
NC_SIZEHINT_DEFAULT= 0

"""
 * In nc__enddef(), align to the chunk size.
 """
NC_ALIGN_CHUNK= ((-1))

"""
 * 'size' argument to ncdimdef for an unlimited dimension
 """
NC_UNLIMITED= 0L

"""
 * attribute id to put/get a global attribute
 """
NC_GLOBAL= -1

"""
 * These maximums are enforced by the interface, to facilitate writing
 * applications and utilities.  However, nothing is statically allocated to
 * these sizes internally.
 """
NC_MAX_DIMS	=1024
NC_MAX_ATTRS	=8192
NC_MAX_VARS	=8192	 
NC_MAX_NAME	=256	 
NC_MAX_VAR_DIMS	=NC_MAX_DIMS


FILL_BYTE=	NC_FILL_BYTE
FILL_CHAR	=NC_FILL_CHAR
FILL_SHORT=	NC_FILL_SHORT
FILL_LONG=	NC_FILL_INT
FILL_FLOAT=	NC_FILL_FLOAT
FILL_DOUBLE=	NC_FILL_DOUBLE

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
                        
""" N.B. following must match value in ncx.h """
NC_ERANGE	=(-60)
NC_ENOMEM	=(-61)

NC_EVARSIZE  =   (-62) 
		
NC_EDIMSIZE   =  (-63)
NC_ETRUNC     =  (-64)
NC_EAXISTYPE  =  (-65)

""" Following errors are added for DAP """
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






# MINC
MI_ORIGINAL_TYPE = 0 

MI_EMPTY_STRING = ""

# Error flags
MI_ERROR = -1
MI_NOERROR = 0

# Maximum length of standard attributes
MI_MAX_ATTSTR_LEN = 64

# Number of spatial dimensions
MI_NUM_SPACE_DIMS = 3

""" Maximum number of image dimensions for image conversion """

""" Bert 10-Aug-2004 - MI_MAX_IMGDIMS used to be defined to be MAX_VAR_DIMS,
 * a constant defined in netcdf.h. For many years MAX_VAR_DIMS was 100,
 * but in netCDF 3.5.1 the value was changed to 512.
 * Unfortunately, the definitions of MI_ICV_DIM_SIZE, MI_ICV_DIM_STEP,
 * and MI_ICV_DIM_START assume that MI_MAX_IMGDIMS is less than or
 * equal to 100.  To avoid changing the MINC API, we have to define
 * MI_MAX_IMGDIMS to 100 here.  Otherwise the miicv_inqdbl() function
 * will return bogus values for these ICV properties.
 """
MI_MAX_IMGDIMS = 100

""" NetCDF standard attributes """
MIunits    =   "units"
MIlong_name  = "long_name"
MIvalid_range ="valid_range"
MIvalid_max   ="valid_max"
MIvalid_min   ="valid_min"
MI_FillValue  ="_FillValue"
MItitle       ="title"
MIhistory     ="history"

""" General variable attributes """
MIvartype  ="vartype"
MIvarid    ="varid"
MIsigntype ="signtype"
MIparent   ="parent"
MIchildren ="children"
MIcomments ="comments"
MIversion  ="version"

""" General attribute constants """
"""    Prefix for identifying a variable attribute pointer """
MI_VARATT_POINTER_PREFIX ="--->"
"""    Separator for elements of MIchildren """
MI_CHILD_SEPARATOR ="\n"
"""    MIvartype values """
MI_GROUP     ="group________"
MI_DIMENSION ="dimension____"
MI_DIM_WIDTH ="dim-width____"
MI_VARATT    ="var_attribute"
"""    MIvarid value """
MI_STDVAR ="MINC standard variable"
"""    MIsigntype values """
MI_SIGNED   ="signed__"
MI_UNSIGNED ="unsigned"
"""    MIversion value """
MI_VERSION_1_0 ="MINC Version    1.0"
MI_CURRENT_VERSION = MI_VERSION_1_0
""" Generally useful values for boolean attributes """
MI_TRUE  ="true_"
MI_FALSE ="false"

""" Dimension names and names of associated variables """
MIxspace      =     "xspace"
MIyspace       =    "yspace"
MIzspace        =   "zspace"
MItime           =  "time"
MItfrequency      = "tfrequency"
MIxfrequency      = "xfrequency"
MIyfrequency   =    "yfrequency"
MIzfrequency    =   "zfrequency"
MIvector_dimension= "vector_dimension"
MIxspace_width =    "xspace-width"
MIyspace_width  =   "yspace-width"
MIzspace_width   =  "zspace-width"
MItime_width     =  "time-width"
MItfrequency_width= "tfrequency-width"
MIxfrequency_width= "xfrequency-width"
MIyfrequency_width= "yfrequency-width"
MIzfrequency_width= "zfrequency-width"

""" Dimension variable attribute names """
""" For dimension variables (MIspacing is also for dimension width vars) """
MIspacing   =        "spacing"
MIstep       =       "step"
MIstart       =      "start"
MIspacetype    =     "spacetype"
MIalignment     =    "alignment"
MIdirection_cosines= "direction_cosines"
""" For dimension width variables """
MIwidth     =        "width"
MIfiltertype =       "filtertype"

""" Dimension attribute constants """
"""    MIgridtype values """
MI_REGULAR =  "regular__"
MI_IRREGULAR ="irregular"
"""    MIspacetype values """
MI_NATIVE  =  "native____"
MI_TALAIRACH ="talairach_"
MI_CALLOSAL = "callosal__"
"""    MIalignment values """
MI_START = "start_"
MI_CENTRE= "centre"
MI_END   = "end___"
MI_CENTER= MI_CENTRE
"""    MIfiltertype values """
MI_SQUARE   =  "square____"
MI_GAUSSIAN  = "gaussian__"
MI_TRIANGULAR ="triangular"

""" The root variable """
MIrootvariable ="rootvariable"

""" The image variable and its attributes """
MIimage  =  "image"
MIimagemax= "image-max"
MIimagemin= "image-min"
MIcomplete= "complete"

""" The patient variable and its attributes """
MIpatient    =    "patient"
MIfull_name   =   "full_name"
MIother_names  =  "other_names"
MIidentification= "identification"
MIother_ids  =    "other_ids"
MIbirthdate   =   "birthdate"
MIsex       =     "sex"
MIage        =    "age"
MIweight      =   "weight"
MIsize         =  "size"
MIaddress       = "address"
MIinsurance_id =  "insurance_id"

""" Patient attribute constants """
MI_MALE  = "male__"
MI_FEMALE= "female"
MI_OTHER = "other_"

""" The study variable and its attributes """
MIstudy        =       "study"
MIstart_time    =      "start_time"
MIstart_year     =     "start_year"
MIstart_month     =    "start_month"
MIstart_day    =       "start_day"
MIstart_hour    =      "start_hour"
MIstart_minute   =     "start_minute"
MIstart_seconds  =     "start_seconds"
MImodality        =    "modality"
MImanufacturer     =   "manufacturer"
MIdevice_model  =      "device_model"
MIinstitution    =     "institution"
MIdepartment      =    "department"
MIstation_id    =      "station_id"
MIreferring_physician= "referring_physician"
MIattending_physician= "attending_physician"
MIradiologist    =     "radiologist"
MIoperator      =      "operator"
MIadmitting_diagnosis ="admitting_diagnosis"
MIprocedure     =      "procedure"
MIstudy_id       =     "study_id"

""" Study attribute constants """
MI_PET =  "PET__"
MI_SPECT= "SPECT"
MI_GAMMA= "GAMMA"
MI_MRI  = "MRI__"
MI_MRS =  "MRS__"
MI_MRA =  "MRA__"
MI_CT  =  "CT___"
MI_DSA =  "DSA__"
MI_DR  =  "DR___"
MI_LABEL= "label"

""" The acquisition variable and its attributes """
MIacquisition =          "acquisition"
MIprotocol     =         "protocol"
MIscanning_sequence=     "scanning_sequence"
MIrepetition_time   =    "repetition_time"
MIecho_time     =        "echo_time"
MIinversion_time =       "inversion_time"
MInum_averages    =      "num_averages"
MIimaging_frequency=     "imaging_frequency"
MIimaged_nucleus    =    "imaged_nucleus"
MIradionuclide       =   "radionuclide"
MIcontrast_agent      =  "contrast_agent"
MIradionuclide_halflife= "radionuclide_halflife"
MItracer           =     "tracer"
MIinjection_time    =    "injection_time"
MIinjection_year     =   "injection_year"
MIinjection_month     =  "injection_month"
MIinjection_day        = "injection_day"
MIinjection_hour  =      "injection_hour"
MIinjection_minute =     "injection_minute"
MIinjection_seconds =    "injection_seconds"
MIinjection_length   =   "injection_length"
MIinjection_dose      =  "injection_dose"
MIdose_units         =   "dose_units"
MIinjection_volume    =  "injection_volume"
MIinjection_route      = "injection_route"

""" Constants for image conversion variable (icv) properties """
""" Maximum number of icv's allowed """
#MI_MAX_NUM_ICV = MAX_NC_OPEN
""" Default max and min for normalization """
MI_DEFAULT_MAX =1.0
MI_DEFAULT_MIN = 0.0
""" For converting data type """
MI_ICV_TYPE     =        1
MI_ICV_SIGN      =       2
MI_ICV_DO_RANGE   =      3
MI_ICV_VALID_MAX   =     4
MI_ICV_VALID_MIN    =    5
""" For doing normalization """
MI_ICV_DO_NORM    =      6
MI_ICV_USER_NORM   =     7
MI_ICV_IMAGE_MAX    =    8
MI_ICV_IMAGE_MIN     =   9
""" Values actually used in normalization - read-only """
MI_ICV_NORM_MAX    =    10
MI_ICV_NORM_MIN    =    11
""" For doing dimension conversions """
MI_ICV_DO_DIM_CONV  =   12
""" For converting vector fields to scalar """
MI_ICV_DO_SCALAR    =   13
""" For flipping axis direction """
MI_ICV_XDIM_DIR   =     14
MI_ICV_YDIM_DIR   =     15
MI_ICV_ZDIM_DIR   =     16
""" For changing size of first two dimensions (excluding MIvector_dimension) """
MI_ICV_ADIM_SIZE   =    17
MI_ICV_BDIM_SIZE    =   18
MI_ICV_KEEP_ASPECT   =  19
""" The pixel size and location of first two dimensions (these are readonly) """
MI_ICV_ADIM_STEP  =     20
MI_ICV_BDIM_STEP   =    21
MI_ICV_ADIM_START   =   22
MI_ICV_BDIM_START   =   23
""" Number of image dimensions for dimension conversion """
MI_ICV_NUM_IMGDIMS  =   24
""" Number of dimensions of image variable taking into account vector/scalar
   data (read-only property) """
MI_ICV_NUM_DIMS     =   25
""" Id of file and image variable (read-only properties) """
MI_ICV_CDFID     =      26
MI_ICV_VARID      =     27
""" Names of MIimagemax and MIimagemin variables """
MI_ICV_MAXVAR      =    28
MI_ICV_MINVAR       =   29
""" For setting input values to a specified fillvalue """
MI_ICV_DO_FILLVALUE  =  30
MI_ICV_FILLVALUE      = 31
""" Image dimension properties. For each dimension, add the dimension 
   number (counting from fastest to slowest). """
MI_ICV_DIM_SIZE   =     1000
MI_ICV_DIM_STEP    =    1100
MI_ICV_DIM_START    =   1200

""" Constants that can be used as values for the above properties. """
""" Possible values for MI_ICV_?DIM_DIR """
MI_ICV_POSITIVE    =     1
MI_ICV_NEGATIVE    =   (-1)
MI_ICV_ANYDIR      =     0
""" Possible value for MI_ICV_?DIM_SIZE """
MI_ICV_ANYSIZE   =     (-1)

""" Error codes.
   Note that they must not conflict with NetCDF error codes since
   they are stored in the same global variable. """
MI_ERR_NONNUMERIC =      1331 
""" Non-numeric type """
MI_ERR_NONCHAR     =     1332 
""" Non-character type """
MI_ERR_NONSCALAR    =    1333 
""" Non-scalar attribute """
MI_ERR_BADOP         =   1334 
""" Bad operation for MI_varaccess """
MI_ERR_NOTPOINTER  =     1335 
""" Attribute is not a pointer """
MI_ERR_BAD_STDVAR   =    1336 
""" Not a standard variable """
MI_ERR_BADSUFFIX     =   1337  
""" Bad dimension width suffix """
MI_ERR_NOICV          =  1338 
""" Out of icv slots """
MI_ERR_BADICV          = 1339 
""" Illegal icv identifier """
MI_ERR_BADPROP     =     1340 
""" Unknown icv property """
MI_ERR_ICVATTACHED  =    1341 
""" Tried to modify attached icv """
MI_ERR_TOOFEWDIMS    =   1342 
""" Too few dimensions to be an image """
MI_ERR_ICVNOTATTACHED =  1343 
""" Tried to access an unattached icv """
MI_ERR_DIMSIZE      =    1344 
""" Dimensions differ in size """
MI_ERR_ICV_INVCOORDS =   1345 
""" Invalid icv coordinates """
MI_ERR_WRONGNDIMS    =   1346 
""" Too many dimensions for a dim var """
MI_ERR_BADMATCH      =   1347 
""" Variables do not match for copy """
MI_ERR_MAXMIN_DIMS   =   1348 
""" Imagemax/min variables vary over image dimensions """
MI_ERR_UNCOMPRESS    =   1349 
""" Not able to uncompress file """

