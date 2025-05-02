from netCDF cimport nc_type
from builtins cimport *

cdef extern from "minc.h":

	cdef enum:
		MI_ORIGINAL_TYPE = 0
		NC_UNSPECIFIED = MI_ORIGINAL_TYPE
	
	
	cdef char *MI_EMPTY_STRING = ""
	
	cdef enum:

		# Error flags
		MI_ERROR = -1
		MI_NOERROR = 0

		# Maximum length of standard attributes
		MI_MAX_ATTSTR_LEN = 64

		# Number of spatial dimensions
		MI_NUM_SPACE_DIMS = 3

#		""" Maximum number of image dimensions for image conversion """

#		""" Bert 10-Aug-2004 - MI_MAX_IMGDIMS used to be defined to be MAX_VAR_DIMS,
#		 * a constant defined in netcdf.h. For many years MAX_VAR_DIMS was 100,
#		 * but in netCDF 3.5.1 the value was changed to 512.
#		 * Unfortunately, the definitions of MI_ICV_DIM_SIZE, MI_ICV_DIM_STEP,
#		 * and MI_ICV_DIM_START assume that MI_MAX_IMGDIMS is less than or
#		 * equal to 100.  To avoid changing the MINC API, we have to define
#		 * MI_MAX_IMGDIMS to 100 here.  Otherwise the miicv_inqdbl() function
#		 * will return bogus values for these ICV properties.
#		 """
		MI_MAX_IMGDIMS = 100

#		""" NetCDF standard attributes """


	cdef char *MIunits    =   "units"
	cdef char *MIlong_name  = "long_name"
	cdef char *MIvalid_range ="valid_range"
	cdef char *MIvalid_max   ="valid_max"
	cdef char *MIvalid_min   ="valid_min"
	cdef char *MI_FillValue  ="_FillValue"
	cdef char *MItitle       ="title"
	cdef char *MIhistory     ="history"

#		""" General variable attributes """
	cdef char *MIvartype  ="vartype"
	cdef char *MIvarid    ="varid"
	cdef char *MIsigntype ="signtype"
	cdef char *MIparent   ="parent"
	cdef char *MIchildren ="children"
	cdef char *MIcomments ="comments"
	cdef char *MIversion  ="version"

#		""" General attribute constants """
#		"""    Prefix for identifying a variable attribute pointer """
	cdef char *MI_VARATT_POINTER_PREFIX ="--->"
#		"""    Separator for elements of MIchildren """
	cdef char *MI_CHILD_SEPARATOR ="\n"
#		"""    MIvartype values """
	cdef char *MI_GROUP     ="group________"
	cdef char *MI_DIMENSION ="dimension____"
	cdef char *MI_DIM_WIDTH ="dim-width____"
	cdef char *MI_VARATT    ="var_attribute"
#		"""    MIvarid value """
	cdef char *MI_STDVAR ="MINC standard variable"
#		"""    MIsigntype values """
	cdef char *MI_SIGNED   ="signed__"
	cdef char *MI_UNSIGNED ="unsigned"
#		"""    MIversion value """
	cdef char *MI_VERSION_1_0 ="MINC Version    1.0"
	cdef char *MI_CURRENT_VERSION = MI_VERSION_1_0
#	cdef char *""" Generally useful values for boolean attributes """
	cdef char *MI_TRUE  ="true_"
	cdef char *MI_FALSE ="false"

#		""" Dimension names and names of associated variables """
	cdef char *MIxspace      =     "xspace"
	cdef char *MIyspace       =    "yspace"
	cdef char *MIzspace        =   "zspace"
	cdef char *MItime           =  "time"
	cdef char *MItfrequency      = "tfrequency"
	cdef char *MIxfrequency      = "xfrequency"
	cdef char *MIyfrequency   =    "yfrequency"
	cdef char *MIzfrequency    =   "zfrequency"
	cdef char *MIvector_dimension= "vector_dimension"
	cdef char *MIxspace_width =    "xspace-width"
	cdef char *MIyspace_width  =   "yspace-width"
	cdef char *MIzspace_width   =  "zspace-width"
	cdef char *MItime_width     =  "time-width"
	cdef char *MItfrequency_width= "tfrequency-width"
	cdef char *MIxfrequency_width= "xfrequency-width"
	cdef char *MIyfrequency_width= "yfrequency-width"
	cdef char *MIzfrequency_width= "zfrequency-width"

#		""" Dimension variable attribute names """
#		""" For dimension variables (MIspacing is also for dimension width vars) """
	cdef char *MIspacing   =        "spacing"
	cdef char *MIstep       =       "step"
	cdef char *MIstart       =      "start"
	cdef char *MIspacetype    =     "spacetype"
	cdef char *MIalignment     =    "alignment"
	cdef char *MIdirection_cosines= "direction_cosines"
#		""" For dimension width variables """
	cdef char *MIwidth     =        "width"
	cdef char *MIfiltertype =       "filtertype"

#		""" Dimension attribute constants """
#		"""    MIgridtype values """
	cdef char *MI_REGULAR =  "regular__"
	cdef char *MI_IRREGULAR ="irregular"
#		"""    MIspacetype values """
	cdef char *MI_NATIVE  =  "native____"
	cdef char *MI_TALAIRACH ="talairach_"
	cdef char *MI_CALLOSAL = "callosal__"
#		"""    MIalignment values """
	cdef char *MI_START = "start_"
	cdef char *MI_CENTRE= "centre"
	cdef char *MI_END   = "end___"
	cdef char *MI_CENTER= MI_CENTRE
#		"""    MIfiltertype values """
	cdef char *MI_SQUARE   =  "square____"
	cdef char *MI_GAUSSIAN  = "gaussian__"
	cdef char *MI_TRIANGULAR ="triangular"

#		""" The root variable """
	cdef char *MIrootvariable ="rootvariable"

#		""" The image variable and its attributes """
	cdef char *MIimage  =  "image"
	cdef char *MIimagemax= "image-max"
	cdef char *MIimagemin= "image-min"
	cdef char *MIcomplete= "complete"

#		""" The patient variable and its attributes """
	cdef char *MIpatient    =    "patient"
	cdef char *MIfull_name   =   "full_name"
	cdef char *MIother_names  =  "other_names"
	cdef char *MIidentification= "identification"
	cdef char *MIother_ids  =    "other_ids"
	cdef char *MIbirthdate   =   "birthdate"
	cdef char *MIsex       =     "sex"
	cdef char *MIage        =    "age"
	cdef char *MIweight      =   "weight"
	cdef char *MIsize         =  "size"
	cdef char *MIaddress       = "address"
	cdef char *MIinsurance_id =  "insurance_id"

#		""" Patient attribute constants """
	cdef char *MI_MALE  = "male__"
	cdef char *MI_FEMALE= "female"
	cdef char *MI_OTHER = "other_"

#		""" The study variable and its attributes """
	cdef char *MIstudy        =       "study"
	cdef char *MIstart_time    =      "start_time"
	cdef char *MIstart_year     =     "start_year"
	cdef char *MIstart_month     =    "start_month"
	cdef char *MIstart_day    =       "start_day"
	cdef char *MIstart_hour    =      "start_hour"
	cdef char *MIstart_minute   =     "start_minute"
	cdef char *MIstart_seconds  =     "start_seconds"
	cdef char *MImodality        =    "modality"
	cdef char *MImanufacturer     =   "manufacturer"
	cdef char *MIdevice_model  =      "device_model"
	cdef char *MIinstitution    =     "institution"
	cdef char *MIdepartment      =    "department"
	cdef char *MIstation_id    =      "station_id"
	cdef char *MIreferring_physician= "referring_physician"
	cdef char *MIattending_physician= "attending_physician"
	cdef char *MIradiologist    =     "radiologist"
	cdef char *MIoperator      =      "operator"
	cdef char *MIadmitting_diagnosis ="admitting_diagnosis"
	cdef char *MIprocedure     =      "procedure"
	cdef char *MIstudy_id       =     "study_id"

#		""" Study attribute constants """
	cdef char *MI_PET =  "PET__"
	cdef char *MI_SPECT= "SPECT"
	cdef char *MI_GAMMA= "GAMMA"
	cdef char *MI_MRI  = "MRI__"
	cdef char *MI_MRS =  "MRS__"
	cdef char *MI_MRA =  "MRA__"
	cdef char *MI_CT  =  "CT___"
	cdef char *MI_DSA =  "DSA__"
	cdef char *MI_DR  =  "DR___"
	cdef char *MI_LABEL= "label"

#		""" The acquisition variable and its attributes """
	cdef char *MIacquisition =          "acquisition"
	cdef char *MIprotocol     =         "protocol"
	cdef char *MIscanning_sequence=     "scanning_sequence"
	cdef char *MIrepetition_time   =    "repetition_time"
	cdef char *MIecho_time     =        "echo_time"
	cdef char *MIinversion_time =       "inversion_time"
	cdef char *MInum_averages    =      "num_averages"
	cdef char *MIimaging_frequency=     "imaging_frequency"
	cdef char *MIimaged_nucleus    =    "imaged_nucleus"
	cdef char *MIradionuclide       =   "radionuclide"
	cdef char *MIcontrast_agent      =  "contrast_agent"
	cdef char *MIradionuclide_halflife= "radionuclide_halflife"
	cdef char *MItracer           =     "tracer"
	cdef char *MIinjection_time    =    "injection_time"
	cdef char *MIinjection_year     =   "injection_year"
	cdef char *MIinjection_month     =  "injection_month"
	cdef char *MIinjection_day        = "injection_day"
	cdef char *MIinjection_hour  =      "injection_hour"
	cdef char *MIinjection_minute =     "injection_minute"
	cdef char *MIinjection_seconds =    "injection_seconds"
	cdef char *MIinjection_length   =   "injection_length"
	cdef char *MIinjection_dose      =  "injection_dose"
	cdef char *MIdose_units         =   "dose_units"
	cdef char *MIinjection_volume    =  "injection_volume"
	cdef char *MIinjection_route      = "injection_route"


#		""" Constants for image conversion variable (icv) properties """
#		""" Maximum number of icv's allowed """
		#MI_MAX_NUM_ICV = MAX_NC_OPEN
#		""" Default max and min for normalization """
	cdef int MI_DEFAULT_MAX =1.0
	cdef int MI_DEFAULT_MIN = 0.0

	cdef enum:
#		""" For converting data type """
		MI_ICV_TYPE     =        1
		MI_ICV_SIGN      =       2
		MI_ICV_DO_RANGE   =      3
		MI_ICV_VALID_MAX   =     4
		MI_ICV_VALID_MIN    =    5
#		""" For doing normalization """
		MI_ICV_DO_NORM    =      6
		MI_ICV_USER_NORM   =     7
		MI_ICV_IMAGE_MAX    =    8
		MI_ICV_IMAGE_MIN     =   9
#		""" Values actually used in normalization - read-only """
		MI_ICV_NORM_MAX    =    10
		MI_ICV_NORM_MIN    =    11
#		""" For doing dimension conversions """
		MI_ICV_DO_DIM_CONV  =   12
#		""" For converting vector fields to scalar """
		MI_ICV_DO_SCALAR    =   13
#		""" For flipping axis direction """
		MI_ICV_XDIM_DIR   =     14
		MI_ICV_YDIM_DIR   =     15
		MI_ICV_ZDIM_DIR   =     16
#		""" For changing size of first two dimensions (excluding MIvector_dimension) """
		MI_ICV_ADIM_SIZE   =    17
		MI_ICV_BDIM_SIZE    =   18
		MI_ICV_KEEP_ASPECT   =  19
#		""" The pixel size and location of first two dimensions (these are readonly) """
		MI_ICV_ADIM_STEP  =     20
		MI_ICV_BDIM_STEP   =    21
		MI_ICV_ADIM_START   =   22
		MI_ICV_BDIM_START   =   23
#		""" Number of image dimensions for dimension conversion """
		MI_ICV_NUM_IMGDIMS  =   24
#		""" Number of dimensions of image variable taking into account vector/scalar
#		   data (read-only property) """
		MI_ICV_NUM_DIMS     =   25
#		""" Id of file and image variable (read-only properties) """
		MI_ICV_CDFID     =      26
		MI_ICV_VARID      =     27
#		""" Names of MIimagemax and MIimagemin variables """
		MI_ICV_MAXVAR      =    28
		MI_ICV_MINVAR       =   29
#		""" For setting input values to a specified fillvalue """
		MI_ICV_DO_FILLVALUE  =  30
		MI_ICV_FILLVALUE      = 31
#		""" Image dimension properties. For each dimension, add the dimension 
#		   number (counting from fastest to slowest). """
		MI_ICV_DIM_SIZE   =     1000
		MI_ICV_DIM_STEP    =    1100
		MI_ICV_DIM_START    =   1200

#		""" Constants that can be used as values for the above properties. """
#		""" Possible values for MI_ICV_?DIM_DIR """
		MI_ICV_POSITIVE    =     1
		MI_ICV_NEGATIVE    =   (-1)
		MI_ICV_ANYDIR      =     0
#		""" Possible value for MI_ICV_?DIM_SIZE """
		MI_ICV_ANYSIZE   =     (-1)

#		""" Error codes.
#		   Note that they must not conflict with NetCDF error codes since
#		   they are stored in the same global variable. """
		MI_ERR_NONNUMERIC =      1331 
#		""" Non-numeric type """
		MI_ERR_NONCHAR     =     1332 
#		""" Non-character type """
		MI_ERR_NONSCALAR    =    1333 
#		""" Non-scalar attribute """
		MI_ERR_BADOP         =   1334 
#		""" Bad operation for MI_varaccess """
		MI_ERR_NOTPOINTER  =     1335 
#		""" Attribute is not a pointer """
		MI_ERR_BAD_STDVAR   =    1336 
#		""" Not a standard variable """
		MI_ERR_BADSUFFIX     =   1337  
#		""" Bad dimension width suffix """
		MI_ERR_NOICV          =  1338 
#		""" Out of icv slots """
		MI_ERR_BADICV          = 1339 
#		""" Illegal icv identifier """
		MI_ERR_BADPROP     =     1340 
#		""" Unknown icv property """
		MI_ERR_ICVATTACHED  =    1341 
#		""" Tried to modify attached icv """
		MI_ERR_TOOFEWDIMS    =   1342 
#		""" Too few dimensions to be an image """
		MI_ERR_ICVNOTATTACHED =  1343 
#		""" Tried to access an unattached icv """
		MI_ERR_DIMSIZE      =    1344 
#		""" Dimensions differ in size """
		MI_ERR_ICV_INVCOORDS =   1345 
#		""" Invalid icv coordinates """
		MI_ERR_WRONGNDIMS    =   1346 
#		""" Too many dimensions for a dim var """
		MI_ERR_BADMATCH      =   1347 
#		""" Variables do not match for copy """
		MI_ERR_MAXMIN_DIMS   =   1348 
#		""" Imagemax/min variables vary over image dimensions """
		MI_ERR_UNCOMPRESS    =   1349 
#		""" Not able to uncompress file """
	
	
	
	char *miexpand_file(char *path, char *tempfile, int header_only, int *created_tempfile)	
	int miopen(char *path, int mode)
	int micreate(char *path, int cmode)
	int miclose(int cdfid)
	int miattget_with_sign(int cdfid, int varid, char *name, 
		                   char *insign, nc_type datatype, char *outsign,
		                   int max_length, void *value, int *att_length)
	int miattget(int cdfid, int varid, char *name, nc_type datatype,
		                    int max_length, void *value, int *att_length)
	int miattget1(int cdfid, int varid, char *name, nc_type datatype,
		                     void *value)
	char *miattgetstr(int cdfid, int varid, char *name, 
		                         int maxlen, char *value)
	int miattputint(int cdfid, int varid, char *name, int value)
	int miattputdbl(int cdfid, int varid, char *name, double value)
	int miattputstr(int cdfid, int varid, char *name, char *value)
	int mivarget(int cdfid, int varid, long start[], long count[],
		                    nc_type datatype, char *sign, void *values)
	int mivarget1(int cdfid, int varid, long mindex[],
	                     nc_type datatype, char *sign, void *value)
	int mivarput(int cdfid, int varid, long start[], long count[],
	                    nc_type datatype, char *sign, void *values)
	int mivarput1(int cdfid, int varid, long mindex[],
	                     nc_type datatype, char *sign, void *value)
	long *miset_coords(int nvals, long value, long coords[])
	long *mitranslate_coords(int cdfid, 
	                                int invar,  long incoords[],
	                                int outvar, long outcoords[])
	int micopy_all_atts(int incdfid, int invarid, 
	                           int outcdfid, int outvarid)
	int micopy_var_def(int incdfid, int invarid, int outcdfid)
	int micopy_var_values(int incdfid, int invarid, 
	                             int outcdfid, int outvarid)
	int micopy_all_var_defs(int incdfid, int outcdfid, int nexclude,
	                               int excluded_vars[])
	int micopy_all_var_values(int incdfid, int outcdfid, int nexclude,
	                                 int excluded_vars[])
	char *micreate_tempfile()

	# From minc_convenience.c
	int miget_datatype(int cdfid, int imgid, 
	                          nc_type *datatype, int *is_signed)
	int miget_default_range(nc_type datatype, int is_signed, 
	                               double default_range[])
	int miget_valid_range(int cdfid, int imgid, double valid_range[])
	int miset_valid_range(int cdfid, int imgid, double valid_range[])
	int miget_image_range(int cdfid, double image_range[])
	int mivar_exists(int cdfid, char *varname)
	int miattput_pointer(int cdfid, int varid, char *name, int ptrvarid)
	int miattget_pointer(int cdfid, int varid, char *name)
	int miadd_child(int cdfid, int parent_varid, int child_varid)
	int micreate_std_variable(int cdfid, char *name, nc_type datatype, 
	                                 int ndims, int dim[])
	int micreate_group_variable(int cdfid, char *name)
	const_char_ptr miget_version()
	int miappend_history(int fd, const_char_ptr tm_stamp)

	int miicv_create()
	int miicv_free(int icvid)
	int miicv_setdbl(int icvid, int icv_property, double value)
	int miicv_setint(int icvid, int icv_property, int value)
	int miicv_setlong(int icvid, int icv_property, long value)
	int miicv_setstr(int icvid, int icv_property, char *value)
	int miicv_inqdbl(int icvid, int icv_property, double *value)
	int miicv_inqint(int icvid, int icv_property, int *value)
	int miicv_inqlong(int icvid, int icv_property, long *value)
	int miicv_inqstr(int icvid, int icv_property, char *value)
	int miicv_ndattach(int icvid, int cdfid, int varid)
	int miicv_detach(int icvid)
	int miicv_get(int icvid, long start[], long count[], void *values)
	int miicv_put(int icvid, long start[], long count[], void *values)
	
	int miicv_attach(int icvid, int cdfid, int varid)
