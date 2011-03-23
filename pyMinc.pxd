cdef extern from *:
	ctypedef char* const_char_ptr "const char*"

cdef extern from "string.h":
	void* memcpy(void *s1, void *s2, size_t n)

cdef extern from "stdlib.h":
	void free(void *ptr)


cdef extern from "netcdf.h":

	ctypedef int nc_type
	int ncvarid(int ncid, const_char_ptr name)
	int ncvarinq(int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp, int *dimidsp, int *nattsp)
	int ncdiminq(int ncid, int dimid, char *name, long *lenp)	


cdef extern from "minc.h":
	
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

