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

	cdef enum:
		MI_ORIGINAL_TYPE = 0
		NC_UNSPECIFIED = MI_ORIGINAL_TYPE
	
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
	


cdef extern from "volume_io.h":
	
	# Constants (these are #defines in the C code)
	cdef enum:
		VIO_MAX_DIMENSIONS   =  5
		VIO_N_DIMENSIONS 	= 3
		VIO_X = 0
		VIO_Y = 1
		VIO_Z = 2
		FALSE = 0
		TRUE = 1
		
	
	# Basic types
	ctypedef char *VIO_STR
	ctypedef int VIO_BOOL
	ctypedef double VIO_Real
	ctypedef signed char VIO_SCHAR
	ctypedef unsigned char VIO_UCHAR
	ctypedef enum VIO_Status:
		OK,
		ERROR,
		INTERNAL_ERROR,
		END_OF_FILE,
		QUIT


	ctypedef  enum  VIO_Data_types:
		NO_DATA_TYPE
		UNSIGNED_BYTE
		SIGNED_BYTE
		UNSIGNED_SHORT
		SIGNED_SHORT
		UNSIGNED_INT
		SIGNED_INT
		FLOAT
		DOUBLE
		MAX_DATA_TYPE
		
		
	ctypedef  struct minc_output_options :
		VIO_Real global_image_range[2]
		VIO_STR  dimension_names[VIO_MAX_DIMENSIONS]
		VIO_BOOL use_starts_set
		VIO_BOOL use_volume_starts_and_steps
		
	
	# Geometry.h
	ctypedef  float   VIO_Point_coord_type

	ctypedef  struct VIO_Point:
		VIO_Point_coord_type   coords[VIO_N_DIMENSIONS]

	ctypedef  struct VIO_Vector:
		VIO_Point_coord_type   coords[VIO_N_DIMENSIONS]


	# transforms.h
	ctypedef  double  VIO_Transform_elem_type
	
	ctypedef  void   (*VIO_User_transform_function)( void  *user_data,
												VIO_Real  x,
												VIO_Real  y,
												VIO_Real  z,
												VIO_Real  *x_trans,
												VIO_Real  *y_trans,
												VIO_Real  *z_trans )
	

	ctypedef  struct VIO_Transform_2d:
		VIO_Transform_elem_type    m2d[2][3]

	ctypedef  struct VIO_Transform:
		VIO_Transform_elem_type    m[4][4]

		
	ctypedef enum VIO_Transform_types:
		LINEAR, THIN_PLATE_SPLINE, USER_TRANSFORM,
		CONCATENATED_TRANSFORM, GRID_TRANSFORM
		

	ctypedef struct VIO_General_transform:
		VIO_Transform_types         type
		VIO_BOOL                    inverse_flag

		#// --- linear transform 

		VIO_Transform               *linear_transform
		VIO_Transform               *inverse_linear_transform

		#/* --- non-linear transform */

		int                         n_points
		int                         n_dimensions
		VIO_Real                    **points
		VIO_Real                    **displacements   #/* n_points + n_dim + 1 by */
		                                             #/* n_dim */

		#/* --- grid transform */

		void                        *displacement_volume

		#/* --- user_defined */

		void                        *user_data
		size_t                      size_user_data
		VIO_User_transform_function     user_transform_function
		VIO_User_transform_function     user_inverse_transform_function

		#/* --- concatenated transform */

		int                         n_transforms
		VIO_General_transform    *transforms


	# multidim.h
	ctypedef  struct VIO_multidim_array:
		int                     n_dimensions
		int                     sizes[VIO_MAX_DIMENSIONS]
		VIO_Data_types          data_type
		void                    *data
		
		
	int  get_type_size(
	    VIO_Data_types   type )
		
		
			
	# volume_cache.h
	ctypedef  struct  VIO_cache_block_struct:
		int                         block_index
		VIO_SCHAR                modified_flag
		VIO_multidim_array              array
		VIO_cache_block_struct  *prev_used
		VIO_cache_block_struct  *next_used
		VIO_cache_block_struct  **prev_hash
		VIO_cache_block_struct  *next_hash


	ctypedef  struct VIO_cache_lookup_struct:
		int       block_index_offset
		int       block_offset


	ctypedef struct VIO_volume_cache_struct:
		int                         n_dimensions
		int                         file_offset[VIO_MAX_DIMENSIONS]
		VIO_STR                     input_filename

		VIO_STR                     output_filename
		nc_type                     file_nc_data_type
		VIO_BOOL                    file_signed_flag
		VIO_Real                    file_voxel_min
		VIO_Real                    file_voxel_max
		VIO_STR                     original_filename
		VIO_STR                     history
		minc_output_options         options

		VIO_BOOL                    writing_to_temp_file
		int                         total_block_size
		int                         block_sizes[VIO_MAX_DIMENSIONS]
		int                         blocks_per_dim[VIO_MAX_DIMENSIONS]
		VIO_BOOL                    output_file_is_open
		VIO_BOOL                    must_read_blocks_before_use
		void                        *minc_file
		int                         n_blocks
		int                         max_cache_bytes
		int                         max_blocks
		int                         hash_table_size
		VIO_cache_block_struct      *head
		VIO_cache_block_struct      *tail
		VIO_cache_block_struct      **hash_table

		VIO_cache_lookup_struct     *lookup[VIO_MAX_DIMENSIONS]
		VIO_cache_block_struct      *previous_block
		int                         previous_block_index

		VIO_BOOL                    debugging_on
		int                         n_accesses
		int                         output_every
		int                         n_hits
		int                         n_prev_hits
		
		

	# volume.h

		
	ctypedef  struct volume_struct:
		VIO_BOOL                is_cached_volume
		VIO_volume_cache_struct cache

		VIO_multidim_array      array

		VIO_STR                 dimension_names[VIO_MAX_DIMENSIONS]
		int                     spatial_axes[VIO_N_DIMENSIONS]
		nc_type                 nc_data_type
		VIO_BOOL                signed_flag
		VIO_BOOL                is_rgba_data

		VIO_Real                voxel_min
		VIO_Real                voxel_max
		VIO_BOOL                real_range_set
		VIO_Real                real_value_scale
		VIO_Real                real_value_translation

		VIO_Real                separations[VIO_MAX_DIMENSIONS]
		VIO_Real                starts[VIO_MAX_DIMENSIONS]
		VIO_Real                direction_cosines[VIO_MAX_DIMENSIONS][VIO_N_DIMENSIONS]

		VIO_BOOL                voxel_to_world_transform_uptodate
		VIO_General_transform   voxel_to_world_transform

		VIO_STR                 coordinate_system_name

		VIO_Real               *irregular_starts[VIO_MAX_DIMENSIONS]
		VIO_Real               *irregular_widths[VIO_MAX_DIMENSIONS]


	ctypedef  volume_struct  *VIO_Volume
		

	ctypedef  struct minc_input_options:
		VIO_BOOL    promote_invalid_to_zero_flag
		VIO_BOOL    convert_vector_to_scalar_flag
		VIO_BOOL    convert_vector_to_colour_flag
		int         dimension_size_for_colour_data
		int         max_dimension_size_for_colour_data
		int         rgba_indices[4]
		double      user_real_range[2]
	

	void ALLOC(void *ptr, int N)
	
	char** get_default_dim_names(int    n_dimensions )


	VIO_Status  input_volume(
	    char                 filename[],
	    int                  n_dimensions,
	    char                 *dim_names[],
	    nc_type              volume_nc_data_type,
	    VIO_BOOL              volume_signed_flag,
	    VIO_Real                 volume_voxel_min,
	    VIO_Real                 volume_voxel_max,
	    VIO_BOOL              create_volume_flag,
	    VIO_Volume               *volume,
	    minc_input_options   *options )
	
	
	VIO_Status  output_volume(
	    char                filename[],
	    nc_type               file_nc_data_type,
	    VIO_BOOL               file_signed_flag,
	    VIO_Real                  file_voxel_min,
	    VIO_Real                  file_voxel_max,
	    VIO_Volume                volume,
	    char                history[],
	    minc_output_options   *options )
	

	VIO_Volume   create_volume(
	    int         n_dimensions,
	    char*      dimension_names[],
	    nc_type     nc_data_type,
	    VIO_BOOL     signed_flag,
	    VIO_Real        voxel_min,
	    VIO_Real        voxel_max )
	

	void  alloc_volume_data(
	    VIO_Volume   volume )

	void  set_volume_type(
	    VIO_Volume       volume,
	    nc_type      nc_data_type,
	    VIO_BOOL      signed_flag,
	    VIO_Real         voxel_min,
	    VIO_Real         voxel_max )

	nc_type  get_volume_nc_data_type(
	    VIO_Volume       volume,
	    VIO_BOOL      *signed_flag )

	VIO_Data_types  get_volume_data_type(
	    VIO_Volume       volume )


	VIO_BOOL  volume_is_alloced(
	    VIO_Volume   volume )

	void  free_volume_data(
	    VIO_Volume   volume )

	void  delete_volume(
	    VIO_Volume   volume )

	int  get_volume_n_dimensions(
	    VIO_Volume   volume )

	void  get_volume_sizes(
	    VIO_Volume   volume,
	    int      sizes[] )

	void  set_volume_sizes(
	    VIO_Volume       volume,
	    int          sizes[] )

	unsigned int  get_volume_total_n_voxels(
	    VIO_Volume    volume )
	
	void  assign_voxel_to_world_transform(
	    VIO_Volume             volume,
	    VIO_General_transform  *transform )
	
	void  compute_world_transform(
	    int                 spatial_axes[VIO_N_DIMENSIONS],
	    VIO_Real                separations[],
	    VIO_Real                direction_cosines[][VIO_N_DIMENSIONS],
	    VIO_Real                starts[],
	    VIO_General_transform   *world_transform )


	void  check_recompute_world_transform(
	    VIO_Volume  volume )

	void  set_voxel_to_world_transform(
	    VIO_Volume             volume,
	    VIO_General_transform  *transform )

	VIO_General_transform  *get_voxel_to_world_transform(
	    VIO_Volume   volume )

	char** get_volume_dimension_names(
	    VIO_Volume   volume )

	void  get_volume_separations(
	    VIO_Volume   volume,
	    VIO_Real     separations[] )

	void  set_volume_separations(
	    VIO_Volume   volume,
	    VIO_Real     separations[] )

	void  set_volume_starts(
	    VIO_Volume  volume,
	    VIO_Real    starts[] )

	void  get_volume_starts(
	    VIO_Volume  volume,
	    VIO_Real    starts[] )


	void  set_volume_direction_unit_cosine(
	    VIO_Volume   volume,
	    int      axis,
	    VIO_Real     dir[] )


	void  set_volume_direction_cosine(
	    VIO_Volume   volume,
	    int      axis,
	    VIO_Real     dir[] )

	void  get_volume_direction_cosine(
	    VIO_Volume   volume,
	    int      axis,
	    VIO_Real     dir[] )

	void  set_volume_translation(
	    VIO_Volume  volume,
	    VIO_Real    voxel[],
	    VIO_Real    world_space_voxel_maps_to[] )


	void  get_volume_translation(
	    VIO_Volume  volume,
	    VIO_Real    voxel[],
	    VIO_Real    world_space_voxel_maps_to[] )


	void  convert_voxel_to_world(
	    VIO_Volume   volume,
	    VIO_Real     voxel[],
	    VIO_Real     *x_world,
	    VIO_Real     *y_world,
	    VIO_Real     *z_world )


	void  convert_world_to_voxel(
	    VIO_Volume   volume,
	    VIO_Real     x_world,
	    VIO_Real     y_world,
	    VIO_Real     z_world,
	    VIO_Real     voxel[] )
	
	
	void GET_VOXEL_PTR(void* src, VIO_Volume volume, int x, int y, int z, int t, int v)

		
#cdef extern from "cUtils.h":
#	void* getDataPtr(VIO_Volume volume, int x, int y, int z, int t, int v)
	
