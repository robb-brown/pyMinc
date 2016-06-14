from netCDF cimport nc_type

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
		
		
	void  delete_general_transform(VIO_General_transform   *transform )	
	
	int  get_n_concated_transforms(VIO_General_transform   *transform )
	
	void  general_transform_point(
		VIO_General_transform   *transform,
		VIO_Real                x,
		VIO_Real                y,
		VIO_Real                z,
		VIO_Real                *x_transformed,
		VIO_Real                *y_transformed,
		VIO_Real                *z_transformed )

	void  general_inverse_transform_point(
		VIO_General_transform   *transform,
		VIO_Real                x,
		VIO_Real                y,
		VIO_Real                z,
		VIO_Real                *x_transformed,
		VIO_Real                *y_transformed,
		VIO_Real                *z_transformed )

	# void  transform_or_invert_point(
	# 	VIO_General_transform   *transform,
	# 	VIO_BOOL             inverse_flag,
	# 	VIO_Real                x,
	# 	VIO_Real                y,
	# 	VIO_Real                z,
	# 	VIO_Real                *x_transformed,
	# 	VIO_Real                *y_transformed,
	# 	VIO_Real                *z_transformed )
		
		
	void  concat_general_transforms(
	    VIO_General_transform   *first,
	    VIO_General_transform   *second,
	    VIO_General_transform   *result )
			
		
		
		
	# vol_io_prototypes.h
	VIO_Status  output_transform_file(
	    char*              filename,
	    char*              comments,
	    VIO_General_transform   *transform )

	VIO_Status  input_transform_file(
	    char*              filename,
	    VIO_General_transform   *transform )
	


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
	
	void FREE(void *ptr)
	
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
