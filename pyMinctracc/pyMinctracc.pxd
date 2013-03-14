from volume_io cimport VIO_Volume, VIO_General_transform, VIO_N_DIMENSIONS, VIO_Real

cdef extern from "point_vector.h":

	ctypedef  struct PointR:
		VIO_Real   coords[VIO_N_DIMENSIONS]

	ctypedef  struct VectorR:
		VIO_Real   coords[VIO_N_DIMENSIONS]
		
		
cdef extern from "constants.h":
	
	ctypedef enum:
		TRANS_LSQ         =1
		TRANS_LSQ3        =2
		TRANS_LSQ6        =3
		TRANS_LSQ7        =4
		TRANS_LSQ9        =5
		TRANS_LSQ10       =6
		TRANS_LSQ12       =7
		TRANS_PAT         =8
		TRANS_NONLIN      =9
		TRANS_IDENT       =10
		TRANS_ROT 			= 0
		TRANS_QUAT =  1
		
		NONLIN_XCORR     =     0
		NONLIN_DIFF      =     1
		NONLIN_LABEL      =    2
		NONLIN_CHAMFER     =   3
		NONLIN_OPTICALFLOW  =  4
		NONLIN_CORRCOEFF    =  5
		NONLIN_SQDIFF       =  6

		OPT_SIMPLEX     =  0
		OPT_BFGS	=	1
		
	

	

cdef extern from "arg_data.h":

	ctypedef  Arg_Data_struct Arg_Data

	ctypedef float (*Objective_Function) (VIO_Volume, VIO_Volume, VIO_Volume, VIO_Volume, Arg_Data*)

	ctypedef void (*Transform_Function)(PointR *result, VIO_General_transform *trans_data, PointR *coordinate)

	ctypedef int (*Interpolating_Function)(VIO_Volume volume, PointR *coord, double *result)

	ctypedef enum Interpolating_Type:
		TRILINEAR,
		TRICUBIC,
		N_NEIGHBOUR,

	ctypedef enum Objective_Type: 
		XCORR, 
		ZSCORE, 
		SSC, 
		VR, 
		MUTUAL_INFORMATION, 
		NORMALIZED_MUTUAL_INFORMATION


	ctypedef struct Program_Flags:
		int verbose
		int debug

	ctypedef struct Transform_Flags:
		int estimate_center
		int estimate_scale
		int estimate_trans
		int estimate_rots
		int estimate_quats

	ctypedef struct Program_Filenames:
		char *data
		char *model
		char *mask_data
		char *mask_model
		char *output_trans
		char *measure_file
		char *matlab_file

	ctypedef struct Feature_volumes:
		int number_of_features
		VIO_Volume *data
		VIO_Volume *model
		VIO_Volume *data_mask
		VIO_Volume *model_mask
		char **data_name
		char **model_name
		char **mask_data_name
		char **mask_model_name
		char *obj_func
		VIO_Real *weight
		VIO_Real *thresh_data
		VIO_Real *thresh_model

	ctypedef struct Program_Transformation:
		int use_identity
		int use_default
		int use_magnitude
		VIO_Real max_def_magnitude
		int use_simplex
		int use_super
		int use_local_smoothing
		int use_local_isotropic
		char *file_name
		char *file_contents
		long buffer_length
		VIO_General_transform *transformation
		VIO_General_transform *orig_transformation
		int transform_type
		double center[3]         
		double scales[3]
		double shears[3]
		double translations[3]
		double quaternions[4]
		double rotations[3]
		double weights[12]       
		int invert_mapping_flag
		int rotation_type    

	ctypedef struct Arg_Data_struct:
		Program_Filenames      filenames
		Program_Flags          flags  
		Program_Transformation trans_info
		Feature_volumes        features  

		Interpolating_Function interpolant
		Interpolating_Type interpolant_type
		Objective_Function     obj_function
		Objective_Type    obj_function_type
		int                    optimize_type
		int                    force_lattice
		double                 step[3] 
		double                 lattice_width[3]
		double                 start[3]
		int                    count[3]
		VectorR                directions[3]
		int                    smallest_vol 
		Transform_Flags        trans_flags  
		double                 threshold[2]
		double                 speckle  
		int                    groups       
		int                    blur_pdf

			


cdef extern from "libminctracc.h":

	void initializeArgs(Arg_Data *args)
	VIO_General_transform* minctracc( VIO_Volume source, VIO_Volume target, VIO_Volume sourceMask, VIO_Volume targetMask, VIO_General_transform *initialXFM, int iterations, Arg_Data *args)



