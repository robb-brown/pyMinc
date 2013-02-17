from volume_io cimport VIO_Volume, VIO_General_transform

cdef extern from "libminctracc.h":

	int minctracc ( VIO_Volume source, VIO_Volume target, VIO_Volume sourceMask, VIO_Volume targetMask, VIO_General_transform *initialXFM, VIO_General_transform *outputXFM, int argc, char* argv[] )



