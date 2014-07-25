cdef extern from *:
	ctypedef char* const_char_ptr "const char*"
	

cdef extern from "string.h":
	void* memcpy(void *s1, void *s2, size_t n)

cdef extern from "stdlib.h":
	void *malloc(size_t size)
	void free(void *ptr)
