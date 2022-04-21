#include "stdio.h"
#include "ISO_Fortran_binding.h"

void get_stride(CFI_cdesc_t *a, int *inca) {


	printf("sm %ld\n",a->dim[0].sm);
	printf("elem_len %ld\n",a->elem_len);
	
	*inca = (int) (a->dim[0].sm / a->elem_len);


}