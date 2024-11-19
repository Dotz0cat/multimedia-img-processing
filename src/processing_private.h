#ifndef PROCESSING_PRIVATE_H
#define PROCESSING_PRIVATE_H

#include "processing.h"

#include <math.h>

struct edit {
	
	size_t loc;

	double new_value;
};

static void apply_edit(double *image_arr, struct edit *edit_to_make);

#endif /* PROCESSING_PRIVATE_H */