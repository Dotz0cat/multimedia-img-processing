#ifndef PROCESSING_PRIVATE_H
#define PROCESSING_PRIVATE_H

#include "processing.h"

#include <math.h>
#include <stdbit.h>

struct edit {
	
	size_t loc;

	double new_value;
};

static void apply_edit(double *restrict image_arr, struct edit *restrict edit_to_make);

static void noise_reduction_lchuv(struct jpeg_img_data *restrict img);
static void noise_reduction_rgb(struct jpeg_img_data *restrict img);

static void noise_detection_lchuv(struct jpeg_img_data *restrict img);
static void noise_detection_rgb(struct jpeg_img_data *restrict img);

#endif /* PROCESSING_PRIVATE_H */