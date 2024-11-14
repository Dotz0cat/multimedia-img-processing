#include "processing.h"

void apply_edit(double *image_arr, struct edit *edit_to_make) {

	if (image_arr == NULL) return;
	if (edit_to_make == NULL) return;

	image_arr[edit_to_make->loc] = edit_to_make->new_value;

	return;
	
}

void noise_reduction(double *image, size_t length) {
	return;
}
