#ifndef COLOR_H
#define COLOR_H

#include <stdlib.h>
#include <math.h>

#include <openblas64/cblas.h>
#include <jpeglib.h>

#include "image.h"

double * convert_to_lchuv_space(struct jpeg_img_data *restrict img);
uint8_t * rgb_space_from_lchuv(double *lchuv, size_t length);

#endif /* COLOR_H */
