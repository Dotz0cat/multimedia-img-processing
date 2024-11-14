#ifndef COLOR_PRIVATE_H
#define COLOR_PRIVATE_H

#include "color.h"

#define CIE_EPLSION (216.0 / 24389.0)
#define CIE_KAPPA (24389.0 / 27.0)

static const double D65[] = {0.95047, 1.00, 1.08883};
static const double D50[] = {0.9642, 1.00, 0.8521};

static double * lchuv_space_from_ycbcr(struct jpeg_img_data *img);
static double * lchuv_space_from_grayscale(struct jpeg_img_data *img);
static double * lchuv_space_from_rgb(struct jpeg_img_data *img);

static void xyz_to_luv(double *xyz, const double ref_white[3], size_t length);
static void luv_to_lchuv(double *luv, size_t length);
static void rgb_to_xyz(const uint8_t *rgb_data, double *xyz, size_t length);
static void lchuv_to_luv(double *lchuv, size_t length);
static void luv_to_xyz(double *luv, const double ref_white[3], size_t length);
static void xyz_to_rgb(double *xyz, uint8_t *rgb, size_t length);

#endif /* COLOR_PRIVATE_H */
