#ifndef PROCESSING_H
#define PROCESSING_H

#include <stdlib.h>

#include "image.h"

void noise_reduction(double *image, size_t length);
void noise_detection(double *image, size_t length);

double * calculate_avg_lchuv(const double *image, size_t length);
double * calculate_avg_rgb(const uint8_t *image, size_t length);

double * calculate_mse_lchuv(const double *image, const double *average, size_t length);
double * calculate_mse_rgb(const uint8_t *image, const double *average, size_t length);

double * find_max_lchuv(const double *image, size_t length);
uint8_t * find_max_rgb(const uint8_t *image, size_t length);

double * psnr_lchuv(const double *max, const double *mse);
double * psnr_rgb(const uint8_t *max, const double *mse);

double * snr(const double *mean, const double *mse);

#endif /* PROCESSING_H */
