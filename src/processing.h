#ifndef PROCESSING_H
#define PROCESSING_H

#include <stdlib.h>

#include "image.h"

void process_image(struct jpeg_img_data *restrict img);

//void noise_reduction(double *restrict image, size_t length);
void noise_detection(struct jpeg_img_data *restrict img);

double * calculate_avg_lchuv(const double *restrict image, size_t length);
double * calculate_avg_rgb(const uint8_t *restrict image, size_t length);

double * calculate_mse_lchuv(const double *restrict image, const double *restrict average, size_t length);
double * calculate_mse_rgb(const uint8_t *restrict image, const double *restrict average, size_t length);

double * find_max_lchuv(const double *restrict image, size_t length);
uint8_t * find_max_rgb(const uint8_t *restrict image, size_t length);

double * psnr_lchuv(const double *restrict max, const double *restrict mse);
double * psnr_rgb(const uint8_t *restrict max, const double *restrict mse);

double * snr(const double *restrict mean, const double *restrict mse);

#endif /* PROCESSING_H */
