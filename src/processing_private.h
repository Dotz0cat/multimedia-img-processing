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

static inline void compute_mask(uint8_t *restrict mask, const size_t i, const size_t j, const size_t height, const size_t width);

static inline void select_neighbors_double(const uint8_t mask, double *restrict neighbors, const double *restrict image, const size_t width, const size_t i, const size_t j);
static inline void select_neighbors_uint8_t(const uint8_t mask, double *restrict neighbors, const uint8_t *restrict image, const size_t width, const size_t i, const size_t j);

static inline void add_inital_center_to_stat_lchuv(double *restrict stat, const double *restrict image, const size_t width, const size_t i, const size_t j);
static inline void calc_8_neighbor_average_lchuv(const uint8_t mask, const double *restrict neighbors, double *restrict average);
static inline void calc_stdev_lchuv(const uint8_t mask, const double *restrict neighbors, const double *restrict average, double *restrict stdev);

static inline void add_inital_center_to_stat_rgb(double *restrict stat, const uint8_t *restrict image, const size_t width, const size_t i, const size_t j);
static inline void calc_8_neighbor_average_rgb(const uint8_t mask, const double *restrict neighbors, double *restrict average);
static inline void calc_stdev_rgb(const uint8_t mask, const double *restrict neighbors, const double *restrict average, double *restrict stdev);

#endif /* PROCESSING_PRIVATE_H */