#include "processing_private.h"

static void apply_edit(double *image_arr, struct edit *edit_to_make) {

	if (image_arr == NULL) return;
	if (edit_to_make == NULL) return;

	image_arr[edit_to_make->loc] = edit_to_make->new_value;

	return;
	
}

void noise_reduction(double *image, size_t length) {
	return;
}

void noise_detection(double *image, size_t length) {
	return;
}

double * calculate_avg_lchuv(const double *image, size_t length) {
	double *average = malloc(sizeof(double) * 3);

	double indvidual_length = (double) (length / 3UL);

	average[0] = 0.0;
	average[1] = 0.0;
	average[2] = 0.0;

	#pragma omp parallel for shared(image) reduction(+:average[0:3])
	for (size_t i = 0UL; i < length; i += 3UL) {
		average[0] += image[i + 0];
		average[1] += image[i + 1];
		average[2] += image[i + 2];
	}

	average[0] = average[0] / indvidual_length;
	average[1] = average[1] / indvidual_length;
	average[2] = average[2] / indvidual_length;

	return average;
}

double * calculate_avg_rgb(const uint8_t *image, size_t length) {
	uint64_t sum[3] = {0, 0, 0};

	double indvidual_length = (double) (length / 3UL);

	double *average = malloc(sizeof(double) * 3);

	#pragma omp parallel for shared(image) reduction(+:average[0:3])
	for (size_t i = 0UL; i < length; i += 3UL) {
		sum[0] += image[i + 0];
		sum[1] += image[i + 1];
		sum[2] += image[i + 2];
	}

	average[0] = (double) sum[0] / indvidual_length;
	average[1] = (double) sum[1] / indvidual_length;
	average[2] = (double) sum[2] / indvidual_length;

	return average;
}

double * calculate_mse_lchuv(const double *image, const double *average, size_t length) {

	double *mse = malloc(sizeof(double) * 3);

	double indvidual_length = (double) (length / 3UL);

	mse[0] = 0.0;
	mse[1] = 0.0;
	mse[2] = 0.0;

	#pragma omp parallel for shared(image, average) reduction(+:mse[0:3])
	for (size_t i = 0UL; i < length; i += 3UL) {
		mse[0] += pow(average[0] - image[i + 0], 2);
		mse[1] += pow(average[1] - image[i + 1], 2);
		mse[2] += pow(average[2] - image[i + 2], 2);
	}

	mse[0] /= indvidual_length;
	mse[1] /= indvidual_length;
	mse[2] /= indvidual_length;

	return mse;
}

double * calculate_mse_rgb(const uint8_t *image, const double *average, size_t length) {

	double *mse = malloc(sizeof(double) * 3);

	double indvidual_length = (double) (length / 3UL);

	mse[0] = 0.0;
	mse[1] = 0.0;
	mse[2] = 0.0;

	#pragma omp parallel for shared(image, average) reduction(+:mse[0:3])
	for (size_t i = 0UL; i < length; i += 3UL) {
		mse[0] += pow(average[0] - (double) image[i + 0], 2);
		mse[1] += pow(average[1] - (double) image[i + 1], 2);
		mse[2] += pow(average[2] - (double) image[i + 2], 2);
	}

	mse[0] /= indvidual_length;
	mse[1] /= indvidual_length;
	mse[2] /= indvidual_length;

	return mse;
}

double * find_max_lchuv(const double *image, size_t length) {
	double *max = malloc(sizeof(double) * 3);

	#pragma omp parallel for shared(image) reduction(max:max[0:3])
	for (size_t i = 0UL; i < length; i += 3UL) {
		max[0] = image[i + 0] > max[0] ? image[i + 0] : max[0];
		max[1] = image[i + 1] > max[1] ? image[i + 1] : max[1];
		max[2] = image[i + 2] > max[2] ? image[i + 2] : max[2];
	}

	return max;
}

uint8_t * find_max_rgb(const uint8_t *image, size_t length) {
	uint8_t *max = malloc(sizeof(uint8_t) * 3);

	#pragma omp parallel for shared(image) reduction(max:max[0:3])
	for (size_t i = 0UL; i < length; i += 3UL) {
		max[0] = image[i + 0] > max[0] ? image[i + 0] : max[0];
		max[1] = image[i + 1] > max[1] ? image[i + 1] : max[1];
		max[2] = image[i + 2] > max[2] ? image[i + 2] : max[2];
	}

	return max;
}

double * psnr_lchuv(const double *max, const double *mse) {
	double *psnr = malloc(sizeof(double) * 3);

	psnr[0] = 20.0 * log10(max[0]) - 10.0 * log10(mse[0]);
	psnr[1] = 20.0 * log10(max[1]) - 10.0 * log10(mse[1]);
	psnr[2] = 20.0 * log10(max[2]) - 10.0 * log10(mse[2]);

	return psnr;
}

double * psnr_rgb(const uint8_t *max, const double *mse) {
	double *psnr = malloc(sizeof(double) * 3);

	psnr[0] = 20.0 * log10((double) max[0]) - 10.0 * log10(mse[0]);
	psnr[1] = 20.0 * log10((double) max[1]) - 10.0 * log10(mse[1]);
	psnr[2] = 20.0 * log10((double) max[2]) - 10.0 * log10(mse[2]);

	return psnr;
}

double * snr(const double *mean, const double *mse) {
	double *snr = malloc(sizeof(double) * 3);

	snr[0] = 20.0 * log10(mean[0]) - 20.0 * log10(mse[0]);
	snr[1] = 20.0 * log10(mean[1]) - 20.0 * log10(mse[1]);
	snr[2] = 20.0 * log10(mean[2]) - 20.0 * log10(mse[2]);

	return snr;
}
