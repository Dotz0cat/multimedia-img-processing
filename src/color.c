#include "color_private.h"

double * convert_to_lchuv_space(struct jpeg_img_data *restrict img) {
	switch (img->color_space) {
		case JCS_YCbCr:
			return lchuv_space_from_ycbcr(img);
			break;
		case JCS_GRAYSCALE:
			return lchuv_space_from_grayscale(img);
			break;
		case JCS_RGB:
		default:
			return lchuv_space_from_rgb(img);
			break;
	}

	return NULL;
}

static double * lchuv_space_from_ycbcr(struct jpeg_img_data *restrict img) {
	// struct YCrCb *yuv = img->raw_data;

	// size_t length = (img->height * img->width);
	// struct rgb *rgb_from_y = malloc(sizeof(struct rgb) * length);

	// #pragma omp parallel for
	// for (size_t i = 0u; i < length; i++) {
	// 	rgb_from_y[i].r = 77 * yuv[i].Y + 150 * yuv[i].cr + 29* yuv[i];
	// 	rgb_from_y[i].g = yuv[i] + yuv[i] + yuv[i];
	// 	rgb_from_y[i].b = yuv[i] + yuv[i] + yuv[i];
	// }

	// void *img_data = malloc(sizeof(struct XYZ) * length);

	// bt601rgb_to_xyz(rgb_from_y, img_data, length);

	// free(rgb_from_y);

	// xyz_to_luv(img_data, D65, length);
	// luv_to_lchuv(img_data, length);

	// return img_data;

	return NULL;
}

static double * lchuv_space_from_grayscale(struct jpeg_img_data *restrict img) {
	// size_t length = (img->height * img->width) * 3;

	// uint8_t *img_rgb = malloc(length);

	// double *img_data = malloc(sizeof(double) * length);

	// #pragma omp parallel for
	// for (size_t i = 0u; i < length; i += 3) {
	// 	img_rgb[i + 0] = img->raw_data[i + 0];
	// 	img_rgb[i + 1] = img->raw_data[i + 0];
	// 	img_rgb[i + 2] = img->raw_data[i + 0];
	// }

	// bt601rgb_to_xyz(img_rgb, img_data, length);

	// xyz_to_luv(img_data, D65, length);
	// luv_to_lchuv(img_data, length);

	// return img_data;

	return NULL;
}

static double * lchuv_space_from_rgb(struct jpeg_img_data *restrict img) {
	size_t length = (img->height * img->width) * 3UL;

	double *restrict img_data = aligned_alloc(64, sizeof(double) * length);

	rgb_to_xyz(img->raw_data, img_data, length);

	xyz_to_luv(img_data, D65, length);
	//xyz_to_luv(img_data, D50, length);

	luv_to_lchuv(img_data, length);

	return img_data;
}

static void xyz_to_luv(double *restrict xyz, const double ref_white[3], size_t length) {
	const double ref_u = ( 4.0 * ref_white[0] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );
	const double ref_v = ( 9.0 * ref_white[1] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );

	#pragma omp parallel for shared(ref_u, ref_v, xyz)
	for (size_t i = 0UL; i < length; i += 3UL) {

		double y_r = xyz[i + 0UL] / ref_white[0];
		
		//double L = 0.0;
		double L = (y_r > CIE_EPLSION) ? 116.0 * cbrt(y_r) - 16.0 : CIE_KAPPA * y_r;

		// if (y_r > CIE_EPLSION) {
		// 	L = 116.0 * cbrt(y_r) - 16.0;
		// }
		// else {
		// 	L = CIE_KAPPA * y_r;
		// }

		double u_prime = ( 4.0 * xyz[i + 0UL] ) / ( xyz[i + 0UL] + 15.0 * xyz[i + 1UL] + 3.0 * xyz[i + 2UL] );
		double v_prime = ( 9.0 * xyz[i + 1UL] ) / ( xyz[i + 0UL] + 15.0 * xyz[i + 1UL] + 3.0 * xyz[i + 2UL] );

		double u = 13.0 * L * (u_prime - ref_u);
		double v = 13.0 * L * (v_prime - ref_v);

		xyz[i + 0UL] = L;
		xyz[i + 1UL] = u;
		xyz[i + 2UL] = v;
	}
}

#ifndef M_PI
#define M_PI 4 * atan2(1,1)
#endif /* M_PI */

static void luv_to_lchuv(double *restrict luv, size_t length) {
	//L componet stays the same

	#pragma omp parallel for shared(luv)
	for (size_t i = 0UL; i < length; i += 3UL) {
		double C = sqrt(pow(luv[i + 1UL], 2.0) + pow(luv[i + 2UL], 2.0));
		double H = atan2(luv[i + 2UL], luv[i + 1UL]);
		
		if (H <= 0.0) H = H + (2.0 * M_PI);

		luv[i + 1UL] = C;
		luv[i + 2UL] = H;
	}
}

static void rgb_to_xyz(const uint8_t *restrict rgb_data, double *restrict xyz, size_t length) {

	//initalize the xyz array with the rgb values
	#pragma omp parallel for
	for (size_t i = 0UL; i < length; i += 3UL) {
		xyz[i + 0UL] = rgb_data[i + 0UL] / 255.0;
		xyz[i + 1UL] = rgb_data[i + 1UL] / 255.0;
		xyz[i + 2UL] = rgb_data[i + 2UL] / 255.0;

		xyz[i + 0UL] = xyz[i + 0UL] <= 0.04045 ? xyz[i + 0UL] / 12.92 : pow((xyz[i + 0UL] + 0.055) / 1.055, 2.4);
		xyz[i + 1UL] = xyz[i + 1UL] <= 0.04045 ? xyz[i + 1UL] / 12.92 : pow((xyz[i + 1UL] + 0.055) / 1.055, 2.4);
		xyz[i + 2UL] = xyz[i + 2UL] <= 0.04045 ? xyz[i + 2UL] / 12.92 : pow((xyz[i + 2UL] + 0.055) / 1.055, 2.4);

		// xyz[i + 0] = xyz[i + 0] <= 0.08 ? 100.0 * xyz[i + 0] / CIE_KAPPA : pow((xyz[i + 0] + 0.16) / 1.16, 3);
		// xyz[i + 1] = xyz[i + 1] <= 0.08 ? 100.0 * xyz[i + 1] / CIE_KAPPA : pow((xyz[i + 1] + 0.16) / 1.16, 3);
		// xyz[i + 2] = xyz[i + 2] <= 0.08 ? 100.0 * xyz[i + 2] / CIE_KAPPA : pow((xyz[i + 2] + 0.16) / 1.16, 3);
	}

	// sRGB
	static const double conversion_matrix[] = {
		0.4124564, 0.3575761, 0.1804375,
		0.2126729, 0.7151522, 0.0721750,
		0.0193339, 0.1191920, 0.9503041
	};

	static const double conversion_matrix_t[] = {
		0.4124564, 0.2126729, 0.0193339,
		0.3575761, 0.7151522, 0.1191920,
		0.1804375, 0.0721750, 0.9503041
	};

	
	//double *rgb_value;
	double xyz_output_loc[3] = {0.0, 0.0, 0.0};

	#pragma omp parallel for private(xyz_output_loc) shared(conversion_matrix, xyz)
	for (size_t i = 0u; i < length; i += 3u) {
		//rgb_value = &xyz[i];

		xyz_output_loc[0] = conversion_matrix[0] * xyz[i + 0] + conversion_matrix[1] * xyz[i + 1] + conversion_matrix[2] * xyz[i + 2];
		xyz_output_loc[1] = conversion_matrix[3] * xyz[i + 0] + conversion_matrix[4] * xyz[i + 1] + conversion_matrix[5] * xyz[i + 2];
		xyz_output_loc[2] = conversion_matrix[6] * xyz[i + 0] + conversion_matrix[7] * xyz[i + 1] + conversion_matrix[8] * xyz[i + 2];

		//cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, conversion_matrix, 3, (xyz + i), 1, 0.0, xyz_output_loc, 1);
		//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0, conversion_matrix, 3, rgb_value, 1, 0.0, xyz_output_loc, 1);

		xyz[i + 0] = xyz_output_loc[0];
		xyz[i + 1] = xyz_output_loc[1];
		xyz[i + 2] = xyz_output_loc[2];
	}

	// double *extra_xyz = aligned_alloc(64, sizeof(double) * length);

	// cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, (length / 3UL), 3, 1.0, conversion_matrix_t, 3, xyz, 3, 0.0, extra_xyz, (length / 3UL));

	// memcpy(xyz, extra_xyz, length);

	// free(extra_xyz);
	
}

uint8_t * rgb_space_from_lchuv(double *restrict img_data, size_t length) {
	uint8_t *restrict rgb_data = malloc(length);

	//double *img_data = lchuv;

	lchuv_to_luv(img_data, length);

	luv_to_xyz(img_data, D65, length);
	//luv_to_xyz(img_data, D50, length);

	xyz_to_rgb(img_data, rgb_data, length);

	return rgb_data;
}

static void lchuv_to_luv(double *restrict lchuv, size_t length) {

	#pragma omp parallel for
	for (size_t i = 0UL; i < length; i += 3UL) {
		//L stays the same
		double u = lchuv[i + 1UL] * cos(lchuv[i + 2UL]);
		double v = lchuv[i + 1UL] * sin(lchuv[i + 2UL]);

		lchuv[i + 1UL] = u;
		lchuv[i + 2UL] = v;
	}
}

static void luv_to_xyz(double *restrict luv, const double ref_white[3], size_t length) {
	const double ref_u = ( 4.0 * ref_white[0] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );
	const double ref_v = ( 9.0 * ref_white[1] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );

	#pragma omp parallel for
	for (size_t i = 0UL; i < length; i += 3UL) {
		double Y = (luv[i + 0UL] > (CIE_KAPPA * CIE_EPLSION)) ? pow( ( luv[i + 0UL] + 16.0 ) / 116.0, 3.0) : luv[i + 0UL] / CIE_KAPPA;
		//double Y = 0.0;

		// if (luv[i + 0UL] > (CIE_KAPPA * CIE_EPLSION)) {
		// 	Y = pow( (luv[i + 0UL] + 16.0) / 116.0, 3.0);
		// }
		// else {
		// 	Y = luv[i + 0UL] / CIE_KAPPA;
		// }

		double a = ( 1.0 / 3.0 ) * ( ( 52.0 * luv[i + 0UL] ) / ( luv[i + 1UL] + 13.0 * luv[i + 0UL] * ref_u ) - 1.0 );
		double b = -5.0 * Y;
		double c = -( 1.0 / 3.0 );
		double d = Y * ( ( 39.0 * luv[i + 0UL] ) / ( luv[i + 2UL] + 13.0 * luv[i + 0UL] * ref_v ) - 5.0 );

		double X = ( d - b ) / ( a - c );
		double Z = X * a + b;

		// X = X < 0.0 ? 0.0 : X > 1.0 ? 1.0 : X;
		// Y = Y < 0.0 ? 0.0 : Y > 1.0 ? 1.0 : Y;
		// Z = Z < 0.0 ? 0.0 : Z > 1.0 ? 1.0 : Z;

		luv[i + 0UL] = X;
		luv[i + 1UL] = Y;
		luv[i + 2UL] = Z;
	}
}

static void xyz_to_rgb(double *restrict xyz, uint8_t *restrict rgb, size_t length) {

	//pal/secam / bt601 D65
	// const double conversion_matrix[] = {
	// 	 3.0628971, -1.3931791, -0.4757517,
	// 	-0.9692660,  1.8760108,  0.0415560,
	// 	 0.0678775, -0.2288548,  1.0693490
	// };

	//pal/secam / bt601 D50
	// const double conversion_matrix[] = {
	// 	 2.9603944, -1.4678519, -0.4685105,
	// 	-0.9787684,  1.9161415,  0.0334540,
	// 	 0.0844874, -0.2545973,  1.4216174
	// };

	// sRGB D65
	const double conversion_matrix[] = {
		 3.2404542, -1.5371385, -0.4985314,
		-0.9692660,  1.8760108,  0.0415560,
		 0.0556434, -0.2040259,  1.0572252
	};

	// const double conversion_matrix[] = {
	// 	 3.1338561, -1.6168667, -0.4906146,
	// 	-0.9787684,  1.9161415,  0.0334540,
	// 	 0.0719453, -0.2289914,  1.4052427
 // 	};

	// const double conversion_matrix[] = {
	// 	-1.28873e-6, 10.3463, -62.0777,
	// 	-0.699122, 0.885554, 5.1735,
	// 	6.92711, -25.3994, 131.615
	// };

	// const double conversion_matrix[] = {
	// 	 2.0413690, -0.5649464, -0.3446944,
	// 	-0.9692660,  1.8760108,  0.0415560,
	// 	 0.0134474, -0.1183897,  1.0154096
	// };

	
	double rgb_value[3] = {0.0, 0.0, 0.0};
	//double *xyz_value;

	#pragma omp parallel for private(rgb_value) shared(rgb, conversion_matrix)
	for (size_t i = 0UL; i < length; i += 3UL) {

		//xyz_value = &xyz[i];

		//cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, conversion_matrix, 3, (xyz + i), 1, 0, rgb_value, 1);
		//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0, conversion_matrix, 3, xyz_value, 1, 0.0, rgb_value, 1);

		rgb_value[0] = conversion_matrix[0] * xyz[i + 0] + conversion_matrix[1] * xyz[i + 1] + conversion_matrix[2] * xyz[i + 2];
		rgb_value[1] = conversion_matrix[3] * xyz[i + 0] + conversion_matrix[4] * xyz[i + 1] + conversion_matrix[5] * xyz[i + 2];
		rgb_value[2] = conversion_matrix[6] * xyz[i + 0] + conversion_matrix[7] * xyz[i + 1] + conversion_matrix[8] * xyz[i + 2];

		rgb_value[0] = rgb_value[0] <= 0.0031308 ? rgb_value[0] * 12.92 : 1.055 * ( pow(rgb_value[0], 1.0 / 2.4) ) - 0.055;
		rgb_value[1] = rgb_value[1] <= 0.0031308 ? rgb_value[1] * 12.92 : 1.055 * ( pow(rgb_value[1], 1.0 / 2.4) ) - 0.055;
		rgb_value[2] = rgb_value[2] <= 0.0031308 ? rgb_value[2] * 12.92 : 1.055 * ( pow(rgb_value[2], 1.0 / 2.4) ) - 0.055;

		// rgb_value[0] = rgb_value[0] <= CIE_EPLSION ? ( rgb_value[0] * CIE_KAPPA ) / 100.0 : 1.16 * cbrt(rgb_value[0]) - 0.16;
		// rgb_value[1] = rgb_value[1] <= CIE_EPLSION ? ( rgb_value[1] * CIE_KAPPA ) / 100.0 : 1.16 * cbrt(rgb_value[1]) - 0.16;
		// rgb_value[2] = rgb_value[2] <= CIE_EPLSION ? ( rgb_value[2] * CIE_KAPPA ) / 100.0 : 1.16 * cbrt(rgb_value[2]) - 0.16;

		rgb_value[0] = rint(rgb_value[0] * 255.0);
		rgb_value[1] = rint(rgb_value[1] * 255.0);
		rgb_value[2] = rint(rgb_value[2] * 255.0);

		// rgb[i + 0] = (uint8_t) fmin(fmax(rgb_value[0], 0.0), 255.0);
		// rgb[i + 1] = (uint8_t) fmin(fmax(rgb_value[1], 0.0), 255.0);
		// rgb[i + 2] = (uint8_t) fmin(fmax(rgb_value[2], 0.0), 255.0);

		rgb[i + 0UL] = (uint8_t) (rgb_value[0] <= 0.0) ? 0U : (rgb_value[0] >= 255.0) ? 255U : rgb_value[0];
		rgb[i + 1UL] = (uint8_t) (rgb_value[1] <= 0.0) ? 0U : (rgb_value[1] >= 255.0) ? 255U : rgb_value[1];
		rgb[i + 2UL] = (uint8_t) (rgb_value[2] <= 0.0) ? 0U : (rgb_value[2] >= 255.0) ? 255U : rgb_value[2];
	}
}
