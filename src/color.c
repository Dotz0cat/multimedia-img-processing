#include "color_private.h"

double * convert_to_lchuv_space(struct jpeg_img_data *img) {
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

static double * lchuv_space_from_ycbcr(struct jpeg_img_data *img) {
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

static double * lchuv_space_from_grayscale(struct jpeg_img_data *img) {
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

static double * lchuv_space_from_rgb(struct jpeg_img_data *img) {
	size_t length = (img->height * img->width) * 3;

	double *img_data = aligned_alloc(64, sizeof(double) * length);

	rgb_to_xyz(img->raw_data, img_data, length);

	xyz_to_luv(img_data, D65, length);
	//xyz_to_luv(img_data, D50, length);

	luv_to_lchuv(img_data, length);

	return img_data;
}

static void xyz_to_luv(double *xyz, const double ref_white[3], size_t length) {
	const double ref_u = ( 4.0 * ref_white[0] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );
	const double ref_v = ( 9.0 * ref_white[1] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );

	#pragma omp parallel for shared(ref_u, ref_v, xyz)
	for (size_t i = 0u; i < length; i += 3u) {

		double y_r = xyz[i + 0] / ref_white[0];
		double L = 0.0;

		if (y_r > CIE_EPLSION) {
			L = 116.0 * cbrt(y_r) - 16.0;
		}
		else {
			L = CIE_KAPPA * y_r;
		}

		double u_prime = ( 4.0 * xyz[i + 0] ) / ( xyz[i + 0] + 15.0 * xyz[i + 1] + 3.0 * xyz[i + 2] );
		double v_prime = ( 9.0 * xyz[i + 1] ) / ( xyz[i + 0] + 15.0 * xyz[i + 1] + 3.0 * xyz[i + 2] );

		double u = 13.0 * L * (u_prime - ref_u);
		double v = 13.0 * L * (v_prime - ref_v);

		xyz[i + 0] = L;
		xyz[i + 1] = u;
		xyz[i + 2] = v;
	}
}

#ifndef M_PI
#define M_PI 4 * atan2(1,1)
#endif /* M_PI */

static void luv_to_lchuv(double *luv, size_t length) {
	//L componet stays the same

	#pragma omp parallel for shared(luv)
	for (size_t i = 0u; i < length; i += 3u) {
		double C = sqrt(pow(luv[i + 1], 2.0) + pow(luv[i + 2], 2.0));
		double H = atan2(luv[i + 2], luv[i + 1]);
		
		if (H <= 0.0) H = H + (2 * M_PI);

		luv[i + 1] = C;
		luv[i + 2] = H;
	}
}

static void rgb_to_xyz(const uint8_t *rgb_data, double *xyz, size_t length) {

	//initalize the xyz array with the rgb values
	#pragma omp parallel for
	for (size_t i = 0u; i < length; i += 3u) {
		xyz[i + 0] = rgb_data[i + 0] / 255.0;
		xyz[i + 1] = rgb_data[i + 1] / 255.0;
		xyz[i + 2] = rgb_data[i + 2] / 255.0;

		xyz[i + 0] = xyz[i + 0] <= 0.04045 ? xyz[i + 0] / 12.92 : pow((xyz[i + 0] + 0.055) / 1.055, 2.4);
		xyz[i + 1] = xyz[i + 1] <= 0.04045 ? xyz[i + 1] / 12.92 : pow((xyz[i + 1] + 0.055) / 1.055, 2.4);
		xyz[i + 2] = xyz[i + 2] <= 0.04045 ? xyz[i + 2] / 12.92 : pow((xyz[i + 2] + 0.055) / 1.055, 2.4);

		// xyz[i + 0] = xyz[i + 0] <= 0.08 ? 100.0 * xyz[i + 0] / CIE_KAPPA : pow((xyz[i + 0] + 0.16) / 1.16, 3);
		// xyz[i + 1] = xyz[i + 1] <= 0.08 ? 100.0 * xyz[i + 1] / CIE_KAPPA : pow((xyz[i + 1] + 0.16) / 1.16, 3);
		// xyz[i + 2] = xyz[i + 2] <= 0.08 ? 100.0 * xyz[i + 2] / CIE_KAPPA : pow((xyz[i + 2] + 0.16) / 1.16, 3);
	}

	//pal/secam / bt601 D65
	// const double conversion_matrix[] = {
	// 	0.4306190, 0.3415419, 0.1783091,
	// 	0.2220379, 0.7066384, 0.0713236,
	// 	0.0201853, 0.1295504, 0.9390944
	// };

	// const double conversion_matrix[] = {
	// 	0.4552773,  0.3675500,  0.1413926,
	// 	0.2323025,  0.7077956,  0.0599019,
	// 	0.0145457,  0.1049154,  0.7057489
	// };

	// sRGB
	const double conversion_matrix[] = {
		0.4124564, 0.3575761, 0.1804375,
		0.2126729, 0.7151522, 0.0721750,
		0.0193339, 0.1191920, 0.9503041
	};

	// const double conversion_matrix[] = {
	// 	0.4360747,  0.3850649,  0.1430804,
	//  	0.2225045,  0.7168786,  0.0606169,
	//  	0.0139322,  0.0971045,  0.7141733
 // 	};

	// const double conversion_matrix[] = {
	// 	0.41238656, 0.35759149, 0.18045049,
	// 	0.21263682, 0.71518298, 0.0721802,
	// 	0.01933062, 0.11919716, 0.01203003
	// };

	//adobe rgb
	// const double conversion_matrix[] = {
	// 	0.5767309, 0.1855540, 0.1881852,
	// 	0.2973769, 0.6273491, 0.0752741,
	// 	0.0270343, 0.0706872, 0.9911085
 // 	};

	
	double *rgb_value;
	double xyz_output_loc[3] = {0.0, 0.0, 0.0};
	// 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	#pragma omp parallel for private(rgb_value, xyz_output_loc) shared(conversion_matrix, xyz)
	for (size_t i = 0u; i < length; i += 3u) {
		rgb_value = &xyz[i];

		//cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, conversion_matrix, 3, rgb_value, 1, 0.0, xyz_output_loc, 1);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0, conversion_matrix, 3, rgb_value, 1, 0.0, xyz_output_loc, 1);

		xyz[i + 0] = xyz_output_loc[0] < 0.0 ? 0.0 : xyz_output_loc[0] > 1.0 ? 1.0 : xyz_output_loc[0];
		xyz[i + 1] = xyz_output_loc[1] < 0.0 ? 0.0 : xyz_output_loc[1] > 1.0 ? 1.0 : xyz_output_loc[1];
		xyz[i + 2] = xyz_output_loc[2] < 0.0 ? 0.0 : xyz_output_loc[2] > 1.0 ? 1.0 : xyz_output_loc[2];
	// 	// xyz[i + 3] = xyz_output_loc[3] < 0.0 ? 0.0 : xyz_output_loc[3] > 1.0 ? 1.0 : xyz_output_loc[3];
	// 	// xyz[i + 4] = xyz_output_loc[4] < 0.0 ? 0.0 : xyz_output_loc[4] > 1.0 ? 1.0 : xyz_output_loc[4];
	// 	// xyz[i + 5] = xyz_output_loc[5] < 0.0 ? 0.0 : xyz_output_loc[5] > 1.0 ? 1.0 : xyz_output_loc[5];
	// 	// xyz[i + 6] = xyz_output_loc[6] < 0.0 ? 0.0 : xyz_output_loc[6] > 1.0 ? 1.0 : xyz_output_loc[6];
	// 	// xyz[i + 7] = xyz_output_loc[7] < 0.0 ? 0.0 : xyz_output_loc[7] > 1.0 ? 1.0 : xyz_output_loc[7];
	// 	// xyz[i + 8] = xyz_output_loc[8] < 0.0 ? 0.0 : xyz_output_loc[8] > 1.0 ? 1.0 : xyz_output_loc[8];
	}

	//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, (length / 3), 3, 1.0, conversion_matrix, 3, xyz, (length / 3), 0.0, xyz, (length / 3));
	
}

uint8_t * rgb_space_from_lchuv(double *lchuv, size_t length) {
	uint8_t *rgb_data = malloc(length);

	double *img_data = lchuv;

	lchuv_to_luv(img_data, length);

	luv_to_xyz(img_data, D65, length);
	//luv_to_xyz(img_data, D50, length);

	xyz_to_rgb(img_data, rgb_data, length);

	return rgb_data;
}

static void lchuv_to_luv(double *lchuv, size_t length) {

	#pragma omp parallel for
	for (size_t i = 0u; i < length; i += 3u) {
		//L stays the same
		double u = lchuv[i + 1] * cos(lchuv[i + 2]);
		double v = lchuv[i + 1] * sin(lchuv[i + 2]);

		lchuv[i + 1] = u;
		lchuv[i + 2] = v;
	}
}

static void luv_to_xyz(double *luv, const double ref_white[3], size_t length) {
	const double ref_u = ( 4.0 * ref_white[0] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );
	const double ref_v = ( 9.0 * ref_white[1] ) / ( ref_white[0] + 15.0 * ref_white[1] + 3.0 * ref_white[2] );

	#pragma omp parallel for
	for (size_t i = 0u; i < length; i += 3u) {
		double Y = 0.0;
		if (luv[i + 0] > (CIE_KAPPA * CIE_EPLSION)) {
			Y = pow( (luv[i + 0] + 16.0) / 116.0, 3.0);
		}
		else {
			Y = luv[i + 0] / CIE_KAPPA;
		}

		double a = ( 1.0 / 3.0 ) * ( ( 52.0 * luv[i + 0] ) / ( luv[i + 1] + 13.0 * luv[i + 0] * ref_u ) - 1.0 );
		double b = -5.0 * Y;
		double c = -( 1.0 / 3.0 );
		double d = Y * ( ( 39.0 * luv[i + 0] ) / ( luv[i + 2] + 13.0 * luv[i + 0] * ref_v ) - 5.0 );

		double X = ( d - b ) / ( a - c );
		double Z = X * a + b;

		X = X < 0.0 ? 0.0 : X > 1.0 ? 1.0 : X;
		Y = Y < 0.0 ? 0.0 : Y > 1.0 ? 1.0 : Y;
		Z = Z < 0.0 ? 0.0 : Z > 1.0 ? 1.0 : Z;

		luv[i + 0] = X;
		luv[i + 1] = Y;
		luv[i + 2] = Z;
	}
}

static void xyz_to_rgb(double *xyz, uint8_t *rgb, size_t length) {

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
	double *xyz_value;

	#pragma omp parallel for private(rgb_value, xyz_value) shared(rgb, conversion_matrix)
	for (size_t i = 0u; i < length; i += 3u) {

		xyz_value = &xyz[i];

		//cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, conversion_matrix, 3, xyz_value, 1, 0, rgb_value, 1);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0, conversion_matrix, 3, xyz_value, 1, 0.0, rgb_value, 1);

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

		rgb[i + 0] = (uint8_t) (rgb_value[0] <= 0.0) ? 0 : (rgb_value[0] >= 255.0) ? 255 : rgb_value[0];
		rgb[i + 1] = (uint8_t) (rgb_value[1] <= 0.0) ? 0 : (rgb_value[1] >= 255.0) ? 255 : rgb_value[1];
		rgb[i + 2] = (uint8_t) (rgb_value[2] <= 0.0) ? 0 : (rgb_value[2] >= 255.0) ? 255 : rgb_value[2];
	}
}
