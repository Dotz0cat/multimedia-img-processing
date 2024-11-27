#include "processing_private.h"

void process_image(struct jpeg_img_data *restrict img) {
    if (!img->commands) return;

    // Mutual exclusivity: Noise detection and noise reduction cannot both be selected
    if ((img->commands & PROCESS_NOISE_DETECTION) && (img->commands & PROCESS_NOISE_REDUCTION)) {
        return;
    }

    // Dependency checks and processing in order
    if (img->commands & PROCESS_AVG) {
        if (img->commands & PROCESS_COLOR_SPACE_RGB) {
            img->avg = calculate_avg_rgb(img->raw_data, (img->height * img->width * 3));
        } else {
            img->avg = calculate_avg_lchuv(img->lchuv_data, (img->height * img->width * 3));
        }
    }

    if (img->commands & PROCESS_MSE) {
        if (img->avg) {
            if (img->commands & PROCESS_COLOR_SPACE_RGB) {
                img->mse = calculate_mse_rgb(img->raw_data, img->avg, (img->height * img->width * 3));
            } else {
                img->mse = calculate_mse_lchuv(img->lchuv_data, img->avg, (img->height * img->width * 3));
            }
        }
    }

    if (img->commands & PROCESS_MAX) {
        if (img->commands & PROCESS_COLOR_SPACE_RGB) {
            img->max = find_max_rgb(img->raw_data, (img->height * img->width * 3));
        } else {
            img->max = find_max_lchuv(img->lchuv_data, (img->height * img->width * 3));
        }
    }

    if (img->commands & PROCESS_PSNR) {
        if (img->max && img->mse) {
            if (img->commands & PROCESS_COLOR_SPACE_RGB) {
                img->psnr = psnr_rgb(img->max, img->mse);
            } else {
                img->psnr = psnr_lchuv(img->max, img->mse);
            }
        }
    }

    if (img->commands & PROCESS_SNR) {
        if (img->avg && img->mse) {
            img->snr = snr(img->avg, img->mse);
        }
    }

    if (img->commands & PROCESS_NOISE_DETECTION) {
        // Perform noise detection as a standalone operation
        noise_detection(img);
    }

    if (img->commands & PROCESS_NOISE_REDUCTION) {
        // Perform noise reduction as a standalone operation
        if (img->commands & PROCESS_COLOR_SPACE_RGB)
            noise_reduction_rgb(img);
        else
            noise_reduction_lchuv(img);
    }
}

static void apply_edit(double *restrict image_arr, struct edit *restrict edit_to_make) {

    if (image_arr == NULL) return;
    if (edit_to_make == NULL) return;

    image_arr[edit_to_make->loc] = edit_to_make->new_value;

    return;
    
}

static void noise_reduction_lchuv(struct jpeg_img_data *restrict img) {

    double *image = img->lchuv_data;
    double *new_image = aligned_alloc(64, sizeof(double) * ((img->height * img->width) * 3));

    if (new_image == NULL) {
        fprintf(stderr, "I broke\n");
        abort();
    }

    double neighbors[24] = { 0 };
    double average[4] = { 0 };
    double stdev[4] = { 0 };
    uint8_t mask = 0u;

    #pragma omp parallel for shared(image, new_image) firstprivate(neighbors, mask, average, stdev) collapse(2)
    for (size_t i = 0UL; i < img->height - 1; i++) {
        for (size_t j = 0UL; j < img->width - 1; j++) {
            //#pragma omp task
            {
                compute_mask(&mask, i, j, img->height, img->width);
            }

            //#pragma omp task
            {
                select_neighbors_double(mask, neighbors, image, img->width, i, j);
            }

            //#pragma omp task
            {
                add_inital_center_to_stat_lchuv(average, image, img->width, i, j);
                calc_8_neighbor_average_lchuv(mask, neighbors, average);
                add_inital_center_to_stat_lchuv(stdev, image, img->width, i, j);
                calc_stdev_lchuv(mask, neighbors, average, stdev);
            }

            //#pragma omp task
            {
                int loc = ( i * img->width + j ) * 3;
                new_image[loc + 0] = (fabs(image[loc + 0] - average[0]) / stdev[0] > 2.0) ? average[0] : image[loc + 0];
                new_image[loc + 1] = (fabs(image[loc + 1] - average[1]) / stdev[1] > 2.0) ? average[1] : image[loc + 1];
                new_image[loc + 2] = (fabs(image[loc + 2] - average[2]) / stdev[2] > 2.0) ? average[2] : image[loc + 2];
            }
        }
    }

    free(img->lchuv_data);
    img->lchuv_data = new_image;

    return; 
}

static void noise_reduction_rgb(struct jpeg_img_data *restrict img) {

    uint8_t *image = img->raw_data;
    uint8_t *new_image = aligned_alloc(64, (img->height * img->width) * 3);

    double neighbors[24] = { 0 };
    double average[3] = { 0 };
    double stdev[3] = { 0 };
    uint8_t mask = 0;

    #pragma omp parallel for shared(image, new_image) firstprivate(neighbors, mask) collapse(2)
    for (size_t i = 0UL; i < img->height - 1; i++) {
        for (size_t j = 0UL; j < img->width - 1; j++) {
            //#pragma omp task
            {
                compute_mask(&mask, i, j, img->height, img->width);
            }

            //#pragma omp task
            {
                select_neighbors_uint8_t(mask, neighbors, image, img->width, i, j);
            }

            //#pragma omp task
            {
                add_inital_center_to_stat_rgb(average, image, img->width, i, j);
                calc_8_neighbor_average_rgb(mask, neighbors, average);
                add_inital_center_to_stat_rgb(stdev, image, img->width, i, j);
                calc_stdev_rgb(mask, neighbors, average, stdev);
            }

            //#pragma omp task
            {
                                int loc = ( i * img->width + j ) * 3;
                new_image[loc + 0] = (fabs(image[loc + 0] - average[0]) / stdev[0] > 2.0) ? average[0] : image[loc + 0];
                new_image[loc + 1] = (fabs(image[loc + 1] - average[1]) / stdev[1] > 2.0) ? average[1] : image[loc + 1];
                new_image[loc + 2] = (fabs(image[loc + 2] - average[2]) / stdev[2] > 2.0) ? average[2] : image[loc + 2];
            }
        }
    }

    free(img->raw_data);
    img->raw_data = new_image;

    return; 
}

void noise_detection(struct jpeg_img_data *restrict img) {
    if (img->commands & PROCESS_COLOR_SPACE_RGB)
        noise_detection_rgb(img);
    else
        noise_detection_lchuv(img);

    return;
}

static void noise_detection_lchuv(struct jpeg_img_data *restrict img) {

    double *image = img->lchuv_data;
    double *new_image = aligned_alloc(64, sizeof(double) * ((img->height * img->width) * 3));

    if (new_image == NULL) {
        fprintf(stderr, "I broke\n");
        abort();
    }

    double neighbors[24] = { 0 };
    double average[4] = { 0 };
    double stdev[4] = { 0 }; 
    uint8_t mask = 0u;

    #pragma omp parallel for shared(image, new_image) firstprivate(neighbors, mask, average, stdev) collapse(2)
    for (size_t i = 0UL; i < img->height - 1; i++) {
        for (size_t j = 0UL; j < img->width - 1; j++) {
            //#pragma omp task
            {
                compute_mask(&mask, i, j, img->height, img->width);
            }

            //#pragma omp task
            {
                select_neighbors_double(mask, neighbors, image, img->width, i, j);
            }

            //#pragma omp task
            {
                add_inital_center_to_stat_lchuv(average, image, img->width, i, j);
                calc_8_neighbor_average_lchuv(mask, neighbors, average);
                add_inital_center_to_stat_lchuv(stdev, image, img->width, i, j);
                calc_stdev_lchuv(mask, neighbors, average, stdev);
            }

            //#pragma omp task
            {
                int loc = ( i * img->width + j ) * 3;
                new_image[loc + 0] = (fabs(image[loc + 0] - average[0]) / stdev[0] > 2.0) ? 0.0 : 100.0;
                new_image[loc + 1] = (fabs(image[loc + 1] - average[1]) / stdev[1] > 2.0) ? 0.0 : 100.0;
                new_image[loc + 2] = (fabs(image[loc + 2] - average[2]) / stdev[2] > 2.0) ? 0.0 : average[2];
            }
        }
    }

    free(img->lchuv_data);
    img->lchuv_data = new_image;

    return; 
}

static void noise_detection_rgb(struct jpeg_img_data *restrict img) {

    uint8_t *image = img->raw_data;
    uint8_t *new_image = aligned_alloc(64, (img->height * img->width) * 3);

    double neighbors[24] = { 0 };
    double average[3] = { 0 };
    double stdev[3] = { 0 };
    uint8_t mask = 0;

    #pragma omp parallel for shared(image, new_image) firstprivate(neighbors, mask) collapse(2)
    for (size_t i = 0UL; i < img->height - 1; i++) {
        for (size_t j = 0UL; j < img->width - 1; j++) {
            //#pragma omp task
            {
                compute_mask(&mask, i, j, img->height, img->width);
            }

            //#pragma omp task
            {
                select_neighbors_uint8_t(mask, neighbors, image, img->width, i, j);
            }

            //#pragma omp task
            {
                add_inital_center_to_stat_rgb(average, image, img->width, i, j);
                calc_8_neighbor_average_rgb(mask, neighbors, average);
                add_inital_center_to_stat_rgb(stdev, image, img->width, i, j);
                calc_stdev_rgb(mask, neighbors, average, stdev);
            }

            //#pragma omp task
            {
                int loc = ( i * img->width + j ) * 3;
                new_image[loc + 0] = (fabs(image[loc + 0] - average[0]) / stdev[0] > 2.0) ? 0u : 255u;
                new_image[loc + 1] = (fabs(image[loc + 1] - average[1]) / stdev[1] > 2.0) ? 0u : 255u;
                new_image[loc + 2] = (fabs(image[loc + 2] - average[2]) / stdev[2] > 2.0) ? 0u : 255u;
            }
        }
    }

    free(img->raw_data);
    img->raw_data = new_image;

    return; 
}

double * calculate_avg_lchuv(const double *restrict image, size_t length) {
    double *restrict average = malloc(sizeof(double) * 3);

    double angle_avg[2] = {0}; 

    double indvidual_length = (double) (length / 3UL);

    average[0] = 0.0;
    average[1] = 0.0;
    average[2] = 0.0;

    #pragma omp parallel for shared(image) reduction(+:average[0:3])
    for (size_t i = 0UL; i < length; i += 3UL) {
        average[0] += image[i + 0];
        average[1] += image[i + 1];
        //average[2] += image[i + 2];
        angle_avg[0] += sin(image[i + 2]);
        angle_avg[1] += cos(image[i + 2]);
    }

    average[0] = average[0] / indvidual_length;
    average[1] = average[1] / indvidual_length;

    average[2] = atan2(angle_avg[0], angle_avg[1]);

    return average;
}

double * calculate_avg_rgb(const uint8_t *restrict image, size_t length) {
    uint64_t sum[3] = {0, 0, 0};

    double indvidual_length = (double) (length / 3UL);

    double *restrict average = malloc(sizeof(double) * 3);

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

double * calculate_mse_lchuv(const double *restrict image, const double *restrict average, size_t length) {

    double *restrict mse = malloc(sizeof(double) * 3);

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

double * calculate_mse_rgb(const uint8_t *restrict image, const double *restrict average, size_t length) {

    double *restrict mse = malloc(sizeof(double) * 3);

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

double * find_max_lchuv(const double *restrict image, size_t length) {
    double *restrict max = malloc(sizeof(double) * 3);

    #pragma omp parallel for shared(image) reduction(max:max[0:3])
    for (size_t i = 0UL; i < length; i += 3UL) {
        max[0] = image[i + 0] > max[0] ? image[i + 0] : max[0];
        max[1] = image[i + 1] > max[1] ? image[i + 1] : max[1];
        max[2] = image[i + 2] > max[2] ? image[i + 2] : max[2];
    }

    return max;
}

uint8_t * find_max_rgb(const uint8_t *restrict image, size_t length) {
    uint8_t *restrict max = malloc(sizeof(uint8_t) * 3);

    #pragma omp parallel for shared(image) reduction(max:max[0:3])
    for (size_t i = 0UL; i < length; i += 3UL) {
        max[0] = image[i + 0] > max[0] ? image[i + 0] : max[0];
        max[1] = image[i + 1] > max[1] ? image[i + 1] : max[1];
        max[2] = image[i + 2] > max[2] ? image[i + 2] : max[2];
    }

    return max;
}

double * psnr_lchuv(const double *restrict max, const double *restrict mse) {
    double *restrict psnr = malloc(sizeof(double) * 3);

    psnr[0] = 20.0 * log10(max[0]) - 10.0 * log10(mse[0]);
    psnr[1] = 20.0 * log10(max[1]) - 10.0 * log10(mse[1]);
    psnr[2] = 20.0 * log10(max[2]) - 10.0 * log10(mse[2]);

    return psnr;
}

double * psnr_rgb(const uint8_t *restrict max, const double *restrict mse) {
    double *restrict psnr = malloc(sizeof(double) * 3);

    psnr[0] = 20.0 * log10((double) max[0]) - 10.0 * log10(mse[0]);
    psnr[1] = 20.0 * log10((double) max[1]) - 10.0 * log10(mse[1]);
    psnr[2] = 20.0 * log10((double) max[2]) - 10.0 * log10(mse[2]);

    return psnr;
}

double * snr(const double *restrict mean, const double *restrict mse) {
    double *restrict snr = malloc(sizeof(double) * 3);

    snr[0] = 20.0 * log10(mean[0]) - 20.0 * log10(mse[0]);
    snr[1] = 20.0 * log10(mean[1]) - 20.0 * log10(mse[1]);
    snr[2] = 20.0 * log10(mean[2]) - 20.0 * log10(mse[2]);

    return snr;
}

static inline void compute_mask(uint8_t *restrict mask, const size_t i, const size_t j, const size_t height, const size_t width) {
    *mask |= (i > 0 && j > 0) ? 1u << 0u : 0u;
    *mask |= (i > 0) ? 1u << 1u : 0u;
    *mask |= (i > 0 && j < width - 1) ? 1u << 2u : 0u;
    *mask |= (j > 0) ? 1u << 3u : 0u;
    *mask |= (j < width -1) ? 1u << 4u : 0u;
    *mask |= (i < height - 1 && j > 0) ? 1u << 5u : 0u;
    *mask |= (i < height - 1) ? 1u << 6u : 0u;
    *mask |= (i < height - 1 && j < width - 1) ? 1u << 7u : 0u;
}

static inline void select_neighbors_double(const uint8_t mask, double *restrict neighbors, const double *restrict image, const size_t width, const size_t i, const size_t j) {
    if (mask & (1u << 0u)) {
        int loc = ( ( i - 1 ) * width + ( j - 1 ) ) * 3;
        neighbors[0]  = image[loc + 0];
        neighbors[1]  = image[loc + 1];
        neighbors[2]  = image[loc + 2];
    }

    if (mask & (1u << 1u)) {
        int loc = ( ( i - 1 ) * width + ( j + 0 ) ) * 3;
        neighbors[3]  = image[loc + 0];
        neighbors[4]  = image[loc + 1];
        neighbors[5]  = image[loc + 2];
    }

    if (mask & (1u << 2u)) {
        int loc = ( ( i - 1 ) * width + ( j + 1 ) ) * 3;
        neighbors[6]  = image[loc + 0];
        neighbors[7]  = image[loc + 1];
        neighbors[8]  = image[loc + 2];
    }

    if (mask & (1u << 3u)) {
        int loc = ( ( i + 0 ) * width + ( j - 1 ) ) * 3;
        neighbors[9]  = image[loc + 0];
        neighbors[10] = image[loc + 1];
        neighbors[11] = image[loc + 2];
    }

    if (mask & (1u << 4u)) {
        int loc = ( ( i + 0 ) * width + ( j + 1 ) ) * 3;
        neighbors[12] = image[loc + 0];
        neighbors[13] = image[loc + 1];
        neighbors[14] = image[loc + 2];
    }

    if (mask & (1u << 5u)) {
        int loc = ( ( i + 1 ) * width + ( j - 1 ) ) * 3;
        neighbors[15] = image[loc + 0];
        neighbors[16] = image[loc + 1];
        neighbors[17] = image[loc + 2];
    }

    if (mask & (1u << 6u)) {
        int loc = ( ( i + 1 ) * width + ( j + 0 ) ) * 3;
        neighbors[18] = image[loc + 0];
        neighbors[19] = image[loc + 1];
        neighbors[20] = image[loc + 2];
    }

    if (mask & (1u << 7u)) {
        int loc = ( ( i + 1 ) * width + ( j + 1 ) ) * 3;
        neighbors[21] = image[loc + 0];
        neighbors[22] = image[loc + 1];
        neighbors[23] = image[loc + 2];
    }
}

static inline void select_neighbors_uint8_t(const uint8_t mask, double *restrict neighbors, const uint8_t *restrict image, const size_t width, const size_t i, const size_t j) {
    if (mask & (1u << 0u)) {
        int loc = ( ( i - 1 ) * width + ( j - 1 ) ) * 3;
        neighbors[0]  = image[loc + 0];
        neighbors[1]  = image[loc + 1];
        neighbors[2]  = image[loc + 2];
    }

    if (mask & (1u << 1u)) {
        int loc = ( ( i - 1 ) * width + ( j + 0 ) ) * 3;
        neighbors[3]  = image[loc + 0];
        neighbors[4]  = image[loc + 1];
        neighbors[5]  = image[loc + 2];
    }

    if (mask & (1u << 2u)) {
        int loc = ( ( i - 1 ) * width + ( j + 1 ) ) * 3;
        neighbors[6]  = image[loc + 0];
        neighbors[7]  = image[loc + 1];
        neighbors[8]  = image[loc + 2];
    }

    if (mask & (1u << 3u)) {
        int loc = ( ( i + 0 ) * width + ( j - 1 ) ) * 3;
        neighbors[9]  = image[loc + 0];
        neighbors[10] = image[loc + 1];
        neighbors[11] = image[loc + 2];
    }

    if (mask & (1u << 4u)) {
        int loc = ( ( i + 0 ) * width + ( j + 1 ) ) * 3;
        neighbors[12] = image[loc + 0];
        neighbors[13] = image[loc + 1];
        neighbors[14] = image[loc + 2];
    }

    if (mask & (1u << 5u)) {
        int loc = ( ( i + 1 ) * width + ( j - 1 ) ) * 3;
        neighbors[15] = image[loc + 0];
        neighbors[16] = image[loc + 1];
        neighbors[17] = image[loc + 2];
    }

    if (mask & (1u << 6u)) {
        int loc = ( ( i + 1 ) * width + ( j + 0 ) ) * 3;
        neighbors[18] = image[loc + 0];
        neighbors[19] = image[loc + 1];
        neighbors[20] = image[loc + 2];
    }

    if (mask & (1u << 7u)) {
        int loc = ( ( i + 1 ) * width + ( j + 1 ) ) * 3;
        neighbors[21] = image[loc + 0];
        neighbors[22] = image[loc + 1];
        neighbors[23] = image[loc + 2];
    }
}

static inline void add_inital_center_to_stat_lchuv(double *restrict stat, const double *restrict image, const size_t width, const size_t i, const size_t j) {
    int loc = ( i * width + j ) * 3;
    stat[0] += image[loc + 0];
    stat[1] += image[loc + 1];
    stat[2] += sin(image[loc + 2]);
    stat[3] += cos(image[loc + 2]);
}

static inline void calc_8_neighbor_average_lchuv(const uint8_t mask, const double *restrict neighbors, double *restrict average) {
    int count = 1;
    double temp;
    /*
         1   2   4
         8   X  16
        32  64 128
    */
    switch (mask) {
        /*
            1 1 1
            1 1 1
            1 1 1
        */
        case 255u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[6] + neighbors[9]  + neighbors[12] + neighbors[15] + neighbors[18] + neighbors[21];
            average[1] += neighbors[1] + neighbors[4] + neighbors[7] + neighbors[10] + neighbors[13] + neighbors[16] + neighbors[19] + neighbors[22];
            average[2] += sin(neighbors[2]) + sin(neighbors[5]) + sin(neighbors[8]) + sin(neighbors[11]) + sin(neighbors[14]) + sin(neighbors[17]) + sin(neighbors[20]) + sin(neighbors[23]);
            average[3] += cos(neighbors[2]) + cos(neighbors[5]) + cos(neighbors[8]) + cos(neighbors[11]) + cos(neighbors[14]) + cos(neighbors[17]) + cos(neighbors[20]) + cos(neighbors[23]);
            //                     1              2              4               8              16              32              64             128
            count += 8;
            break;
        /*
            0 0 0
            1 1 1
            1 1 1
        */
        case 248u:
            average[0] += neighbors[9]  + neighbors[12] + neighbors[15] + neighbors[18] + neighbors[21];
            average[1] += neighbors[10] + neighbors[13] + neighbors[16] + neighbors[19] + neighbors[22];
            average[2] += sin(neighbors[11]) + sin(neighbors[14]) + sin(neighbors[17]) + sin(neighbors[20]) + sin(neighbors[23]);
            average[3] += cos(neighbors[11]) + cos(neighbors[14]) + cos(neighbors[17]) + cos(neighbors[20]) + cos(neighbors[23]);
            //                      8              16              32              64             128
            count += 5;
            break;
        /*
            0 0 0
            0 1 1
            0 1 1
        */
        case 208u:
            average[0] += neighbors[12] + neighbors[18] + neighbors[21];
            average[1] += neighbors[13] + neighbors[19] + neighbors[22];
            average[2] += sin(neighbors[14]) + sin(neighbors[20]) + sin(neighbors[23]);
            average[3] += cos(neighbors[14]) + cos(neighbors[20]) + cos(neighbors[23]);
            //                     16              64             128
            count += 3;
            break;

        /*
            0 1 1
            0 1 1
            0 1 1
        */
        case 214u:
            average[0] += neighbors[3] + neighbors[6] + neighbors[12] + neighbors[18] + neighbors[21];
            average[1] += neighbors[4] + neighbors[7] + neighbors[13] + neighbors[19] + neighbors[22];
            average[2] += sin(neighbors[5]) + sin(neighbors[8]) + sin(neighbors[14]) + sin(neighbors[20]) + sin(neighbors[23]);
            average[3] += cos(neighbors[5]) + cos(neighbors[8]) + cos(neighbors[14]) + cos(neighbors[20]) + cos(neighbors[23]);
            //                     2              4              16              64             128
            count += 5;
            break;

        /*
            0 1 1
            0 1 1
            0 0 0
        */
        case 22u:
            average[0] += neighbors[3] + neighbors[6] + neighbors[12];
            average[1] += neighbors[4] + neighbors[7] + neighbors[13];
            average[2] += sin(neighbors[5]) + sin(neighbors[8]) + sin(neighbors[14]);
            average[3] += cos(neighbors[5]) + cos(neighbors[8]) + cos(neighbors[14]);
            //                     2              4              16 
            count += 3;
            break;

        /*
            1 1 1
            1 1 1
            0 0 0
        */
        case 31u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[6] + neighbors[9]  + neighbors[12];
            average[1] += neighbors[1] + neighbors[4] + neighbors[7] + neighbors[10] + neighbors[13];
            average[2] += sin(neighbors[2]) + sin(neighbors[5]) + sin(neighbors[8]) + sin(neighbors[11]) + sin(neighbors[14]);
            average[3] += cos(neighbors[2]) + cos(neighbors[5]) + cos(neighbors[8]) + cos(neighbors[11]) + cos(neighbors[14]);
            //                     1              2              4               8              16 
            count += 5;
            break;

        /*
            1 1 0
            1 1 0
            0 0 0
        */
        case 11u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[9] ;
            average[1] += neighbors[1] + neighbors[4] + neighbors[10];
            average[2] += sin(neighbors[2]) + sin(neighbors[5]) + sin(neighbors[11]);
            average[3] += cos(neighbors[2]) + cos(neighbors[5]) + cos(neighbors[11]);
            //                     1              2               8 
            count += 3;
            break;

        /*
            1 1 0
            1 1 0
            1 1 0
        */
        case 107u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[9]  + neighbors[15] + neighbors[18];
            average[1] += neighbors[1] + neighbors[4] + neighbors[10] + neighbors[16] + neighbors[19];
            average[2] += sin(neighbors[2]) + sin(neighbors[5]) + sin(neighbors[11]) + sin(neighbors[17]) + sin(neighbors[20]);
            average[3] += cos(neighbors[2]) + cos(neighbors[5]) + cos(neighbors[11]) + cos(neighbors[17]) + cos(neighbors[20]);
            //                     1              2               8              32              64 
            count += 5;
            break;

        /*
            0 0 0
            1 1 0
            1 1 0
        */
        case 104u:
            average[0] += neighbors[9]  + neighbors[15] + neighbors[18];
            average[1] += neighbors[10] + neighbors[16] + neighbors[19];
            average[2] += sin(neighbors[11]) + sin(neighbors[17]) + sin(neighbors[20]);
            average[3] += cos(neighbors[11]) + cos(neighbors[17]) + cos(neighbors[20]);
            //                      8              32              64 
            count += 3;
            break;
    }

    average[0] /= count;
    average[1] /= count;
    temp = atan2(average[2], average[3]);
    average[3] = sqrt(1 - pow(average[2], 2) + pow(average[3], 2));
    average[2] = temp;
}

static inline void calc_stdev_lchuv(const uint8_t mask, const double *restrict neighbors, const double *restrict average, double *restrict stdev) {
    int count = 1;
    switch (mask) {
    /*
        1 1 1
        1 1 1
        1 1 1
    */
    case 255u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);
        count += 8;
        break;
    /*
        0 0 0
        1 1 1
        1 1 1
    */
    case 248u:
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);
        count += 5;
        break;
    /*
        0 0 0
        0 1 1
        0 1 1
    */
    case 208u:
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);
        count += 3;
        break;

    /*
        0 1 1
        0 1 1
        0 1 1
    */
    case 214u:
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);
        count += 5;
        break;

    /*
        0 1 1
        0 1 1
        0 0 0
    */
    case 22u:
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        count += 3;
        break;

    /*
        1 1 1
        1 1 1
        0 0 0
    */
    case 31u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        count += 5;
        break;

    /*
        1 1 0
        1 1 0
        0 0 0
    */
    case 11u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);
        count += 3;
        break;

    /*
        1 1 0
        1 1 0
        1 1 0
    */
    case 107u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        count += 5;
        break;

    /*
        0 0 0
        1 1 0
        1 1 0
    */
    case 104u:
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        count += 3;
        break;
    }

    stdev[0] /= count;
    stdev[1] /= count;

    stdev[0] = sqrt(stdev[0]);
    stdev[1] = sqrt(stdev[1]);
    stdev[2] = asin(average[3]) * ( (1 + (2.0 / sqrt(3.0) - 1) * pow(average[3], 3) ) );
}

static inline void add_inital_center_to_stat_rgb(double *restrict stat, const uint8_t *restrict image, const size_t width, const size_t i, const size_t j) {
    int loc = ( i * width + j ) * 3;
    stat[0] += image[loc + 0];
    stat[1] += image[loc + 1];
    stat[2] += image[loc + 2];
}

static inline void calc_8_neighbor_average_rgb(const uint8_t mask, const double *restrict neighbors, double *restrict average) {
    int count = 1;
    /*
         1   2   4
         8   X  16
        32  64 128
    */
    switch (mask) {
        /*
            1 1 1
            1 1 1
            1 1 1
        */
        case 255u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[6] + neighbors[9]  + neighbors[12] + neighbors[15] + neighbors[18] + neighbors[21];
            average[1] += neighbors[1] + neighbors[4] + neighbors[7] + neighbors[10] + neighbors[13] + neighbors[16] + neighbors[19] + neighbors[22];
            average[2] += neighbors[2] + neighbors[5] + neighbors[8] + neighbors[11] + neighbors[14] + neighbors[17] + neighbors[20] + neighbors[23];
            //                     1              2              4               8              16              32              64             128
            count += 8;
            break;
        /*
            0 0 0
            1 1 1
            1 1 1
        */
        case 248u:
            average[0] += neighbors[9]  + neighbors[12] + neighbors[15] + neighbors[18] + neighbors[21];
            average[1] += neighbors[10] + neighbors[13] + neighbors[16] + neighbors[19] + neighbors[22];
            average[2] += neighbors[11] + neighbors[14] + neighbors[17] + neighbors[20] + neighbors[23];
            //                      8              16              32              64             128
            count += 5;
            break;
        /*
            0 0 0
            0 1 1
            0 1 1
        */
        case 208u:
            average[0] += neighbors[12] + neighbors[18] + neighbors[21];
            average[1] += neighbors[13] + neighbors[19] + neighbors[22];
            average[2] += neighbors[14] + neighbors[20] + neighbors[23];
            //                     16              64             128
            count += 4;
            break;

        /*
            0 1 1
            0 1 1
            0 1 1
        */
        case 214u:
            average[0] += neighbors[3] + neighbors[6] + neighbors[12] + neighbors[18] + neighbors[21];
            average[1] += neighbors[4] + neighbors[7] + neighbors[13] + neighbors[19] + neighbors[22];
            average[2] += neighbors[5] + neighbors[8] + neighbors[14] + neighbors[20] + neighbors[23];
            //                     2              4              16              64             128
            count += 5;
            break;

        /*
            0 1 1
            0 1 1
            0 0 0
        */
        case 22u:
            average[0] += neighbors[3] + neighbors[6] + neighbors[12];
            average[1] += neighbors[4] + neighbors[7] + neighbors[13];
            average[2] += neighbors[5] + neighbors[8] + neighbors[14];
            //                     2              4              16 
            count += 3;
            break;

        /*
            1 1 1
            1 1 1
            0 0 0
        */
        case 31u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[6] + neighbors[9]  + neighbors[12];
            average[1] += neighbors[1] + neighbors[4] + neighbors[7] + neighbors[10] + neighbors[13];
            average[2] += neighbors[2] + neighbors[5] + neighbors[8] + neighbors[11] + neighbors[14];
            //                     1              2              4               8              16 
            count += 5;
            break;

        /*
            1 1 0
            1 1 0
            0 0 0
        */
        case 11u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[9] ;
            average[1] += neighbors[1] + neighbors[4] + neighbors[10];
            average[2] += neighbors[2] + neighbors[5] + neighbors[11];
            //                     1              2               8 
            count += 3;
            break;

        /*
            1 1 0
            1 1 0
            1 1 0
        */
        case 107u:
            average[0] += neighbors[0] + neighbors[3] + neighbors[9]  + neighbors[15] + neighbors[18];
            average[1] += neighbors[1] + neighbors[4] + neighbors[10] + neighbors[16] + neighbors[19];
            average[2] += neighbors[2] + neighbors[5] + neighbors[11] + neighbors[17] + neighbors[20];
            //                     1              2               8              32              64 
            count += 5;
            break;

        /*
            0 0 0
            1 1 0
            1 1 0
        */
        case 104u:
            average[0] += neighbors[9]  + neighbors[15] + neighbors[18];
            average[1] += neighbors[10] + neighbors[16] + neighbors[19];
            average[2] += neighbors[11] + neighbors[17] + neighbors[20];
            //                      8              32              64 
            count += 3;
            break;
    }

    average[0] /= count;
    average[1] /= count;
    average[2] /= count;
}

static inline void calc_stdev_rgb(const uint8_t mask, const double *restrict neighbors, const double *restrict average, double *restrict stdev) {
    int count = 1;
    switch (mask) {
    /*
        1 1 1
        1 1 1
        1 1 1
    */
    case 255u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);

        stdev[2] += pow(neighbors[2]  - average[2], 2);
        stdev[2] += pow(neighbors[5]  - average[2], 2);
        stdev[2] += pow(neighbors[8]  - average[2], 2);
        stdev[2] += pow(neighbors[11] - average[2], 2);
        stdev[2] += pow(neighbors[14] - average[2], 2);
        stdev[2] += pow(neighbors[17] - average[2], 2);
        stdev[2] += pow(neighbors[20] - average[2], 2);
        stdev[2] += pow(neighbors[23] - average[2], 2);

        count += 8;
        break;
    /*
        0 0 0
        1 1 1
        1 1 1
    */
    case 248u:
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);

        stdev[2] += pow(neighbors[11] - average[2], 2);
        stdev[2] += pow(neighbors[14] - average[2], 2);
        stdev[2] += pow(neighbors[17] - average[2], 2);
        stdev[2] += pow(neighbors[20] - average[2], 2);
        stdev[2] += pow(neighbors[23] - average[2], 2);

        count += 5;
        break;
    /*
        0 0 0
        0 1 1
        0 1 1
    */
    case 208u:
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);

        stdev[2] += pow(neighbors[14] - average[2], 2);
        stdev[2] += pow(neighbors[20] - average[2], 2);
        stdev[2] += pow(neighbors[23] - average[2], 2);

        count += 3;
        break;

    /*
        0 1 1
        0 1 1
        0 1 1
    */
    case 214u:
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        stdev[0] += pow(neighbors[21] - average[0], 2);
        
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);
        stdev[1] += pow(neighbors[22] - average[1], 2);

        stdev[2] += pow(neighbors[5]  - average[2], 2);
        stdev[2] += pow(neighbors[8]  - average[2], 2);
        stdev[2] += pow(neighbors[14] - average[2], 2);
        stdev[2] += pow(neighbors[20] - average[2], 2);
        stdev[2] += pow(neighbors[23] - average[2], 2);

        count += 5;
        break;

    /*
        0 1 1
        0 1 1
        0 0 0
    */
    case 22u:
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);

        stdev[2] += pow(neighbors[5]  - average[2], 2);
        stdev[2] += pow(neighbors[8]  - average[2], 2);
        stdev[2] += pow(neighbors[14] - average[2], 2);

        count += 3;
        break;

    /*
        1 1 1
        1 1 1
        0 0 0
    */
    case 31u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[6]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[12] - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[7]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[13] - average[1], 2);

        stdev[2] += pow(neighbors[2]  - average[2], 2);
        stdev[2] += pow(neighbors[5]  - average[2], 2);
        stdev[2] += pow(neighbors[8]  - average[2], 2);
        stdev[2] += pow(neighbors[11] - average[2], 2);
        stdev[2] += pow(neighbors[14] - average[2], 2);

        count += 5;
        break;

    /*
        1 1 0
        1 1 0
        0 0 0
    */
    case 11u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);

        stdev[2] += pow(neighbors[2]  - average[2], 2);
        stdev[2] += pow(neighbors[5]  - average[2], 2);
        stdev[2] += pow(neighbors[11] - average[2], 2);
        count += 3;
        break;

    /*
        1 1 0
        1 1 0
        1 1 0
    */
    case 107u:
        stdev[0] += pow(neighbors[0]  - average[0], 2);
        stdev[0] += pow(neighbors[3]  - average[0], 2);
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        
        stdev[1] += pow(neighbors[1]  - average[1], 2);
        stdev[1] += pow(neighbors[4]  - average[1], 2);
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);

        stdev[2] += pow(neighbors[2]  - average[2], 2);
        stdev[2] += pow(neighbors[5]  - average[2], 2);
        stdev[2] += pow(neighbors[11] - average[2], 2);
        stdev[2] += pow(neighbors[17] - average[2], 2);
        stdev[2] += pow(neighbors[20] - average[2], 2);
        count += 5;
        break;

    /*
        0 0 0
        1 1 0
        1 1 0
    */
    case 104u:
        stdev[0] += pow(neighbors[9]  - average[0], 2);
        stdev[0] += pow(neighbors[15] - average[0], 2);
        stdev[0] += pow(neighbors[18] - average[0], 2);
        
        stdev[1] += pow(neighbors[10] - average[1], 2);
        stdev[1] += pow(neighbors[16] - average[1], 2);
        stdev[1] += pow(neighbors[19] - average[1], 2);

        stdev[2] += pow(neighbors[11] - average[2], 2);
        stdev[2] += pow(neighbors[17] - average[2], 2);
        stdev[2] += pow(neighbors[20] - average[2], 2);
        count += 3;
        break;
    }

    stdev[0] /= count;
    stdev[1] /= count;
    stdev[2] /= count;

    stdev[0] = sqrt(stdev[0]);
    stdev[1] = sqrt(stdev[1]);
    stdev[2] = sqrt(stdev[2]);
}
