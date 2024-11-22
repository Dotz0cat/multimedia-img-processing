#ifndef IMAGE_H
#define IMAGE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <jpeglib.h>

struct marker_carrier {
    size_t marker_length;
    uint8_t *marker_data;
    uint8_t marker_number;
};

enum command_map {
    PROCESS_COLOR_SPACE_RGB     = 1u << 0u,
    PROCESS_NOISE_DETECTION     = 1u << 1u,
    PROCESS_NOISE_REDUCTION     = 1u << 2u,
    PROCESS_AVG                 = 1u << 3u,
    PROCESS_MSE                 = 1u << 4u,
    PROCESS_MAX                 = 1u << 5u,
    PROCESS_PSNR                = 1u << 6u,
    PROCESS_SNR                 = 1u << 7u
};

struct jpeg_img_data {
    unsigned int num_components;

    uint64_t height;
    uint64_t width;

    J_COLOR_SPACE color_space;
    
    double gamma;

    unsigned int num_markers;
    struct marker_carrier *markers;

    uint8_t *raw_data;
    double *lchuv_data;

    enum command_map commands;

    double *avg;
    double *mse;
    void *max;
    double *psnr;
    double *snr;
};

struct jpeg_file_list {
    struct jpeg_img_data *data;

    struct jpeg_file_list *next;
};

static const int JPEG_MARKERS[] = {
    JPEG_APP0 + 0,
    JPEG_APP0 + 1,
    JPEG_APP0 + 2,
    JPEG_APP0 + 3,
    JPEG_APP0 + 4,
    JPEG_APP0 + 5,
    JPEG_APP0 + 6,
    JPEG_APP0 + 7,
    JPEG_APP0 + 8,
    JPEG_APP0 + 9,
    JPEG_APP0 + 10,
    JPEG_APP0 + 11,
    JPEG_APP0 + 12,
    JPEG_APP0 + 13,
    JPEG_APP0 + 14,
    JPEG_APP0 + 15,
    JPEG_COM
};

struct jpeg_img_data * read_image(FILE *fp);
void write_image(struct jpeg_img_data *image, const char *file_name);

void free_jpeg_img_data(struct jpeg_img_data *data);
void free_marker_carrier(struct marker_carrier *marker);
void free_jpeg_file_list(struct jpeg_file_list *list);

#endif /* IMAGE_H */
