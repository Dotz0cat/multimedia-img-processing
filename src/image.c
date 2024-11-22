#include "image.h"

struct jpeg_img_data * read_image(FILE *fp) {
    struct jpeg_img_data *image = malloc(sizeof(struct jpeg_img_data));

    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr err;

    size_t total_img_bytes;
    uint8_t *buffer;

    uint8_t *row_ptr[1];

    cinfo.err = jpeg_std_error(&err);
    jpeg_create_decompress(&cinfo);

    jpeg_stdio_src(&cinfo, fp);

    for (unsigned int i = 0u; i < sizeof(JPEG_MARKERS) / sizeof(JPEG_MARKERS[0]); i++) {
        jpeg_save_markers(&cinfo, JPEG_MARKERS[i], 0xFFFF);
    }

    jpeg_read_header(&cinfo, 1);

    jpeg_start_decompress(&cinfo);

    image->num_components = cinfo.num_components;
    image->width = cinfo.output_width;
    image->height = cinfo.output_height;
    image->color_space = cinfo.out_color_space;
    image->gamma = cinfo.output_gamma;
    cinfo.CCIR601_sampling = 0;

    total_img_bytes = image->width * image->height * image->num_components;

    buffer = malloc(sizeof(uint8_t) * total_img_bytes);

    image->raw_data = buffer;

    while (cinfo.output_scanline < cinfo.output_height) {
        row_ptr[0] = &buffer[cinfo.num_components * cinfo.output_width * cinfo.output_scanline];
        jpeg_read_scanlines(&cinfo, row_ptr, 1);
    }

    unsigned int j = 0u;
    struct jpeg_marker_struct *marker_list = cinfo.marker_list;
    while (marker_list) {
        j++;
        marker_list = marker_list->next;
    }

    image->num_markers = j;
    image->markers = malloc(sizeof(struct marker_carrier) * j);

    j = 0u;
    marker_list = cinfo.marker_list;
    while (marker_list) {
        size_t size = marker_list->data_length;
        image->markers[j].marker_length = size;
        image->markers[j].marker_number = marker_list->marker;
        uint8_t *marker_data = malloc(size);
        if (!marker_data) {
            fprintf(stderr, "Could not get memory to save marker\n");
            abort();
        }
        memcpy(marker_data, marker_list->data, size);
        image->markers[j++].marker_data = marker_data;
        marker_list = marker_list->next;
    }

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(fp);

    image->commands = 0;
    image->lchuv_data = NULL;
    image->avg = NULL;
    image->mse = NULL;
    image->max = NULL;
    image->psnr = NULL;
    image->snr = NULL;

    return image;
}

void write_image(struct jpeg_img_data *image, const char *file_name) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr err;

    uint8_t *row_ptr[1];

    FILE *fp;

    if ((fp = fopen(file_name, "wb")) == NULL) {
            fprintf(stderr, "Could not open %s for writing\n", file_name);
            return;
    }

    jpeg_create_compress(&cinfo);
    
    cinfo.err = jpeg_std_error(&err);

    cinfo.image_width = image->width;
    cinfo.image_height = image->height;
    cinfo.input_components = image->num_components;
    cinfo.in_color_space = image->color_space;
    cinfo.input_gamma = image->gamma;

    jpeg_stdio_dest(&cinfo, fp);

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 96, 1);
    cinfo.write_JFIF_header = 0;

    jpeg_start_compress(&cinfo, 1);

    for (unsigned int i = 0u; i < image->num_markers; i++) {
        jpeg_write_marker(&cinfo, image->markers[i].marker_number, image->markers[i].marker_data, image->markers[i].marker_length);
    }

    while(cinfo.next_scanline < cinfo.image_height) {
        row_ptr[0] = &(image->raw_data[cinfo.next_scanline * (image->width * image->num_components)]);
        jpeg_write_scanlines(&cinfo, row_ptr, 1);
    }

    jpeg_finish_compress(&cinfo);

    fclose(fp);

    jpeg_destroy_compress(&cinfo);
}

void free_jpeg_img_data(struct jpeg_img_data *data) {
    for (unsigned int i = 0u; i < data->num_markers; i++) {
        free_marker_carrier(&data->markers[i]);
    }
    free(data->markers);
    free(data->raw_data);

    if (data->lchuv_data != NULL) free(data->lchuv_data);
    if (data->avg != NULL) free(data->avg);
    if (data->mse != NULL) free(data->mse);
    if (data->max != NULL) free(data->max);
    if (data->psnr != NULL) free(data->psnr);
    if (data->snr != NULL) free(data->snr);

    free(data);
    data = NULL;
}

void free_marker_carrier(struct marker_carrier *marker) {
    free(marker->marker_data);
    //free(marker);
    //marker = NULL;
}

void free_jpeg_file_list(struct jpeg_file_list *list) {
    struct jpeg_file_list *list1 = list;

    while(list1) {
        free_jpeg_img_data(list1->data);
        struct jpeg_file_list *temp = list1;
        list1 = list1->next;
        free(temp);
        temp = NULL;
    }

    list = NULL;
}
