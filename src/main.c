#include "main.h"

int main(int argc, char **argv) {
    
    if (argc <= 1) {
        fprintf(stderr, "Usage: %s file ...\n", argv[0]);
        return 1;
    }
    
    struct jpeg_file_list *file_list_head = NULL;
    struct jpeg_file_list **file_list = &file_list_head;


    for (int i = 1; i < argc; i++) {
        (*file_list) = malloc(sizeof(struct jpeg_file_list));
        (*file_list)->data = NULL;
        (*file_list)->next = NULL;
        
        FILE *fp;

        if ((fp = fopen(argv[i], "rb")) == NULL) {
            fprintf(stderr, "Could not open %s\n", argv[i]);
            return 1;
        }
        
        (*file_list)->data = read_image(fp);
        file_list = &(*file_list)->next;
    }

    if (!file_list_head) {
        fprintf(stderr, "An issue with memory");
        return 1;
    }
    
    //do stuff here
    file_list = &file_list_head;

    for (int i = 1; i < argc; i++) {
        (*file_list)->data->lchuv_data = convert_to_lchuv_space((*file_list)->data);
        (*file_list)->data->commands = PROCESS_NOISE_REDUCTION | PROCESS_AVG | PROCESS_MSE;
        process_image((*file_list)->data);
        uint8_t *rgb = rgb_space_from_lchuv((*file_list)->data->lchuv_data, ((*file_list)->data->width * (*file_list)->data->height) * 3);

        free((*file_list)->data->raw_data);
        (*file_list)->data->raw_data = rgb;

        file_list = &(*file_list)->next;
    }
   
    file_list = &file_list_head;

    for (int i = 1; i < argc; i++) {
        size_t len = snprintf(NULL, 0, "%s_output.jpg", argv[i]);
        char *file_name = malloc(sizeof(char) * len + 1u);
        if (!file_name) {
            fprintf(stderr, "I couldn't get the memory I needed\n Good luck!\n");
            return 1;
        }
        snprintf(file_name, len + 1u, "%s_output.jpg", argv[i]);
        
        write_image((*file_list)->data, file_name);

        free(file_name);

        file_list = &(*file_list)->next;
    }

    free_jpeg_file_list(file_list_head);

    return 0;
}
