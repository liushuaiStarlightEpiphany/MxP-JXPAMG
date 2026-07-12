#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int level, fine_size, ahi_level_num_fpts, ahi_level_num_ei_fpts;
    double ahi_level_ratio_node, ahi_level_ratio_F;

    if (argc > 1) {
        level = atoi(argv[1]);
        fine_size = argc > 2 ? atoi(argv[2]) : 1000;
        ahi_level_num_fpts = argc > 3 ? atoi(argv[3]) : 500;
        ahi_level_num_ei_fpts = argc > 4 ? atoi(argv[4]) : 100;
        ahi_level_ratio_node = argc > 5 ? atof(argv[5]) : 25.5;
        ahi_level_ratio_F = argc > 6 ? atof(argv[6]) : 12.3;
    } else {
        level = 0;
        fine_size = 900;
        ahi_level_num_fpts = 450;
        ahi_level_num_ei_fpts = 120;
        ahi_level_ratio_node = 30.5;
        ahi_level_ratio_F = 15.2;
    }

    printf("[AHI-STAT] Level %d: "
           "N = %d, N_F = %d, N_EI = %d, "
           "rEI,node = %.4f %%, rEI,F = %.4f %%\n",
           level, fine_size, ahi_level_num_fpts, ahi_level_num_ei_fpts,
           ahi_level_ratio_node, ahi_level_ratio_F);

    return 0;
}
