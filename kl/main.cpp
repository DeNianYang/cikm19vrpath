#include "kernighanlin.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage: %s graph_file size_limit cost_ratio\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    KernighanLin k;
    k.read_graph(argv[1], static_cast<size_t>(std::atoi(argv[2])), std::atof(argv[3]));

    return 0;
}
