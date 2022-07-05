#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#define SPLITS 128

int main(int argc, char **argv)
{
    char *line = malloc(16384);
    size_t linecap = 16384;
    ssize_t linelen;
    FILE *splits[SPLITS];

    for (int n = 0; n < SPLITS; n++)
    {
        char fname[256];
        sprintf(fname, "zr5654_10_S6_R1_%03d.fastq", n + 1);

        splits[n] = fopen(fname, "w+");
    }

    int cur = 0;
    while (1)
    {
        int slot = cur++ % SPLITS;

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        fwrite(line, linelen, 1, splits[slot]);

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        fwrite(line, linelen, 1, splits[slot]);

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        fwrite(line, linelen, 1, splits[slot]);

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        fwrite(line, linelen, 1, splits[slot]);
    }

    free(line);

    for (int n = 0; n < SPLITS; n++)
    {
        fflush(splits[n]);
        fclose(splits[n]);
    }
}