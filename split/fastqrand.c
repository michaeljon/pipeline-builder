#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#define GAPSIZE 128

int main(int argc, char **argv)
{
    char *line = malloc(16384);
    size_t linecap = 16384;
    ssize_t linelen;

    char *fname = argc == 0 ? "zr5654_4_S1_sample_R1.fastq" : argv[1];
    FILE *fp = fopen(fname, "w+");

    int cur = 0;
    while (1)
    {
        int slot = cur++ % GAPSIZE;

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        if (slot == 0)
            fwrite(line, linelen, 1, fp);

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        if (slot == 0)
            fwrite(line, linelen, 1, fp);

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        if (slot == 0)
            fwrite(line, linelen, 1, fp);

        if ((linelen = getline(&line, &linecap, stdin)) <= 0)
            break;
        if (slot == 0)
            fwrite(line, linelen, 1, fp);
    }

    free(line);

    fflush(fp);
    fclose(fp);
}