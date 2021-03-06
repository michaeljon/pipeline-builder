#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define EXPERIMENTAL_VCF

typedef unsigned char byte;

// human
#define REF_HG38_P13 1UL
#define REF_HG38_P14 2UL

// coronovirus
#define HCOV_229E 10UL
#define HCOV_HKU1 11UL
#define HCOV_NL63 12UL
#define HCOV_OC63 13UL
#define SARS_COV_2 15UL

typedef struct _chromSize
{
    char *chrom;
    long length;
} chromSize;

const chromSize sizes[] = {
    {.chrom = "1", .length = 248956422},
    {.chrom = "2", .length = 242193529},
    {.chrom = "3", .length = 198295559},
    {.chrom = "4", .length = 190214555},
    {.chrom = "5", .length = 181538259},
    {.chrom = "6", .length = 170805979},
    {.chrom = "7", .length = 159345973},
    {.chrom = "8", .length = 145138636},
    {.chrom = "9", .length = 138394717},
    {.chrom = "10", .length = 133797422},
    {.chrom = "11", .length = 135086622},
    {.chrom = "12", .length = 133275309},
    {.chrom = "13", .length = 114364328},
    {.chrom = "14", .length = 107043718},
    {.chrom = "15", .length = 101991189},
    {.chrom = "16", .length = 90338345},
    {.chrom = "17", .length = 83257441},
    {.chrom = "18", .length = 80373285},
    {.chrom = "19", .length = 58617616},
    {.chrom = "20", .length = 64444167},
    {.chrom = "21", .length = 46709983},
    {.chrom = "22", .length = 50818468},
    {.chrom = "X", .length = 156040895},
    {.chrom = "Y", .length = 57227415}};

const char neucleotides[] = {'A', 'C', 'T', 'G'};
const int ncount = sizeof(neucleotides) / sizeof(neucleotides[0]);

const double variantChance = 0.999;
const double indelChance = 0.18;

static inline int randomRange(int lower, int upper)
{
    return (rand() % (upper - lower + 1)) + lower;
}

static inline double randomPercent()
{
    return ((rand() % (100 - 0 + 1)) + 0) / 100.0;
}

static inline char randomRead()
{
    return neucleotides[randomRange(0, ncount - 1)];
}

static inline char *generateIndel()
{
    int length = randomRange(2, 16);
    static char indel[16];

    for (int n = 0; n < length; n++)
    {
        indel[n] = randomRead();
    }

    indel[length - 1] = 0;
    return indel;
}

static inline unsigned long chromToKey(char *chr)
{
    if (*chr == 'X')
        return 22UL;
    if (*chr == 'Y')
        return 23UL;
    return (unsigned long)atol(chr);
}

static inline void writeSNP(unsigned long referenceId, char *chr, long locus, char r, char a)
{
#ifdef EXPERIMENTAL_VCF
    const unsigned long key =
        referenceId << 56 |
        chromToKey(chr) << 32 |
        locus;

    printf("%016lx\tchr%s\t%ld\t.\t%c\t%c\t0.0\tDP=0\n", key, chr, locus, r, a);
#else
    printf("chr%s\t%ld\t.\t%c\t%c\t0.0\tDP=0\n", chr, locus, r, a);
#endif
}

static inline void writeInsert(unsigned long referenceId, char *chr, long locus, char r, char *a)
{
#ifdef EXPERIMENTAL_VCF
    const unsigned long key =
        referenceId << 56 |
        chromToKey(chr) << 32 |
        locus;

    printf("%016lx\tchr%s\t%ld\t.\t%c\t%s\t0.0\tINDEL;DP=0\n", key, chr, locus, r, a);
#else
    printf("chr%s\t%ld\t.\t%c\t%s\t0.0\tINDEL;DP=0\n", chr, locus, r, a);
#endif
}

static inline void writeDelete(unsigned long referenceId, char *chr, long locus, char *r, char a)
{
#ifdef EXPERIMENTAL_VCF
    const unsigned long key =
        referenceId << 56 |
        chromToKey(chr) << 32 |
        locus;

    printf("%016lx\tchr%s\t%ld\t.\t%s\t%c\t0.0\tINDEL;DP=0\n", key, chr, locus, r, a);
#else
    printf("chr%s\t%ld\t.\t%s\t%c\t0.0\tINDEL;DP=0\n", chr, locus, r, a);
#endif
}

int variantCount = 0;
int indelCount = 0;

void generateVariants(char *chr, long size)
{
    for (long locus = 1; locus <= size; locus++)
    {
        if (randomPercent() > variantChance)
        {
            variantCount++;

            if (randomPercent() > indelChance)
            {
                // is snp
                writeSNP(REF_HG38_P13, chr, locus, randomRead(), randomRead());
            }
            else
            {
                // is indel
                char *indel = generateIndel();

                indelCount++;
                if (randomPercent() > 0.5)
                {
                    // insertion
                    writeInsert(REF_HG38_P13, chr, locus, indel[0], indel);
                }
                else
                {
                    // deletion
                    writeDelete(REF_HG38_P13, chr, locus, indel, indel[0]);
                }

                locus += strlen(indel);
            }
        }
    }
}

int main(int argc, char **argv)
{
    char buffer[BUFSIZ];
    setbuf(stdout, buffer);

    srand(time(0));

    for (int n = 0; n < sizeof(sizes) / sizeof(sizes[0]); n++)
    {
        generateVariants(sizes[n].chrom, sizes[n].length);
    }

    printf("Number of variants %d\n", variantCount);
    printf("Number of indels %d\n", indelCount);
}