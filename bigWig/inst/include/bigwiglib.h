#ifndef BIGWIGLIB_H
#define BIGWIGLIB_H

typedef struct bbiFile bigWig_t;

int is_bigwig(char * filename);

bigWig_t * bigwig_load(const char * filename, const char * udc_dir);
void bigwig_free(bigWig_t * bw);

long bigwig_chrom_size(bigWig_t * bw, const char * chrom);
int bigwig_valid_chrom(bigWig_t * bw, const char * chrom);

double * bigwig_readf(bigWig_t * bw, const char * chrom, int start, int end, int step, int abs, int * out_length, int * out_is_blank);

int * bigwig_readi(bigWig_t * bw, const char * chrom, int start, int end, int step, int abs, int * out_length, int * out_is_blank);

#endif
