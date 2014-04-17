#ifndef BW_BASE_H
#define BW_BASE_H

typedef struct bbiFile bigWig_t;

int bw_has_chrom(bigWig_t * bw, const char * chromName);
long bw_chrom_size(bigWig_t * bw, const char * chromName);

#endif
