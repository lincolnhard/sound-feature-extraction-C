#ifndef FFT_H
#define FFT_H

#include <stdio.h> // for printf
#include <stdlib.h> // for posix_memalign
#include <string.h> // for memset
#include <math.h> // for sqrt
#include <limits.h> //for SHRT_MAX
#include <time.h> // for clock_gettime
#include <wchar.h>
#define FEATURE_LEN 34

typedef struct
    {
    double re;
    double im;
    }complex;

typedef struct
    {
    int fs;                             // sampling rate
    int signalsize;                     // input signal length (long-term)
    double **filterbanks;               // the triangular filterbank for MFCC computation
    int *chroma_nummat;                 // the chroma matrices (first half: nChroma, second half: nFreqsPerChroma[nChroma])
    double *features;                   // features (output)
    double *signal_prev_freq_domain;
    double signal_prev_freq_domain_sum;
    }feature_extraction_struct;

void feature_extraction_release
    (
    void
    );

double *feature_extraction_main
    (
    short *signal,
    int frame_idx
    );

void feature_extraction_init
    (
    const int FS,
    const int SIGLEN
    );

#endif // FFT_H
