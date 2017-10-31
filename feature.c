#include "feature.h"

#define VERBOSE
static feature_extraction_struct featext_wksp;

double stddev
    (
    const double *src,
    const int LEN
    )
{
    double mean = 0.0;
    double devsum = 0.0;
    int i = 0;
    for(i = 0; i < LEN; ++i)
        {
        mean += src[i];
        }
    mean /= LEN;
    for(i = 0; i < LEN; ++i)
        {
        double dev = src[i] - mean;
        devsum += (dev * dev);
        }
    return sqrt(devsum / (LEN - 0));
}

static void zcr
    (
    const double *src,
    double *dst
    )
{
    const int DIFFLEN = featext_wksp.signalsize - 1;
    int i = 0;
    int sum = 0;
    for (i = 0; i < DIFFLEN; ++i)
        {
        if(signbit(src[i+1]) - signbit(src[i]))
            {
            ++sum;
            }
        }
    *dst = (double)sum / DIFFLEN;
}

static void energy
    (
    const double *src,
    double *dst
    )
{
    const int LEN = featext_wksp.signalsize;
    double sum = 0;
    int i = 0;
    for (i = 0; i < LEN; ++i)
        {
        sum += (src[i] * src[i]);
        }
    *dst = sum / LEN;
}

static void energy_entropy
    (
    const double *src,
    const double ENERGY,
    double *dst
    )
{
    const int LEN = featext_wksp.signalsize;
    const int NUM_BLOCKS = 10;
    const double EPS = 0.00000001;
    const double ENERGYSUM_RECIPROCAL = 1.0 / (ENERGY * LEN + EPS);
    const int BLOCK_LEN = LEN / NUM_BLOCKS;
    double H = 0.0;
    int i = 0;
    int j = 0;
    for (j = 0; j < NUM_BLOCKS; ++j)
        {
        double block = 0.0;
        const double *srcptr = src + j*BLOCK_LEN;
        for (i = 0; i < BLOCK_LEN; ++i)
            {
            block += srcptr[i] * srcptr[i];
            }
        block *= ENERGYSUM_RECIPROCAL;
        H += (-block * log2(block + EPS));
        }
    *dst = H;
}

static void dft_normalize
    (
    const double *src,
    double *dst,
    double *dstsum,
    double *dstsqsum
    )
{
    const int LEN = featext_wksp.signalsize;
    const int FLEN = LEN >> 1;
    const double FLEN_RECIPROCAL= 1.0 / FLEN;
    double sum = 0.0;
    double sqsum = 0.0;
    int i = 0;
    int j = 0;
    for(i = 0; i < FLEN; ++i)
        {
        double re = 0.0;
        double im = 0.0;
        for(j = 0; j < LEN; ++j)
            {
            //TODO: sin, cos can build table or init first
            re += (src[j] * cos(-2 * M_PI * j * i / LEN));
            im += (src[j] * sin(-2 * M_PI * j * i / LEN));
            }
        dst[i] = sqrt(re*re + im*im) * FLEN_RECIPROCAL;
        sum += dst[i];
        sqsum += (dst[i] * dst[i]);
        }
    *dstsum = sum;
    *dstsqsum = sqsum;
}

static void spectral_centroid_spread
    (
    const double *fsrc,
    const double FSRCSUM,
    double *dst
    )
{
    const int LEN = featext_wksp.signalsize;
    const int FS = featext_wksp.fs;
    const int SEGMENT = FS / LEN;
    const int FLEN = LEN >> 1;
    const double EPS = 0.00000001;
    double num = 0.0;
    double den = FSRCSUM + EPS;
    double centroid = 0.0;
    double spread = 0.0;
    int i = 0;
    for(i = 0; i < FLEN; ++i)
        {
        num += fsrc[i] * (i + 1) * SEGMENT;
        }
    centroid = num / den;
    for(i = 0; i < FLEN; ++i)
        {
        double temp = (i + 1) * SEGMENT - centroid;
        spread += temp * temp * fsrc[i];
        }
    spread = sqrt(spread / den);
    dst[0] = centroid / (FS >> 1);
    dst[1] = spread / (FS >> 1);
}

static void spectral_entropy
    (
    const double *fsrc,
    const double FSRC_SQSUM,
    double *dst
    )
{
    const int FLEN = featext_wksp.signalsize >> 1;
    const int NUM_BLOCKS = 10;
    const double EPS = 0.00000001;
    const int BLOCK_LEN = FLEN / NUM_BLOCKS;
    const double ENERGY_RECIPROCAL = 1.0 / (FSRC_SQSUM + EPS);
    double H = 0.0;
    int i = 0;
    int j = 0;
    for (j = 0; j < NUM_BLOCKS; ++j)
        {
        double block = 0.0;
        const double *fsrcptr = fsrc + j * BLOCK_LEN;
        for (i = 0; i < BLOCK_LEN; ++i)
            {
            block += fsrcptr[i] * fsrcptr[i];
            }
        block *= ENERGY_RECIPROCAL;
        H += (-block * log2(block + EPS));
        }
    *dst = H;
}

static void spectral_flux
    (
    const double *fsrc,
    const double FSRCSUM,
    const double *fsrc_prev,
    const double FSRCSUM_PREV,
    double *dst
    )
{
    const int FLEN = featext_wksp.signalsize >> 1;
    const double S_RECIPROCAL = 1.0 / FSRCSUM;
    const double SPREV_RECIPROCAL = 1.0 / FSRCSUM_PREV;
    double sum = 0.0;
    int i = 0;
    for(i = 0; i < FLEN; ++i)
        {
        double temp = fsrc[i] * S_RECIPROCAL - fsrc_prev[i] * SPREV_RECIPROCAL;
        sum += temp * temp;
        }
    *dst = sum;
}

static void spectral_rolloff
    (
    const double *fsrc,
    const double FSRC_SQSUM,
    double *dst
    )
{
    const int FLEN = featext_wksp.signalsize >> 1;
    const double PERCENT = 0.9;
    const double TH = PERCENT * FSRC_SQSUM;
    double cumsum = 0.0;
    int rollidx = 0;
    int i = 0;
    for(i = 0; i < FLEN; ++i)
        {
        cumsum += (fsrc[i] * fsrc[i]);
        if(cumsum > TH)
            {
            rollidx = i;
            break;
            }
        }
    *dst = (double)rollidx / FLEN;
}

static void mfcc
    (
    const double *fsrc,
    double *dst
    )
{
    const int NUM_LIN_FILTER = 13; //number of Mel-Frequency Cepstral Coefficients
    const int NUM_LOG_FILTER = 27;
    const int NUM_TOTAL_FILTER = NUM_LIN_FILTER + NUM_LOG_FILTER;
    const int FLEN = featext_wksp.signalsize >> 1;
    double **fbank = featext_wksp.filterbanks;
    double *mspec = (double *)malloc(NUM_TOTAL_FILTER * sizeof(double));
    const double PI_F = M_PI / NUM_TOTAL_FILTER;
    const double SF1 = sqrt(0.5); // TODO: pre-compute first
    const double SF2 = sqrt(2.0 / NUM_TOTAL_FILTER); // TODO: pre-compute first
    int i = 0;
    int j = 0;
    for(i = 0; i < NUM_TOTAL_FILTER; ++i)
        {
        double dot = 0.0;
        for(j = 0; j < FLEN; ++j)
            {
            dot += fsrc[j] * fbank[i][j];
            }
        mspec[i] = log10(dot + 0.00000001);
        }
    // compute DCT
    for(i = 0; i < NUM_LIN_FILTER; ++i)
        {
        double sum = 0;
        for(j = 0; j < NUM_TOTAL_FILTER; ++j)
            {
            sum += mspec[j] * cos(PI_F * ( j + 0.5 ) * i );
            }
        if(i == 0)
            {
            sum *= SF1;
            }
        sum *= SF2;
        dst[i] = sum;
        }
    free(mspec);
}

static void chroma
    (
    const double *fsrc,
    const double FSRC_SQSUM,
    double *dst
    )
{
    const int NUM_FEATURE = 12;
    const int FLEN = featext_wksp.signalsize >> 1;
    int *chromaidx = featext_wksp.chroma_nummat;
    int *reciprocalidx = featext_wksp.chroma_nummat + FLEN;
    const int FITDIM = (int)ceil((double)FLEN / NUM_FEATURE);
    double *C = (double*)calloc(FITDIM * NUM_FEATURE, sizeof(double));
    const double FSRC_SQSUM_RECIPROCAL = 1.0 / FSRC_SQSUM;
    int i = 0;
    int j = 0;
    for(i = 0; i < FLEN; ++i)
        {
        C[chromaidx[i]] = fsrc[i] * fsrc[i];
        }
    for(i = 0; i < FLEN; ++i)
        {
        C[i] /= reciprocalidx[i];
        }
    for(i = 0; i < NUM_FEATURE; ++i)
        {
        double temp = 0.0;
        for(j = 0; j < FITDIM; ++j)
            {
            temp += C[NUM_FEATURE * j + i];
            }
        dst[i] = temp * FSRC_SQSUM_RECIPROCAL;
        }
    dst[NUM_FEATURE] = stddev(dst, NUM_FEATURE);
    free(C);
}

static double ** mfcc_filter_banks_init // TODO: should build fbank table manually for faster index
    (
    void
    )
{
    const int NUM_LIN_FILTER = 13; //number of Mel-Frequency Cepstral Coefficients
    const int NUM_LOG_FILTER = 27;
    const int NUM_TOTAL_FILTER = NUM_LOG_FILTER + NUM_LIN_FILTER;
    const int FREQSIZE = featext_wksp.signalsize >> 1;
    const int FS = featext_wksp.fs;
    const double LOWFREQ = 133.33;
    const double LINSC = 200 / 3.0;
    const double LOGSC = 1.0711703;
    double *freqs = (double *)malloc((NUM_TOTAL_FILTER + 2) * sizeof(double));
    double **fbank = (double **)malloc(NUM_TOTAL_FILTER * sizeof(double*));
    int i = 0;
    int j = 0;
    for(i = 0; i < NUM_LIN_FILTER; ++i)
        {
        freqs[i] = LOWFREQ + i * LINSC;
        }
    for(i = 0; i < NUM_LOG_FILTER + 2; ++i)
        {
        freqs[i + NUM_LIN_FILTER] = freqs[NUM_LIN_FILTER - 1] * pow(LOGSC, i + 1);
        }
    for(i = 0; i < NUM_TOTAL_FILTER; ++i)
        {
        fbank[i] = (double *)calloc(FREQSIZE, sizeof(double));
        }
    for(i = 0; i < NUM_TOTAL_FILTER; ++i)
        {
        double height = 2.0 / (freqs[2 + i] - freqs[i]);
        double low_triangle_freq = freqs[i];
        double mid_triangle_freq = freqs[i + 1];
        double high_triangle_freq = freqs[i + 2];
        int lidstart = (int)floor(low_triangle_freq * FREQSIZE / FS) + 1;
        int lidend= (int)floor(mid_triangle_freq * FREQSIZE / FS) + 1;
        double lslope = height / (mid_triangle_freq - low_triangle_freq);
        int ridstart = (int)floor(mid_triangle_freq * FREQSIZE / FS) + 1;
        int ridend= (int)floor(high_triangle_freq * FREQSIZE / FS) + 1;
        double rslope = height / (high_triangle_freq - mid_triangle_freq);
        for(j = lidstart; j < lidend; ++j)
            {
            fbank[i][j] = lslope * ((j * FS / FREQSIZE) - low_triangle_freq);
            }
        for(j = ridstart; j < ridend; ++j)
            {
            fbank[i][j] = rslope * (high_triangle_freq - (j * FS / FREQSIZE));
            }
        }
    free(freqs);
    return fbank;
}

static int *chroma_matrices_init // TODO: should build table for faster initialization (if we can decide appopriate signal length)
    (
    void
    )
{
    const int SIGNALSIZE = featext_wksp.signalsize;
    const int FREQSIZE = SIGNALSIZE >> 1;
    const int FS = featext_wksp.fs;
    const int FSEG = FS / SIGNALSIZE;
    int *chroma_mat = (int *)malloc(SIGNALSIZE * sizeof(int));
    int *num_freqs_per_chroma = (int *)malloc(FREQSIZE * sizeof(int));
    int i = 0;
    for(i = 0; i < FREQSIZE; ++i)
        {
        chroma_mat[i] = (int)round(12.0 * log2((i + 1) * FSEG / 27.50));
        if(chroma_mat[i] < 0)
            {
            chroma_mat[i] += FREQSIZE;
            }
        }
    i = 0;
    while(i < FREQSIZE)
        {
        int j = i + 1;
        int accum = 1;
        while(chroma_mat[j] == chroma_mat[i])
            {
            ++accum;
            ++j;
            }
        // TODO: need check, maybe unportable, sizeof(wchar_t) should equal to 32
        wmemset(num_freqs_per_chroma + i, accum, accum);
        i = j;
        }
    for(i = 0; i < FREQSIZE; ++i)
        {
        chroma_mat[i + FREQSIZE] = num_freqs_per_chroma[chroma_mat[i]];
        }
    free(num_freqs_per_chroma);
    return chroma_mat;
}

static void normalize
    (
    const short *src,
    double *dst
    )
{
    const int SIGNALSIZE = featext_wksp.signalsize;
    const double SHRT_RANGE = -SHRT_MIN;
    double signal_dc = 0.0;
    double signal_max = 0.0;
    int i = 0;
    for(i = 0; i < SIGNALSIZE; ++i)
        {
        double absdst = 0.0;
        dst[i] = src[i] / SHRT_RANGE;
        signal_dc += dst[i];
        absdst = fabs(dst[i]);
        if(absdst > signal_max)
            {
            signal_max = absdst;
            }
        }
    signal_dc /= SIGNALSIZE;
    for(i = 0; i < SIGNALSIZE; ++i)
        {
        dst[i] = (dst[i] - signal_dc) / signal_max;
        }
#ifdef VERBOSE
    printf("mean:%f, max:%f\n", signal_dc, signal_max);
#endif
}

void feature_extraction_release
    (
    void
    )
{
    const int NUM_LIN_FILTER = 13;
    const int NUM_LOG_FILTER = 27;
    const int NUM_TOTAL_FILTER = NUM_LOG_FILTER + NUM_LIN_FILTER;
    int i = 0;
    free(featext_wksp.features);
    free(featext_wksp.chroma_nummat);
    for(i = 0; i < NUM_TOTAL_FILTER; ++i)
        {
        free(featext_wksp.filterbanks[i]);
        }
    free(featext_wksp.filterbanks);
    free(featext_wksp.signal_prev_freq_domain);
}

double *feature_extraction_main
    (
    short *signal,
    int frame_idx
    )
{
    const int SIGNALSIZE = featext_wksp.signalsize;
    const int FREQSIZE = SIGNALSIZE >> 1;
    double *feature = featext_wksp.features;
    double *x = (double *)malloc(SIGNALSIZE * sizeof(double));      // normalized time signal
    double *fx = (double *)calloc(FREQSIZE, sizeof(double));        // normalized fft magnitude
    double *fx_prev = featext_wksp.signal_prev_freq_domain;
    double fx_sum = 0.0;
    double fx_sum_prev = featext_wksp.signal_prev_freq_domain_sum;
    double fx_squaresum = 0.0;

    normalize(signal, x);

    // time-domain features
    zcr(x, &(feature[0]));
    energy(x, &(feature[1]));
    energy_entropy(x, feature[1], &(feature[2]));
    // to frequency-domain
    dft_normalize(x, fx, &fx_sum, &fx_squaresum);
    // frequency-domain features
    spectral_centroid_spread(fx, fx_sum, &(feature[3]));
    spectral_entropy(fx, fx_squaresum, &(feature[5]));
    if(frame_idx == 0)
        {
        feature[6] = 0.0; // TODO: should remove if() in the while loop
        }
    else
        {
        spectral_flux(fx, fx_sum, fx_prev, fx_sum_prev, &(feature[6]));
        }
    spectral_rolloff(fx, fx_squaresum, &(feature[7]));
    mfcc(fx, &(feature[8]));
    chroma(fx, fx_squaresum, &(feature[21]));

#ifdef VERBOSE
    int i = 0;
    for(i = 0; i < FEATURE_LEN; ++i)
        {
        printf("[%d]:%f ", i, feature[i]);
        }
    printf("\n");
#endif

    featext_wksp.signal_prev_freq_domain_sum = fx_sum;
    memcpy(fx_prev, fx, FREQSIZE * sizeof(double));

    free(x);
    free(fx);
    return feature;
}

void feature_extraction_init
    (
    const int FS,
    const int SIGLEN
    )
{
    featext_wksp.fs = FS;
    featext_wksp.signalsize = SIGLEN;
    featext_wksp.filterbanks = mfcc_filter_banks_init();
    featext_wksp.chroma_nummat = chroma_matrices_init();
    featext_wksp.features = (double *)malloc(FEATURE_LEN * sizeof(double));
    featext_wksp.signal_prev_freq_domain = (double *)malloc((SIGLEN >> 1) * sizeof(double));
    featext_wksp.signal_prev_freq_domain_sum = 0.0;
}
