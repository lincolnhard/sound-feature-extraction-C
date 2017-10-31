#include <sndfile.h>
#include "feature.h"

static void print_usage
    (
    void
    )
{
    printf("\nUsage : ./alarm-detector [input sound file] [winsize (secs)] [stepsize (secs)]\n");
}

int main
    (
    int argc,
    char *argv[]
    )
{
    char *infilename = NULL;
    SNDFILE	*infile = NULL;
    SF_INFO	sfinfo;
    double win = 0.0;
    double step = 0.0;
    int winsize = 0;
    int stepsize = 0;
    int frameidx = 0;
    int readcount = 0;
    short *clip = NULL;

    if(argc != 4)
        {
        print_usage();
        return EXIT_FAILURE;
        }

    memset(&sfinfo, 0, sizeof(sfinfo));
    infilename = argv[1];

    if((infile = sf_open(infilename, SFM_READ, &sfinfo)) == NULL)
        {
        printf ("Not able to open input sound file %s.\n", infilename);
        puts (sf_strerror (NULL));
        return EXIT_FAILURE;
        }

    win = atof(argv[2]);
    step = atof(argv[3]);
    winsize = (int)(win * sfinfo.samplerate);
    stepsize = (int)(step * sfinfo.samplerate);
    clip = (short *)malloc(winsize * sizeof(short));
    feature_extraction_init(sfinfo.samplerate, winsize);
    while(1)
        {
        sf_seek(infile, frameidx * stepsize, SEEK_SET);
        readcount = sf_readf_short(infile, clip, winsize);
        if(readcount != winsize)
            {
            break;
            }
        printf("-%d-\n", frameidx);

        feature_extraction_main(clip, frameidx);

        ++frameidx;
        }
    feature_extraction_release();

    sf_close(infile);
    free(clip);
    return 0 ;
}
