#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include <likwid.h>

extern void sweep(float *, float *, float *, float, float, float, float, float, int, int, int);

uint64_t
get_time()
{
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_sec * 1000000 + now.tv_usec;
}

void
init(float *grid, int X, int Y, int Z)
{
    uint64_t x, y, z;
#pragma omp parallel for schedule(static)
    for (z=0; z<Z; ++z)
        for (y=0; y<Y; ++y)
            for (x=0; x<X; ++x)
                if (y == 0 || x == 0)
                    grid[z*X*Y+y*X+x] = 1.0;
                else
                    grid[z*X*Y+y*X+x] = 0.0;
}

int
main(int argc, char *argv[])
{
    uint64_t data_set_size, X, Y, Z;
    X = Y = Z = 1400; // 1400^3*4 = ~10.x GB per grid

    if (argc > 1)
        X = Y = Z = atoi(argv[1]);

    fprintf(stderr, "grid size: %dx%dx%d\n", X, Y, Z);

    float *U, *V, *ROC;
    if ((U = (float *)_mm_malloc(X*Y*Z*sizeof(float), 64)) == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    if ((V = (float *)_mm_malloc(X*Y*Z*sizeof(float), 64)) == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    if ((ROC = (float *)_mm_malloc(X*Y*Z*sizeof(float), 64)) == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }


    // init grids
    init(U, X, Y, Z);
    init(V, X, Y, Z);
    init(ROC, X, Y, Z);

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_START("stencil");

    int T;
    float runtime=0.0;
    uint64_t start = get_time();
    for (T=0; runtime < 300.0; ++T) {

        sweep(U, V, ROC, 1.0, 1.0, 1.0, 1.0, 1.0, X, Y, Z);

        runtime = ((float)get_time() - (float)start) / 1e6;
    }

    LIKWID_MARKER_STOP("stencil");
    LIKWID_MARKER_CLOSE;

    // calculate and output MUp/s
    double MUps = ((double)(Z-8)*(X-8)*(Y-8)*T)/runtime/1e6;
    printf("size: %dx%dx%d MLUP/s: %f iterations: %d\n", X, Y, Z, MUps, T);

    // free memory
    _mm_free(U);
    _mm_free(V);
    _mm_free(ROC);

    return 0;
}
