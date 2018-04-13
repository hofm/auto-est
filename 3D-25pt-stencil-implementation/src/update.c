#include <stdint.h>

void
sweep(float *U, float *V, float *ROC, float c0, float c1, float c2, float c3, float c4, uint64_t X, uint64_t Y, uint64_t Z)
{
    uint64_t x, y, z;
    float lap, c0t, c1t, c2t, c3t, c4t;
#pragma omp parallel for schedule(static)
    for (z=4; z<(Z-4); ++z) {
        for (y=4; y<(Y-4); ++y) {
//#pragma vector aligned
            for (x=4; x<(X-4); ++x) {
                c4t  = V[(z  )*Y*X+(y  )*X+(x-4)];
                c3t  = V[(z  )*Y*X+(y  )*X+(x-3)];
                c2t  = V[(z  )*Y*X+(y  )*X+(x-2)];
                c1t  = V[(z  )*Y*X+(y  )*X+(x-1)];
                c0t  = V[(z  )*Y*X+(y  )*X+(x  )];
                c1t += V[(z  )*Y*X+(y  )*X+(x+1)];
                c2t += V[(z  )*Y*X+(y  )*X+(x+2)];
                c3t += V[(z  )*Y*X+(y  )*X+(x+3)];
                c4t += V[(z  )*Y*X+(y  )*X+(x+4)];

                c4t += V[(z  )*Y*X+(y-4)*X+(x  )];
                c3t += V[(z  )*Y*X+(y-3)*X+(x  )];
                c2t += V[(z  )*Y*X+(y-2)*X+(x  )];
                c1t += V[(z  )*Y*X+(y-1)*X+(x  )];
                c1t += V[(z  )*Y*X+(y+1)*X+(x  )];
                c2t += V[(z  )*Y*X+(y+2)*X+(x  )];
                c3t += V[(z  )*Y*X+(y+3)*X+(x  )];
                c4t += V[(z  )*Y*X+(y+4)*X+(x  )];

                c4t += V[(z-4)*Y*X+(y  )*X+(x  )];
                c3t += V[(z-3)*Y*X+(y  )*X+(x  )];
                c2t += V[(z-2)*Y*X+(y  )*X+(x  )];
                c1t += V[(z-1)*Y*X+(y  )*X+(x  )];
                c1t += V[(z+1)*Y*X+(y  )*X+(x  )];
                c2t += V[(z+2)*Y*X+(y  )*X+(x  )];
                c3t += V[(z+3)*Y*X+(y  )*X+(x  )];
                c4t += V[(z+4)*Y*X+(y  )*X+(x  )];

                lap = c0*c0t + c1*c1t + c2*c2t + c3*c3t + c4*c4t;

                U[z*Y*X+y*X+x] = 2.0f * V[z*Y*X+y*X+x] - U[z*Y*X+y*X+x] + ROC[z*Y*X+y*X+x] *  lap;
            }
        }
    }
}
