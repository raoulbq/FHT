#ifndef PARAMETERS_H
#define PARAMETERS_H

// Recursion parameters
#define PI 3.1415926535897932385
#define BIGN 128
#define LITN (BIGN*5)

#define BIGC sqrt((double)2*BIGN + 1)

#define ALPHA (2.0 / BIGC)
#define BETA 0
#define GAMMA (-1.0)

#define D0 ((double)1.0 / pow(PI, (double)1.0/4.0))

// Define al, bl, cl, ul, vl and wl
#define AL(l) (double)sqrt((double)2.0 / (l+1.0))
#define BL(l) (double)0
#define CL(l) (double)((-1.0) * sqrt((double)l / (l+1.0)))

#define UL(l) (double)(AL(l) / ALPHA)
#define VL(l) (double)(BL(l) + (-1) * (UL(l) * BETA))
#define WL(l) (double)((-1) * GAMMA * UL(l))

#endif
