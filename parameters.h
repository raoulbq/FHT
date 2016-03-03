#ifndef PARAMETERS_H
#define PARAMETERS_H

#define PI 3.1415926535897932385

// Number of basis functions
#define BIGN 128

// Number of data points
#define LITN (BIGN*5)

// Location of the turning point
#define BIGC sqrt(2.0*BIGN + 1.0)

// General recursion formula
// P_{n+1}(x) = (a_{n} x + b_{n}) P_{n}(x) - c_{n} P_{n-1}(x)

// Chebyshev recursion:
// T_{-1}(x)  = 0
// T_{0}(x)   = 1
// T_{n+1}(x) = 2 x T_{n}(x) - T_{n-1}(x)
// Evaluation on [-1,1] instead of [-BIGC, BIGC]
#define ALPHA (2.0 / BIGC)
#define BETA 0.0
#define GAMMA (-1.0)

// Hermite recursion:
// H_{-1}(x)  = 0
// H_{0}(x)   = D0
// H_{n+1}(x) = \sqrt{2}{n+1} x H_{n}(x) - \sqrt{n}{n+1} H_{n-1}(x)
#define D0 (1.0 / pow(PI, 1.0/4.0))
#define AL(l) sqrt(2.0 / (l + 1.0))
#define BL(l) 0.0
#define CL(l) (-sqrt(l / (l + 1.0)))

// Define ul, vl and wl
#define UL(l) (AL(l) / ALPHA)
#define VL(l) (BL(l) - BETA * UL(l))
#define WL(l) (- GAMMA * UL(l))

#endif
