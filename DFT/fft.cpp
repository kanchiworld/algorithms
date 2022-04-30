/*
 * Compile as g++ -std=gnu++11 fft.cpp
 * std::tuple is available only in c++11
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <map>
#include <vector>
#include <tuple>

struct Complex
{
    double real;
    double imag;
};

// These functions can be rewritten by operator overloading
struct Complex complexMultiply(struct Complex c1, struct Complex c2)
{
    struct Complex retVal = { 0, 0 };
    retVal.real = c1.real * c2.real - c1.imag * c2.imag;
    retVal.imag = c1.real * c2.imag + c2.real * c1.imag;

    return retVal;
}

struct Complex complexAdd(struct Complex c1, struct Complex c2)
{
    struct Complex retVal = { 0, 0 };
    retVal.real = c1.real + c2.real ;
    retVal.imag = c1.imag + c2.imag ;

    return retVal;
}

struct Complex kth_nthRootsOfUnity(unsigned int k, unsigned int n)
{
// if n < 1 return error
// if k > n-1 return error

    struct Complex retVal = { 0, 0};
    retVal.real = cos((2* M_PI / n)*k);
    retVal.imag = -1 * sin((2* M_PI / n )*k);

    return retVal;
}

// Canonical implementation of DFT. Very rudimentary. O(n**2) complexity. Creates a DFT matrix and uses the definition of DFT to compute DFT vector
void DFT_Simple(unsigned int numSamples, struct Complex *timeDomainSignals, struct Complex *freqDomainSignals)
{
    struct Complex *nthRoU = (struct Complex *) malloc(numSamples * sizeof(struct Complex));
    for (unsigned int idx=0; idx <= numSamples-1; idx++)
    {
        nthRoU[idx] = kth_nthRootsOfUnity(idx, numSamples);
        //printf("debug: nthRoU[%d] = { %f, %f }\n", idx, nthRoU[idx].real, nthRoU[idx].imag);
    }

    for (unsigned int k=0; k <= numSamples-1; k++)
    {
        freqDomainSignals[k].real = 0; freqDomainSignals[k].imag = 0 ;
        for (unsigned int n=0; n<=numSamples-1; n++)
        {
        freqDomainSignals[k] = complexAdd( freqDomainSignals[k], complexMultiply( nthRoU[(k*n)%numSamples], timeDomainSignals[n] ));
        //printf("debug: [ %d %d ] ", (k*n)%numSamples, n);
        }
    }
}

struct Complex DFT_Row(unsigned int numSamples, unsigned int offset, unsigned int step,  struct Complex *timeDomainSignals, int row)
{
    struct Complex *nthRoU = (struct Complex *) malloc(numSamples * sizeof(struct Complex));
    for (unsigned int idx=0; idx <= numSamples-1; idx++)
    {
        nthRoU[idx] = kth_nthRootsOfUnity(idx, numSamples);
    }

    struct Complex   freqDomainSignal = { 0, 0 } ;
    for (unsigned int n=0; n<=numSamples-1; n++)
    {
        freqDomainSignal = complexAdd( freqDomainSignal, complexMultiply( nthRoU[(row*n)%numSamples], timeDomainSignals[offset + n*step] ));
    }
    return freqDomainSignal;
}

typedef std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> ffttuple;
struct Complex FFT_Row(unsigned int numSamples, unsigned int offset, unsigned int step,  struct Complex *timeDomainSignals, unsigned int row)
{
//    if numSamples is ! power of 2 return error
//    printf("debug: FFT_Row, numSamples = %d, offset = %d, step = %d [", numSamples, offset, step);
//    for (int idx=0; idx<=numSamples-1; idx++) printf("%d,", offset+idx*step);
//    printf("]\n");

    static std::map<ffttuple, struct Complex > precompute;

    struct Complex retVal = { 0, 0 };

    if (numSamples == 2)
    {
    retVal = DFT_Row(numSamples, offset, step, timeDomainSignals, row);
    }
    else
    {
    struct Complex fftEven = { 0, 0 };
    struct Complex fftOdd = {0, 0};
    if ( row < numSamples/2 )
    {
        fftEven = FFT_Row( numSamples/2, offset, 2*step,  timeDomainSignals, row);
        ffttuple t1{numSamples/2, offset, 2*step, row};
        precompute[t1] = fftEven;

        fftOdd = FFT_Row( numSamples/2, offset+step, 2*step ,  timeDomainSignals, row);
        ffttuple t2{numSamples/2, offset+step, 2*step, row};
        precompute[t2] = fftOdd;
    }
    else
    {
        ffttuple t1{numSamples/2, offset, 2*step, row-numSamples/2};
        ffttuple t2{numSamples/2, offset+step, 2*step, row-numSamples/2};
        fftEven = precompute[t1];
        fftOdd = precompute[t2];
    }

    retVal = complexAdd( fftEven, complexMultiply(kth_nthRootsOfUnity(row, numSamples), fftOdd ));
    }
    return retVal;
}



int main(void)
{
#define MAX_N 1024*64
    struct Complex timeDomainSignals[MAX_N];
    struct Complex freqDomainSignals1[MAX_N];
    struct Complex freqDomainSignals2[MAX_N];

    int N = 0;
    std::cout << "specify array dimension (MUST be power of 2)" << std::endl;
    std::cin >> N;

    double sample = 0.0;
    for(int i = 0; i <= N-1; i++)
    {
    std::cout << "specify element number: " << i << std::endl;
    std::cin >> sample;
    timeDomainSignals[i] = {sample, 0};
    }

    DFT_Simple(N, timeDomainSignals, freqDomainSignals1);
    for (int idx=0; idx <= N-1; idx++)
    {
        printf("freqDomainSignals1[%d] = { %f, %f }\n", idx, freqDomainSignals1[idx].real, freqDomainSignals1[idx].imag);
    }
    printf("\n");

    for (unsigned int idx=0; idx<=N-1; idx++)
    {
        freqDomainSignals2[idx] = FFT_Row( N, 0, 1, timeDomainSignals, idx);
        printf("freqDomainSignals2[%d] = { %f, %f }\n", idx, freqDomainSignals2[idx].real, freqDomainSignals2[idx].imag);
    }

    return 0;
}


