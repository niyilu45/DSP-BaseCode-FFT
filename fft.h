#ifndef FFT_H
#define FFT_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include "ComplexNum.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

typedef struct ComplexNum{
	double re;
	double im;
}ComplexNum;

typedef struct OmegaLibNode* OmegaLibList;
typedef struct OmegaLibNode{
	OmegaLibList Next;
	int Len;
	ComplexNum* Omega;
}OmegaLibNode;


int CalFFTLen(int InputLen);
ComplexNum* FFT(ComplexNum* Input, int InputLen);
ComplexNum* IFFT(ComplexNum* Input, int InputLen);
void ClearOmegaLib(void); // After FFT, should clear lib.



static ComplexNum CalOmega(int i, int N);
static ComplexNum ComplexAddInFFT(ComplexNum A, ComplexNum B);
static ComplexNum ComplexSubInFFT(ComplexNum A, ComplexNum B);
static ComplexNum ComplexMulInFFT(ComplexNum A, ComplexNum B);
static ComplexNum* FFTDataFlip(ComplexNum* Input, int InputLen);
static unsigned int BitFlipInt(register unsigned int val, int BitWidth);
static int CalFFTStateNum(int InputLen);

// Omega Lib
static OmegaLibList FindOmegaLib(int FFTLen, OmegaLibList Head);
static OmegaLibList FindAndInsertOmegaLib(int FFTLen, OmegaLibList Head);
static OmegaLibList InsertOmegaLib(int FFTLen, OmegaLibList L);
static ComplexNum* CalcOmegaLib(int FFTLen);
static OmegaLibList DestroyOmegaLib(OmegaLibList L);
#endif
