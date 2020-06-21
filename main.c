#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"
static const double CoeffLib[75] = {
0.000449632392674000,	0.00201762209861100,	0.00347151717132100,	
0.00445799229265400,	0.00465235358163400,	0.00383452697156200,	
0.00195681033836200,	-0.000812373082690000,	-0.00408199174157700,	
-0.00727682517446500,	-0.00972036497659300,	-0.0107526452539830,	
-0.00986496581225800,	-0.00682743828878100,	-0.00178322708657800,	
0.00471382655278100,	0.0117337888426980,		0.0180813319417600,	
0.0224742520685990,		0.0237616546120790,		0.0211456776448470,	
0.0143844278169160,		0.00392032744109200,	-0.00909200084337300,	
-0.0228742456014840,	-0.0352265516543490,	-0.0438207845632670,	
-0.0465450405666470,	-0.0418477862504380,	-0.0290263026288240,	
-0.00840954554753000,	0.0186003313777910,		0.0496427151033210,	
0.0816601889554240,		0.111274036392512,		0.135221358061849,	
0.150790311277688,		0.156188706604522,		0.150790311277688,	
0.135221358061849,		0.111274036392512,		0.0816601889554240,		
0.0496427151033210,		0.0186003313777910,		-0.00840954554753000,	
-0.0290263026288240,	-0.0418477862504380,	-0.0465450405666470,	
-0.0438207845632670,	-0.0352265516543490,	-0.0228742456014840,	
-0.00909200084337300,	0.00392032744109200,	0.0143844278169160,	
0.0211456776448470,		0.0237616546120790,		0.0224742520685990,	
0.0180813319417600,		0.0117337888426980,		0.00471382655278100,	
-0.00178322708657800,	-0.00682743828878100,	-0.00986496581225800,	
-0.0107526452539830,	-0.00972036497659300,	-0.00727682517446500,	
-0.00408199174157700,	-0.000812373082690000,	0.00195681033836200,	
0.00383452697156200,	0.00465235358163400,	0.00445799229265400,	
0.00347151717132100,	0.00201762209861100,	0.000449632392674000
};

int main()
{
	int i;
	int FFTLen = 2048;
	ComplexNum* Coeff;

	Coeff = (ComplexNum*)calloc(FFTLen, sizeof(ComplexNum));
	for(i=0;i<75;i++){
		Coeff[i].re = CoeffLib[i];
	}

	ComplexNum* FFTOut;
	ComplexNum* IFFTOut;

	FFTOut = FFT(Coeff, FFTLen);
	IFFTOut = IFFT(FFTOut, FFTLen);


	printf("输出FFT后的结果\n");
	for(i=0;i<20;i++){
		printf("%f+i%f\n",FFTOut[i].re, FFTOut[i].im);
	}
	printf("输出FFT后的结果\n");
	for(i=0;i<20;i++){
		printf("%f+i%f\n",IFFTOut[i].re, IFFTOut[i].im);
	}

	free(Coeff);
	free(FFTOut);
	free(IFFTOut);
	return 0;
}

