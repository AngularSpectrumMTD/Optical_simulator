#include"WaveFront.h"

WaveFront& WaveFront::SetQuadraticPhase(double f)
{
	double k = 2*PI/_lambda;
	double x, y;
	double phase;
	double amp;
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; ++j)
		for (i = 0; i < _nx; ++i)
		{
			x = itox(i); y = jtoy(j);
			phase = -k * (x * x + y * y) / (2 * f);
			amp = GetAmplitude(i,j);
			SetPixel(i, j, complex<double>(amp * cos(phase), amp * sin(phase)));
		}
	return *this;
}
WaveFront& WaveFront::MultiplyQuadraticPhase(double f)
{
	double k = 2 * PI / _lambda;
	double x, y;
	double phase;
	int i, j;
	complex<double> val;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; ++j)
		for (i = 0; i < _nx; ++i)
		{
			x = itox(i); y = jtoy(j);
			phase = -k * (x * x + y * y) / (2 * f);
			val = GetPixel(i, j);
			SetPixel(i, j, val * complex<double>(cos(phase), sin(phase)));
		}
	return *this;
}