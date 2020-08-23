#include"../include/WaveFront.h"
using namespace std;
WaveFront& WaveFront::SetQuadraticPhase(const double f)
{
	double k = 2*PI/w_lambda;
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			double x = itox(i); double y = jtoy(j);
			double phase = -k * (x * x + y * y) / (2 * f);
			double amp = GetAmplitude(i, j);
			SetPixel(i, j, complex<double>(amp * cos(phase), amp * sin(phase)));
		}
	}
	return *this;
}
WaveFront& WaveFront::MultiplyQuadraticPhase(const double f)
{
	double k = 2 * PI / w_lambda;
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			double x = itox(i); double y = jtoy(j);
			double phase = -k * (x * x + y * y) / (2 * f);
			complex<double> val = GetPixel(i, j);
			SetPixel(i, j, val * complex<double>(cos(phase), sin(phase)));
		}
	}
	return *this;
}