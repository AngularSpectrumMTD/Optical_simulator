#include"WaveFront.h"

// set basic aperture
void WaveFront::SetCirc(double r)
{
	double x = 0, y = 0;
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < _nx; i++)
		for (j = 0; j < _ny; j++)
		{
			x = itox(i);
			y = jtoy(j);

			if (x * x + y * y <= r*r)
				SetPixel(i, j, complex<double>(1.0, 0.0));
			else
				SetPixel(i, j, complex<double>(0.0, 0.0));
		}
}
void WaveFront::SetRect(double wx, double wy)
{
	double x = 0, y = 0;
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < _nx; i++)
		for (j = 0; j < _ny; j++)
		{
			x = itox(i);
			y = jtoy(j);

			if (abs(x) <= wx / 2 && abs(y) <= wy / 2)
				SetPixel(i, j, complex<double>(1.0, 0.0));
			else
				SetPixel(i, j, complex<double>(0.0, 0.0));
		}
}
void WaveFront::SetGaussian(double r, double n)
{
	double x = 0, y = 0;
	int i, j;
	double power;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < _nx; i++)
		for (j = 0; j < _ny; j++)
		{
			x = itox(i);
			y = jtoy(j);
			power = pow(sqrt(x * x + y * y) / r, n);
			SetPixel(i, j, complex<double>(exp(-power), 0.0));
		}
}