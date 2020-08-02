#include"..\include\WaveFront.h"
using namespace std;
// set basic aperture
void WaveFront::GenerateCirc(double r)
{
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < w_nx; i++)
	{
		for (j = 0; j < w_ny; j++)
		{
			double x = itox(i);
			double y = jtoy(j);

			if (x * x + y * y <= r * r)
			{
				SetPixel(i, j, complex<double>(1.0, 0.0));
			}
			else
			{
				SetPixel(i, j, complex<double>(0.0, 0.0));
			}
		}
	}
}
void WaveFront::GenerateRect(double wx, double wy)
{
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < w_nx; i++)
	{
		for (j = 0; j < w_ny; j++)
		{
			double x = itox(i);
			double y = jtoy(j);

			if (abs(x) <= wx / 2 && abs(y) <= wy / 2)
			{ 
				SetPixel(i, j, complex<double>(1.0, 0.0)); 
			}
			else
			{
				SetPixel(i, j, complex<double>(0.0, 0.0));
			}
		}
	}
		
}
void WaveFront::GenerateGaussian(double r, double n)
{
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < w_nx; i++)
	{
		for (j = 0; j < w_ny; j++)
		{
			double x = itox(i);
			double y = jtoy(j);
			double power = pow(sqrt(x * x + y * y) / r, n);
			SetPixel(i, j, complex<double>(exp(-power), 0.0));
		}
	}	
}