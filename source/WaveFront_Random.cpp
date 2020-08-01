#include"..\include\WaveFront.h"
using namespace std;
WaveFront& WaveFront::ModRandomphase()
{
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			double phase = static_cast<double>(rand()) / RAND_MAX * 2 * PI - PI;
			double amp = GetAmplitude(i, j);
			SetPixel(i, j, complex<double>(amp * cos(phase), amp * sin(phase)));
		}
	}

	return *this;
}
double WaveFront::getrandom(double min, double max)
{
	random_device rnd;// generate seed
	mt19937 mt(rnd());// 32bitMT
	uniform_real_distribution<> randvul(min, max);
	return randvul(mt);
}