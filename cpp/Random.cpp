#include"WaveFront.h"

WaveFront& WaveFront::ModRandomphase()
{
	int i, j;
	double amp = 0;
	double phase = 0;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; ++j)
		for (i = 0; i < _nx; ++i)
		{
			phase = (float)((double)rand() / RAND_MAX * 2 * PI - PI);
			amp = GetAmplitude(i, j);
			SetPixel(i, j, complex<double>(amp * cos(phase), amp * sin(phase)));
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