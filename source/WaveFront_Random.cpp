#include"../include/WaveFront.h"
using namespace std;

std::random_device WaveFront::w_rnd;
std::uniform_real_distribution<> WaveFront::w_randvul(-1.0, 1.0);
std::mt19937 WaveFront::w_mt(w_rnd());

WaveFront& WaveFront::ModRandomphase()
{
	QueryPerformanceFrequency(&w_freq);
	QueryPerformanceCounter(&w_start);
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			//double phase = static_cast<double>(rand()) / RAND_MAX * 2 * PI - PI;
			//double phase = static_cast<double>(getrandom(-1,1)) * PI;
			double phase = static_cast<double>(getrandomMinusOneToOne()) * PI;
			double amp = GetAmplitude(i, j);
			SetPixel(i, j, complex<double>(amp * cos(phase), amp * sin(phase)));
		}
	}
	QueryPerformanceCounter(&w_end);
	w_time_random += getdeltatime();
	return *this;
}
double WaveFront::getrandom(double min, double max)
{
	//random_device rnd;// generate seed
	//mt19937 mt(w_rnd());// 32bitMT
	uniform_real_distribution<> randvul(min, max);
	//return randvul(mt);

	return randvul(w_mt);
}
double WaveFront::getrandomMinusOneToOne()
{
	//mt19937 mt(w_rnd());// 32bitMT
	//return w_randvul(mt);

	return w_randvul(w_mt);
}