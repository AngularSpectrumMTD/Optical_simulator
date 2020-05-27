//#define STB_IMAGE_IMPLEMENTATION
//#define STB_IMAGE_WRITE_IMPLEMENTATION
#include"WaveFront.h"

// basic function
unsigned int WaveFront::xtoi(double x)const
{
	return unsigned int(x / _px + _nx / 2);
}
unsigned int WaveFront::ytoj(double y)const
{
	return unsigned int(y / _py + _ny / 2);
}
double WaveFront::itox(unsigned int i)const
{
	return double(((double)i - (double)_nx / 2.0f) * _px);
}
double WaveFront::jtoy(unsigned int j)const
{
	return double(((double)j - (double)_ny / 2.0f) * _py);
}
complex<double> WaveFront::GetPixel(unsigned int i, unsigned int j)const
{
	return _data[idxij(i, j)];
}
double WaveFront::GetReal(unsigned int i, unsigned int j)const
{
	double ret;
	ret = GetPixel(i, j).real();
	return ret;
}
double WaveFront::GetImage(unsigned int i, unsigned int j)const
{
	double ret;
	ret = GetPixel(i, j).imag();
	return ret;
}
double WaveFront::GetAmplitude(unsigned int i, unsigned int j)const
{
	double ret;
	ret = abs(GetPixel(i, j));
	return ret;
}
double WaveFront::GetIntensity(unsigned int i, unsigned int j)const
{
	double ret;
	ret = GetAmplitude(i, j) * GetAmplitude(i, j);
	return ret;
}
double WaveFront::GetPhase(unsigned int i, unsigned int j)const
{
	double ret;
	ret = arg(GetPixel(i, j));
	return ret;
}
double WaveFront::GetEnergy()const
{
	double ret = 0;
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; ++j)
		for (i = 0; i < _nx; ++i)
			ret += GetIntensity(i, j);
	return ret;
}
double WaveFront::GetMaxAmplitude()const
{
	double ret = 0;
	int i, j;
	vector<double> amplitude;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; ++j)
		for (i = 0; i < _nx; ++i)
			amplitude.push_back(GetAmplitude(i, j));

	vector<double>::iterator ite = max_element(amplitude.begin(),amplitude.end());
	size_t index = distance(amplitude.begin(),ite);

	ret = amplitude[index];

	return ret;
}
void WaveFront::Clear()
{
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; ++j)
		for (i = 0; i < _nx; ++i)
			SetPixel(i, j, complex<double>(0.0, 0.0));
}
void WaveFront::DispMat(mat3 mat)const
{
	printf("\n%f %f %f\n%f %f %f\n%f %f %f"
		, mat.getCol0().getX(), mat.getCol1().getX(), mat.getCol2().getX()
		, mat.getCol0().getY(), mat.getCol1().getY(), mat.getCol2().getY()
		, mat.getCol0().getZ(), mat.getCol1().getZ(), mat.getCol2().getZ());
}
void WaveFront::DispVec(vec3 vec)const
{
	printf("\n%f \n%f \n%f\n"
		, vec.getX(),vec.getY(),vec.getZ());
}
unsigned int WaveFront::nearPow2(int n)// return most near power of 2
{
	if (n <= 0) return 0;

	if ((n & (n - 1)) == 0) return (unsigned int)n;

	unsigned int ret = 1;
	while (n > 0) { ret <<= 1; n >>= 1; }
	return ret;
}




