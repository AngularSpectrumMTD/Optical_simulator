//#define STB_IMAGE_IMPLEMENTATION
//#define STB_IMAGE_WRITE_IMPLEMENTATION
#include"WaveFront.h"
using namespace std;
// basic function
unsigned int WaveFront::xtoi(double x)const
{
	return unsigned int(x / w_px + w_nx / 2);
}
unsigned int WaveFront::ytoj(double y)const
{
	return unsigned int(y / w_py + w_ny / 2);
}
double WaveFront::itox(unsigned int i)const
{
	return double((static_cast<double>(i) - static_cast<double>(w_nx) / 2.0f)* w_px);
}
double WaveFront::jtoy(unsigned int j)const
{
	return double((static_cast<double>(j) - static_cast<double>(w_ny) / 2.0f)* w_py);
}
complex<double> WaveFront::GetPixel(unsigned int i, unsigned int j)const
{
	return w_data[idxij(i, j)];
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
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			ret += GetIntensity(i, j);
		}
	}
	return ret;
}
double WaveFront::GetMaxAmplitude()const
{
	int i, j;
	vector<double> amplitude(w_nx * w_ny);
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			amplitude[idxij(i, j)] = GetAmplitude(i, j);
			//amplitude.push_back(GetAmplitude(i, j));
		}
	}

	vector<double>::iterator ite = max_element(amplitude.begin(), amplitude.end());
	return *ite;
}
void WaveFront::Clear()
{
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, complex<double>(0.0, 0.0));
		}
	}
}
void WaveFront::DispMat(mat3& mat)const
{
	printf("\n%f %f %f\n%f %f %f\n%f %f %f"
		, mat.getCol0().getX(), mat.getCol1().getX(), mat.getCol2().getX()
		, mat.getCol0().getY(), mat.getCol1().getY(), mat.getCol2().getY()
		, mat.getCol0().getZ(), mat.getCol1().getZ(), mat.getCol2().getZ());
}
void WaveFront::DispVec(vec3& vec)const
{
	printf("\n%f \n%f \n%f\n"
		, vec.getX(), vec.getY(), vec.getZ());
}
void WaveFront::DispParam()const
{
	printf("\nNx: %d\nNy: %d\nPx: %lf\nPy: %lf\nLambda: %lf",
		w_nx, w_ny, w_px, w_py, w_lambda);
}
unsigned int WaveFront::nearPow2(int n)// return most near power of 2
{
	if (n <= 0) { return 0; }

	if ((n & (n - 1)) == 0) { return (unsigned int)n; }

	unsigned int ret = 1;
	while (n > 0) { ret <<= 1; n >>= 1; }
	return ret;
}
WaveFront& WaveFront::AllSet(double val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, std::complex<double>(val, 0.0));
		}
	}
	return *this;
}
WaveFront& WaveFront::MultiplyPlaneWave(double u, double v, double phase)
{
	if (!((u >= -1.0 && u <= +1.0) && (v >= -1.0 && v <= +1.0) && (u * u + v * v <= 1.0)))
	{
		printf(">>ERROR: Absolute value of element of direction must be fitted in 0 to 1...\n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}

	double k = 2 * PI / w_lambda;

	double	phase0 = k * (u * GetOrigin().getX() + v * GetOrigin().getY());

	int i, j;

#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < GetNy(); j++)
		for (i = 0; i < GetNx(); i++)
		{
			double x = (double)itox(i);
			double y = (double)jtoy(j);
			double p = k * (u * x + v * y) + phase0 + phase;
			complex<double> c = complex<double>(cos(p), sin(p));
			this->SetPixel(i, j, this->GetPixel(i, j) * c);
		}

	return *this;
}
WaveFront& WaveFront::Add(const WaveFront& source)
{
	WaveFront& frame = *this;
	//matching coordinate
	double x = (source.itox(0) + source.GetOrigin().getX()) - (frame.itox(0) + frame.GetOrigin().getX());
	double y = (source.jtoy(0) + source.GetOrigin().getY()) - (frame.jtoy(0) + frame.GetOrigin().getY());

	//converting from real coordinate value to integer index
	int ii = int(x / source.GetPx() + 0.5);
	int jj = int(y / source.GetPy() + 0.5);

	//source's coordinate + (ii, jj) = frame's coordinate
	if ((0 + ii) >= frame.GetNx() || (0 + jj) >= frame.GetNy())
	{
		//this->Clear();
		return *this;	//No overlap
	}
	if ((source.GetNx() + ii) < 0 || (source.GetNy() + jj) < 0)
	{
		//this->Clear();
		return *this;	//No overlap
	}

	int is0 = 0, is1 = source.GetNx();
	int js0 = 0, js1 = source.GetNy();

	if ((0 + ii) < 0)
	{
		is0 = 0 - ii;
	}
	if ((0 + jj) < 0)
	{
		js0 = 0 - jj;
	}
	if ((source.GetNx() + ii) > frame.GetNx())
	{
		is1 = frame.GetNx() - ii;
	}
	if ((source.GetNy() + jj) > frame.GetNy())
	{
		js1 = frame.GetNy() - jj;
	}

	// adding field from source to frame
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (int j = js0; j < js1; j++)
	{
		for (int i = is0; i < is1; i++)
		{
			frame.SetPixel(i + ii, j + jj, frame.GetPixel(i + ii, j + jj) + source.GetPixel(i, j));
		}
	}

	return *this;
}
WaveFront& WaveFront::TransformforBrainImage()
{
	WaveFront tmp(*this);
	for (int i = 0; i < w_nx; i++)
	{
		for (int j = 0; j < w_ny; j++)
		{
			SetPixel(i, j, tmp.GetPixel(w_nx - i - 1, w_ny - j - 1));
		}
	}
	return *this;
}




