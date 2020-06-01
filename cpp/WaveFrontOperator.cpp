#include"WaveFront.h"

WaveFront& WaveFront::operator+=(double val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) + val);
		}
	}
	return *this;
}

WaveFront& WaveFront::operator-=(double val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) - val);
		}
	}
	return *this;
}

WaveFront& WaveFront::operator*=(double val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) * val);
		}
	}
	return *this;
}

WaveFront& WaveFront::operator/=(double val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) / val);
		}
	}
	return *this;
}

WaveFront& WaveFront::operator+=(const WaveFront& val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) + val.GetPixel(i, j));
		}
	}
	return *this;
}

WaveFront& WaveFront::operator-=(const WaveFront& val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) - val.GetPixel(i, j));
		}
	}
	return *this;
}

WaveFront& WaveFront::operator*=(const WaveFront& val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) * val.GetPixel(i, j));
		}
	}
	return *this;
}

WaveFront& WaveFront::operator/=(const WaveFront& val)
{
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, GetPixel(i, j) / val.GetPixel(i, j));
		}
	}
	return *this;
}

WaveFront& WaveFront::operator =(const WaveFront& val)
{
	if (GetNx() == val.GetNx() && GetNy() == val.GetNy())
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < w_ny; ++j)
		{
			for (i = 0; i < w_nx; ++i)
			{
				SetPixel(i, j, val.GetPixel(i, j));
			}
		}
	}
	return *this;
}