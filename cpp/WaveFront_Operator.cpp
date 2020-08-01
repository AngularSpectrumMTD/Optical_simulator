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
	if (!(w_nx == val.w_nx && w_ny == val.w_ny && w_px == val.w_px && w_py == val.w_py && w_lambda == val.w_lambda))
	{
		printf(">>ERROR: operator +=: all parameters must be same\n");
		printf(">>Process is terminated forcibly...\n");
		system("pause");
		exit(0);
	}
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
	if (!(w_nx == val.w_nx && w_ny == val.w_ny && w_px == val.w_px && w_py == val.w_py && w_lambda == val.w_lambda))
	{
		printf(">>ERROR: operator +=: all parameters must be same\n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
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
	if (!(w_nx == val.w_nx && w_ny == val.w_ny && w_px == val.w_px && w_py == val.w_py && w_lambda == val.w_lambda))
	{
		printf(">>ERROR: operator +=: all parameters must be same\n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
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
	if (!(w_nx == val.w_nx && w_ny == val.w_ny && w_px == val.w_px && w_py == val.w_py && w_lambda == val.w_lambda))
	{
		printf(">>ERROR: operator +=: all parameters must be same\n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
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
	if (!(w_nx == val.w_nx && w_ny == val.w_ny && w_px == val.w_px && w_py == val.w_py && w_lambda == val.w_lambda))
	{
		printf(">>ERROR: operator +=: all parameters must be same\n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
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