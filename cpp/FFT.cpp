#include"WaveFront.h"

// function related to fft
void WaveFront::swap()
{
	int i, j;
	int x, y;
	WaveFront tmp(_nx, _ny, _px, _py, _lambda);
#pragma omp parallel for private(x, y) num_threads(omp_get_num_threads())
	for (y = 0; y < _ny; y++)
	{
		if (_ny / 2 <= y && y < _ny)
			j = (y - _ny / 2);
		else
			j = (y + _ny / 2);
		for (x = 0; x < _nx; x++)
		{
			if (_nx / 2 <= x && x < _nx)
				i = (x - _nx / 2);
			else
				i = (x + _nx / 2);
			tmp.SetPixel(i, j, this->GetPixel(x, y));
		}
	}
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; j++)
		for (i = 0; i < _nx; i++)
			this->SetPixel(i, j, tmp.GetPixel(i, j));
}
int inv2pow(int n) {
	int m = 0;
	m = (int)(log((double)n) / log(2.0));
	return m;
}
void WaveFront::fft1D(unique_ptr <complex<double>[]>& x, int n, int func)
{
	int pow_2;
	int mainloop;
	int mx, i, j, k;
	double phase_int, arg;
	complex<double> w, w_tmp, x_tmp1, x_tmp2;
	pow_2 = inv2pow(n);
	mx = n;
	phase_int = 2.0 * PI / (double)n;

	for (mainloop = 0; mainloop < pow_2; mainloop++)
	{
		int switchloop;
		int mx_1;
		mx_1 = mx - 1;
		mx /= 2;
		arg = 0.0;
		for (switchloop = 0; switchloop < mx; switchloop++)
		{
			w = complex<double>(cos(-func * arg), sin(-func * arg));
			arg += phase_int;
			for (i = mx_1; i < n; i += (mx_1 + 1))
			{
				int j1, j2;
				j1 = i - mx_1 + switchloop;
				j2 = j1 + mx;
				x_tmp1 = x[j1] + x[j2];
				x_tmp2 = x[j1] - x[j2];
				x[j1] = x_tmp1;
				x[j2] = x_tmp2 * w;
			}
		}
		phase_int *= 2.0;
	}
	j = 0;
	for (i = 0; i < n - 1; i++)
	{
		complex<double> x_tmp;
		if (i < j)
		{
			x_tmp = x[i];
			x[i] = x[j];
			x[j] = x_tmp;
		}
		k = n / 2;
		while (k <= j)
		{
			j = j - k;
			k /= 2;
		}
		j = j + k;
	}
}
void WaveFront::fft2D(int func)
{
	this->pitchtrans();
	this->swap();
	int i, j;
	unique_ptr <complex<double>[]> tempx;
	tempx.reset(new complex<double>[_nx]);

#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < _ny; j++)
	{
		for (i = 0; i < _nx; i++)
			tempx[i] = this->GetPixel(i, j);
		fft1D(tempx, _nx, func);
		for (i = 0; i < _nx; i++)
			SetPixel(i, j, tempx[i]);
	}
	unique_ptr <complex<double>[]> tempy;
	tempy.reset(new complex<double>[_ny]);
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < _nx; i++)
	{
		for (j = 0; j < _ny; j++)
			tempy[j] = this->GetPixel(i, j);
		fft1D(tempy, _ny, func);
		for (j = 0; j < _ny; j++)
			SetPixel(i, j, tempy[j]);
	}
	this->swap();
}