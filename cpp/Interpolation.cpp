#include"WaveFront.h"

// interpolation method
complex<double> WaveFront::GetInterpolatedValueNEAREST_NEIGHBOR(double u, double v)const
{
	// caluculate index as double format
	double tx = u / _px + _nx / 2;
	double ty = v / _py + _ny / 2;

	// caluculate index as integer format
	int i1 = (int)ceil(tx);
	int j1 = (int)ceil(ty);

	if (i1 < 0)
		return complex<double>(0.0, 0.0);
	if (i1 >= _nx)
		return complex<double>(0.0, 0.0);
	if (j1 < 0)
		return complex<double>(0.0, 0.0);
	if (j1 >= _ny)
		return complex<double>(0.0, 0.0);

	complex<double> ret = GetPixel(i1, j1);

	return ret;
}
complex<double> WaveFront::GetInterpolatedValueBILINEAR(double u, double v)const
{
	// caluculate index as double format
	double tx = u / _px + _nx / 2;
	double ty = v / _py + _ny / 2;

	// caluculate index as integer format
	int i1 = (int)floor(tx);
	int j1 = (int)floor(ty);

	if (i1 < 0 || i1 > _nx - 1 || j1 < 0 || j1 > _ny - 1)
		return complex<double>(0.0, 0.0);

	int i2, j2;
	if (i1 == _nx - 1)
		i2 = i1;
	else
		i2 = i1 + 1;

	if (j1 == _ny - 1)
		j2 = j1;
	else
		j2 = j1 + 1;

	double tt = (tx - i1);
	double uu = (ty - j1);
	double t1 = 1.0 - tt;
	double u1 = 1.0 - uu;

	complex<double> ret = 
		  GetPixel(i2, j2) * (tt * uu)
		+ GetPixel(i1, j2) * (uu * t1)
		+ GetPixel(i1, j1) * (t1 * u1)
		+ GetPixel(i2, j1) * (u1 * tt);

	return ret;
}
static void Cubic4(double bufx[8], double x1)
{
	double t1 = x1 - 1.0;
	double t2 = x1 * t1;
	double t3 = t2 - 1.0;

	double p[4] = { -t1, t1, -x1, x1 };
	double q[4] = { t2,  t3, t3,  t2 };
	int i;
//#pragma omp simd
	for (i = 0; i < 4; i++) {
		bufx[i] = p[i] * q[i];
	}
	return;
}
complex<double> WaveFront::GetInterpolatedValueBICUBIC(double u, double v)const
{
	// caluculate index as double format
	double tx = u / _px + _nx / 2;
	double ty = v / _py + _ny / 2;

	int i, j;
	int tx_int = (int)tx;
	int ty_int = (int)ty;
	double tx_dec = tx - tx_int;
	double ty_dec = ty - ty_int;

	int is = 0;
	int js = 0;
	int in = 4;
	int jn = 4;
	int iflag = 1;
	int i_bufpos = tx_int - 1;
	int j_bufpos = ty_int - 1;


	if (i_bufpos < 0) {
		is = -i_bufpos;
		iflag = 0;
	}
	if (j_bufpos < 0) {
		js = -j_bufpos;
	}
	if (in > _nx - i_bufpos) {
		in = _nx - i_bufpos;
		iflag = 0;
	}
	if (jn > _ny - j_bufpos) {
		jn = _ny - j_bufpos;
	}
	if (is >= in || js >= jn) {
		return complex<double>(0.0, 0.0);
	}

	double bufx[4];
	double bufy[4];
	double dbuf[4][8] = { 0 };

	for (j = js; j < jn; j++) {
		for (i = is; i < in; i++) {
			dbuf[j][i * 2] = (GetPixel(i + i_bufpos, j + j_bufpos + js).real());
			dbuf[j][i * 2 + 1] = (GetPixel(i + i_bufpos, j + j_bufpos + js).imag());
		}
	}

	Cubic4(bufx, tx_dec);
	Cubic4(bufy, ty_dec);

	double c1r, c1i;
	double c2r = 0, c2i = 0;

	for (j = 0; j < 4; j++) {
		c1r = 0;
		c1i = 0;
		for (i = 0; i < 4; i++) {
			c1r += dbuf[j][i * 2] * bufx[i];
			c1i += dbuf[j][i * 2 + 1] * bufx[i];
		}
		c2r += c1r * bufy[j];
		c2i += c1i * bufy[j];
	}
	return complex<double>(c2r, c2i);
}
complex<double> WaveFront::GetInterpolatedValue(double u, double v, Interpol interp)const
{
	complex<double> ret;
	switch (interp)
	{
	case NEAREST_NEIGHBOR:
		ret = GetInterpolatedValueNEAREST_NEIGHBOR(u, v);
		return ret;
	case BILINEAR:
		ret = GetInterpolatedValueBILINEAR(u, v);
		return ret;
	case BICUBIC:
		ret = GetInterpolatedValueBICUBIC(u, v);
		return ret;
	}
}