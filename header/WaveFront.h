#ifndef WAVEFRONT_H
#define WAVEFRONT_H

//#include"Image.h"

#include<stdio.h>
#include<math.h>
#include<complex>
#include<memory>
#include<iostream>
#include <algorithm>
#include <functional>
#include<random>

#include <fstream>
#include <sstream>

#include"Image.h"

#ifndef OUT
typedef enum OUT
{
	REAL,
	IMAGE,
	AMPLITUDE,
	INTENSITY,
	PHASE
}Outputformat;
#endif

#ifndef INTPOL
typedef enum INTPOL
{
	NEAREST_NEIGHBOR,
	BILINEAR,
	BICUBIC
}Interpol;
#endif

#ifndef AXIS
typedef enum AXIS
{
	X_AXIS,
	Y_AXIS
}Axis;
#endif

class WaveFront {

	int _nx = 0;
	int _ny = 0;
	double _px = 1e-6;
	double _py = 1e-6;
	double _lambda = 633e-9;
	vec3 _origin = vec3{0.0,0.0,0.0};
	bool _real = true;
	vec3 _normal = vec3{0.0,0.0,1.0};
	unique_ptr<complex<double>[]> _data;
public:
	bool ispow2(unsigned int n) { return (n != 0) && (n & (n - 1)) == 0; }
	WaveFront() {};
	WaveFront(int nx, int ny, double px, double py, double lambda)
		:_nx(nx), _ny(ny), _px(px), _py(py), _lambda(lambda) {

		double expX = log((double)nx) / log(2.0);
		double expY = log((double)ny) / log(2.0);

		if (ispow2(nx) == false)
		_nx = nearPow2(pow(2.0, expX));
		if (ispow2(ny) == false)
		_ny = nearPow2(pow(2.0, expY));

		_data.reset(new complex<double>[_nx * _ny]);
		//_normal = vec3{ 0.0, 0.0, 0.0 };
		_normal = vec3{ 0.0, 0.0, 1.0 };
	};
	WaveFront(const WaveFront& a) :_nx(a._nx), _ny(a._ny), _px(a._px), _py(a._py), _lambda(a._lambda) {
		_origin = a._origin;
		_real = a._real;
		_normal = a._normal;
		_data.reset(new complex<double>[_nx * _ny]);
		for (int i = 0; i < _nx; i++)
			for (int j = 0; j < _ny; j++)
				_data[idxij(i, j)] = a._data[idxij(i, j)];
	};
	~WaveFront() {};
	// getter
	int GetNx() const{ return _nx; }
	int GetNy() const { return _ny; }
	int GetN() const { return _nx * _ny; }
	double GetPx() const { return _px; }
	double GetPy() const { return _py; }
	double GetLambda() const { return _lambda; }
	vec3 GetOrigin() const { return _origin; }
	bool GetSpace() const { return _real; }
	vec3 GetNormal() const { return _normal; }
	// setter
	void SetNx(unsigned int nx) { 
		if (ispow2(nx) == false)
		{
			double expX = log((double)nx) / log(2.0);
			_nx = nearPow2(pow(2.0, expX));
		}
		else
			_nx = nx;
	}
	void SetNy(unsigned int ny) { 
		if (ispow2(ny) == false)
		{
			double expY = log((double)ny) / log(2.0);
			_ny = nearPow2(pow(2.0, expY));
		}
		else
			_ny = ny;
	}
	void SetPx(double px) { _px = px; }
	void SetPy(double py) { _py = py; }
	void SetLambda(double lambda) { _lambda = lambda; }
	void SetSpace(bool flug) { _real = flug; }
	void SetPixel(unsigned int i, unsigned int j, complex<double> val) { _data[idxij(i, j)] = val; }
	void SetReal(unsigned int i, unsigned int j, double val) { 
		double amp = GetAmplitude(i, j);
		double image = sqrt(amp*amp - val*val);
		_data[idxij(i, j)] = complex<double>(val, image); 
	}
	void SetImage(unsigned int i, unsigned int j, double val) { 
		double amp = GetAmplitude(i, j);
		double real = sqrt(amp * amp - val * val);
		_data[idxij(i, j)] = complex<double>(real, val);
	}
	void SetOrigin(vec3 vec) { _origin = vec; }
	void SetNormal(vec3 vec) { _normal = vec; }
	// basic function
	int idxij(int i, int j) const { return j * _nx + i; }// return index
	void CopyParam(WaveFront& wave) {
		_px = wave._px;
		_py = wave._py;
		_lambda = wave._lambda;
		_origin = wave._origin;
		_real = wave._real;
		_normal = wave._normal;
	};
	void Init() { _data.reset(new complex<double>[_nx * _ny]); };
	unsigned int xtoi(double x)const;// convert distance along to x-axis to index i
	unsigned int ytoj(double y)const;// convert distance along to y-axis to index j 
	double itox(unsigned int i)const;// convert index i to distance along to x-axis
	double jtoy(unsigned int j)const;// convert index j to distance along to y-axis
	complex<double> GetPixel(unsigned int i, unsigned int j)const;
	double GetReal(unsigned int i, unsigned int j)const;
	double GetImage(unsigned int i, unsigned int j)const;
	double GetAmplitude(unsigned int i, unsigned int j)const;
	double GetIntensity(unsigned int i, unsigned int j)const;
	double GetPhase(unsigned int i, unsigned int j)const;
	double GetEnergy()const;
	double GetMaxAmplitude()const;
	void Clear();// reset all complex value as 0
	void DispMat(mat3 mat)const;
	void DispVec(vec3 vec)const;
	unsigned int nearPow2(int n);// return most near power of 2
	// function related to save as image and load real part from image
	void SaveBmp(const char* filename, Outputformat type);// save as image
	void Normalize();// normalize
	WaveFront& LoadBmp(const char* filename);// load bmp
	// set basic aperture
	void SetCirc(double r);// circular aperture
	void SetRect(double wx, double wy);// rectangular aperture
	void SetGaussian(double r, double n);// gaussian beam
	// interpolation method
	complex<double> GetInterpolatedValueNEAREST_NEIGHBOR(double u, double v)const;// nearest neighbor
	complex<double> GetInterpolatedValueBILINEAR(double u, double v)const;// bi-linear
	complex<double> GetInterpolatedValueBICUBIC(double u, double v)const;// bi-cubic
	complex<double> GetInterpolatedValue(double u, double v, Interpol interp)const;// return interpolated value by use of chosed method
	// function related to fft
	void swap();// swap orthant
	void fft1D(unique_ptr <complex<double>[]>& x, int n, int func);// 1D fft
	void fft2D(int func);// 2D fft
	// function related to angular spectrum method
	void pitchtrans() { _px = 1.0 / (_px * _nx); _py = 1.0 / (_py * _ny); };// convert sampling interval
	void generateFRF(double distance);// generate frequential function
	void bandlimit(double uband, double vband);// bandlimit by use of chosed frequency
	void AsmProp(double R);// execute calculation of optical diffraction by use of angular spectrum method
	void AsmPropInFourierSpace(double R);
	void Embed();// embed distribution by 4 times
	void Extract();// extruct distribution by 4 times
	void ExactAsmProp(double R);// exacterASM
	// function related to rotational transformation of wave-field
	double GetWpow2(double u, double v);// return power of w-element of spatial frequency(distinguish evanescent region)
	mat3 GetRotMat(vec3 v);// get rotational matrix from 2 vectors v and normal vector of this field
	mat3 GetMat(vec3 v1, vec3 v2);// get matrix from 2 vectors
	mat3 GetRotMat2(vec3 v1, vec3 v2);
	void RotInFourierSpace(WaveFront& reference, Interpol interp);// rotation in fourier space
	void TiltedAsmProp(WaveFront& reference, Interpol interp);//execute calculation of optical diffraction between non-parallel planes
	// function related to random
	WaveFront& ModRandomphase();// randomize phase
	double getrandom(double min, double max);// generate random
	// function related to save as csv
	void SaveAsCsv(const char* fname, Axis axis, int ij);
	// QuadraticPhase
	WaveFront& SetQuadraticPhase(double f);
	WaveFront& MultiplyQuadraticPhase(double f);
	// function related to complex distribution
	void SaveAsWaveFront(const char* filename);
	WaveFront& LoadAsWaveFront(const char* filename);

	double GetAmplitudeMax() const
	{
		vector<double> amp;
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				amp.push_back(GetAmplitude(i, j));
		
		vector<double>::iterator ite = max_element(amp.begin(), amp.end());
		return *ite;
	}

	WaveFront& AllSet(double val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, complex<double>(val, 0.0));
		return *this;
	}

	WaveFront& operator+=(double val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i, j) + val);
		return *this;
	}

	WaveFront& operator-=(double val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i, j) - val);
		return *this;
	}

	WaveFront& operator*=(double val)
	{ 
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i,j)*val);
		return *this;
	}

	WaveFront& operator/=(double val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i, j) / val);
		return *this;
	}

	WaveFront& operator+=(const WaveFront& val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i, j) + val.GetPixel(i, j));
		return *this;
	}

	WaveFront& operator-=(const WaveFront& val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i, j) - val.GetPixel(i, j));
		return *this;
	}

	WaveFront& operator*=(const WaveFront &val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i, j) * val.GetPixel(i, j));
		return *this;
	}

	WaveFront& operator/=(const WaveFront& val)
	{
		int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (j = 0; j < _ny; ++j)
			for (i = 0; i < _nx; ++i)
				SetPixel(i, j, GetPixel(i, j) / val.GetPixel(i, j));
		return *this;
	}

	WaveFront& operator =(const WaveFront& val)
	{
		if (GetNx() == val.GetNx() && GetNy() == val.GetNy())
		{
			int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
			for (j = 0; j < _ny; ++j)
				for (i = 0; i < _nx; ++i)
					SetPixel(i, j, val.GetPixel(i, j));
		}
		return *this;
	}
};

#endif