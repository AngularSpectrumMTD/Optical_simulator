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
enum Out
{
	REAL,
	IMAGE,
	AMPLITUDE,
	INTENSITY,
	PHASE
};
#endif

#ifndef INTPOL
enum Interp
{
	NEAREST_NEIGHBOR,
	BILINEAR,
	BICUBIC
};
#endif

#ifndef AXIS
enum Axis
{
	X_AXIS,
	Y_AXIS
};
#endif
#include<windows.h>
class WaveFront {
	int w_nx = 0;
	int w_ny = 0;
	double w_px = 1e-6;
	double w_py = 1e-6;
	double w_lambda = 633e-9;
	vec3 w_origin = vec3{0.0,0.0,0.0};
	bool w_real = true;
	vec3 w_normal = vec3{0.0,0.0,1.0};
	std::unique_ptr<std::complex<double>[]> w_data;
	static std::random_device w_rnd;

	double w_time_fft = 0;
	double w_time_interpol = 0;
	double w_time_random = 0;

	LARGE_INTEGER w_freq;
	LARGE_INTEGER w_start;
	LARGE_INTEGER w_end;
	double getdeltatime() { return (w_end.QuadPart - w_start.QuadPart) * 1000.0 / w_freq.QuadPart; }

public:
	void dispTotalTime() { printf("FFT: %lf [ms] INTERPOL: %lf [ms] RANDOM: %lf [ms]", w_time_fft, w_time_interpol, w_time_random); }
	bool ispow2(unsigned int n) { return (n != 0) && (n & (n - 1)) == 0; }
	WaveFront() {};
	WaveFront(int nx, int ny, double px, double py, double lambda);
	WaveFront(const WaveFront& a);
	~WaveFront() {};
	// getter
	int GetNx() const{ return w_nx; }
	int GetNy() const { return w_ny; }
	int GetN() const { return w_nx * w_ny; }
	double GetPx() const { return w_px; }
	double GetPy() const { return w_py; }
	double GetLambda() const { return w_lambda; }
	vec3 GetOrigin() const { return w_origin; }
	bool GetSpace() const { return w_real; }
	vec3 GetNormal() const { return w_normal; }
	double GetWidth() const { return w_px * w_nx; }
	double GetHeight() const { return w_py * w_ny; }
	// setter
	void SetNx(const unsigned int nx) { 
		if (ispow2(nx) == false)
		{
			double expX = log(static_cast<double>(nx)) / log(2.0);
			w_nx = nearPow2(pow(2.0, expX));
		}
		else
		{
			w_nx = nx;
		}
	}
	void SetNy(const unsigned int ny) {
		if (ispow2(ny) == false)
		{
			double expY = log(static_cast<double>(ny)) / log(2.0);
			w_ny = nearPow2(pow(2.0, expY));
		}
		else
		{
			w_ny = ny;
		}
	}
	void SetPx(const double px) { w_px = px; }
	void SetPy(const double py) { w_py = py; }
	void SetLambda(const double lambda) { w_lambda = lambda; }
	void SetSpace(const bool flug) { w_real = flug; }
	void SetPixel(const unsigned int i, const unsigned int j, const std::complex<double> val) { w_data[idxij(i, j)] = val; }
	void SetReal(const unsigned int i, const unsigned int j, const double val) {
		double amp = GetAmplitude(i, j);
		double image = sqrt(amp*amp - val*val);
		w_data[idxij(i, j)] = std::complex<double>(val, image);
	}
	void SetImage(const unsigned int i, const unsigned int j, const double val) {
		double amp = GetAmplitude(i, j);
		double real = sqrt(amp * amp - val * val);
		w_data[idxij(i, j)] = std::complex<double>(real, val);
	}
	void SetOrigin(const const vec3 &vec) { w_origin = vec; }
	void SetNormal(const const vec3 &vec) { w_normal = vec; }
	// basic function
	int idxij(const int i, const int j) const { return j * w_nx + i; }// return index
	void CopyParam(const WaveFront& wave) {
		w_px = wave.w_px;
		w_py = wave.w_py;
		w_lambda = wave.w_lambda;
		w_origin = wave.w_origin;
		w_real = wave.w_real;
		w_normal = wave.w_normal;
	};
	void Init() { 
		//w_data.reset(new std::complex<double>[w_nx * w_ny]);

		std::complex<double>* mempool = (std::complex<double>*)malloc(sizeof(std::complex<double>) * w_nx * w_ny);

		if (mempool == NULL)
		{
			printf(">>ERROR: Required memory is too large. Memory allocation is mistaken...\n");
			printf(">>Process is terminated forcibly...\n");
			exit(0);
		}

		w_data.reset(mempool);
	};

	//BasicFunc
	unsigned int xtoi(double x)const;// convert distance along to x-axis to index i
	unsigned int ytoj(double y)const;// convert distance along to y-axis to index j 
	double itox(unsigned int i)const;// convert index i to distance along to x-axis
	double jtoy(unsigned int j)const;// convert index j to distance along to y-axis
	std::complex<double> GetPixel(unsigned int i, unsigned int j)const;
	double GetReal(unsigned int i, unsigned int j)const;
	double GetImage(unsigned int i, unsigned int j)const;
	double GetAmplitude(unsigned int i, unsigned int j)const;
	double GetIntensity(unsigned int i, unsigned int j)const;
	double GetPhase(unsigned int i, unsigned int j)const;
	double ComputeEnergy()const;
	double ComputeMaxAmplitude()const;
	void Clear();// reset all complex value as 0
	void DispMat(mat3 &mat)const;
	void DispVec(vec3 &vec)const;
	void DispParam()const;
	unsigned int nearPow2(int n);// return most near power of 2
	WaveFront& AllSet(double val);
	WaveFront& MultiplyPlaneWave(double u, double v, double phase = 0.0);
	WaveFront& Add(const WaveFront &source);
	WaveFront& TransformforBrainImage();

	// SaveLoadImage (function related to save as image and load real part from image)
	void SaveBmp(const char* filename, Out type);// save as image
	void Normalize();// normalize
	WaveFront& LoadBmp(const char* filename);// load bmp

	// Aperture (set basic aperture)
	void GenerateCirc(double r);// circular aperture
	void GenerateRect(double wx, double wy);// rectangular aperture
	void GenerateGaussian(double r, double n);// gaussian beam

	// Interpolation (interpolation method)
	std::complex<double> GetInterpolatedValueNEAREST_NEIGHBOR(double u, double v)const;// nearest neighbor
	std::complex<double> GetInterpolatedValueBILINEAR(double u, double v)const;// bi-linear
	std::complex<double> GetInterpolatedValueBICUBIC(double u, double v)const;// bi-cubic
	std::complex<double> GetInterpolatedValue(double u, double v, Interp interp)const;// return interpolated value by use of chosed method
	
	// FFT (function related to fft)
	void swap();// swap orthant
	void fft1D(std::unique_ptr <std::complex<double>[]>& x, int n, int func);// 1D fft
	void fft2D(int func);// 2D fft

	// function related to angular spectrum method
	void pitchtrans() { w_px = 1.0 / (w_px * w_nx); w_py = 1.0 / (w_py * w_ny); };// convert sampling interval
	// Propagation
	void generateFRF(double distance);// generate frequential function
	void bandlimit(double uband, double vband);// bandlimit by use of chosed frequency
	void AsmProp(const double R);// execute calculation of optical diffraction by use of angular spectrum method
	WaveFront& ShiftedAsmProp(const WaveFront& source);// off-axis propagation 
	WaveFront& ShiftedAsmPropAdd(const WaveFront& source);
	WaveFront& ShiftedAsmPropEx(const WaveFront& source);
	WaveFront& ShiftedAsmPropAddEx(const WaveFront& source);
	void AsmPropInFourierSpace(const double R);
	void Embed();// embed distribution by 4 times
	void Extract();// extruct distribution by 4 times
	void ExactAsmProp(const double R);// exacterASM
	double GetWpow2(double u, double v);// return power of w-element of spatial frequency(distinguish evanescent region)
	mat3 GetRotMat(const vec3 &v) const;// get rotational matrix from 2 vectors v and normal vector of this field
	mat3 GetRotMat(const vec3& global, const vec3& local) const;
	mat3 L2G() const;
	vec3 GetUnitVector_alongX() const;
	vec3 GetUnitVector_alongY() const;
	void RotInFourierSpace(WaveFront& reference, Interp interp, vec3 *carrier = nullptr);// rotation in fourier space
	void TiltedAsmProp(WaveFront& reference, Interp interp, vec3* carrier = nullptr);//execute calculation of optical diffraction between non-parallel planes
	
	// Random (function related to random)
	WaveFront& ModRandomphase();// randomize phase
	double getrandom(double min, double max);// generate random

	// SaveCsv (function related to save as csv)
	void SaveAsCsv(const char* fname, Axis axis, int ij);

	// QuadraticPhase
	WaveFront& SetQuadraticPhase(const double f);
	WaveFront& MultiplyQuadraticPhase(const double f);

	// SaveLoad (function related to complex distribution)
	void SaveAsWaveFront(const char* filename);
	WaveFront& LoadAsWaveFront(const char* filename);
	
	// Operator
	WaveFront& operator+=(double val);
	WaveFront& operator-=(double val);
	WaveFront& operator*=(double val);
	WaveFront& operator/=(double val);
	WaveFront& operator+=(const WaveFront& val);
	WaveFront& operator-=(const WaveFront& val);
	WaveFront& operator*=(const WaveFront& val);
	WaveFront& operator/=(const WaveFront& val);
	WaveFront& operator =(const WaveFront& val);
};

#endif