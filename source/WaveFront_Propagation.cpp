#include"../include/WaveFront.h"
using namespace std;
// function related to angular spectrum method
void WaveFront::generateFRF(double distance)
{
	double k = 1 / w_lambda / w_lambda;
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_max_threads())
	for (i = 0; i < w_nx; i++)
	{
		for (j = 0; j < w_ny; j++)
		{
			double kx = itox(i);
			double ky = jtoy(j);
			double w = k - kx * kx - ky * ky;
			if (w > 0)
			{
				double phase = 2 * PI * sqrt(w) * distance;
				this->SetPixel(i, j, complex<double>(cos(phase), sin(phase)));
			}
			else//evanescent element
			{
				this->SetPixel(i, j, complex<double>(1.0, 0.0));
			}
		}
	}
}
void WaveFront::bandlimit(double uband, double vband)
{
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_max_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			if (abs(itox(i)) <= uband / 2 && abs(jtoy(j)) <= vband / 2)
			{
				continue;
			}
			else
			{
				SetPixel(i, j, complex<double>(0.0, 0.0));
			}
		}
	}
}

void WaveFront::AsmProp(const double R)
{
	fft2D(-1);
	AsmPropInFourierSpace(R);
	fft2D(1);
	*this /= GetN();
}
void WaveFront::AsmPropInFourierSpace(const double R)
{
	WaveFront h(w_nx, w_ny, w_px, w_py, w_lambda);
	h.generateFRF(R);
	double uband = 2 / (sqrt((2 * h.w_px * R) * (2 * h.w_px * R) + 1) * h.w_lambda);
	double vband = 2 / (sqrt((2 * h.w_py * R) * (2 * h.w_py * R) + 1) * h.w_lambda);
	h.bandlimit(uband, vband);// bandlimit
	int x, y;
#pragma omp parallel for private(x, y) num_threads(omp_get_max_threads())
	for (y = 0; y < w_ny; ++y)
	{
		for (x = 0; x < w_nx; ++x)
		{
			SetPixel(x, y, GetPixel(x, y) * h.GetPixel(x, y));
		}
	}
	SetOrigin(GetOrigin() + GetNormal() * R);
}
void WaveFront::Embed()
{
	WaveFront tmp(*this);// generate same distribution
	w_nx *= 2;
	w_ny *= 2;
	Init();
	Clear();
	int i, j;
#pragma omp parallel for num_threads(omp_get_max_threads())
	for (j = 0; j < tmp.w_ny; ++j)
	{
		for (i = 0; i < tmp.w_nx; ++i)
		{
			SetPixel(i + w_nx / 4, j + w_ny / 4, tmp.GetPixel(i, j));
		}
	}
}
void WaveFront::Extract()
{
	WaveFront tmp(*this);// generate same distribution
	w_nx /= 2;
	w_ny /= 2;
	Init();
	Clear();
	int i, j;
#pragma omp parallel for num_threads(omp_get_max_threads())
	for (j = 0; j < w_ny; ++j)
	{
		for (i = 0; i < w_nx; ++i)
		{
			SetPixel(i, j, tmp.GetPixel(i + tmp.w_nx / 4, j + tmp.w_ny / 4));
		}
	}
}
void WaveFront::ExactAsmProp(const double R)
{
	Embed();
	AsmProp(R);
	Extract();
}
// function related to rotational transformation of wave-field
double WaveFront::GetWpow2(double u, double v)
{
	return 1 / w_lambda / w_lambda - u * u - v * v;
}
mat3 WaveFront::GetRotMat(const vec3& v) const
{
	mat3 ret = mat3::identity();
	if ((v.getX() == w_normal.getX()) && (v.getY() == w_normal.getY()) && (v.getY() == w_normal.getY()))
	{
		return ret;
	}
	else
	{
		vec3 normalw = normalize(w_normal), normalv = normalize(v);
		double angle = acos(dot(normalw, normalv));
		vec3 axis = normalize(cross(normalw, normalv));
		mat3 ret = mat3::rotation(angle, axis);
		return ret;
	}
}
mat3 WaveFront::GetRotMat(const vec3& global, const vec3& local) const
{
	mat3 ret = mat3::identity();
	vec3 nglobal = normalize(global);
	vec3 nlocal = normalize(local);

	if (dot(nglobal, nlocal) > 0.999)//if surface normal vector's z_elem nearly 1, matrix is invalid, so return identity matrix
	{
		return ret;
	}
	double angle = acos(dot(global, local));
	vec3 axis = normalize(cross(global, local));
	ret = mat3::rotation(angle, axis);
	return ret;
}
mat3 WaveFront::L2G() const
{
	mat3 ret = GetRotMat(GetNormal(), vec3{ 0.0,0.0,1.0 });
	return ret;
}
vec3 WaveFront::GetUnitVector_alongX() const
{
	vec3 vec{ 1.0,0.0,0.0 };
	mat3 mat = L2G();
	return mat * vec;
}
vec3 WaveFront::GetUnitVector_alongY() const
{
	vec3 vec{ 0.0,1.0,0.0 };
	mat3 mat = L2G();
	return mat * vec;
}

void WaveFront::RotInFourierSpace(WaveFront& source, Interp interp, vec3* carrier)
{
	QueryPerformanceFrequency(&w_freq);
	QueryPerformanceCounter(&w_start);
	WaveFront& reference = *this;
	mat3 rot = reference.GetRotMat(source.w_normal);//get rotation matrix

	mat3 invrot = transpose(rot);//get inverce matrix of it
	double invL1 = 1 / w_lambda;
	vec3 source0{ 0.0,0.0,(float)invL1 };
	source0 = rot * source0;//carrier frequency
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_max_threads())
	for (j = 0; j < reference.w_ny; ++j)
	{
		for (i = 0; i < reference.w_nx; ++i)
		{
			//frequency in shifted fourier space
			vec3 ref = vec3{ static_cast<float>(reference.itox(i)),static_cast<float>(reference.jtoy(j)),static_cast<float>(invL1) };
			//w element of ref
			double w_ref = reference.GetWpow2(ref.getX(), ref.getY());
			if (w_ref < 0)//if frequency is evanescent, must be ignored
			{
				reference.SetPixel(i, j, complex<double>(0.0, 0.0));
				continue;
			}
			//frequency in reference fourier space
			vec3 shift = ref + source0;
			//w element of shift
			double w_shift = reference.GetWpow2(shift.getX(), shift.getY());
			if (w_shift < 0)//if frequency is evanescent, must be ignored
			{
				reference.SetPixel(i, j, complex<double>(0.0, 0.0));
				continue;
			}
			shift.setZ(static_cast<float>(sqrt(w_shift)));
			//frequency in source fourier space
			vec3 sourcefreq = invrot * shift;
			reference.SetPixel(i, j, source.GetInterpolatedValue(sourcefreq.getX(), sourcefreq.getY(), interp));
		}
	}
	//reference *= static_cast<float>(sqrt(GetEnergy() / reference.GetEnergy()));

	if (carrier != nullptr)
	{
		*carrier = source0;
	}
	QueryPerformanceCounter(&w_end);
	w_time_interpol += getdeltatime();
}
void WaveFront::TiltedAsmProp(WaveFront& source, Interp interp, vec3* carrier)
{
	WaveFront& reference = *this;
	source.fft2D(-1);
	reference.pitchtrans();// convert sampling interval from real to frequential

	reference.RotInFourierSpace(source, interp, carrier);// rotational transform

	reference.fft2D(1);
	reference /= GetN();
}
void CalcBandLimitFreq(double* high, double* low, double xy, double z, double duv)
{
	double positive = xy + 1.0 / (2.0 * duv);
	double negative = xy - 1.0 / (2.0 * duv);
	double freq_limitP = 1.0 / sqrt(1.0 / positive / positive *z * z + 1);
	double freq_limitN = 1.0 / sqrt(1.0 / negative / negative * z * z + 1);

	if (xy > 1.0 / (2 * duv))
	{
		*low = +freq_limitN;
		*high = +freq_limitP;
	}
	else if (xy < -1.0 / (2 * duv))
	{
		*low = -freq_limitN;
		*high = -freq_limitP;
	}
	else
	{
		*low = -freq_limitN;
		*high = +freq_limitP;
	}
}
WaveFront& WaveFront::ShiftedAsmPropAdd(const WaveFront& source)
{
	WaveFront& reference = *this;

	if (source.GetLambda() != reference.GetLambda()
		|| source.GetNx() != reference.GetNx()
		|| source.GetNy() != reference.GetNy()
		)
	{
		printf(">>ERROR: all parameters of reference must be equal to source's ones \n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
	if ((source.GetNormal().getX() != reference.GetNormal().getX())
		|| (source.GetNormal().getY() != reference.GetNormal().getY())
		|| (source.GetNormal().getZ() != reference.GetNormal().getZ()))
	{
		printf(">>ERROR: normalvector of reference must be equal to source's one \n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}

	WaveFront result = *this;
	result.ShiftedAsmProp(source);
	(*this).Add(result);

	return *this;
}
WaveFront& WaveFront::ShiftedAsmProp(const WaveFront& source)
{
	WaveFront& reference = *this;//alias of "this" object
	if (source.GetLambda() != reference.GetLambda()
		|| source.GetNx() != reference.GetNx()
		|| source.GetNy() != reference.GetNy()
		|| source.GetPx() != reference.GetPx()
		|| source.GetPy() != reference.GetPy()
		)
	{
		printf(">>ERROR: all parameters of reference must be equal to source's ones \n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}

	if ((source.GetNormal().getX() != reference.GetNormal().getX())
		|| (source.GetNormal().getY() != reference.GetNormal().getY())
		|| (source.GetNormal().getZ() != reference.GetNormal().getZ()))
	{
		printf(">>ERROR: normalvector of reference must be equal to source's one \n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}

	Plane refPlane(reference.GetNormal(), reference.GetOrigin());
	double z = refPlane.GetDistance(source.GetOrigin());
	vec3 reference_origin = reference.GetOrigin();
	vec3 trans = reference.GetOrigin() - source.GetOrigin();

	double x0 = trans.getX();
	double y0 = trans.getY();
	//printf("[x0 %f  y0 %f]",x0,y0);

	*this = source;//after this line, "this" object is copy of source
	this->SetOrigin(reference_origin);

	//calculate band width for band-limit
	double Sx = this->GetWidth();
	double Sy = this->GetHeight();
	double band_u = 1.0 / (2.0 * Sx);
	double band_v = 1.0 / (2.0 * Sy);//Ç±ÇÍÇ»ãCÇ™Ç∑ÇÈ(2èúéZ)

	double u_low, u_high;
	CalcBandLimitFreq(&u_high, &u_low, x0, z, band_u);
	u_high /=source.GetLambda();
	u_low /= source.GetLambda();

	double v_low, v_high;
	CalcBandLimitFreq(&v_high, &v_low, y0, z, band_v);
	v_high /= source.GetLambda();
	v_low /= source.GetLambda();

	double u_s = 1.0 / source.GetPx() / 2.0;
	double v_s = 1.0 / source.GetPy() / 2.0;

	if (+u_s < u_low || -u_s > u_high || +v_s < v_low || -v_s > v_high)
	{
		this->Clear();
		printf("<<SKIP>>");
		return *this;
	}

	WaveFront& result = *this;
	result.Embed();
	result.fft2D(-1);

	double invLambdaPow = 1.0 / result.GetLambda() / result.GetLambda();
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_max_threads())
	for (j = 0; j < result.w_ny; ++j)
	{
		for (i = 0; i < result.w_nx; ++i)
		{
			double u = result.itox(i);
			double v = result.jtoy(j);

			if (u > u_low && u < u_high && v > v_low && v < v_high)
			{
				double w = sqrt(invLambdaPow - u * u - v * v);
				double phase = 2 * PI * (w * z + u * x0 + v * y0);
				complex<double> val = complex<double>(cos(phase), sin(phase));
				result.SetPixel(i, j, result.GetPixel(i, j) * val);
			}
			else
			{
				result.SetPixel(i, j, complex<double>(0.0, 0.0));
			}
		}
	}
	//result.SaveBmp("maltipledSpectrum.bmp",INTENSITY);//debug
	//printf("[stop]");
	//system("pause");
	result.fft2D(+1);
	result *= 1.0 / result.GetN();
	result.Extract();

	return *this;
}
WaveFront& WaveFront::ShiftedAsmPropGeneralPurpose(const WaveFront& source)
{
	this->Clear();
	return this->ShiftedAsmPropAddGeneralPurpose(source);
}
WaveFront& WaveFront::ShiftedAsmPropAddGeneralPurpose(const WaveFront& source) //both object mustbe in real space 
{
	WaveFront& reference = *this;
	if (!source.GetSpace() || !reference.GetSpace())
	{
		printf(">>ERROR: this function can be used for wave front only \n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
	if (source.GetLambda() != reference.GetLambda())
	{
		printf(">>ERROR: wavelength of reference must be equal to source's one \n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
	if ((source.GetNormal().getX() != reference.GetNormal().getX())
		|| (source.GetNormal().getY() != reference.GetNormal().getY())
		|| (source.GetNormal().getZ() != reference.GetNormal().getZ()))
	{
		printf(">>ERROR: normalvector of reference must be equal to source's one \n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}

	bool islargeRefX = source.GetNx() < reference.GetNx();
	bool islargeRefY = source.GetNy() < reference.GetNy();

	int reference_divx, reference_divy, source_divx, source_divy;

	if (islargeRefX)
	{
		reference_divx = reference.GetNx() / source.GetNx();
		source_divx = 1;
	}
	else
	{
		reference_divx = 1;
		source_divx = source.GetNx() / reference.GetNx();
	}

	if (islargeRefY)
	{
		reference_divy = reference.GetNy() / source.GetNy();
		source_divy = 1;
	}
	else
	{
		reference_divy = 1;
		source_divy = source.GetNy() / reference.GetNy();
	}
	// in case of source size and refence size is equal
	if (reference_divx == 1 && reference_divy == 1 && source_divx == 1 && source_divy == 1)
	{
		printf("\n<<patern>>1\n");
		reference.ShiftedAsmPropAdd(source);
	}

	// in case of source is smaller than reference
	else if (source_divx == 1 && source_divy == 1)
	{
		printf("\n<<patern>>2\n");
		WaveFront sub_reference(reference.GetNx() / reference_divx, reference.GetNy() / reference_divy);
		sub_reference.CopyParam(reference);
		int ri, rj;
		for (rj = 0; rj < reference_divy; rj++)
		{
			for (ri = 0; ri < reference_divx; ri++)
			{
				sub_reference.SetOrigin(reference.GetOrigin()
					+ reference.GetUnitVector_alongX() * sub_reference.GetWidth() * ((double)(1 - reference_divx) / 2.0 + ri)
					+ reference.GetUnitVector_alongY() * sub_reference.GetHeight() * ((double)(1 - reference_divy) / 2.0 + rj));

				sub_reference.ShiftedAsmProp(source);
				sub_reference.SetOrigin(-sub_reference.GetOrigin());//sentence for bug fix
				reference.Add(sub_reference);
			}
		}
	}
	// in case of reference is smaller than source
	else if (reference_divx == 1 && reference_divy == 1)
	{
		printf("\n<<patern>>3\n");
		WaveFront sub_source(source.GetNx() / source_divx, source.GetNy() / source_divy);
		sub_source.CopyParam(source);

		int si, sj;
		for (sj = 0; sj < source_divy; sj++)
		{
			for (si = 0; si < source_divx; si++)
			{
				sub_source.SetOrigin(source.GetOrigin()
					+ source.GetUnitVector_alongX() * sub_source.GetWidth() * ((double)(1 - source_divx) / 2.0 + si)
					+ source.GetUnitVector_alongY() * sub_source.GetHeight() * ((double)(1 - source_divy) / 2.0 + sj));
				{
					int nx = sub_source.GetNx(), ny = sub_source.GetNy();
					int i, j;
#					pragma omp parallel for private(i, j) num_threads(omp_get_max_threads())
					for (j = 0; j < ny; j++)
					{
						for (i = 0; i < nx; i++)
						{
							sub_source.SetPixel(i, j, source.GetPixel(si * nx + i, sj * ny + j));
						}
					}
				}
				sub_source.SetOrigin(-sub_source.GetOrigin());//sentence for bug fix
				reference.ShiftedAsmPropAdd(sub_source);
			}
		}
	}
	//other case
	else
	{
		printf("\n<<patern>>4\n");
		WaveFront sub_source(source.GetNx() / source_divx, source.GetNy() / source_divy);
		sub_source.CopyParam(source);

		WaveFront sub_reference(reference.GetNx() / reference_divx, reference.GetNy() / reference_divy);
		sub_reference.CopyParam(reference);

		int si, sj;

		for (sj = 0; sj < source_divy; sj++)
		{
			for (si = 0; si < source_divx; si++)
			{
				sub_source.SetOrigin(source.GetOrigin()
					+ source.GetUnitVector_alongX() * sub_source.GetWidth() * ((double)(1 - source_divx) / 2.0 + si)
					+ source.GetUnitVector_alongY() * sub_source.GetHeight() * ((double)(1 - source_divy) / 2.0 + sj));
				//copy to sub_source
				{
					int nx = sub_source.GetNx(), ny = sub_source.GetNy();
					int i, j;
#					pragma omp parallel for private(i,j) num_threads(omp_get_max_threads())
					for (j = 0; j < ny; j++)
					{
						for (i = 0; i < nx; i++)
						{
							sub_source.SetPixel(i, j, source.GetPixel(si * nx + i, sj * ny + j));
						}
					}
				}

				int ri, rj;
				for (rj = 0; rj < reference_divy; rj++)
				{
					for (ri = 0; ri < reference_divx; ri++)
					{
						sub_reference.SetOrigin(reference.GetOrigin()
							+ reference.GetUnitVector_alongX() * sub_reference.GetWidth() * ((double)(1 - reference_divx) / 2.0 + ri)
							+ reference.GetUnitVector_alongY() * sub_reference.GetHeight() * ((double)(1 - reference_divy) / 2.0 + rj));
						
						sub_source.SetOrigin(-sub_source.GetOrigin());//sentence for bug fix
						sub_reference.ShiftedAsmProp(sub_source);
						sub_reference.SetOrigin(-sub_reference.GetOrigin());//sentence for bug fix
						reference.Add(sub_reference);
					}
				}
			}


		}
	}
	return *this;
}