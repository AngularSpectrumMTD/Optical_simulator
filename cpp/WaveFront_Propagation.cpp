#include"WaveFront.h"
using namespace std;
// function related to angular spectrum method
void WaveFront::generateFRF(double distance)
{
	double k = 1 / w_lambda / w_lambda;
	pitchtrans();
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
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
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
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
	pitchtrans();
	WaveFront h(w_nx, w_ny, w_px, w_py, w_lambda);
	h.generateFRF(R);
	double uband = 2 / (sqrt((2 * h.w_px * R) * (2 * h.w_px * R) + 1) * h.w_lambda);
	double vband = 2 / (sqrt((2 * h.w_py * R) * (2 * h.w_py * R) + 1) * h.w_lambda);
	h.bandlimit(uband, vband);// bandlimit
	pitchtrans();
	int x, y;
#pragma omp parallel for private(x, y) numw_threads(omp_get_num_threads())
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
	w_data.reset(new complex<double>[w_nx * w_ny]);
	Clear();
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
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
	w_data.reset(new complex<double>[w_nx * w_ny]);
	Clear();
	int i, j;
#pragma omp parallel for num_threads(omp_get_num_threads())
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

		/*double angle = acos(dot(w_normal, v));
		vec3 axis = normalize(cross(w_normal, v));
		mat3 ret = mat3::rotation(angle, axis);
		return ret;*/
	}
}
mat3 WaveFront::GetRotMat(const vec3& global, const vec3& local) const
{
	mat3 ret = mat3::identity();
	vec3 nglobal = normalize(global);
	vec3 nlocal = normalize(local);
	//if (global.getZ() > 0.999 || local.getZ() > 0.999)//if surface normal vector's z_elem nearly 1, matrix is invalid, so return identity matrix
	//{
	//	printf("returned");
	//	return ret;
	//}
	if (dot(nglobal,nlocal) > 0.999)//if surface normal vector's z_elem nearly 1, matrix is invalid, so return identity matrix
	{
		//printf("returned");
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
//mat3 WaveFront::GetRotMat2(vec3 v1,vec3 v2)
//{
//	v1 = normalize(v1);
//	v2 = normalize(v2);
//	double cost = dot(v1,v2);
//	v1 = cross(v2,v1);
//	double sint = length(v1);
//	double cost1 = 1.0 - cost;
//
//	double m00, m10, m20, m01, m11, m21, m02, m12, m22;
//
//	m00 = cost + v1.getX() * v1.getX() * cost1;
//	m10 = v1.getX() * v1.getY() * cost1 - v1.getZ() * sint;
//	m20 = v1.getX() * v1.getZ() * cost1 + v1.getY() * sint;
//
//	m01 = v1.getY() * v1.getX() * cost1 + v1.getZ() * sint;
//	m11 = cost + v1.getY() * v1.getY() * cost1;
//	m21 = v1.getY() * v1.getZ() * cost1 - v1.getX() * sint;
//
//	m02 = v1.getZ() * v1.getX() * cost1 - v1.getY() * sint;
//	m12 = v1.getZ() * v1.getY() * cost1 + v1.getX() * sint;
//	m22 = cost + v1.getZ() * v1.getZ() * cost1;
//
//	return mat3(
//		vec3(m00, m01, m02),
//		vec3(m10, m11, m12),
//		vec3(m20, m21, m22)
//	);
//}
//mat3 WaveFront::GetMat(vec3 v1, vec3 v2)
//{
//	return mat3(
//		vec3(v1.getX() * v2.getX(), v1.getX() * v2.getY(), v1.getX() * v2.getZ()),
//		vec3(v1.getY() * v2.getX(), v1.getY() * v2.getY(), v1.getY() * v2.getZ()),
//		vec3(v1.getZ() * v2.getX(), v1.getZ() * v2.getY(), v1.getZ() * v2.getZ())
//	);
//}
//void WaveFront::RotInFourierSpace(WaveFront& reference, Interp interp, vec3 * carrier)
//{
//	mat3 rot = reference.GetRotMat(w_normal);//get rotation matrix
//	mat3 invrot = transpose(rot);//get inverce matrix of it
//	double invL1 = 1 / w_lambda;
//	vec3 source0{ 0.0,0.0,(float)invL1 };
//	source0 = rot * source0;//carrier frequency
//	int i, j;
//#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
//	for (j = 0; j < reference.w_ny; ++j)
//	{
//		for (i = 0; i < reference.w_nx; ++i)
//		{
//			//frequency in shifted fourier space
//			vec3 ref = vec3{ static_cast<float>(reference.itox(i)),static_cast<float>(reference.jtoy(j)),static_cast<float>(invL1) };
//			//w element of ref
//			double w_ref = reference.GetWpow2(ref.getX(), ref.getY());
//			if (w_ref < 0)//if frequency is evanescent, must be ignored
//			{
//				reference.SetPixel(i, j, complex<double>(0.0, 0.0));
//				continue;
//			}
//			//frequency in reference fourier space
//			vec3 shift = ref + source0;
//			//w element of shift
//			double w_shift = reference.GetWpow2(shift.getX(), shift.getY());
//			if (w_shift < 0)//if frequency is evanescent, must be ignored
//			{
//				reference.SetPixel(i, j, complex<double>(0.0, 0.0));
//				continue;
//			}
//			shift.setZ(static_cast<float>(sqrt(w_shift)));
//			//frequency in source fourier space
//			vec3 source = invrot * shift;
//			reference.SetPixel(i, j, GetInterpolatedValue(source.getX(), source.getY(), interp));
//		}
//	}
//	//reference *= static_cast<float>(sqrt(GetEnergy() / reference.GetEnergy()));
//
//	if (carrier != nullptr)
//	{
//		*carrier = source0;
//	}
//}
//void WaveFront::TiltedAsmProp(WaveFront& reference, Interp interp, vec3* carrier)
//{
//	this->fft2D(-1);
//	reference.pitchtrans();// convert sampling interval from real to frequential
//	this->RotInFourierSpace(reference, interp, carrier);// rotational transform
//	reference.fft2D(1);
//	reference /= GetN();
//}
void WaveFront::RotInFourierSpace(WaveFront& source, Interp interp, vec3* carrier)
{
	WaveFront& reference = *this;
	mat3 rot = reference.GetRotMat(source.w_normal);//get rotation matrix
	/*printf("\nRotation Matrix =");
	DispMat(rot);*/
	mat3 invrot = transpose(rot);//get inverce matrix of it
	double invL1 = 1 / w_lambda;
	vec3 source0{ 0.0,0.0,(float)invL1 };
	source0 = rot * source0;//carrier frequency
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
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
}
void WaveFront::TiltedAsmProp(WaveFront& source, Interp interp, vec3* carrier)
{
	WaveFront& reference = *this;
	source.fft2D(-1);
	reference.pitchtrans();// convert sampling interval from real to frequential
	/*printf("\nreference vec");
	vec3 a = GetNormal();
	DispVec(a);
	vec3 b = source.GetNormal();
	printf("\nsource vec");
	DispVec(b);
	source.SaveBmp("before.bmp",INTENSITY);*/
	reference.RotInFourierSpace(source, interp, carrier);// rotational transform
	//reference.SaveBmp("after.bmp", INTENSITY);
	reference.fft2D(1);
	reference /= GetN();
}
void CalcBandLimitFreq(double* high, double* low, double xy, double z, double duv)
{
	double positive = xy + 1.0 / (2.0 * duv);
	double negative = xy - 1.0 / (2.0 * duv);
	double freq_limitP = 1.0 / sqrt(z * z / positive / positive + 1);
	double freq_limitN = 1.0 / sqrt(z * z / negative / negative + 1);

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
	*this += result;

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

	vec3 trans = reference.GetOrigin() - source.GetOrigin();
	float z = trans.getZ();
	float x0 = trans.getX();
	float y0 = trans.getY();
	vec3 reference_origin = reference.GetOrigin();

	//printf("\n----------------------------");
	//printf("\n source <max samp %f>", source.GetMaxAmplitude());
	*this = source;//after this line, "this" object is copy of source
	this->SetOrigin(reference_origin);

	//calculate band width for band-limit
	double Sx = this->GetPx() * this->GetNx();
	double Sy = this->GetPy() * this->GetNy();
	double band_u = 1.0 / Sx;
	double band_v = 1.0 / Sy;

	double u_low, u_high;
	CalcBandLimitFreq(&u_high, &u_low, x0, z, band_u);
	u_high /= this->GetLambda();
	u_low /= this->GetLambda();

	double v_low, v_high;
	CalcBandLimitFreq(&v_high, &v_low, y0, z, band_v);
	v_high /= this->GetLambda();
	v_low /= this->GetLambda();

	//printf("\nU:(%lf ,%lf) V:(%lf ,%lf)\n",u_low,u_high,v_low,v_high);

	double u_s = 1.0 / source.GetPx() / 2.0;
	double v_s = 1.0 / source.GetPy() / 2.0;

	if (+u_s < u_low || -u_s > u_high || +v_s < v_low || -v_s > v_high)
	{
		this->Clear();
		return *this;
	}

	//printf("\n this <max samp %f>", this->GetMaxAmplitude());
	WaveFront& result = *this;

	result.Embed();

	//printf("\n result before <max samp %f>", result.GetMaxAmplitude());

	result.fft2D(-1);

	double invLambdaPow = 1.0 / result.GetLambda() / result.GetLambda();
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < result.w_ny; ++j)
	{
		for (i = 0; i < result.w_nx; ++i)
		{
			double u = result.itox(i);
			double v = result.jtoy(j);

			if (u > u_low&& u < u_high && v > v_low&& v < v_high)
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

	result.fft2D(+1);
	result *= 1.0 / result.GetN();
	result.Extract();

	//printf("\n result after <max samp %f>",result.GetMaxAmplitude());

	return *this;
}
WaveFront& WaveFront::ShiftedAsmPropEx(const WaveFront& source)
{
	this->Clear();
	return this->ShiftedAsmPropAddEx(source);
}
WaveFront& WaveFront::ShiftedAsmPropAddEx(const WaveFront& source) //both object mustbe in real space 
{
	WaveFront& reference = *this;
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
		printf(">>ERROR: wavelength of reference must be equal to source's one \n");
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
	if (reference_divx * reference_divy * source_divx * source_divy == 1)
	{
		printf("<1>");
		reference.ShiftedAsmPropAdd(source);
	}

	// in case of source is smaller than reference
	else if (source_divx * source_divy == 1)
	{
		printf("<2>");
		WaveFront sub_reference(reference.GetNx() / reference_divx, reference.GetNx() / reference_divx,
			reference.GetPx(), reference.GetPy(), reference.GetLambda());
		int ri, rj;
		for (rj = 0; rj < reference_divy; rj++)
		{
			for (ri = 0; ri < reference_divx; ri++)
			{
				//printf("(%d,%d)",ri,rj);
				sub_reference.SetOrigin(reference.GetOrigin()
					+ reference.GetUnitVector_alongX() * sub_reference.GetNx() * sub_reference.GetPx() * ((double)(1 - reference_divx) / 2.0 + ri)
					+ reference.GetUnitVector_alongX() * sub_reference.GetNy() * sub_reference.GetPy() * ((double)(1 - reference_divy) / 2.0 + rj));

				sub_reference.ShiftedAsmProp(source);
				reference.Add(sub_reference);
			}
		}
	}

	// in case of reference is smaller than source
	else if (reference_divx * reference_divy == 1)
	{
		printf("<3>");
		WaveFront sub_source(source.GetNx() / source_divx, source.GetNx() / source_divx,
			source.GetPx(), source.GetPy(), source.GetLambda());

		int si, sj;
		for (sj = 0; sj < source_divy; sj++)
		{
			for (si = 0; si < source_divx; si++)
			{
				sub_source.SetOrigin(source.GetOrigin()
					+ source.GetUnitVector_alongX() * sub_source.GetNx() * sub_source.GetPx() * ((double)(1 - source_divx) / 2.0 + si)
					+ source.GetUnitVector_alongX() * sub_source.GetNy() * sub_source.GetPy() * ((double)(1 - source_divy) / 2.0 + sj));
				{
					int nx = sub_source.GetNx(), ny = sub_source.GetNy();
					int i, j;
#					pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
					for (j = 0; j < ny; j++)
					{
						for (i = 0; i < nx; i++)
						{
							sub_source.SetPixel(i, j, source.GetPixel(si * nx + i, sj * ny + j));
						}
					}
				}
				reference.ShiftedAsmPropAdd(sub_source);
			}
		}

	}
	//other case
	else
	{
		printf("<4>");
		WaveFront sub_source(source.GetNx() / source_divx, source.GetNx() / source_divx,
			source.GetPx(), source.GetPy(), source.GetLambda());

		WaveFront sub_reference(reference.GetNx() / reference_divx, reference.GetNx() / reference_divx,
			reference.GetPx(), reference.GetPy(), reference.GetLambda());

		int si, sj;

		for (sj = 0; sj < source_divy; sj++)
		{
			for (si = 0; si < source_divx; si++)
			{
				sub_source.SetOrigin(source.GetOrigin()
					+ source.GetUnitVector_alongX() * sub_source.GetNx() * sub_source.GetPx() * ((double)(1 - source_divx) / 2.0 + si)
					+ source.GetUnitVector_alongX() * sub_source.GetNy() * sub_source.GetPy() * ((double)(1 - source_divy) / 2.0 + sj));
				//copy to sub_source
				{
					int nx = sub_source.GetNx(), ny = sub_source.GetNy();
					int i, j;
#					pragma omp parallel for private(i,j) num_threads(wfl::GetNumThreads())
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
							+ reference.GetUnitVector_alongX() * sub_reference.GetNx() * sub_reference.GetPx() * ((double)(1 - reference_divx) / 2.0 + ri)
							+ reference.GetUnitVector_alongX() * sub_reference.GetNy() * sub_reference.GetPy() * ((double)(1 - reference_divy) / 2.0 + rj));
						sub_reference.ShiftedAsmProp(sub_source);

						reference.Add(sub_reference);
					}
				}
			}
		}
	}

	return *this;
}