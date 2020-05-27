#include"WaveFront.h"
#include"Model.h"
#include <numeric>

void MODEL::RotInFourierSpaceForward(const WaveFront& origin, WaveFront& reference, vec3* c, Interpol interp)
{
	reference.Clear();
	mat3 rot = reference.GetRotMat(origin.GetNormal());//get rotation matrix
	mat3 invrot = transpose(rot);//get inverce matrix of it
	double invL1 = 1 / origin.GetLambda();
	vec3 source0{ 0.0,0.0,(float)invL1 };
	source0 = rot * source0;//carrier frequency
	vec3 ref;//frequency in shifted fourier space
	vec3 shift;//frequency in reference fourier space
	vec3 source;//frequency in source fourier space
	double w_ref;//w element of ref
	double w_shift;//w element of shift
	double w_source;//w element of source
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < reference.GetNy(); ++j)
		for (i = 0; i < reference.GetNx(); ++i)
		{
			ref = vec3{ (float)reference.itox(i),(float)reference.jtoy(j),(float)invL1 };
			w_ref = reference.GetWpow2(ref.getX(), ref.getY());
			if (w_ref < 0.0)//if frequency is evanescent, must be ignored
				continue;

			shift = ref + source0;
			w_shift = reference.GetWpow2(shift.getX(), shift.getY());
			if (w_shift == 0.0)//if frequency is evanescent, must be ignored
				continue;

			shift.setZ((float)sqrt(w_shift));
			source = invrot * shift;
			w_source = reference.GetWpow2(source.getX(), source.getY());
			if (w_source < 0.0)//remove the backface element
				continue;

			reference.SetPixel(i, j, origin.GetInterpolatedValue(source.getX(), source.getY(), interp));
		}
	if (c != nullptr)
		*c = source0;
}
void MODEL::RotInFourierSpaceBackward(const WaveFront& origin, WaveFront& reference, vec3& c, Interpol interp)
{
	reference.Clear();
	mat3 rot = reference.GetRotMat(origin.GetNormal());//get rotation matrix
	mat3 invrot = transpose(rot);//get inverce matrix of it
	double invL1 = 1 / origin.GetLambda();
	vec3 ref;//frequency in shifted fourier space
	vec3 shift;//frequency in reference fourier space
	vec3 source;//frequency in source fourier space
	double w_ref;//w element of ref
	double w_shift;//w element of shift
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < reference.GetNy(); ++j)
		for (i = 0; i < reference.GetNx(); ++i)
		{
			ref = vec3{ (float)reference.itox(i),(float)reference.jtoy(j),(float)invL1 };
			w_ref = reference.GetWpow2(ref.getX(), ref.getY());
			if (w_ref < 0.0)//if frequency is evanescent, must be ignored
				continue;

			ref.setZ((float)sqrt(w_ref));
			shift = invrot * ref;
			w_shift = reference.GetWpow2(shift.getX(), shift.getY());
			if (w_shift < 0.0)//if frequency is evanescent, must be ignored
				continue;

			shift.setZ((float)sqrt(w_shift));
			source = shift - c;
			reference.SetPixel(i, j, origin.GetInterpolatedValue(source.getX(), source.getY(), interp));
		}
}
mat3 MODEL::RotMatFromG2L(vec3 global, vec3 local)
{
	double angle = acos(dot(global, local));
	vec3 axis = normalize(cross(global, local));
	mat3 ret = mat3::rotation(angle, axis);
	return ret;
}
mat3 MODEL::G2L()
{
	mat3 ret = mat3::identity();
	vec3 unitn = normalize(_currentpolygon._surfacenormal);
	if (unitn.getZ() > 0.995)//if surface normal vector's z_elem nearly 1, matrix is invalid, so return identity matrix
		return ret;
	ret = RotMatFromG2L(unitn, vec3{ 0.0,0.0,1.0 });
	return ret;
}
BoundingBox MODEL::GetDiffractionRect(double targetZ)
{
	double z;
	z = _currentpolygon._center.getZ();

	double Z = abs(z - targetZ);

	double bb0xp, bb0yp;//bbox param of p0
	double bb1xp, bb1yp;//bbox param of p1
	double bb2xp, bb2yp;//bbox param of p2

	double phasex = _lambda / 2 / _px, phasey = _lambda / 2 / _py;

	double diffx = Z * phasex / sqrt(1 - phasex * phasex);
	double diffy = Z * phasey / sqrt(1 - phasey * phasey);

	//calcurate x side
	bb0xp = _currentpolygon._vertex[0].getX() + diffx;
	bb1xp = _currentpolygon._vertex[1].getX() + diffx;
	bb2xp = _currentpolygon._vertex[2].getX() + diffx;
	vector<double> xp{ bb0xp, bb1xp ,bb2xp };
	vector<double>::iterator xmaxiter = max_element(xp.begin(), xp.end());
	vector<double>::iterator xminiter = min_element(xp.begin(), xp.end());
	float xpp = *xmaxiter;
	float xnn = *xminiter - 2 * diffx;

	//calcurate y side
	bb0yp = _currentpolygon._vertex[0].getY() + diffy;
	bb1yp = _currentpolygon._vertex[1].getY() + diffy;
	bb2yp = _currentpolygon._vertex[2].getY() + diffy;
	vector<double> yp{ bb0yp, bb1yp ,bb2yp };
	vector<double>::iterator ymaxiter = max_element(yp.begin(), yp.end());
	vector<double>::iterator yminiter = min_element(yp.begin(), yp.end());
	float ypp = *ymaxiter;
	float ynn = *yminiter - 2 * diffy;

	BoundingBox bb;
	bb._min = vec3{ xnn, ynn, (float)targetZ };
	bb._max = vec3{ xpp, ypp, (float)targetZ };
	bb._center = (bb._min + bb._max) / 2.0;

	return bb;
}
void MODEL::SetCurrentPolygon(depthList dpl)
{
	_currentpolygon._objidx = dpl.objidx;
	_currentpolygon._faceidx = dpl.faceidx;
	int n = dpl.objidx, m = dpl.faceidx;
	int index0, index1, index2;
	index0 = (*this)._Object[n]._Triangle[m]._Index[0];
	index1 = (*this)._Object[n]._Triangle[m]._Index[1];
	index2 = (*this)._Object[n]._Triangle[m]._Index[2];

	//setting vertex
	_currentpolygon._vertex[0] = (*this)._Object[n]._Vertex[index0]._Coord;
	_currentpolygon._vertex[1] = (*this)._Object[n]._Vertex[index1]._Coord;
	_currentpolygon._vertex[2] = (*this)._Object[n]._Vertex[index2]._Coord;
	//setting vertexNV
	_currentpolygon._vertexnormal[0] = (*this)._Object[n]._Vertex[index0]._VertexNV;
	_currentpolygon._vertexnormal[1] = (*this)._Object[n]._Vertex[index1]._VertexNV;
	_currentpolygon._vertexnormal[2] = (*this)._Object[n]._Vertex[index2]._VertexNV;
	//setting surfaceNV
	_currentpolygon._surfacenormal = (*this)._Object[n]._Triangle[m]._SurfaceNV;
	//setting center
	_currentpolygon._center = (*this)._Object[n]._Triangle[m]._center;
	//setting material ID
	_currentpolygon._MaterialID = (*this)._Object[n]._Triangle[m]._MaterialID;
	//setting uv map
	_currentpolygon._uv[0] = (*this)._Object[n]._Triangle[m]._uv[0];
	_currentpolygon._uv[1] = (*this)._Object[n]._Triangle[m]._uv[1];
	_currentpolygon._uv[2] = (*this)._Object[n]._Triangle[m]._uv[2];

}
void ClipSubfield(WaveFront& sub, WaveFront& frame, bool dir = true)
{
	//matching coordinate
	double x = (sub.itox(0) + sub.GetOrigin().getX()) - (frame.itox(0) + frame.GetOrigin().getX());
	double y = (sub.jtoy(0) + sub.GetOrigin().getY()) - (frame.jtoy(0) + frame.GetOrigin().getY());

	//converting from real coordinate value to integer index
	int ii = int(x / sub.GetPx() + 0.5);
	int jj = int(y / sub.GetPy() + 0.5);

	//sub's coordinate + (ii, jj) = frame's coordinate
	if ((0 + ii) >= frame.GetNx() || (0 + jj) >= frame.GetNy())	return;	//No overlap
	if ((sub.GetNx() + ii) < 0 || (sub.GetNy() + jj) < 0) return;		//No overlap

	int is0 = 0, is1 = sub.GetNx();
	int js0 = 0, js1 = sub.GetNy();

	if ((0 + ii) < 0)	is0 = 0 - ii;
	if ((0 + jj) < 0)	js0 = 0 - jj;
	if ((sub.GetNx() + ii) > frame.GetNx())		is1 = frame.GetNx() - ii;
	if ((sub.GetNy() + jj) > frame.GetNy())		js1 = frame.GetNy() - jj;

	if (dir)
	{
		// clipping field from frame to sub
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (int j = js0; j < js1; j++)
			for (int i = is0; i < is1; i++)
				sub.SetPixel(i, j, frame.GetPixel(i + ii, j + jj));
	}
	else
	{
		// adding field from sub to frame
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (int j = js0; j < js1; j++)
			for (int i = is0; i < is1; i++)
				frame.SetPixel(i + ii, j + jj, frame.GetPixel(i + ii, j + jj) + sub.GetPixel(i, j));
	}
}
void MODEL::AddFieldToMFB(WaveFront& mfb)
{
	double z_mfb = mfb.GetOrigin().getZ();//z pos of mfb

	double z_ap = _currentpolygon._center.getZ();//z pos of aperture

	double d = z_ap - z_mfb;//propagate distance

	BoundingBox bbox = GetDiffractionRect(z_mfb);//calclate diff rect

	WaveFront sfb;
	sfb.CopyParam(mfb);

	double gapx, gapy;
	gapx = fmod(_currentpolygon._center.getX(), sfb.GetPx());
	gapy = fmod(_currentpolygon._center.getY(), sfb.GetPy());

	//setting center of polygon as origin of sfb
	sfb.SetOrigin(vec3{ (float)(_currentpolygon._center.getX() - gapx),
		(float)(_currentpolygon._center.getY() - gapy), (float)z_mfb });

	sfb.SetNx(int(bbox.Width() / sfb.GetPx()));
	sfb.SetNy(int(bbox.Height() / sfb.GetPy()));

	if (sfb.GetN() <= 0 || sfb.GetN() > mfb.GetN())
		return;

	sfb.Init(); sfb.Clear();
	//clipping field from mfb to sfb
	ClipSubfield(sfb, mfb, true);

	//-----------------------------Silhouette
	if (_shieldmtd == SILHOUETTE)
	{
		if (d > _lambda)
			sfb.AsmProp(d);//mfb to sfb

		SilhouetteShieldingAddingField(sfb);//execute occlusion processing

		if (d > _lambda)
			sfb.AsmProp(-d);//sfb to mfb
	}
	//-----------------------------

	//-----------------------------Ex
	if (_shieldmtd == EXACT)
	{
		sfb.fft2D(-1);

		if (d > _lambda)
			sfb.AsmPropInFourierSpace(d);//mfb to sfb

		ExShieldingAddingField(sfb);//execute occlusion processing

		if (d > _lambda)
			sfb.AsmPropInFourierSpace(-d);//sfb to mfb

		sfb.fft2D(1);
		sfb /= sfb.GetN();//get field on real space
	}
	//-----------------------------
	//pasting field from sfb to mfb
	ClipSubfield(sfb, mfb, false);
}
void MODEL::AddObjectFieldPersubmodel(WaveFront& mfb, depthListArray& list)
{
	//double frameZ = mfb.GetOrigin().getZ();
	SortByDepth(list);
	int n = 0;
	int div = (int)(list._list.size()/20);
	int count = 5;

	for (auto& dl : list._list)
	{
		//setting polygon must be handle according to list
		SetCurrentPolygon(dl);
		if (PolygonIsVisible())
		{
			double polyZ = _currentpolygon._center.getZ();
			//occlusion processing
			AddFieldToMFB(mfb);

			if (div >= 20 && n % div == 0)
			{
				printf("%d%% ", count);
				count += 5;
			}
		}
		n++;
	}
}
bool MODEL::PolygonIsVisible()
{
	bool visible = true;

	mat3 rot = G2L();
	//calculate using domain in fourier space
	vector<vec3> freq;
	MarkingRectangularPointsInFourierSpace(freq);
	BoundingBox originfspace = GetBoundingBox(freq);
	double originsurface = originfspace.Width() * originfspace.Height();
	mul(freq, rot);
	BoundingBox fspace = GetBoundingBox(freq);
	double transedsurface = fspace.Width() * fspace.Height();
	if (transedsurface / originsurface < 0.5)
		visible = false;

	return visible;
}
void MODEL::CalcCenterOfModel(depthListArray& model)
{
	vector<vec3> vec;
	vec3 vec0, vec1, vec2;
	WaveFront a;
	for (auto& mm : model._list)
	{
		vec0 = _Object[mm.objidx]._Vertex[_Object[mm.objidx]._Triangle[mm.faceidx]._Index[0]]._Coord;
		vec1 = _Object[mm.objidx]._Vertex[_Object[mm.objidx]._Triangle[mm.faceidx]._Index[1]]._Coord;
		vec2 = _Object[mm.objidx]._Vertex[_Object[mm.objidx]._Triangle[mm.faceidx]._Index[2]]._Coord;
		vec.push_back(vec0);
		vec.push_back(vec1);
		vec.push_back(vec2);
	}
	BoundingBox bbdx = GetBoundingBox(vec);
	model._modelcenter = bbdx._center;
}
void MODEL::SetUp(mat3& rot)
{
	AccommodatePolygonInBB();
	CalcModelCenter();
	_center = _bbox._center;
	*this += _center;
	*this *= rot;
	CalcSurfaceNV();
	CalcVertexNV();
	GenDepthList();
	CalcPolygonCenter();
}
void MODEL::AddObjectField(WaveFront& mfb, unsigned int div, mat3& rot, bool exact, bool back)
{
	switch (_shieldmtd)
	{
	case SILHOUETTE:
		printf("SILHOUETTE MTD>>");
		break;
	case EXACT:
		printf("EXACT MTD>>");
	}
	_px = mfb.GetPx();
	_py = mfb.GetPy();
	_lambda = mfb.GetLambda();
	vec3 originpos = mfb.GetOrigin();

	SetUp(rot);
	vector<depthListArray> model;
	depthListArray tmp = _depthlistArray;
	double backpos = _bbox._min.getZ();
	double depth = _bbox.Depth();

	for (int i = 0; i < div - 1; i++)
	{
		depthListArray ddiv;
		//devide the list tmp as front, ddiv as back
		DivideByDepth(tmp, ddiv, backpos + depth / div * (i + 1));
		if (ddiv._list.size() != 0)//back list insert
			model.push_back(ddiv);
	}
	model.push_back(tmp);
	for (auto& dla : model)
		CalcCenterOfModel(dla);

	double d0 = (double)model[0]._modelcenter.getZ() - (double)originpos.getZ();

	if (back)// move the pos of frame from first pos to center of first submodel
	{
		if (exact)
			mfb.ExactAsmProp(d0);
		else
			mfb.AsmProp(d0);
	}
	else
	{
		vec3 neworigin = mfb.GetOrigin();
		neworigin.setZ(model[0]._modelcenter.getZ());
		mfb.SetOrigin(neworigin);
		mfb.Clear();
	}

	for (int n = 0; n < model.size(); n++)
	{
		AddObjectFieldPersubmodel(mfb, model[n]);
		if (n < model.size() - 1)
		{
			double dn = (double)model[n + 1]._modelcenter.getZ() - (double)originpos.getZ();
			if (exact)
				mfb.ExactAsmProp(dn);
			else
				mfb.AsmProp(dn);
		}
	}
	printf("\nCalculation Fnish");
}
void MODEL::ExShieldingAddingField(WaveFront& pfb)
{
	//saving sampling interval in real space
	double pfbpx = 1.0 / pfb.GetNx() / pfb.GetPx();
	double pfbpy = 1.0 / pfb.GetNy() / pfb.GetPy();
	//saving polygon infomation
	CurrentPolygon grobal = _currentpolygon;
	CurrentPolygon local = grobal;
	//getting rotation matrix
	mat3 rot = G2L();
	//calculate using domain in fourier space
	vector<vec3> freq;
	MarkingRectangularPointsInFourierSpace(freq);
	mul(freq, rot);
	BoundingBox fspace = GetBoundingBox(freq);
	float u = (fspace._max.getX() + fspace._min.getX()) / 2,
		v = (fspace._max.getY() + fspace._min.getY()) / 2, w;
	vec3 shiftsource{ u, v, (float)sqrt(1 / _lambda / _lambda - u * u - v * v) };
	sub(freq, shiftsource);
	fspace = GetBoundingBox(freq);
	//calculate sampling interval on fourier space
	double tfbpx, tfbpy;
	if (-fspace._min.getX() < fspace._max.getX())
		tfbpx = 0.5 / fspace._max.getX();
	else tfbpx = -0.5 / fspace._min.getX();
	if (-fspace._min.getY() < fspace._max.getY())
		tfbpy = 0.5 / fspace._max.getY();
	else tfbpy = -0.5 / fspace._min.getY();
	//transform origin to local coordinate
	vec3 pfborigin{ pfb.GetOrigin() };
	pfb.SetOrigin(vec3{ 0.0, 0.0, 0.0 });
	vector<vec3> localpoint = local.GetVertex();
	sub(localpoint, local._center);
	mul(localpoint, rot);
	local.SetVertex(localpoint);

	WaveFront tfb;
	tfb.SetPx(tfbpx); tfb.SetPy(tfbpy);
	tfb.SetLambda(_lambda);
	unsigned int tfbnx = (1 / fspace.Width() / pfbpx) * pfb.GetNx();
	unsigned int tfbny = (1 / fspace.Height() / pfbpy) * pfb.GetNy();
	tfb.SetNx(tfbnx); tfb.SetNy(tfbny);
	//tfb.SetNx(2*tfbnx); tfb.SetNy(2*tfbny);
    //tfb.SetNx(4 * tfbnx); tfb.SetNy(4 * tfbny);

	tfb.Init();
	tfb.pitchtrans();//handle as fourier spectrum
	//calculate polygon field
	vec3 carrier;
    tfb.SetNormal(rot * pfb.GetNormal());
	tfb.SetOrigin(pfb.GetOrigin());
	//--------------------------------------execute on polygon surface
	RotInFourierSpaceForward(pfb, tfb, &carrier, BICUBIC);//pfb to tfb
	tfb.fft2D(1);
	tfb /= tfb.GetN();
	MultiplyAperture(tfb, local);
	//tfb.SaveBmp("äJå˚èÊéZå„.bmp",INTENSITY);//Ç®Ç©ÇµÇ≠ÇÕÇ»Ç¢
	if (_surface)
	{
		WaveFront surface(tfb);
		surface.SetOrigin(tfb.GetOrigin());
		surface.Clear();
		GeneratePolygonField(surface, local);
		tfb += surface;
	}
	tfb.fft2D(-1);
	//----------------------------------------------------------------
	RotInFourierSpaceBackward(tfb, pfb, carrier, BICUBIC);//tfb to pfb

	pfb.SetOrigin(pfborigin);
}
void MODEL::SilhouetteShieldingAddingField(WaveFront& pfb)
{
	//saving sampling interval in real space
	double pfbpx = pfb.GetPx();
	double pfbpy = pfb.GetPy();
	//saving polygon infomation
	CurrentPolygon grobal = _currentpolygon;
	CurrentPolygon local = grobal;
	//shielding on parallel plane
	grobal._vertex[0].setZ(_currentpolygon._center.getZ());
	grobal._vertex[1].setZ(_currentpolygon._center.getZ());
	grobal._vertex[2].setZ(_currentpolygon._center.getZ());
	vector<vec3> grobalpoint = grobal.GetVertex();
	sub(grobalpoint, grobal._center);
	grobal.SetVertex(grobalpoint);
	MultiplyAperture(pfb, grobal);
	//getting rotation matrix
	mat3 rot = G2L();
	//calculate using domain in fourier space
	vector<vec3> freq;
	MarkingRectangularPointsInFourierSpace(freq);
	mul(freq, rot);
	BoundingBox fspace = GetBoundingBox(freq);
	float u = (fspace._max.getX() + fspace._min.getX()) / 2,
		v = (fspace._max.getY() + fspace._min.getY()) / 2, w;
	vec3 shiftsource{ u, v, (float)sqrt(1 / _lambda / _lambda - u * u - v * v) };
	sub(freq, shiftsource);
	fspace = GetBoundingBox(freq);
	//calculate sampling interval on fourier space
	double tfbpx, tfbpy;
	if (-fspace._min.getX() < fspace._max.getX())
		tfbpx = 0.5 / fspace._max.getX();
	else tfbpx = -0.5 / fspace._min.getX();
	if (-fspace._min.getY() < fspace._max.getY())
		tfbpy = 0.5 / fspace._max.getY();
	else tfbpy = -0.5 / fspace._min.getY();
	//transform origin to local coordinate
	vec3 pfborigin{ pfb.GetOrigin() };
	pfb.SetOrigin(vec3{ 0.0, 0.0, 0.0 });
	vector<vec3> localpoint = local.GetVertex();
	sub(localpoint, local._center);
	mul(localpoint, rot);
	local.SetVertex(localpoint);

	WaveFront tfb;
	tfb.SetPx(tfbpx); tfb.SetPy(tfbpy);
	tfb.SetLambda(_lambda);
	unsigned int tfbnx = (1 / fspace.Width() / pfbpx) * pfb.GetNx();
	unsigned int tfbny = (1 / fspace.Height() / pfbpy) * pfb.GetNy();
	tfb.SetNx(tfbnx); tfb.SetNy(tfbny);

	tfb.Init();
	//calculate polygon field
	vec3 carrier{0.0,0.0,(float)(1/_lambda)};
	carrier = transpose(rot) * carrier;//for inverse trans, take transpose
	tfb.SetNormal(rot * pfb.GetNormal());
	tfb.Clear();//must do
	//--------------------------------------execute on polygon surface
	if (_surface)
		GeneratePolygonField(tfb, local);

	tfb.fft2D(-1);
	//----------------------------------------------------------------
	WaveFront field(pfb);
	field.pitchtrans();
	RotInFourierSpaceBackward(tfb, field, carrier, BICUBIC);//tfb to pfb
	field.fft2D(1);
	field /= field.GetN();
	pfb += field;
	
	pfb.SetOrigin(pfborigin);
}
void MODEL::GeneratePolygonField(WaveFront& field, CurrentPolygon& local)
{
	CurrentPolygon dummyl(local);
	Shading(field, dummyl);
	Mapping(field, dummyl);
}
void MODEL::Shading(WaveFront& field, CurrentPolygon& local)
{
	switch (_shader)
	{
	case FLAT:
		FlatShading(field, local);
		break;
	case SMOOTH:
		SmoothShading(field, local);
		break;
	}
}
void MODEL::FlatShading(WaveFront& field, CurrentPolygon& local)
{
	double brt;
	_currentpolygon._surfacenormal = normalize(_currentpolygon._surfacenormal);
	double coef = dot(-_EMV, _currentpolygon._surfacenormal);
	if (coef < 0)
		coef = 0;
	int matindex = _Object[_currentpolygon._objidx]._Triangle[_currentpolygon._faceidx]._MaterialID;
	float emvcol = (_Material[matindex]._Color._r + _Material[matindex]._Color._r + _Material[matindex]._Color._r) / 3.0f;
	brt = (coef + emvcol) / (1 + emvcol);

	GetCorrectedAmplitude(field, brt);

	PaintTriangle(field, local, brt);
	field.ModRandomphase();
}
double MODEL::GetCorrectedAmplitude(WaveFront& tfb, double brt)
{
	double samplingdensity = _px * _py / (tfb.GetPx() * tfb.GetPy());
	double theta = acos(dot(tfb.GetNormal(), vec3{ 0.0, 0.0, 1.0 }));
	if (theta >= PI / 2.0)
		theta = PI / 2.0;
	return sqrt(brt / samplingdensity * ((cos(theta) + _gamma) / (1 + _gamma)));
}
double sign(vec3 p0, vec3 p1, vec3 p2)
{
	return (p0.getX() - p2.getX()) * (p1.getY() - p2.getY()) - (p1.getX() - p2.getX()) * (p0.getY() - p2.getY());
}
bool IsInTriangle(vec3 p, CurrentPolygon& polyL)
{
	bool b0, b1, b2;
	b0 = sign(p, polyL._vertex[0], polyL._vertex[1]) < 0.0;
	b1 = sign(p, polyL._vertex[1], polyL._vertex[2]) < 0.0;
	b2 = sign(p, polyL._vertex[2], polyL._vertex[0]) < 0.0;

	return ((b0 == b1) && (b1 == b2));
}
void MODEL::SmoothShading(WaveFront& field, CurrentPolygon& polyL)
{
	int i, j;

	vec3 p{ 0.0, 0.0, 0.0 };
	vec3 p0 = polyL._vertex[0];
	vec3 p1 = polyL._vertex[1];
	vec3 p2 = polyL._vertex[2];
	double aa = 0, bb = 0;

	vec3 weightedvec{ 0.0, 0.0, 0.0 };

	int N = polyL._objidx;
	int M = polyL._faceidx;
	int U = 0, V = 0;

	vec3 o{ 0.0, 0.0, 1.0 };

	complex<double> val{ 0.0, 0.0 };

	vec3 e1 = p1 - p0, e2 = p2 - p0;

	vec3 d{ 0.0 ,0.0 ,-1.0 };

	vec3 Éø = cross(d, e2);
	float det = dot(e1, Éø);

	int index0 = 0, index1 = 0, index2 = 0;

	double x, y;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < field.GetNy(); j++)
		for (i = 0; i < field.GetNx(); i++)
		{
			x = field.itox(i);
			y = field.jtoy(j);

			if (IsInTriangle(vec3{ (float)x, (float)y, 0 }, polyL))
			{
				p = vec3{ (float)x, (float)y, 0 };

				vec3 rr = p - p0;

				aa = dot(Éø, rr);
				vec3 É¿ = cross(rr, e1);
				bb = dot(d, É¿);

				aa /= det;
				bb /= det;

				index0 = _Object[N]._Triangle[M]._Index[0];
				index1 = _Object[N]._Triangle[M]._Index[1];
				index2 = _Object[N]._Triangle[M]._Index[2];

				weightedvec = (1 - aa - bb) * _Object[N]._Vertex[index0]._VertexNV
					+ aa * _Object[N]._Vertex[index1]._VertexNV
					+ bb * _Object[N]._Vertex[index2]._VertexNV;

				weightedvec = normalize(weightedvec);

				double coef = dot(-_EMV, weightedvec);
				if (coef < 0)
					coef = 0;
				int matindex = _Object[_currentpolygon._objidx]._Triangle[_currentpolygon._faceidx]._MaterialID;
				float emvcol = (_Material[matindex]._Color._r + _Material[matindex]._Color._r + _Material[matindex]._Color._r) / 3.0f;
				double brt = (coef + emvcol) / (1 + emvcol);

				field.SetPixel(i, j, complex<double>(brt, 0.0));
			}
		}

	field.ModRandomphase();
}
double alpha(vec3 p, vec3 p0, vec3 p1, vec3 p2)
{
	double ret = 0;
	ret = (p1.getY() - p2.getY()) * (p.getX() - p2.getX()) + (p2.getX() - p1.getX()) * (p.getY() - p2.getY());
	ret /= (p1.getY() - p2.getY()) * (p0.getX() - p2.getX()) + (p2.getX() - p1.getX()) * (p0.getY() - p2.getY());
	return ret;
}
double beta(vec3 p, vec3 p0, vec3 p1, vec3 p2)
{
	double ret = 0;
	ret = (p2.getY() - p0.getY()) * (p.getX() - p2.getX()) + (p0.getX() - p2.getX()) * (p.getY() - p2.getY());
	ret /= (p1.getY() - p2.getY()) * (p0.getX() - p2.getX()) + (p2.getX() - p1.getX()) * (p0.getY() - p2.getY());
	return ret;
}
void MODEL::Mapping(WaveFront& field, CurrentPolygon& polyL)
{
	int matindex = _Object[polyL._objidx]._Triangle[polyL._faceidx]._MaterialID;
	if (_Material[matindex]._texexist == false)
		return;

	double r = 0, g = 0, b = 0;
	double gray = 0;
	int i, j;

	vec3 p{ 0.0, 0.0, 0.0 };
	vec3 p0 = polyL._vertex[0];
	vec3 p1 = polyL._vertex[1];
	vec3 p2 = polyL._vertex[2];
	double aa = 0, bb = 0;
	double uu = 0, vv = 0;
	int N = polyL._objidx;
	int M = polyL._faceidx;
	int U = 0, V = 0;

	vec3 o{ 0.0, 0.0, 1.0 };

	complex<double> val{ 0.0, 0.0 };

	vec3 e1 = p1 - p0, e2 = p2 - p0;

	vec3 d{ 0.0 ,0.0 ,-1.0 };

	vec3 Éø = cross(d, e2);
	float det = dot(e1, Éø);

	double x, y;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < field.GetNy(); j++)
		for (i = 0; i < field.GetNx(); i++)
		{
			x = field.itox(i);
			y = field.jtoy(j);

			if (IsInTriangle(vec3{ (float)x, (float)y, 0 }, polyL))
			{
				p = vec3{ (float)x, (float)y, 0 };

				vec3 rr = p - p0;

				aa = dot(Éø, rr);
				vec3 É¿ = cross(rr, e1);
				bb = dot(d, É¿);

				aa /= det;
				bb /= det;

				uu = (1 - aa - bb) * _Object[N]._Triangle[M]._uv[0]._u
					+ aa * _Object[N]._Triangle[M]._uv[1]._u
					+ bb * _Object[N]._Triangle[M]._uv[2]._u;
				vv = (1 - aa - bb) * _Object[N]._Triangle[M]._uv[0]._v
					+ aa * _Object[N]._Triangle[M]._uv[1]._v
					+ bb * _Object[N]._Triangle[M]._uv[2]._v;

				U = (int)(uu * (float)_Material[matindex]._TextureImg->Width());
				V = (int)(vv * (float)_Material[matindex]._TextureImg->Height());

				U = ceil(U % _Material[matindex]._TextureImg->Width());
				V = ceil(V % _Material[matindex]._TextureImg->Height());

				r = _Material[matindex]._TextureImg->Load(U, V).getX();//r
				g = _Material[matindex]._TextureImg->Load(U, V).getY();//g
				b = _Material[matindex]._TextureImg->Load(U, V).getZ();//b

				gray = (r + g + b) / 3;
				gray /= 255.99f;
				gray = sqrt(gray);

				val = field.GetPixel(i, j);
				field.SetPixel(i, j, val * gray);
			}
		}
}
void MODEL::PaintTriangle(WaveFront& tfb, CurrentPolygon& polyL, double amp)
{
	int i, j;
	double x, y;
	//inside: zero clear
	//outside: inverse
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < tfb.GetNy(); j++)
		for (i = 0; i < tfb.GetNx(); i++)
		{
			x = tfb.itox(i);
			y = tfb.jtoy(j);
			if (IsInTriangle(vec3{ (float)x, (float)y, 0 }, polyL))
				tfb.SetPixel(i, j, complex<double>(amp, 0.0));
			else
				tfb.SetPixel(i, j, complex<double>(0.0, 0.0));
		}
}
void MODEL::MultiplyAperture(WaveFront& tfb, CurrentPolygon& polyL)
{
	int i, j;
	double x, y;
	//inside: zero clear
	//outside: inverse
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < tfb.GetNy(); j++)
		for (i = 0; i < tfb.GetNx(); i++)
		{
			x = tfb.itox(i);
			y = tfb.jtoy(j);
			if (IsInTriangle(vec3{ (float)x,(float)y,0 }, polyL))
				tfb.SetPixel(i, j, -tfb.GetPixel(i, j));
			else
				tfb.SetPixel(i, j, complex<double>(0.0, 0.0));
		}
}
void MODEL::MarkingRectangularPointsInFourierSpace(vector<vec3>& vec)
{
	double invlambda2 = 1 / _lambda / _lambda;

	double maxu = 0.5 / _px;
	double maxv = 0.5 / _py;
	double du = (2 * maxu) / 50;
	double dv = (2 * maxv) / 50;

	double u, v, w;
	int i, j;
	for (j = 0; j < 50; j++)
	{
		v = maxv - j * dv;
		for (i = 0; i < 50; i++)
		{
			u = maxu - i * du;
			w = invlambda2 - u * u - v * v;
			if (w < 0)
				continue;
			else
			{
				vec.push_back(vec3{ (float)u, (float)v, (float)sqrt(w) });
			}
		}
	}
}
BoundingBox MODEL::GetBoundingBox(vector<vec3> vec)
{
	vector<float> x;
	vector<float> y;
	vector<float> z;

	for (int i = 0; i < vec.size(); i++)
	{
		x.push_back(vec[i].getX());
		y.push_back(vec[i].getY());
		z.push_back(vec[i].getZ());
	}
	vector<float>::iterator xmaxiter = max_element(x.begin(), x.end());
	vector<float>::iterator xminiter = min_element(x.begin(), x.end());
	
	vector<float>::iterator ymaxiter = max_element(y.begin(), y.end());
	vector<float>::iterator yminiter = min_element(y.begin(), y.end());
	
	vector<float>::iterator zmaxiter = max_element(z.begin(), z.end());
	vector<float>::iterator zminiter = min_element(z.begin(), z.end());
	
	BoundingBox bbox;
	
	bbox._min = vec3{ *xminiter, *yminiter, *zminiter };
	bbox._max = vec3{ *xmaxiter, *ymaxiter, *zmaxiter };
	bbox._center = (bbox._min + bbox._max) / 2;

	return bbox;
}