#include"..\include\WaveFront.h"
#include"..\include\Model.h"
#include <numeric>
using namespace std;
void MODEL::RotInFourierSpaceForward(const WaveFront& origin, WaveFront& reference, vec3* c, Interp interp)
{
	reference.Clear();
	mat3 rot = reference.GetRotMat(origin.GetNormal());//get rotation matrix
	mat3 invrot = transpose(rot);//get inverce matrix of it
	double invL1 = 1 / origin.GetLambda();
	vec3 source0{ 0.0,0.0,static_cast<float>(invL1) };
	source0 = rot * source0;//carrier frequency
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < reference.GetNy(); ++j)
	{
		for (i = 0; i < reference.GetNx(); ++i)
		{
			//frequency in shifted fourier space
			vec3 ref = vec3{ static_cast<float>(reference.itox(i)),static_cast<float>(reference.jtoy(j)),static_cast<float>(invL1) };
			//w element of ref
			double w_ref = reference.GetWpow2(ref.getX(), ref.getY());
			if (w_ref < 0.0)//if frequency is evanescent, must be ignored
			{
				continue;
			}
			//frequency in reference fourier space
			vec3 shift = ref + source0;
			//w element of shift
			double w_shift = reference.GetWpow2(shift.getX(), shift.getY());
			if (w_shift == 0.0)//if frequency is evanescent, must be ignored
			{
				continue;
			}
			shift.setZ(static_cast<float>(sqrt(w_shift)));
			//frequency in source fourier space
			vec3 source = invrot * shift;
			//w element of source
			double w_source = reference.GetWpow2(source.getX(), source.getY());
			if (w_source < 0.0)//remove the backface element
			{
				continue;
			}

			reference.SetPixel(i, j, origin.GetInterpolatedValue(source.getX(), source.getY(), interp));
		}
	}
	if (c != nullptr)
	{
		*c = source0;
	}
}
void MODEL::RotInFourierSpaceBackward(const WaveFront& origin, WaveFront& reference, const vec3& c, Interp interp)
{
	reference.Clear();
	mat3 rot = reference.GetRotMat(origin.GetNormal());//get rotation matrix
	mat3 invrot = transpose(rot);//get inverce matrix of it
	double invL1 = 1 / origin.GetLambda();
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < reference.GetNy(); ++j)
	{
		for (i = 0; i < reference.GetNx(); ++i)
		{
			//frequency in shifted fourier space
			vec3 ref = vec3{ static_cast<float>(reference.itox(i)),static_cast<float>(reference.jtoy(j)),static_cast<float>(invL1) };
			//w element of ref
			double w_ref = reference.GetWpow2(ref.getX(), ref.getY());
			if (w_ref < 0.0)//if frequency is evanescent, must be ignored
			{
				continue;
			}

			ref.setZ(static_cast<float>(sqrt(w_ref)));
			//frequency in reference fourier space
			vec3 shift = invrot * ref;
			//w element of shift
			double w_shift = reference.GetWpow2(shift.getX(), shift.getY());
			if (w_shift < 0.0)//if frequency is evanescent, must be ignored
			{
				continue;
			}

			shift.setZ(static_cast<float>(sqrt(w_shift)));
			//frequency in source fourier space
			vec3 source = shift - c;
			reference.SetPixel(i, j, origin.GetInterpolatedValue(source.getX(), source.getY(), interp));
		}
	}
}
mat3 MODEL::RotMatFromG2L(const vec3 &global, const vec3 &local)
{
	double angle = acos(dot(global, local));
	vec3 axis = normalize(cross(global, local));
	mat3 ret = mat3::rotation(angle, axis);
	return ret;
}
mat3 MODEL::G2L()
{
	mat3 ret = mat3::identity();
	vec3 unitn = normalize(w_currentpolygon.w_surfacenormal);
	if (unitn.getZ() > 0.995)//if surface normal vector's z_elem nearly 1, matrix is invalid, so return identity matrix
	{
		return ret;
	}
	vec3 unit{ 0.0,0.0,1.0 };
	ret = RotMatFromG2L(unitn, unit);
	return ret;
}
BoundingBox MODEL::GetDiffractionRect(double targetZ)
{
	double z;
	z = w_currentpolygon.w_center.getZ();

	double Z = abs(z - targetZ);

	double bb0xp, bb0yp;//bbox param of p0
	double bb1xp, bb1yp;//bbox param of p1
	double bb2xp, bb2yp;//bbox param of p2

	double phasex = w_lambda / 2 / w_px, phasey = w_lambda / 2 / w_py;

	double diffx = Z * phasex / sqrt(1 - phasex * phasex);
	double diffy = Z * phasey / sqrt(1 - phasey * phasey);

	//calcurate x side
	bb0xp = w_currentpolygon.w_vertex[0].getX() + diffx;
	bb1xp = w_currentpolygon.w_vertex[1].getX() + diffx;
	bb2xp = w_currentpolygon.w_vertex[2].getX() + diffx;
	vector<double> xp{ bb0xp, bb1xp ,bb2xp };
	vector<double>::iterator xmaxiter = max_element(xp.begin(), xp.end());
	vector<double>::iterator xminiter = min_element(xp.begin(), xp.end());
	float xpp = *xmaxiter;
	float xnn = *xminiter - 2 * diffx;

	//calcurate y side
	bb0yp = w_currentpolygon.w_vertex[0].getY() + diffy;
	bb1yp = w_currentpolygon.w_vertex[1].getY() + diffy;
	bb2yp = w_currentpolygon.w_vertex[2].getY() + diffy;
	vector<double> yp{ bb0yp, bb1yp ,bb2yp };
	vector<double>::iterator ymaxiter = max_element(yp.begin(), yp.end());
	vector<double>::iterator yminiter = min_element(yp.begin(), yp.end());
	float ypp = *ymaxiter;
	float ynn = *yminiter - 2 * diffy;

	BoundingBox bb;
	bb.w_min = vec3{ xnn, ynn, static_cast<float>(targetZ) };
	bb.w_max = vec3{ xpp, ypp, static_cast<float>(targetZ) };
	bb.w_center = (bb.w_min + bb.w_max) / 2.0;

	return bb;
}
void MODEL::SetCurrentPolygon(depthList dpl)
{
	w_currentpolygon.w_objidx = dpl.objidx;
	w_currentpolygon.w_faceidx = dpl.faceidx;
	int n = dpl.objidx, m = dpl.faceidx;
	int index0, index1, index2;
	index0 = (*this).w_Object[n].w_Triangle[m].w_Index[0];
	index1 = (*this).w_Object[n].w_Triangle[m].w_Index[1];
	index2 = (*this).w_Object[n].w_Triangle[m].w_Index[2];

	//setting vertex
	w_currentpolygon.w_vertex[0] = (*this).w_Object[n].w_Vertex[index0].w_Coord;
	w_currentpolygon.w_vertex[1] = (*this).w_Object[n].w_Vertex[index1].w_Coord;
	w_currentpolygon.w_vertex[2] = (*this).w_Object[n].w_Vertex[index2].w_Coord;
	//setting vertexNV
	w_currentpolygon.w_vertexnormal[0] = (*this).w_Object[n].w_Vertex[index0].w_VertexNV;
	w_currentpolygon.w_vertexnormal[1] = (*this).w_Object[n].w_Vertex[index1].w_VertexNV;
	w_currentpolygon.w_vertexnormal[2] = (*this).w_Object[n].w_Vertex[index2].w_VertexNV;
	//setting surfaceNV
	w_currentpolygon.w_surfacenormal = (*this).w_Object[n].w_Triangle[m].w_SurfaceNV;
	//setting center
	w_currentpolygon.w_center = (*this).w_Object[n].w_Triangle[m].w_center;
	//setting material ID
	w_currentpolygon.w_MaterialID = (*this).w_Object[n].w_Triangle[m].w_MaterialID;
	//setting uv map
	w_currentpolygon.w_uv[0] = (*this).w_Object[n].w_Triangle[m].w_uv[0];
	w_currentpolygon.w_uv[1] = (*this).w_Object[n].w_Triangle[m].w_uv[1];
	w_currentpolygon.w_uv[2] = (*this).w_Object[n].w_Triangle[m].w_uv[2];

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
	if ((0 + ii) >= frame.GetNx() || (0 + jj) >= frame.GetNy())
	{
		return;	//No overlap
	}
	if ((sub.GetNx() + ii) < 0 || (sub.GetNy() + jj) < 0)
	{
		return;		//No overlap
	}

	int is0 = 0, is1 = sub.GetNx();
	int js0 = 0, js1 = sub.GetNy();

	if ((0 + ii) < 0)
	{
		is0 = 0 - ii;
	}
	if ((0 + jj) < 0)
	{
		js0 = 0 - jj;
	}
	if ((sub.GetNx() + ii) > frame.GetNx())
	{
		is1 = frame.GetNx() - ii;
	}
	if ((sub.GetNy() + jj) > frame.GetNy())
	{
		js1 = frame.GetNy() - jj;
	}

	if (dir)
	{
		// clipping field from frame to sub
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (int j = js0; j < js1; j++)
		{
			for (int i = is0; i < is1; i++)
			{
				sub.SetPixel(i, j, frame.GetPixel(i + ii, j + jj));
			}
		}
	}
	else
	{
		// adding field from sub to frame
#pragma omp parallel for num_threads(omp_get_num_threads())
		for (int j = js0; j < js1; j++)
		{
			for (int i = is0; i < is1; i++)
			{
				frame.SetPixel(i + ii, j + jj, frame.GetPixel(i + ii, j + jj) + sub.GetPixel(i, j));
			}
		}
	}
}
void MODEL::AddFieldToMFB(WaveFront& mfb)
{
	double z_mfb = mfb.GetOrigin().getZ();//z pos of mfb

	double z_ap = w_currentpolygon.w_center.getZ();//z pos of aperture

	double d = z_ap - z_mfb;//propagate distance

	BoundingBox bbox = GetDiffractionRect(z_mfb);//calclate diff rect

	WaveFront sfb;
	sfb.CopyParam(mfb);

	double gapx, gapy;
	gapx = fmod(w_currentpolygon.w_center.getX(), sfb.GetPx());
	gapy = fmod(w_currentpolygon.w_center.getY(), sfb.GetPy());

	//setting center of polygon as origin of sfb
	sfb.SetOrigin(vec3{ static_cast<float>(w_currentpolygon.w_center.getX() - gapx),
		static_cast<float>(w_currentpolygon.w_center.getY() - gapy), static_cast<float>(z_mfb) });

	sfb.SetNx(int(bbox.GetWidth() / sfb.GetPx()));
	sfb.SetNy(int(bbox.GetHeight() / sfb.GetPy()));

	if (sfb.GetN() <= 0 || sfb.GetN() > mfb.GetN())
	{
		printf("return->");
		return;
	}

	sfb.Init(); sfb.Clear();
	//clipping field from mfb to sfb
	ClipSubfield(sfb, mfb, true);

	//-----------------------------Silhouette
	if (w_shieldmtd == SILHOUETTE)
	{
		if (d > w_lambda)
		{
			sfb.AsmProp(d);//mfb to sfb
		}

		SilhouetteShieldingAddingField(sfb);//execute occlusion processing

		if (d > w_lambda)
		{
			sfb.AsmProp(-d);//sfb to mfb
		}
	}
	//-----------------------------

	//-----------------------------Ex
	if (w_shieldmtd == EXACT)
	{
		sfb.fft2D(-1);

		if (d > w_lambda)
		{
			sfb.AsmPropInFourierSpace(d);//mfb to sfb
		}

		ExShieldingAddingField(sfb);//execute occlusion processing

		if (d > w_lambda)
		{
			sfb.AsmPropInFourierSpace(-d);//sfb to mfb
		}

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
	int div = static_cast<int>(list.w_list.size()/20);
	int count = 5;

	for (auto& dl : list.w_list)
	{
		//setting polygon must be handle according to list
		SetCurrentPolygon(dl);
		if (PolygonIsVisible())
		{
			double polyZ = w_currentpolygon.w_center.getZ();
			//occlusion processing
			AddFieldToMFB(mfb);
			
			if (div >= 20 && n % div == 0)
			{
				printf("%d%% ", count);
				count += 5;
			}
			else
			{
				printf("(%d)", n);
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
	double originsurface = originfspace.GetWidth() * originfspace.GetHeight();
	mul(freq, rot);
	BoundingBox fspace = GetBoundingBox(freq);
	double transedsurface = fspace.GetWidth() * fspace.GetHeight();
	if (transedsurface / originsurface < 0.5)
	{
		visible = false;
	}

	return visible;
}
void MODEL::CalcCenterOfModel(depthListArray& model)
{
	vector<vec3> vec;
	for (auto& mm : model.w_list)
	{
		vec3 vec0 = w_Object[mm.objidx].w_Vertex[w_Object[mm.objidx].w_Triangle[mm.faceidx].w_Index[0]].w_Coord;
		vec3 vec1 = w_Object[mm.objidx].w_Vertex[w_Object[mm.objidx].w_Triangle[mm.faceidx].w_Index[1]].w_Coord;
		vec3 vec2 = w_Object[mm.objidx].w_Vertex[w_Object[mm.objidx].w_Triangle[mm.faceidx].w_Index[2]].w_Coord;
		vec.push_back(vec0);
		vec.push_back(vec1);
		vec.push_back(vec2);
	}
	BoundingBox bbdx = GetBoundingBox(vec);
	model.w_modelcenter = bbdx.w_center;
}
void MODEL::SetUp(const mat3& rot)
{
	AccommodatePolygonInBB();
	CalcModelCenter();
	w_center = w_bbox.w_center;
	*this += w_center;
	*this *= rot;
	CalcSurfaceNV();
	CalcVertexNV();
	GenDepthList();
	CalcPolygonCenter();
}
void MODEL::AddObjectField(WaveFront& mfb, unsigned int div, const mat3& rot, bool exact, bool back)
{
	switch (w_shieldmtd)
	{
	case SILHOUETTE:
		printf("SILHOUETTE MTD>>");
		break;
	case EXACT:
		printf("EXACT MTD>>");
	}
	w_px = mfb.GetPx();
	w_py = mfb.GetPy();
	w_lambda = mfb.GetLambda();
	vec3 originpos = mfb.GetOrigin();

	SetUp(rot);
	vector<depthListArray> model;
	depthListArray tmp = w_depthlistArray;
	double backpos = w_bbox.w_min.getZ();
	double depth = w_bbox.GetDepth();

	for (int i = 0; i < div - 1; i++)
	{
		depthListArray ddiv;
		//devide the list tmp as front, ddiv as back
		DivideByDepth(tmp, ddiv, backpos + depth / div * (i + 1));
		if (ddiv.w_list.size() != 0)//back list insert
		{
			model.push_back(ddiv);
		}
	}
	model.push_back(tmp);
	for (auto& dla : model)
	{
		CalcCenterOfModel(dla);
	}

	double d0 = (double)model[0].w_modelcenter.getZ() - (double)originpos.getZ();

	if (back)// move the pos of frame from first pos to center of first submodel
	{
		if (exact)
		{
			mfb.ExactAsmProp(d0);
		}
		else
		{
			mfb.AsmProp(d0);
		}
	}
	else
	{
		vec3 neworigin = mfb.GetOrigin();
		neworigin.setZ(model[0].w_modelcenter.getZ());
		mfb.SetOrigin(neworigin);
		mfb.Clear();
	}

	for (int n = 0; n < model.size(); n++)
	{
		AddObjectFieldPersubmodel(mfb, model[n]);
		if (n < model.size() - 1)
		{
			double dn = model[n + 1].w_modelcenter.getZ() - originpos.getZ();
			if (exact)
			{
				mfb.ExactAsmProp(dn);
			}
			else
			{
				mfb.AsmProp(dn);
			}
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
	CurrentPolygon grobal = w_currentpolygon;
	CurrentPolygon local = grobal;
	//getting rotation matrix
	mat3 rot = G2L();
	//calculate using domain in fourier space
	vector<vec3> freq;
	MarkingRectangularPointsInFourierSpace(freq);
	mul(freq, rot);
	BoundingBox fspace = GetBoundingBox(freq);
	float u = (fspace.w_max.getX() + fspace.w_min.getX()) / 2,
		v = (fspace.w_max.getY() + fspace.w_min.getY()) / 2, w;
	vec3 shiftsource{ u, v, static_cast<float>(sqrt(1 / w_lambda / w_lambda - u * u - v * v)) };
	sub(freq, shiftsource);
	fspace = GetBoundingBox(freq);
	//calculate sampling interval on fourier space
	double tfbpx, tfbpy;
	if (-fspace.w_min.getX() < fspace.w_max.getX())
	{
		tfbpx = 0.5 / fspace.w_max.getX();
	}
	else
	{
		tfbpx = -0.5 / fspace.w_min.getX();
	}
	if (-fspace.w_min.getY() < fspace.w_max.getY())
	{
		tfbpy = 0.5 / fspace.w_max.getY();
	}
	else
	{
		tfbpy = -0.5 / fspace.w_min.getY();
	}
	//transform origin to local coordinate
	vec3 pfborigin{ pfb.GetOrigin() };
	pfb.SetOrigin(vec3{ 0.0, 0.0, 0.0 });
	vector<vec3> localpoint = local.GetVertex();
	sub(localpoint, local.w_center);
	mul(localpoint, rot);
	local.SetVertex(localpoint);

	WaveFront tfb;
	tfb.SetPx(tfbpx); tfb.SetPy(tfbpy);
	tfb.SetLambda(w_lambda);
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
	//tfb.SaveBmp("�J����Z��.bmp",INTENSITY);//���������͂Ȃ�
	if (w_surface)
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
	CurrentPolygon grobal = w_currentpolygon;
	CurrentPolygon local = grobal;
	//shielding on parallel plane
	grobal.w_vertex[0].setZ(w_currentpolygon.w_center.getZ());
	grobal.w_vertex[1].setZ(w_currentpolygon.w_center.getZ());
	grobal.w_vertex[2].setZ(w_currentpolygon.w_center.getZ());
	vector<vec3> grobalpoint = grobal.GetVertex();
	sub(grobalpoint, grobal.w_center);
	grobal.SetVertex(grobalpoint);
	MultiplyAperture(pfb, grobal);
	//getting rotation matrix
	mat3 rot = G2L();
	//calculate using domain in fourier space
	vector<vec3> freq;
	MarkingRectangularPointsInFourierSpace(freq);
	mul(freq, rot);
	BoundingBox fspace = GetBoundingBox(freq);
	float u = (fspace.w_max.getX() + fspace.w_min.getX()) / 2,
		v = (fspace.w_max.getY() + fspace.w_min.getY()) / 2, w;
	vec3 shiftsource{ u, v, static_cast<float>(sqrt(1 / w_lambda / w_lambda - u * u - v * v)) };
	sub(freq, shiftsource);
	fspace = GetBoundingBox(freq);
	//calculate sampling interval on fourier space
	double tfbpx, tfbpy;
	if (-fspace.w_min.getX() < fspace.w_max.getX())
	{
		tfbpx = 0.5 / fspace.w_max.getX();
	}
	else
	{
		tfbpx = -0.5 / fspace.w_min.getX();
	}
	if (-fspace.w_min.getY() < fspace.w_max.getY())
	{
		tfbpy = 0.5 / fspace.w_max.getY();
	}
	else
	{
		tfbpy = -0.5 / fspace.w_min.getY();
	}
	//transform origin to local coordinate
	vec3 pfborigin{ pfb.GetOrigin() };
	pfb.SetOrigin(vec3{ 0.0, 0.0, 0.0 });
	vector<vec3> localpoint = local.GetVertex();
	sub(localpoint, local.w_center);
	mul(localpoint, rot);
	local.SetVertex(localpoint);

	WaveFront tfb;
	tfb.SetPx(tfbpx); tfb.SetPy(tfbpy);
	tfb.SetLambda(w_lambda);
	unsigned int tfbnx = (1 / fspace.Width() / pfbpx) * pfb.GetNx();
	unsigned int tfbny = (1 / fspace.Height() / pfbpy) * pfb.GetNy();
	tfb.SetNx(tfbnx); tfb.SetNy(tfbny);

	tfb.Init();
	//calculate polygon field
	vec3 carrier{0.0,0.0,static_cast<float>(1/w_lambda)};
	carrier = transpose(rot) * carrier;//for inverse trans, take transpose
	tfb.SetNormal(rot * pfb.GetNormal());
	tfb.Clear();//must do
	//--------------------------------------execute on polygon surface
	if (w_surface)
	{
		GeneratePolygonField(tfb, local);
	}

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
	switch (w_shader)
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
	w_currentpolygon.w_surfacenormal = normalize(w_currentpolygon.w_surfacenormal);
	double coef = dot(-w_EMV, w_currentpolygon.w_surfacenormal);
	if (coef < 0)
	{
		coef = 0;
	}
	int matindex = w_Object[w_currentpolygon.w_objidx].w_Triangle[w_currentpolygon.w_faceidx].w_MaterialID;
	float emvcol = (w_Material[matindex].w_Color.w_r + w_Material[matindex].w_Color.w_r + w_Material[matindex].w_Color.w_r) / 3.0f;
	brt = (coef + emvcol) / (1 + emvcol);

	GetCorrectedAmplitude(field, brt);

	PaintTriangle(field, local, brt);
	field.ModRandomphase();
}
double MODEL::GetCorrectedAmplitude(WaveFront& tfb, double brt)
{
	double samplingdensity = w_px * w_py / (tfb.GetPx() * tfb.GetPy());
	double theta = acos(dot(tfb.GetNormal(), vec3{ 0.0, 0.0, 1.0 }));
	if (theta >= PI / 2.0)
	{
		theta = PI / 2.0;
	}
	return sqrt(brt / samplingdensity * ((cos(theta) + w_gamma) / (1 + w_gamma)));
}
double sign(vec3 p0, vec3 p1, vec3 p2)
{
	return (p0.getX() - p2.getX()) * (p1.getY() - p2.getY()) - (p1.getX() - p2.getX()) * (p0.getY() - p2.getY());
}
bool IsInTriangle(vec3 p, CurrentPolygon& polyL)
{
	bool b0, b1, b2;
	b0 = sign(p, polyL.w_vertex[0], polyL.w_vertex[1]) < 0.0;
	b1 = sign(p, polyL.w_vertex[1], polyL.w_vertex[2]) < 0.0;
	b2 = sign(p, polyL.w_vertex[2], polyL.w_vertex[0]) < 0.0;

	return ((b0 == b1) && (b1 == b2));
}
void MODEL::SmoothShading(WaveFront& field, CurrentPolygon& polyL)
{
	int i, j;

	vec3 p0 = polyL.w_vertex[0];
	vec3 p1 = polyL.w_vertex[1];
	vec3 p2 = polyL.w_vertex[2];

	vec3 o{ 0.0, 0.0, 1.0 };

	complex<double> val{ 0.0, 0.0 };

	vec3 e1 = p1 - p0, e2 = p2 - p0;

	vec3 d{ 0.0 ,0.0 ,-1.0 };

	vec3 �� = cross(d, e2);
	float det = dot(e1, ��);

#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < field.GetNy(); j++)
	{
		for (i = 0; i < field.GetNx(); i++)
		{
			double x = field.itox(i);
			double y = field.jtoy(j);

			if (IsInTriangle(vec3{ static_cast<float>(x), static_cast<float>(y), 0 }, polyL))
			{
				vec3 p = vec3{ static_cast<float>(x), static_cast<float>(y), 0 };

				vec3 rr = p - p0;

				double aa = dot(��, rr);
				vec3 �� = cross(rr, e1);
				double bb = dot(d, ��);

				aa /= det;
				bb /= det;

				int N = polyL.w_objidx;
				int M = polyL.w_faceidx;

				int index0 = w_Object[N].w_Triangle[M].w_Index[0];
				int index1 = w_Object[N].w_Triangle[M].w_Index[1];
				int index2 = w_Object[N].w_Triangle[M].w_Index[2];

				vec3 weightedvec = (1 - aa - bb) * w_Object[N].w_Vertex[index0].w_VertexNV
					+ aa * w_Object[N].w_Vertex[index1].w_VertexNV
					+ bb * w_Object[N].w_Vertex[index2].w_VertexNV;

				weightedvec = normalize(weightedvec);

				double coef = dot(-w_EMV, weightedvec);
				if (coef < 0)
					coef = 0;
				int matindex = w_Object[w_currentpolygon.w_objidx].w_Triangle[w_currentpolygon.w_faceidx].w_MaterialID;
				float emvcol = (w_Material[matindex].w_Color.w_r + w_Material[matindex].w_Color.w_r + w_Material[matindex].w_Color.w_r) / 3.0f;
				double brt = (coef + emvcol) / (1 + emvcol);

				field.SetPixel(i, j, complex<double>(brt, 0.0));
			}
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
	int matindex = w_Object[polyL.w_objidx].w_Triangle[polyL.w_faceidx].w_MaterialID;
	if (w_Material[matindex].w_texexist == false)
	{
		return;
	}

	int i, j;

	vec3 p0 = polyL.w_vertex[0];
	vec3 p1 = polyL.w_vertex[1];
	vec3 p2 = polyL.w_vertex[2];

	int N = polyL.w_objidx;
	int M = polyL.w_faceidx;

	vec3 o{ 0.0, 0.0, 1.0 };

	vec3 e1 = p1 - p0, e2 = p2 - p0;

	vec3 d{ 0.0 ,0.0 ,-1.0 };

	vec3 �� = cross(d, e2);
	float det = dot(e1, ��);

	double x, y;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < field.GetNy(); j++)
	{
		for (i = 0; i < field.GetNx(); i++)
		{
			x = field.itox(i);
			y = field.jtoy(j);

			if (IsInTriangle(vec3{ static_cast<float>(x), static_cast<float>(y), 0 }, polyL))
			{
				vec3 p = vec3{ static_cast<float>(x), static_cast<float>(y), 0 };

				vec3 rr = p - p0;

				double aa = dot(��, rr);
				vec3 �� = cross(rr, e1);
				double bb = dot(d, ��);

				aa /= det;
				bb /= det;

				double uu = (1 - aa - bb) * w_Object[N].w_Triangle[M].w_uv[0].w_u
					+ aa * w_Object[N].w_Triangle[M].w_uv[1].w_u
					+ bb * w_Object[N].w_Triangle[M].w_uv[2].w_u;
				double vv = (1 - aa - bb) * w_Object[N].w_Triangle[M].w_uv[0].w_v
					+ aa * w_Object[N].w_Triangle[M].w_uv[1].w_v
					+ bb * w_Object[N].w_Triangle[M].w_uv[2].w_v;

				int U = static_cast<int>(uu * w_Material[matindex].w_TextureImg->Width());
				int V = static_cast<int>(vv * w_Material[matindex].w_TextureImg->Height());

				U = ceil(U % w_Material[matindex].w_TextureImg->Width());
				V = ceil(V % w_Material[matindex].w_TextureImg->Height());

				double r = w_Material[matindex].w_TextureImg->Load(U, V).getX();//r
				double g = w_Material[matindex].w_TextureImg->Load(U, V).getY();//g
				double b = w_Material[matindex].w_TextureImg->Load(U, V).getZ();//b

				double gray = (r + g + b) / 3;
				gray /= 255.99f;
				gray = sqrt(gray);

				complex<double> val = field.GetPixel(i, j);
				field.SetPixel(i, j, val * gray);
			}
		}
	}
}
void MODEL::PaintTriangle(WaveFront& tfb, CurrentPolygon& polyL, double amp)
{
	int i, j;
	//inside: zero clear
	//outside: inverse
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < tfb.GetNy(); j++)
	{
		for (i = 0; i < tfb.GetNx(); i++)
		{
			double x = tfb.itox(i);
			double y = tfb.jtoy(j);
			if (IsInTriangle(vec3{ static_cast<float>(x), static_cast<float>(y), 0 }, polyL))
			{
				tfb.SetPixel(i, j, complex<double>(amp, 0.0));
			}
			else
			{
				tfb.SetPixel(i, j, complex<double>(0.0, 0.0));
			}
		}
	}
}
void MODEL::MultiplyAperture(WaveFront& tfb, CurrentPolygon& polyL)
{
	int i, j;
	//inside: zero clear
	//outside: inverse
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (j = 0; j < tfb.GetNy(); j++)
	{
		for (i = 0; i < tfb.GetNx(); i++)
		{
			double x = tfb.itox(i);
			double y = tfb.jtoy(j);
			if (IsInTriangle(vec3{ static_cast<float>(x), static_cast<float>(y), 0 }, polyL))
			{
				tfb.SetPixel(i, j, -tfb.GetPixel(i, j));
			}
			else
			{
				tfb.SetPixel(i, j, complex<double>(0.0, 0.0));
			}
		}
	}
}
void MODEL::MarkingRectangularPointsInFourierSpace(vector<vec3>& vec)
{
	double invlambda2 = 1 / w_lambda / w_lambda;

	double maxu = 0.5 / w_px;
	double maxv = 0.5 / w_py;
	double du = (2 * maxu) / 50;
	double dv = (2 * maxv) / 50;

	int i, j;
	for (j = 0; j < 50; j++)
	{
		double v = maxv - j * dv;
		for (i = 0; i < 50; i++)
		{
			double u = maxu - i * du;
			double w = invlambda2 - u * u - v * v;
			if (w < 0)
			{
				continue;
			}
			else
			{
				vec.push_back(vec3{ static_cast<float>(u), static_cast<float>(v), static_cast<float>(sqrt(w)) });
			}
		}
	}
}
BoundingBox MODEL::GetBoundingBox(vector<vec3> &vec)
{
	BoundingBox bbox;

	vec3 min_value(vec[0]);
	vec3 max_value(vec[0]);

	for (auto& v : vec)
	{
		min_value = minPerElem(min_value, v);
		max_value = maxPerElem(max_value, v);
	}

	bbox.w_min = min_value;
	bbox.w_max = max_value;
	bbox.w_center = (bbox.w_min + bbox.w_max) / 2;

	return bbox;
}