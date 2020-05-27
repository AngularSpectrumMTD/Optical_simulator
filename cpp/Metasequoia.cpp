#include"Model.h"

//spliting string
string MODEL::Split(string* str, char str1, char str2) {
	string::size_type start = str->find(str1);
	string::size_type end = str->rfind(str2);
	return str->substr(start + 1, end - start - 1);
}
//setting material
void MODEL::Material_Set() {
	int MATERIALNUM = 0;
	fscanf_s(_fp, "%d", &MATERIALNUM);//material num
	for (int i = 0; i < MATERIALNUM; i++)
	{
		MATERIAL mtl;
		_Material.push_back(mtl);
		_Material[_Material.size() - 1]._MaterialID = i;
		fscanf_s(_fp, "%s", _buf, 255);//skipping {
		fscanf_s(_fp, "%s", _buf, 255);//sroring name
		_str = _buf;
		_Material[_Material.size() - 1]._MaterialName = Split(&_str, '\"', '\"');//removing ""
		fgets(_buf, 255, _fp);
		char* buf2;
		if ((buf2 = strstr(_buf, "col(")) != NULL) {//storing material color
			sscanf_s(buf2, "col(%f %f %f %f)",
				&_Material[_Material.size() - 1]._Color._r,
				&_Material[_Material.size() - 1]._Color._g,
				&_Material[_Material.size() - 1]._Color._b,
				&_Material[_Material.size() - 1]._Color._a);
		}
		if ((buf2 = strstr(_buf, "dif(")) != NULL) {//diffuse
			sscanf_s(buf2, "dif(%f)", &_Material[_Material.size() - 1]._ReflectionColor._diffuse);
		}
		if ((buf2 = strstr(_buf, "amb(")) != NULL) {//ambient
			sscanf_s(buf2, "amb(%f)", &_Material[_Material.size() - 1]._ReflectionColor._ambient);
		}
		if ((buf2 = strstr(_buf, "emi(")) != NULL) {//emission
			sscanf_s(buf2, "emi(%f)", &_Material[_Material.size() - 1]._ReflectionColor._emission);
		}
		if ((buf2 = strstr(_buf, "spc(")) != NULL) {//specular
			sscanf_s(buf2, "spc(%f)", &_Material[_Material.size() - 1]._ReflectionColor._specular);
		}
		if ((buf2 = strstr(_buf, "power(")) != NULL) {//shiness
			sscanf_s(buf2, "power(%f)", &_Material[_Material.size() - 1]._Power);
		}
		if ((buf2 = strstr(_buf, "tex(")) != NULL) {//tex name
			sscanf_s(buf2, "tex(%[^)])", _buf, 255);
			_Material[_Material.size() - 1]._texexist = true;
			_str = _buf;
			_Material[_Material.size() - 1]._TextureName = Split(&_str, '\"', '\"');//removing ""
			char* ptr = (char*)_Material[_Material.size() - 1]._TextureName.c_str();
			_Material[_Material.size() - 1]._TextureImg = _Material[_Material.size() - 1]._TextureImg->Read_Bmp(ptr);//storing image data
			for (int i = 0; i < _Material[_Material.size() - 1]._TextureImg->Height(); i++) {
				for (int j = 0; j < _Material[_Material.size() - 1]._TextureImg->Width(); j++) {
					_Material[_Material.size() - 1]._TextureImg->Write(j, i,
						_Material[_Material.size() - 1]._TextureImg->Load(j, i).getX()
						, _Material[_Material.size() - 1]._TextureImg->Load(j, i).getY()
						, _Material[_Material.size() - 1]._TextureImg->Load(j, i).getZ());
				}
			}
		}
	}
}
//setting vertex coordinate
void MODEL::Vertex_Set() {
	int Vertex_Max;
	VERTEX V;
	float x, y, z;
	fscanf_s(_fp, "%d", &Vertex_Max);//vertex num

	fscanf_s(_fp, "%s", _buf, 255);//skipping {
	for (int i = 0; i < Vertex_Max; i++) {
		fscanf_s(_fp, "%f %f %f", &x, &y, &z);
		V._Coord.setX(x); V._Coord.setY(y); V._Coord.setZ(z);
		_Object[_Object.size() - 1]._Vertex.push_back(V);//inserting vertex for table
	}
	fscanf_s(_fp, "%s", _buf, 255);//skipping }
	fscanf_s(_fp, "%s", _buf, 255);
}
//setting faceset infomation (because of implement status, user can use only triangular polygon currently)
void MODEL::Face_Set(int OBJECT_index) {
	int Face_Max;
	int Face;
	char* buf2;
	fscanf_s(_fp, "%d", &Face_Max);//face num

	fscanf_s(_fp, "%s", _buf, 255);//skipping {
	for (int i = 0; i < Face_Max; i++) {
		fscanf_s(_fp, "%d", &Face);
		if (Face == 3) {
			TRIANGLE tri;
			_Object[_Object.size() - 1]._Triangle.push_back(tri);
			fgets(_buf, 255, _fp);
			if ((buf2 = strstr(_buf, "V(")) != NULL) {
				sscanf_s(buf2, "V(%d %d %d)",
					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._Index[0],
					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._Index[1],
					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._Index[2]);
			}
			if ((buf2 = strstr(_buf, "M(")) != NULL) {
				sscanf_s(buf2, "M(%d)", &_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._MaterialID);
			}
			if ((buf2 = strstr(_buf, "UV(")) != NULL) {
				sscanf_s(buf2, "UV(%f %f %f %f %f %f)",
					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._uv[0]._u,
					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._uv[0]._v,

					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._uv[1]._u,
					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._uv[1]._v,

					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._uv[2]._u,
					&_Object[_Object.size() - 1]._Triangle[_Object[_Object.size() - 1]._Triangle.size() - 1]._uv[2]._v);
			}
		}
		else if (Face == 4) {
			printf(">>ERROR: User can only use triangular polygon currently through my negligence...\n");
			printf(">>Process is terminated forcibly...\n");
			exit(0);
			QUADRILATERAL quad;
			_Object[_Object.size() - 1]._Quadrilateral.push_back(quad);
			fgets(_buf, 255, _fp);
			if ((buf2 = strstr(_buf, "V(")) != NULL) {
				sscanf_s(buf2, "V(%d %d %d %d)",
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._Index[0],
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._Index[1],
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._Index[2],
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._Index[3]);
			}
			if ((buf2 = strstr(_buf, "M(")) != NULL) {
				sscanf_s(buf2, "M(%d)", &_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._MaterialID);
			}
			if ((buf2 = strstr(_buf, "UV(")) != NULL) {
				sscanf_s(buf2, "UV(%f %f %f %f %f %f %f %f)",
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[0]._u,
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[0]._v,

					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[1]._u,
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[1]._v,

					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[2]._u,
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[2]._v,

					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[3]._u,
					&_Object[_Object.size() - 1]._Quadrilateral[_Object[_Object.size() - 1]._Quadrilateral.size() - 1]._uv[3]._v);
			}
		}
		else
		{
			printf(">>ERROR: User can only use triangular polygon currently through my negligence...\n");
			printf(">>Process is terminated forcibly...\n");
			exit(0);
		}
	}
	fscanf_s(_fp, "%s", _buf, 255);// skipping }
	fscanf_s(_fp, "%s", _buf, 255);
}
//Loading mqo file
bool MODEL::MQO_Load(const char* FileName) {
	int objnum = 0;
	if (fopen_s(&_fp, FileName, "r") != 0) { return false; }
	while (!feof(_fp)) {
		fscanf_s(_fp, "%s", _buf, 255);
		if (!strcmp(_buf, "Material")) { Material_Set(); }
		if (!strcmp(_buf, "Object")) {
			OBJECT obj;
			_Object.push_back(obj);
			fscanf_s(_fp, "%s", _buf, 255);
			_Object[_Object.size() - 1]._Name = _buf;
			while (!feof(_fp)) {
				fscanf_s(_fp, "%s", _buf, 255);
				if (!strcmp(_buf, "vertex")) { Vertex_Set(); }
				if (!strcmp(_buf, "face")) { Face_Set(objnum); }
				if (!strcmp(_buf, "}"))break;
			}
			objnum++;
		}
	}
	fclose(_fp);
	return true;
}

vec3 normaloftriangle(vec3 v0, vec3 v1, vec3 v2)
{
	vec3 v10 = v1 - v0;
	vec3 v20 = v2 - v0;
	vec3 ret = normalize(cross(v20, v10));
	return ret;
}

void MODEL::CalcSurfaceNV() {
	vec3 v0, v1, v2;
	WaveFront a;
	int index0, index1, index2;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
	for (int n = 0; n < (*this)._Object.size(); ++n)
		for (int m = 0; m < (*this)._Object[n]._Triangle.size(); ++m)
		{
			index0 = (*this)._Object[n]._Triangle[m]._Index[0];
			index1 = (*this)._Object[n]._Triangle[m]._Index[1];
			index2 = (*this)._Object[n]._Triangle[m]._Index[2];

			v0 = (*this)._Object[n]._Vertex[index0]._Coord;
			v1 = (*this)._Object[n]._Vertex[index1]._Coord;
			v2 = (*this)._Object[n]._Vertex[index2]._Coord;

			(*this)._Object[n]._Triangle[m]._SurfaceNV = normaloftriangle(v0, v1, v2);
		}
}

void MODEL::CalcVertexNV() {

#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
	for (int n = 0; n < (*this)._Object.size(); ++n)
	{
		for (int m = 0; m < (*this)._Object[n]._Triangle.size(); ++m)
		{
			vec3 current = (*this)._Object[n]._Triangle[m]._SurfaceNV;
			for (int index = 0; index < 3; ++index)
			{
				vec3 vec(0, 0, 0);
				for (int mSERCH = 0; mSERCH < (*this)._Object[n]._Triangle.size(); ++mSERCH)
				{
					for (int indexSERCH = 0; indexSERCH < 3; ++indexSERCH)
					{
						if ((*this)._Object[n]._Triangle[m]._Index[index] == (*this)._Object[n]._Triangle[mSERCH]._Index[indexSERCH])//judging whether polygon's vertex is shared or not 
						{
							if (dot(normalize(current), normalize((*this)._Object[n]._Triangle[mSERCH]._SurfaceNV)) < 0)
								continue;
							else
								vec += (*this)._Object[n]._Triangle[mSERCH]._SurfaceNV;
						}
					}
				}
				normalize(vec);
				(*this)._Object[n]._Vertex[(*this)._Object[n]._Triangle[m]._Index[index]]._VertexNV = vec;
			}
		}
	}

}
void MODEL::GenDepthList()
{
	int N = 0;

	for (int i = 0; i < (*this)._Object.size(); i++)
		N += this->_Object[i]._Triangle.size();

	if (N == 0)
		return;

	_depthlistArray._list.resize(N);

	int i = 0;
	//generating depth list
	for (int obj = 0; obj < (*this)._Object.size(); ++obj)
		for (int face = 0; face < (*this)._Object[obj]._Triangle.size(); ++face)
		{
			_depthlistArray._list[i] = depthList(obj, face, this->_Object[obj]._Triangle[face]._center.getZ());
			i++;
		}
}
//sorting from far to near according to z value of center of polygon
void MODEL::SortByDepth(depthListArray& list)
{
	sort(list._list.begin(), list._list.end());
}
void MODEL::DivideByDepth(depthListArray& front, depthListArray& back, double z)
{
	depthListArray temp(_depthlistArray);
	depthListArray tempfront, tempback;
	depthList list;
	for (int n = 0; n < temp._list.size(); n++)
	{
		list = temp._list[n];
		if ((*this)._Object[list.objidx]._Triangle[list.faceidx]._center.getZ() >= z)
			tempfront._list.push_back(list);
		else
			tempback._list.push_back(list);
	}
	front = tempfront;
	back = tempback;
}

void MODEL::AccommodatePolygonInBB()
{
	vector<vec3> vec;
	vec3 v;
	for (int n = 0; n < _Object.size(); n++)
	{
		for (int m = 0; m < _Object[n]._Vertex.size(); m++)
		{
			v = _Object[n]._Vertex[m]._Coord;
			vec.push_back(v);
		}
	}
	BoundingBox originalbbox = GetBoundingBox(vec);

	double rx = _bbox.Width() / originalbbox.Width(),
		ry = _bbox.Height() / originalbbox.Height(),
		rz = _bbox.Depth() / originalbbox.Depth();
	//calculating ratio use for accomodating

	//accomodating in desire bbox
	for (int n = 0; n < _Object.size(); n++)
	{
		for (int m = 0; m < _Object[n]._Vertex.size(); m++)
		{
			v = _Object[n]._Vertex[m]._Coord;
			switch (_dir)
			{
			case DWIDTH:
			{
				v.setX(v.getX() * rx); v.setY(v.getY() * rx); v.setZ(v.getZ() * rx);
				_Object[n]._Vertex[m]._Coord = v;
				break; }
			case DHEIGHT:
			{
				v.setX(v.getX() * ry); v.setY(v.getY() * ry); v.setZ(v.getZ() * ry);
				_Object[n]._Vertex[m]._Coord = v;
				break; }
			case DDEPTH:
			{
				v.setX(v.getX() * rz); v.setY(v.getY() * rz); v.setZ(v.getZ() * rz);
				_Object[n]._Vertex[m]._Coord = v;
				break; }
			}
		}
	}
}

vec3 MODEL::center(vec3 p0, vec3 p1, vec3 p2)
{
	vector<vec3> vect = { p0, p1, p2 };
	BoundingBox bbox = GetBoundingBox(vect);
	vec3 centerBB = bbox._center, centerG = (p0 + p1 + p2) / 3;
	double BBlength = length(centerBB), Glength = length(centerG);
	
	return centerG;
}

void MODEL::CalcPolygonCenter() {
	vec3 p0, p1, p2;
	int idx0, idx1, idx2;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
	for (int n = 0; n < (*this)._Object.size(); ++n)
		for (int m = 0; m < (*this)._Object[n]._Triangle.size(); ++m)
		{
			idx0 = (*this)._Object[n]._Triangle[m]._Index[0];
			idx1 = (*this)._Object[n]._Triangle[m]._Index[1];
			idx2 = (*this)._Object[n]._Triangle[m]._Index[2];
			p0 = (*this)._Object[n]._Vertex[idx0]._Coord;
			p1 = (*this)._Object[n]._Vertex[idx1]._Coord;
			p2 = (*this)._Object[n]._Vertex[idx2]._Coord;
			(*this)._Object[n]._Triangle[m]._center = center(p0, p1, p2);
		}
}

void MODEL::CalcModelCenter() {
	BoundingBox bb;
	vector<vec3> vec;
	vec3 v;
	int idx0, idx1, idx2;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
	for (int n = 0; n < (*this)._Object.size(); ++n)
		for (int m = 0; m < (*this)._Object[n]._Vertex.size(); ++m)
		{
			v = _Object[n]._Vertex[m]._Coord;
			vec.push_back(v);
		}
	bb = GetBoundingBox(vec);
	_center = bb._center;
}

