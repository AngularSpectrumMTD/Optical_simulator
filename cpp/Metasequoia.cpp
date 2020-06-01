#include"Model.h"
using namespace std;
//spliting string
string MODEL::Split(string* str, char str1, char str2) {
	string::size_type start = str->find(str1);
	string::size_type end = str->rfind(str2);
	return str->substr(start + 1, end - start - 1);
}
//setting material
void MODEL::Material_Set() {
	int MATERIALNUM = 0;
	fscanf_s(w_fp, "%d", &MATERIALNUM);//material num
	for (int i = 0; i < MATERIALNUM; i++)
	{
		MATERIAL mtl;
		w_Material.push_back(mtl);
		w_Material[w_Material.size() - 1].w_MaterialID = i;
		fscanf_s(w_fp, "%s", w_buf, 255);//skipping {
		fscanf_s(w_fp, "%s", w_buf, 255);//sroring name
		w_str = w_buf;
		w_Material[w_Material.size() - 1].w_MaterialName = Split(&w_str, '\"', '\"');//removing ""
		fgets(w_buf, 255, w_fp);
		char* buf2;
		if ((buf2 = strstr(w_buf, "col(")) != NULL) {//storing material color
			sscanf_s(buf2, "col(%f %f %f %f)",
				&w_Material[w_Material.size() - 1].w_Color.w_r,
				&w_Material[w_Material.size() - 1].w_Color.w_g,
				&w_Material[w_Material.size() - 1].w_Color.w_b,
				&w_Material[w_Material.size() - 1].w_Color.w_a);
		}
		if ((buf2 = strstr(w_buf, "dif(")) != NULL) {//diffuse
			sscanf_s(buf2, "dif(%f)", &w_Material[w_Material.size() - 1].w_ReflectionColor.w_diffuse);
		}
		if ((buf2 = strstr(w_buf, "amb(")) != NULL) {//ambient
			sscanf_s(buf2, "amb(%f)", &w_Material[w_Material.size() - 1].w_ReflectionColor.w_ambient);
		}
		if ((buf2 = strstr(w_buf, "emi(")) != NULL) {//emission
			sscanf_s(buf2, "emi(%f)", &w_Material[w_Material.size() - 1].w_ReflectionColor.w_emission);
		}
		if ((buf2 = strstr(w_buf, "spc(")) != NULL) {//specular
			sscanf_s(buf2, "spc(%f)", &w_Material[w_Material.size() - 1].w_ReflectionColor.w_specular);
		}
		if ((buf2 = strstr(w_buf, "power(")) != NULL) {//shiness
			sscanf_s(buf2, "power(%f)", &w_Material[w_Material.size() - 1].w_Power);
		}
		if ((buf2 = strstr(w_buf, "tex(")) != NULL) {//tex name
			sscanf_s(buf2, "tex(%[^)])", w_buf, 255);
			w_Material[w_Material.size() - 1].w_texexist = true;
			w_str = w_buf;
			w_Material[w_Material.size() - 1].w_TextureName = Split(&w_str, '\"', '\"');//removing ""
			char* ptr = (char*)w_Material[w_Material.size() - 1].w_TextureName.c_str();
			w_Material[w_Material.size() - 1].w_TextureImg = w_Material[w_Material.size() - 1].w_TextureImg->Read_Bmp(ptr);//storing image data
			for (int i = 0; i < w_Material[w_Material.size() - 1].w_TextureImg->Height(); i++) {
				for (int j = 0; j < w_Material[w_Material.size() - 1].w_TextureImg->Width(); j++) {
					w_Material[w_Material.size() - 1].w_TextureImg->Write(j, i,
						w_Material[w_Material.size() - 1].w_TextureImg->Load(j, i).getX()
						, w_Material[w_Material.size() - 1].w_TextureImg->Load(j, i).getY()
						, w_Material[w_Material.size() - 1].w_TextureImg->Load(j, i).getZ());
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
	fscanf_s(w_fp, "%d", &Vertex_Max);//vertex num

	fscanf_s(w_fp, "%s", w_buf, 255);//skipping {
	for (int i = 0; i < Vertex_Max; i++) {
		fscanf_s(w_fp, "%f %f %f", &x, &y, &z);
		V.w_Coord.setX(x); V.w_Coord.setY(y); V.w_Coord.setZ(z);
		w_Object[w_Object.size() - 1].w_Vertex.push_back(V);//inserting vertex for table
	}
	fscanf_s(w_fp, "%s", w_buf, 255);//skipping }
	fscanf_s(w_fp, "%s", w_buf, 255);
}
//setting faceset infomation (because of implement status, user can use only triangular polygon currently)
void MODEL::Face_Set(int OBJECT_index) {
	int Face_Max;
	int Face;
	char* buf2;
	fscanf_s(w_fp, "%d", &Face_Max);//face num

	fscanf_s(w_fp, "%s", w_buf, 255);//skipping {
	for (int i = 0; i < Face_Max; i++) {
		fscanf_s(w_fp, "%d", &Face);
		if (Face == 3) {
			TRIANGLE tri;
			w_Object[w_Object.size() - 1].w_Triangle.push_back(tri);
			fgets(w_buf, 255, w_fp);
			if ((buf2 = strstr(w_buf, "V(")) != NULL) {
				sscanf_s(buf2, "V(%d %d %d)",
					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_Index[0],
					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_Index[1],
					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_Index[2]);
			}
			if ((buf2 = strstr(w_buf, "M(")) != NULL) {
				sscanf_s(buf2, "M(%d)", &w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_MaterialID);
			}
			if ((buf2 = strstr(w_buf, "UV(")) != NULL) {
				sscanf_s(buf2, "UV(%f %f %f %f %f %f)",
					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_uv[0].w_u,
					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_uv[0].w_v,

					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_uv[1].w_u,
					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_uv[1].w_v,

					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_uv[2].w_u,
					&w_Object[w_Object.size() - 1].w_Triangle[w_Object[w_Object.size() - 1].w_Triangle.size() - 1].w_uv[2].w_v);
			}
		}
		else if (Face == 4) {
			printf(">>ERROR: User can only use triangular polygon currently through my negligence...\n");
			printf(">>Process is terminated forcibly...\n");
			exit(0);
			QUADRILATERAL quad;
			w_Object[w_Object.size() - 1].w_Quadrilateral.push_back(quad);
			fgets(w_buf, 255, w_fp);
			if ((buf2 = strstr(w_buf, "V(")) != NULL) {
				sscanf_s(buf2, "V(%d %d %d %d)",
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_Index[0],
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_Index[1],
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_Index[2],
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_Index[3]);
			}
			if ((buf2 = strstr(w_buf, "M(")) != NULL) {
				sscanf_s(buf2, "M(%d)", &w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_MaterialID);
			}
			if ((buf2 = strstr(w_buf, "UV(")) != NULL) {
				sscanf_s(buf2, "UV(%f %f %f %f %f %f %f %f)",
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[0].w_u,
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[0].w_v,

					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[1].w_u,
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[1].w_v,

					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[2].w_u,
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[2].w_v,

					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[3].w_u,
					&w_Object[w_Object.size() - 1].w_Quadrilateral[w_Object[w_Object.size() - 1].w_Quadrilateral.size() - 1].w_uv[3].w_v);
			}
		}
		else
		{
			printf(">>ERROR: User can only use triangular polygon currently through my negligence...\n");
			printf(">>Process is terminated forcibly...\n");
			exit(0);
		}
	}
	fscanf_s(w_fp, "%s", w_buf, 255);// skipping }
	fscanf_s(w_fp, "%s", w_buf, 255);
}
//Loading mqo file
bool MODEL::MQO_Load(const char* FileName) {
	int objnum = 0;
	if (fopen_s(&w_fp, FileName, "r") != 0) { return false; }
	while (!feof(w_fp)) {
		fscanf_s(w_fp, "%s", w_buf, 255);
		if (!strcmp(w_buf, "Material")) { Material_Set(); }
		if (!strcmp(w_buf, "Object")) {
			OBJECT obj;
			w_Object.push_back(obj);
			fscanf_s(w_fp, "%s", w_buf, 255);
			w_Object[w_Object.size() - 1].w_Name = w_buf;
			while (!feof(w_fp)) {
				fscanf_s(w_fp, "%s", w_buf, 255);
				if (!strcmp(w_buf, "vertex")) { Vertex_Set(); }
				if (!strcmp(w_buf, "face")) { Face_Set(objnum); }
				if (!strcmp(w_buf, "}"))break;
			}
			objnum++;
		}
	}
	fclose(w_fp);
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
	for (int n = 0; n < (*this).w_Object.size(); ++n)
		for (int m = 0; m < (*this).w_Object[n].w_Triangle.size(); ++m)
		{
			index0 = (*this).w_Object[n].w_Triangle[m].w_Index[0];
			index1 = (*this).w_Object[n].w_Triangle[m].w_Index[1];
			index2 = (*this).w_Object[n].w_Triangle[m].w_Index[2];

			v0 = (*this).w_Object[n].w_Vertex[index0].w_Coord;
			v1 = (*this).w_Object[n].w_Vertex[index1].w_Coord;
			v2 = (*this).w_Object[n].w_Vertex[index2].w_Coord;

			(*this).w_Object[n].w_Triangle[m].w_SurfaceNV = normaloftriangle(v0, v1, v2);
		}
}

void MODEL::CalcVertexNV() {

#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
	for (int n = 0; n < (*this).w_Object.size(); ++n)
	{
		for (int m = 0; m < (*this).w_Object[n].w_Triangle.size(); ++m)
		{
			vec3 current = (*this).w_Object[n].w_Triangle[m].w_SurfaceNV;
			for (int index = 0; index < 3; ++index)
			{
				vec3 vec(0, 0, 0);
				for (int mSERCH = 0; mSERCH < (*this).w_Object[n].w_Triangle.size(); ++mSERCH)
				{
					for (int indexSERCH = 0; indexSERCH < 3; ++indexSERCH)
					{
						if ((*this).w_Object[n].w_Triangle[m].w_Index[index] == (*this).w_Object[n].w_Triangle[mSERCH].w_Index[indexSERCH])//judging whether polygon's vertex is shared or not 
						{
							if (dot(normalize(current), normalize((*this).w_Object[n].w_Triangle[mSERCH].w_SurfaceNV)) < 0)
								continue;
							else
								vec += (*this).w_Object[n].w_Triangle[mSERCH].w_SurfaceNV;
						}
					}
				}
				normalize(vec);
				(*this).w_Object[n].w_Vertex[(*this).w_Object[n].w_Triangle[m].w_Index[index]].w_VertexNV = vec;
			}
		}
	}

}
void MODEL::GenDepthList()
{
	int N = 0;

	for (int i = 0; i < (*this).w_Object.size(); i++)
		N += this->w_Object[i].w_Triangle.size();

	if (N == 0)
		return;

	w_depthlistArray.w_list.resize(N);

	int i = 0;
	//generating depth list
	for (int obj = 0; obj < (*this).w_Object.size(); ++obj)
		for (int face = 0; face < (*this).w_Object[obj].w_Triangle.size(); ++face)
		{
			w_depthlistArray.w_list[i] = depthList(obj, face, this->w_Object[obj].w_Triangle[face].w_center.getZ());
			i++;
		}
}
//sorting from far to near according to z value of center of polygon
void MODEL::SortByDepth(depthListArray& list)
{
	sort(list.w_list.begin(), list.w_list.end());
}
void MODEL::DivideByDepth(depthListArray& front, depthListArray& back, double z)
{
	depthListArray temp(w_depthlistArray);
	depthListArray tempfront, tempback;
	depthList list;
	for (int n = 0; n < temp.w_list.size(); n++)
	{
		list = temp.w_list[n];
		if ((*this).w_Object[list.objidx].w_Triangle[list.faceidx].w_center.getZ() >= z)
			tempfront.w_list.push_back(list);
		else
			tempback.w_list.push_back(list);
	}
	front = tempfront;
	back = tempback;
}

void MODEL::AccommodatePolygonInBB()
{
	vector<vec3> vec;
	vec3 v;
	for (int n = 0; n < w_Object.size(); n++)
	{
		for (int m = 0; m < w_Object[n].w_Vertex.size(); m++)
		{
			v = w_Object[n].w_Vertex[m].w_Coord;
			vec.push_back(v);
		}
	}
	BoundingBox originalbbox = GetBoundingBox(vec);

	double rx = w_bbox.Width() / originalbbox.Width(),
		ry = w_bbox.Height() / originalbbox.Height(),
		rz = w_bbox.Depth() / originalbbox.Depth();
	//calculating ratio use for accomodating

	//accomodating in desire bbox
	for (int n = 0; n < w_Object.size(); n++)
	{
		for (int m = 0; m < w_Object[n].w_Vertex.size(); m++)
		{
			v = w_Object[n].w_Vertex[m].w_Coord;
			switch (w_dir)
			{
			case DWIDTH:
			{
				v.setX(v.getX() * rx); v.setY(v.getY() * rx); v.setZ(v.getZ() * rx);
				w_Object[n].w_Vertex[m].w_Coord = v;
				break; }
			case DHEIGHT:
			{
				v.setX(v.getX() * ry); v.setY(v.getY() * ry); v.setZ(v.getZ() * ry);
				w_Object[n].w_Vertex[m].w_Coord = v;
				break; }
			case DDEPTH:
			{
				v.setX(v.getX() * rz); v.setY(v.getY() * rz); v.setZ(v.getZ() * rz);
				w_Object[n].w_Vertex[m].w_Coord = v;
				break; }
			}
		}
	}
}

vec3 MODEL::center(vec3 &p0, vec3 &p1, vec3 &p2)
{
	vector<vec3> vect = { p0, p1, p2 };
	BoundingBox bbox = GetBoundingBox(vect);
	vec3 centerBB = bbox.w_center, centerG = (p0 + p1 + p2) / 3;
	double BBlength = length(centerBB), Glength = length(centerG);
	
	return centerG;
}

void MODEL::CalcPolygonCenter() {
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
	for (int n = 0; n < (*this).w_Object.size(); ++n)
		for (int m = 0; m < (*this).w_Object[n].w_Triangle.size(); ++m)
		{
			int idx0 = (*this).w_Object[n].w_Triangle[m].w_Index[0];
			int idx1 = (*this).w_Object[n].w_Triangle[m].w_Index[1];
			int idx2 = (*this).w_Object[n].w_Triangle[m].w_Index[2];
			vec3 p0 = (*this).w_Object[n].w_Vertex[idx0].w_Coord;
			vec3 p1 = (*this).w_Object[n].w_Vertex[idx1].w_Coord;
			vec3 p2 = (*this).w_Object[n].w_Vertex[idx2].w_Coord;
			(*this).w_Object[n].w_Triangle[m].w_center = center(p0, p1, p2);
		}
}

void MODEL::CalcModelCenter() {
	BoundingBox bb;
	vector<vec3> vec;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
	for (int n = 0; n < (*this).w_Object.size(); ++n)
		for (int m = 0; m < (*this).w_Object[n].w_Vertex.size(); ++m)
		{
			vec3 v = w_Object[n].w_Vertex[m].w_Coord;
			vec.push_back(v);
		}
	bb = GetBoundingBox(vec);
	w_center = bb.w_center;
}

