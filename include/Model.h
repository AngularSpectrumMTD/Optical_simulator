#ifndef MODELH
#define MODELH
#include<vector>
#include<string>
#include<algorithm>
#include<string>
#include"WaveFront.h"

#ifndef SHADE
	enum Shader
	{
		FLAT,
		SMOOTH
	};
#endif
#ifndef DIRECT
	enum Direction
	{
		DWIDTH,
		DHEIGHT,
		DDEPTH
	};
#endif
#ifndef METHOD
	enum Shield
	{
		SILHOUETTE,
		EXACT
	};
#endif
#ifndef REF4
#define REF4
	struct Reflection4 {
		float w_diffuse;
		float w_ambient;
		float w_emission;
		float w_specular;
	};
#endif
#ifndef COL4
#define COL4
	struct Color4 {
		float w_r;
		float w_g;
		float w_b;
		float w_a;
	};
#endif
#ifndef UVC
#define UVC
	struct UV {
		float w_u;
		float w_v;
	};
#endif
#ifndef TRI
#define TRI
	//triangular face
	struct Triangle {
		int w_MaterialID;//material No.
		int w_Index[3];//index
		UV w_uv[3];//uv info.
		vec3 w_SurfaceNV;//surface normal
		float w_bbox[2][3];
		vec3 w_center;
	};
#endif
#ifndef QUAD
#define QUAD
	//quadrilateral face
	struct Quadrilateral {
		int w_MaterialID;
		int w_Index[4];
		UV w_uv[4];
		vec3 w_SurfaceNV;
	};
#endif
#ifndef VTX
#define VTX
	//vertex
	struct Vertex {
		vec3 w_Coord;//coordinate
		vec3 w_VertexNV;//vertex normal
	};
#endif
#ifndef OBJ
#define OBJ
	//object
	struct Object {
		std::string w_Name;//object name
		std::vector<Vertex> w_Vertex;//vertex data
		std::vector<Triangle> w_Triangle;
		std::vector<Quadrilateral> w_Quadrilateral;
	};
#endif
#ifndef MATE
#define MATE
	//material
	struct Material {
		int w_MaterialID;//ID
		std::string w_MaterialName;//material name
		Color4 w_Color;//color
		Reflection4 w_ReflectionColor;//color of refrection
		float w_Power;//shiness
		std::string w_TextureName;//name of texture
		Image* w_TextureImg;//image of texture
		bool w_texexist = false;
	};
#endif
#ifndef DPTL
#define DPTL
	class depthList
	{
	public:
		int objidx;
		int faceidx;
		double depth = 0;	//z value of coordinate of center of current polygon

		depthList() :objidx(0), faceidx(0), depth(0.0)
		{};

		depthList(int obj, int face, double dep) : objidx(obj), faceidx(face), depth(dep)
		{};

		depthList(const depthList& dl) : objidx(dl.objidx), faceidx(dl.faceidx), depth(dl.depth)
		{};

		bool operator<(const depthList& dl) const
		{
			return this->depth < dl.depth;
		}
	};
#endif
#ifndef DPTLA
#define DPTLA
	class depthListArray
	{
	public:
		std::vector<depthList> w_list;
		vec3 w_modelcenter;
	};
#endif
#ifndef BDB
#define BDB
	class BoundingBox
	{
	public:
		vec3 w_min;
		vec3 w_max;
		vec3 w_center;
		BoundingBox() :w_min(vec3{ 0.0,0.0,0.0 }), 
			w_max(vec3{ 0.0,0.0,0.0 }), 
			w_center(vec3{ 0.0,0.0,0.0 }) {};
		BoundingBox(vec3 min, vec3 max) :w_min(min), w_max(max), w_center( (max + min)/2 ) {};
		BoundingBox(double w, double h, double d, vec3 center):w_center(center)
		{
			w_min = vec3{ static_cast<float>(center.getX() - w / 2), static_cast<float>(center.getY() - h / 2), static_cast<float>(center.getZ() - d / 2) };
			w_max = vec3{ static_cast<float>(center.getX() + w / 2), static_cast<float>(center.getY() + h / 2), static_cast<float>(center.getZ() + d / 2) };
		}
		double GetWidth() const { return w_max.getX() - w_min.getX(); }
		double GetHeight() const { return w_max.getY() - w_min.getY(); }
		double GetDepth() const { return w_max.getZ() - w_min.getZ(); }
	};
#endif
#ifndef CUP
#define CUP
	class CurrentPolygon
	{
	public:
		int w_objidx = 0;
		int w_faceidx = 0;
		vec3 w_vertex[3];
		vec3 w_vertexnormal[3];
		vec3 w_surfacenormal = vec3{0.0,0.0,0.0};
		vec3 w_center = vec3{ 0.0,0.0,0.0 };
		int w_MaterialID = 0;
		UV w_uv[3];

		vec3 w_diffractcenter = vec3{ 0.0,0.0,0.0 };

		CurrentPolygon() {};
		~CurrentPolygon() {};

		void SetVertex(const std::vector<vec3> v)
		{
			if (v.size() != 3)
			{
				printf(">>ERROR: SetVertex: input vector size is not suit!!\n");
				printf(">>Process is terminated forcibly...\n");
				exit(0);
			}
			w_vertex[0] = v[0];
			w_vertex[1] = v[1];
			w_vertex[2] = v[2];
		};

		std::vector<vec3> GetVertex()
		{
			std::vector<vec3> v(3);
			
			for (int i = 0; i < 3; i++)
			{
				v[i] = w_vertex[i];
			}

			return v;
		};
	};
#endif
#ifndef RAY
#define RAY
	class Ray {//to handle line object 
	public:

		Ray() {}

		Ray(const vec3& o, const vec3& dir) : w_origin(o), w_direction(dir) {}

		const vec3& origin() const { return w_origin; }
		const vec3& direction() const { return w_direction; }
		vec3 at(float t) const { return w_origin + t * w_direction; }

	private:
		vec3 w_origin;//start point
		vec3 w_direction;//directional vector
	};
#endif
#ifndef MODELL
#define MODELL
	//model
	class Model {
	protected:
		double w_px = 1e-6;//sampling interval along to x axis of main frame buffer
		double w_py = 1e-6;//along to y axis
		double w_lambda = 633e-9;//wavelength

		bool w_calced_surfaceNV = false;

		FILE* w_fp;//file pointer to handle mqo file
		char w_buf[255];
		std::string w_str;
		depthListArray w_depthlistArray;//array of depthlist(contain polygon index and z value of center of it)
		std::vector<Material> w_Material;
		std::vector<Object> w_Object;
		vec3 w_EMV = vec3{0.0,0.0,-1.0};//vector of ambient light

		Shader w_shader;//shading method

		CurrentPolygon w_currentpolygon;//contaions almost information to calculate polygon field

		Direction w_dir = DWIDTH;//direction means accomodate to bounding box

		vec3 w_center = vec3{0.0,0.0,0.0};//model's center

		bool w_surface = true;//whether to calculate surface function

		vec3 w_subcenter = vec3{ 0.0,0.0,0.0 };//center of submodel

		double w_gamma = 2.2;

		Shield w_shieldmtd = EXACT;//shielding method

		std::unique_ptr<WaveFront> w_randomfield;

		//Metasequoia
		void Vertex_Set();
		void Material_Set();
		void Face_Set(int Object_num);
		std::string Split(std::string* str, char str1, char str2);
		BoundingBox w_bbox;

		double w_time_fft = 0;
		double w_time_interpol = 0;
		double w_time_other = 0;

		LARGE_INTEGER w_freq;
		LARGE_INTEGER w_start;
		LARGE_INTEGER w_end;
		double getdeltatime() { return (w_end.QuadPart - w_start.QuadPart) * 1000.0 / w_freq.QuadPart; }
		void dispTotalTime() { printf("\nFFT: %lf [ms] INTERPOL: %lf [ms] OTHER: %lf [ms]\n", w_time_fft, w_time_interpol, w_time_other); }
	public:
		Model() {
		}
		Model(const Model& model)
		{
			w_fp = model.w_fp;
			strcpy(w_buf, model.w_buf);
			w_str = model.w_str;
			w_depthlistArray = model.w_depthlistArray;
			w_Material = model.w_Material;
			w_Object = model.w_Object;
			w_EMV = model.w_EMV;
			w_bbox = model.w_bbox;
			w_surface = model.w_surface;
			w_px = model.w_px;
			w_py = model.w_py;
			*w_randomfield = *model.w_randomfield;
		}
		Model(const char* FileName, vec3 emv, Shader shade, BoundingBox bb, Direction dir, bool surface): 
			w_shader(shade),w_dir(dir),w_surface(surface) {
			if(MQO_Load(FileName))
              w_bbox = bb;
			else
			{
				printf("Indicated file is not exist under the this directly....");
				system("pause");
			}
			w_randomfield = std::make_unique<WaveFront>(512, 512);
			w_randomfield->AllSet(1.0);
			w_randomfield->ModRandomphase();
		}
		~Model()
		{
		}
		double GetPx() const { return w_px; }
		double GetPy() const { return w_py; }
		double GetLambda() const { return w_lambda; }

		void SetShieldMethod(const Shield mtd) { w_shieldmtd = mtd; }

		FILE* GetFilePointer() const { return w_fp; }
		char GetBuffer() const { return *w_buf; }
		std::string GetString() const { return w_str; }
		depthListArray  GetDepthlist() const { return w_depthlistArray; }
		std::vector<Material> GetMaterial()  const  { return w_Material; }
		std::vector<Object> GetObject()  const  { return w_Object; }
		vec3 GetEmvironment() const { return w_EMV; }
		Shield GetShieldMethod() { return w_shieldmtd; }

		//Metasequoia
		bool MQO_Load(const char* FileName);
		void CalcSurfaceNV();
		void CalcVertexNV();
		void CalcPolygonCenter();
		void CalcModelCenter();
		vec3 center(const vec3& p0, const vec3& p1, const vec3& p2);

		//CalcObjectField
		void GenDepthList();
		void SortByDepth(depthListArray& list);
		void DivideByDepth(depthListArray& front, depthListArray& back, double z);

		void AccommodatePolygonInBB();

		void RotInFourierSpaceForward(const WaveFront& source, WaveFront& reference, vec3* c, Interp interp);
		void RotInFourierSpaceBackward(const WaveFront& source, WaveFront& reference, const vec3& c, Interp interp);
		mat3 RotMatFromG2L(const vec3 &global, const vec3 &local);
		mat3 G2L();
		BoundingBox GetDiffractionRect(const double targetZ);
		void SetCurrentPolygon(const depthList dpl);

		void ClipSubfield(WaveFront& sub, WaveFront& frame, bool dir = true);
		void AddFieldToMFB(WaveFront& mfb);
		void AddObjectFieldPersubmodel(WaveFront& mfb, depthListArray& list);
		void CalcCenterOfModel(depthListArray& model);

		void AddObjectField(WaveFront& mfb, const unsigned int div, const mat3 &rot, bool exact = true, bool back = false);
		void ExShieldingAddingField(WaveFront &pfb);
		void SilhouetteShieldingAddingField(WaveFront &pfb);
		void GeneratePolygonField(WaveFront& field, const CurrentPolygon& polyL);

		bool IsInTriangle(vec3 p, const CurrentPolygon& polyL);
		void PaintTriangle(WaveFront& tfb, const CurrentPolygon &polyL, double amp);
		void MultiplyAperture(WaveFront& tfb, const CurrentPolygon& polyL);
		int LineFunc(int x, int y, int x1, int y1, int x2, int y2)
		{
			return (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1);
		}

		bool PolygonIsVisible();

		void MarkingRectangularPointsInFourierSpace(std::vector<vec3> &vec);
		BoundingBox GetBoundingBox(const std::vector<vec3> &vec);

		void Shading(WaveFront& field, const CurrentPolygon& polyL);
		void FlatShading(WaveFront& field, const CurrentPolygon& polyL);
		double GetCorrectedAmplitude(WaveFront& tfb, double brt);

		void SmoothShading(WaveFront& field, const CurrentPolygon& polyL);
		void Mapping(WaveFront& field, const CurrentPolygon& polyL);

		void SetUp(const mat3 &rot);

		void SetRandomPhase(WaveFront& field);

		//Geometry
		vec3 IntersectPoint(const Ray &ray);

		//Operator
		Model& operator +=(const vec3& vec);
		Model& operator *=(const mat3& mat);
		void mul(std::vector<vec3>& vec, const mat3& mat);
		void sub(std::vector<vec3>& vec, const vec3& vv);
		void fouriermul(std::vector<vec3>& vec, const mat3& mat);
		void fouriersub(std::vector<vec3>& vec, const vec3& vv);
	};
#endif
#endif