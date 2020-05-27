#ifndef MODELH
#define MODELH
#include<vector>
#include<string>
#include<algorithm>
#include<string>
#include"WaveFront.h"

#ifndef VECTORMATH
#define VECTORMATH
#include ".\vectormath\include\vectormath\scalar\cpp/vectormath_aos.h"
#endif

using namespace Vectormath::Aos;

using namespace std;

#ifndef VEC3
#define VEC3
	typedef Vector3 vec3;
#endif

#ifndef COL3
#define COL3
	typedef Vector3 col3;
#endif

#ifndef MAT3
#define MAT3
	typedef Matrix3 mat3;
#endif
#ifndef SHADE
	typedef enum SHADER
	{
		FLAT,
		SMOOTH
	}Shader;
#endif
#ifndef DIRECT
	typedef enum DIRECTION
	{
		DWIDTH,
		DHEIGHT,
		DDEPTH
	}Direction;
#endif
#ifndef METHOD
	typedef enum SHIELD
	{
		SILHOUETTE,
		EXACT
	}Shield;
#endif
#ifndef REF4
#define REF4
	struct Reflection4 {
		float _diffuse;
		float _ambient;
		float _emission;
		float _specular;
	};
#endif
#ifndef COL4
#define COL4
	struct Color4 {
		float _r;
		float _g;
		float _b;
		float _a;
	};
#endif
#ifndef UVC
#define UVC
	struct UV {
		float _u;
		float _v;
	};
#endif
#ifndef TRI
#define TRI
	//triangular face
	struct TRIANGLE {
		int _MaterialID;//material No.
		int _Index[3];//index
		UV _uv[3];//uv info.
		vec3 _SurfaceNV;//surface normal
		float _bbox[2][3];
		vec3 _center;
	};
#endif
#ifndef QUAD
#define QUAD
	//quadrilateral face
	struct QUADRILATERAL {
		int _MaterialID;
		int _Index[4];
		UV _uv[4];
		vec3 _SurfaceNV;
	};
#endif
#ifndef VTX
#define VTX
	//vertex
	struct VERTEX {
		vec3 _Coord;//coordinate
		vec3 _VertexNV;//vertex normal
	};
#endif
#ifndef OBJ
#define OBJ
	//object
	struct OBJECT {
		string _Name;//object name
		vector<VERTEX> _Vertex;//vertex data
		vector<TRIANGLE> _Triangle;
		vector<QUADRILATERAL> _Quadrilateral;
	};
#endif
#ifndef MATE
#define MATE
	//material
	struct MATERIAL {
		int _MaterialID;//ID
		string _MaterialName;//material name
		Color4 _Color;//color
		Reflection4 _ReflectionColor;//color of refrection
		float _Power;//shiness
		string _TextureName;//name of texture
		Image* _TextureImg;//image of texture
		bool _texexist = false;
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
		vector<depthList> _list;
		vec3 _modelcenter;
	};
#endif
#ifndef BDB
#define BDB
	class BoundingBox
	{
	public:
		vec3 _min;
		vec3 _max;
		vec3 _center;
		BoundingBox() :_min(vec3{ 0.0,0.0,0.0 }), 
			_max(vec3{ 0.0,0.0,0.0 }), 
			_center(vec3{ 0.0,0.0,0.0 }) {};
		BoundingBox(vec3 min, vec3 max) :_min(min), _max(max), _center( (max + min)/2 ) {};
		BoundingBox(double w, double h, double d, vec3 center):_center(center)
		{
			_min = vec3{ (float)(center.getX() - w / 2), (float)(center.getY() - h / 2), (float)(center.getZ() - d / 2) };
			_max = vec3{ (float)(center.getX() + w / 2), (float)(center.getY() + h / 2), (float)(center.getZ() + d / 2) };
		}
		double Width() { return _max.getX() - _min.getX(); }
		double Height() { return _max.getY() - _min.getY(); }
		double Depth() { return _max.getZ() - _min.getZ(); }
	};
#endif
#ifndef CUP
#define CUP
	class CurrentPolygon
	{
	public:
		int _objidx = 0;
		int _faceidx = 0;
		vec3 _vertex[3];
		vec3 _vertexnormal[3];
		vec3 _surfacenormal = vec3{0.0,0.0,0.0};
		vec3 _center = vec3{ 0.0,0.0,0.0 };
		int _MaterialID = 0;
		UV _uv[3];

		vec3 _diffractcenter = vec3{ 0.0,0.0,0.0 };

		CurrentPolygon() {};
		~CurrentPolygon() {};

		void SetVertex(const vector<vec3> v) 
		{
			if (v.size() != 3)
			{
				printf("\nSetVertex: input vector size is not suit!!");
				system("pause");
			}
			_vertex[0] = v[0];
			_vertex[1] = v[1];
			_vertex[2] = v[2];
		};

		vector<vec3> GetVertex()
		{
			vector<vec3> v;
			
			for (int i = 0; i < 3; i++)
				v.push_back(_vertex[i]);

			return v;
		};
	};
#endif
#ifndef RAY
#define RAY
	class Ray {//to handle line object 
	public:

		Ray() {}

		Ray(const vec3& o, const vec3& dir) : _origin(o), _direction(dir) {}

		const vec3& origin() const { return _origin; }
		const vec3& direction() const { return _direction; }
		vec3 at(float t) const { return _origin + t * _direction; }

	private:
		vec3 _origin;//start point
		vec3 _direction;//directional vector
	};
#endif
#ifndef MODELL
#define MODELL
	//model
	class MODEL {
	protected:
		double _px;//sampling interval along to x axis of main frame buffer
		double _py;//along to y axis
		double _lambda;//wavelength

		FILE* _fp;//file pointer to handle mqo file
		char _buf[255];
		string _str;
		depthListArray _depthlistArray;//array of depthlist(contain polygon index and z value of center of it)
		vector<MATERIAL> _Material;
		vector<OBJECT> _Object;
		vec3 _EMV = vec3{0.0,0.0,-1.0};//vector of ambient light

		Shader _shader;//shading method

		CurrentPolygon _currentpolygon;//contaions almost information to calculate polygon field

		Direction _dir = DWIDTH;//direction means accomodate to bounding box

		vec3 _center = vec3{0.0,0.0,0.0};//model's center

		bool _surface = true;//whether to calculate surface function

		vec3 _subcenter = vec3{ 0.0,0.0,0.0 };//center of submodel

		double _gamma = 2.2;

		Shield _shieldmtd = EXACT;//shielding method

		void Vertex_Set();//頂点情報セット
		void Material_Set();//マテリアル情報セット
		void Face_Set(int Object_num);//面情報セット
		string Split(string* str, char str1, char str2);//文字列分離
	public:

		BoundingBox _bbox;
		MODEL() {
		}
		MODEL(const MODEL& model)
		{
			_fp = model._fp;
			strcpy(_buf, model._buf);
			_str = model._str;
			_depthlistArray = model._depthlistArray;
			_Material = model._Material;
			_Object = model._Object;
			_EMV = model._EMV;
			_bbox = model._bbox;
			_surface = model._surface;
			_px = model._px;
			_py = model._py;
		}
		MODEL(const char* FileName, vec3 emv, Shader shade, BoundingBox bb, Direction dir, bool surface): _shader(shade),_dir(dir),_surface(surface) {
			if(MQO_Load(FileName))
              _bbox = bb;
			else
			{
				printf("Indicated file is not exist under the this directly....");
				system("pause");
			}
		}
		~MODEL()
		{}
		double GetPx() const { return _px; }
		double GetPy() const { return _py; }
		double GetLambda() const { return _lambda; }

		void SetShieldMethod(Shield mtd) { _shieldmtd = mtd; }

		FILE* FilePointer() { return _fp; }
		char Buffer() { return *_buf; }
		string String() { return _str; }
		depthListArray Depthlist() { return _depthlistArray; }
		vector<MATERIAL> Material() { return _Material; }
		vector<OBJECT> Object() { return _Object; }
		vec3 Emvironment() { return _EMV; }

		bool MQO_Load(const char* FileName);//ロード
		void CalcSurfaceNV();//面法線ベクトルの計算
		void CalcVertexNV();//頂点法線ベクトルの計算
		void CalcPolygonCenter();
		void CalcModelCenter();

		void GenDepthList();
		void SortByDepth(depthListArray& list);
		void DivideByDepth(depthListArray& front, depthListArray& back, double z);

		void AccommodatePolygonInBB();

		vec3 center(vec3 p0, vec3 p1, vec3 p2);
		void RotInFourierSpaceForward(const WaveFront& source, WaveFront& reference, vec3* c, Interpol interp);
		void RotInFourierSpaceBackward(const WaveFront& source, WaveFront& reference, vec3& c, Interpol interp);
		mat3 RotMatFromG2L(vec3 global, vec3 local);
		mat3 G2L();
		BoundingBox GetDiffractionRect(double targetZ);
		void SetCurrentPolygon(depthList dpl);
		void AddFieldToMFB(WaveFront& mfb);
		void AddObjectFieldPersubmodel(WaveFront& mfb, depthListArray& list);
		void CalcCenterOfModel(depthListArray& model);

		void AddObjectField(WaveFront& mfb,unsigned int div, mat3 &rot, bool exact = true, bool back = false);
		void ExShieldingAddingField(WaveFront &pfb);
		void SilhouetteShieldingAddingField(WaveFront &pfb);
		void GeneratePolygonField(WaveFront& field, CurrentPolygon& local);

		void PaintTriangle(WaveFront& tfb, CurrentPolygon &polyL, double amp);
		void MultiplyAperture(WaveFront& tfb, CurrentPolygon& polyL);
		int LineFunc(int x, int y, int x1, int y1, int x2, int y2)
		{
			return (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1);
		}

		bool PolygonIsVisible();

		void MarkingRectangularPointsInFourierSpace(vector<vec3> &vec);
		BoundingBox GetBoundingBox(vector<vec3> vec);

		void Shading(WaveFront& field, CurrentPolygon& local);
		void FlatShading(WaveFront& field, CurrentPolygon& local);
		double GetCorrectedAmplitude(WaveFront& tfb, double brt);

		void SmoothShading(WaveFront& field, CurrentPolygon& local);
		void Mapping(WaveFront& field, CurrentPolygon& local);

		void SetUp(mat3 &rot);

		vec3 IntersectPoint(Ray &ray);

		MODEL& operator +=(vec3 vec)
		{
			vec3 v;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
			for (int n = 0; n < (*this)._Object.size(); ++n)
				for (int m = 0; m < (*this)._Object[n]._Vertex.size(); ++m)
				{
					v = _Object[n]._Vertex[m]._Coord;
					v = v + vec;
					_Object[n]._Vertex[m]._Coord = v;
				}
			return *this;
		}

		MODEL& operator *=(mat3 mat)
		{
			vec3 v;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::thread::hardware_concurrency())
			for (int n = 0; n < (*this)._Object.size(); ++n)
				for (int m = 0; m < (*this)._Object[n]._Vertex.size(); ++m)
				{
					v = _Object[n]._Vertex[m]._Coord;
					v = v - _center;
					v = mat * v;
					v = v + _center;
					_Object[n]._Vertex[m]._Coord = v;
				}
			return *this;
		}

		void mul(vector<vec3> &vec, mat3 mat)
		{
			int i;
			for (auto& v : vec)
				v = mat * v;
		}

		void sub(vector<vec3> &vec, vec3 vv)
		{
			int i;
			for (auto& v : vec)
				v -= vv;
		}

		void fouriermul(vector<vec3> vec, mat3 mat)
		{
			int i;
			double invlambda2 = 1 / _lambda / _lambda;
			vec3 v;
			float w;
			auto itr = vec.begin();
			while (itr != vec.end())
			{
				v = mat * (*itr);//access the element and transform
				w = invlambda2 - v.getX() * v.getX() - v.getY() * v.getY();
				if (w < 0) {
					// return the iterator indicate next element of deleted one
					itr = vec.erase(itr);
				}
				// in case of not delete element, get the next indicator
				else {
					*itr = vec3{v.getX(), v.getY(), sqrt(w)};
					itr++;
				}
			}
		}
		void fouriersub(vector<vec3> vec, vec3 vv)
		{
			int i;
			double invlambda2 = 1 / _lambda / _lambda;
			vec3 v;
			float w;
			auto itr = vec.begin();
			while (itr != vec.end())
			{
				v = (*itr) - vv;//access the element and transform
				w = invlambda2 - v.getX() * v.getX() - v.getY() * v.getY();
				if (w < 0) {
					// return the iterator indicate next element of deleted one
					itr = vec.erase(itr);
				}
				// in case of not delete element, get the next indicator
				else {
					*itr = vec3{ v.getX(), v.getY(), sqrt(w) };
					itr++;
				}
			}
		}
	};
#endif
#endif