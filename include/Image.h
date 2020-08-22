#ifndef IMAGE_H
#define IMAGE_H

#include "stb_image.h"

#include "stb_image_write.h"

#include "BLAS.h"

#include <omp.h>
#include <thread>
#include <string>
#include <vector>

#ifndef VECTORMATH
#define VECTORMATH
//#include "vectormath_aos.h"
#endif

#ifndef PI
#define PI (3.14159265359f)
#endif

#ifndef PI2
#define PI2 (6.28318530718f)
#endif

#ifndef RECIP_PI
#define RECIP_PI (0.31830988618f)
#endif

#ifndef RECIP_PI2
#define RECIP_PI2 (0.15915494f)
#endif

#ifndef LOG2
#define LOG2 (1.442695f)
#endif

#ifndef EPSILON
#define EPSILON (1e-6f)
#endif

#ifndef GAMMA_FACTOR
#define GAMMA_FACTOR (2.2f)
#endif

#ifndef VEC3
#define VEC3
//typedef Vectormath::Aos::Vector3 vec3;
typedef BLASVector vec3;
#endif

#ifndef COL3
#define COL3
//typedef Vectormath::Aos::Vector3 col3;
typedef BLASVector col3;
#endif

#ifndef MAT3
#define MAT3
//typedef Vectormath::Aos::Matrix3 mat3;
typedef BLASMatrix mat3;
#endif

#ifndef FILEHEADERSIZE
#define FILEHEADERSIZE (14)					// header size
#endif

#ifndef INFOHEADERSIZE
#define INFOHEADERSIZE (40)					// info header size
#endif

#ifndef HEADERSIZE 
#define HEADERSIZE (FILEHEADERSIZE+INFOHEADERSIZE)
#endif

//using namespace	std;

#ifndef DEG
#define DEG (PI/180)
#endif

class Image {
public:
	struct RGB {
		unsigned char R;
		unsigned char G;
		unsigned char B;
	};// struct to hundle image as RGB

	Image() :w_nx(1024), w_ny(1024),Image_pix(nullptr) {};

	Image(int W, int H):w_nx(W),w_ny(H) {
		Image_pix.reset(new RGB[w_nx * w_ny]);
	}

	int GetWidth() const { return w_nx; }
	int GetHeight() const { return w_ny; }
	void* GetPixels() const { return Image_pix.get(); }
	int idxij(const int i, const int j) const  { return j * w_nx + i; }

	//Load
	void Write(const int i, const int j, const double R, const double G, const double B, const bool read = true);
	vec3 Load(const int i, const int j);
	Image* Read_Bmp(const char* filename);

private:
	int w_nx;// image width
	int w_ny;// image height
	std::unique_ptr<RGB[]> Image_pix;
};
#endif