#ifndef IMAGE_H
#define IMAGE_H

#include "stb_image.h"

#include "stb_image_write.h"

#include <omp.h>
#include <thread>
#include <string>
#include <vector>

#ifndef VECTORMATH
#define VECTORMATH
#include ".\vectormath\include\vectormath\scalar\cpp/vectormath_aos.h"
#endif

using namespace Vectormath::Aos;

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

#ifndef FILEHEADERSIZE
#define FILEHEADERSIZE (14)					// header size
#endif

#ifndef INFOHEADERSIZE
#define INFOHEADERSIZE (40)					// info header size
#endif

#ifndef HEADERSIZE 
#define HEADERSIZE (FILEHEADERSIZE+INFOHEADERSIZE)
#endif

using namespace	std;

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

	Image() :_nx(1024), _ny(1024),Image_pix(nullptr) {};

	Image(int W, int H):_nx(W),_ny(H) {
		Image_pix.reset(new RGB[_nx * _ny]);
	}

	int Width() const { return _nx; }
	int Height() const { return _ny; }
	void* Pixels() const { return Image_pix.get(); }

	int idxij(int i, int j) { return j * _nx + i; }

	void Write(int i, int j, double R, double G, double B, bool read = true) {

		int Image_index = idxij(i, j);
		if (read == false)
		{
			Image_pix[Image_index].R = static_cast<unsigned char>(R * 255.99f);
			Image_pix[Image_index].G = static_cast<unsigned char>(G * 255.99f);
			Image_pix[Image_index].B = static_cast<unsigned char>(B * 255.99f);
		}

		Image_pix[Image_index].R = (unsigned char)R;
		Image_pix[Image_index].G = (unsigned char)G;
		Image_pix[Image_index].B = (unsigned char)B;

	}// writing data

	vec3 Load(int i, int j) {
		vec3 texcolor(0);
		int Image_index = idxij(i, j);
		texcolor.setX(Image_pix[Image_index].R);
		texcolor.setY(Image_pix[Image_index].G);
		texcolor.setZ(Image_pix[Image_index].B);
		return texcolor;
	}

	Image* Create_Image(int width, int height)
	{
		Image* img;

		if ((img = (Image*)malloc(sizeof(Image))) == NULL) {
			fprintf(stderr, "Allocation error\n");
			return NULL;
		}

		img->_nx = width;
		img->_ny = height;

		return img;
	}

	Image* Read_Bmp(const char* filename) {

		int real_width;					// bite per line
		unsigned int width, height;			// pixel num
		unsigned int color;			// depth of bit
		FILE* fp;
		char header_buf[HEADERSIZE];	// header info
		unsigned char* bmp_line_data;  // box per line
		Image* img;
		std::unique_ptr<Image> image;
		if ((fp = fopen(filename, "rb")) == NULL) {
			fprintf(stderr, "Error: %s could not read.", filename);
			return NULL;
		}
		fread(header_buf, sizeof(unsigned char), HEADERSIZE, fp); // header part
		// distinguish bitmap format
		if (strncmp(header_buf, "BM", 2)) {
			fprintf(stderr, "Error: %s is not Bitmap file.", filename);
			return NULL;
		}
		memcpy(&width, header_buf + 18, sizeof(width)); // get real width
		memcpy(&height, header_buf + 22, sizeof(height)); // get real height
		memcpy(&color, header_buf + 28, sizeof(unsigned int)); // distinguish depth of bit
		// only handle 24 bit
		if (color != 24) {
			fprintf(stderr, "Error: %s is not 24bit color image", filename);
			return NULL;
		}
		// RGB info data byte per line must be magnification of 4byte
		real_width = width * 3 + width % 4;
		// generate buffer to get RGB info per line dynamically
		if ((bmp_line_data = (unsigned char*)malloc(sizeof(unsigned char) * real_width)) == NULL) {
			fprintf(stderr, "Error: Allocation error.\n");
			return NULL;
		}
		// generate buffer to get RGB info per image dynamically
		if ((img = Create_Image(width, height)) == NULL) {
			free(bmp_line_data);
			fclose(fp);
			return NULL;
		}
		img = new Image(width, height);
		// RGB info line up from left-bottom to right and bottom to top
		for (int j = 0; j < img->_ny; j++) {
			fread(bmp_line_data, 1, real_width, fp);
			for (int i = 0; i < img->_nx; i++) {
				img->Write(i, (img->_ny - j - 1), bmp_line_data[i * 3 + 2], bmp_line_data[i * 3 + 1], bmp_line_data[i * 3]);
			}
		}
		free(bmp_line_data);
		fclose(fp);
		return img;
	}

private:
	int _nx;// image width
	int _ny;// image height
	unique_ptr<RGB[]> Image_pix;
};
#endif