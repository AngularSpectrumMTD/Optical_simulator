#include "..\include\Image.h"

void Image::Write(const int i, const int j, const double R, const double G, const double B, const bool read) {

	int Image_index = idxij(i, j);
	if (read == false)
	{
		Image_pix[Image_index].R = static_cast<unsigned char>(R * 255.99f);
		Image_pix[Image_index].G = static_cast<unsigned char>(G * 255.99f);
		Image_pix[Image_index].B = static_cast<unsigned char>(B * 255.99f);
	}

	Image_pix[Image_index].R = static_cast<unsigned char>(R);
	Image_pix[Image_index].G = static_cast<unsigned char>(G);
	Image_pix[Image_index].B = static_cast<unsigned char>(B);

}// writing data

vec3 Image:: Load(const int i, const int j) {
	vec3 texcolor(0);
	int Image_index = idxij(i, j);
	texcolor.setX(Image_pix[Image_index].R);
	texcolor.setY(Image_pix[Image_index].G);
	texcolor.setZ(Image_pix[Image_index].B);
	return texcolor;
}

//Image* Image::Create_Image(int width, int height)
//{
//	Image* img;
//
//	if ((img = (Image*)malloc(sizeof(Image))) == NULL) {
//		fprintf(stderr, "Allocation error\n");
//		return NULL;
//	}
//
//	img->w_nx = width;
//	img->w_ny = height;
//
//	return img;
//}

Image* Image::Read_Bmp(const char* filename) {

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
	if ((img = new (std::nothrow)Image(width, height)) == NULL) {
		free(bmp_line_data);
		fclose(fp);
		return NULL;
	}
	//img = new Image(width, height);
	// RGB info line up from left-bottom to right and bottom to top
	for (int j = 0; j < img->w_ny; j++) {
		fread(bmp_line_data, 1, real_width, fp);
		for (int i = 0; i < img->w_nx; i++) {
			img->Write(i, (img->w_ny - j - 1), bmp_line_data[i * 3 + 2], bmp_line_data[i * 3 + 1], bmp_line_data[i * 3]);
		}
	}
	free(bmp_line_data);
	fclose(fp);
	return img;
}