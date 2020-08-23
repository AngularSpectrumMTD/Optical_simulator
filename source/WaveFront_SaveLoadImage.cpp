#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#endif
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#endif
#include"../include/WaveFront.h"
using namespace std;
// function related to save as image and load real part from image
void WaveFront::SaveBmp(const char* name, Out type)
{
	Image* image = new Image(w_nx, w_ny);
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_max_threads())
	for (i = 0; i < w_nx; i++)
	{
		for (j = 0; j < w_ny; j++)
		{
			double val = 0;
			switch (type)
			{
			case REAL:
				val = 255.99f * (GetReal(i, j) + 1.0f) / 2.0f;
				break;
			case IMAGE:
				val = 255.99f * (GetImage(i, j) + 1.0f) / 2.0f;
				break;
			case AMPLITUDE:
				val = 255.99f * GetAmplitude(i, j);
				break;
			case INTENSITY:
				val = 255.99f * GetIntensity(i, j);
				break;
			case PHASE:
				val = 255.99f * (GetPhase(i, j) / 2.0f / PI + 1.0f / 2.0f);
				break;
			}
			image->Write(i, (w_ny - j - 1), val, val, val, false);// data must be inserted from bottom
		}
	}
	stbi_write_bmp(name, w_nx, w_ny, sizeof(Image::RGB), image->GetPixels());
}
void WaveFront::Normalize()
{
	double max_amp = 0;
	for (int i = 0; i < w_nx; i++)
	{
		for (int j = 0; j < w_ny; j++)
		{
			max_amp = max<double>(max_amp, GetAmplitude(i, j));// seek max value
		}
	}

	for (int i = 0; i < w_nx; i++)
	{
		for (int j = 0; j < w_ny; j++)
		{
			SetPixel(i, j, GetPixel(i, j) / max_amp);// normalization
		}
	}
}
WaveFront& WaveFront::LoadBmp(const char* filename)
{
	int i, j;
	int bpp;//byte per pixel
	unsigned char* pixels = nullptr;
	int loopnx, loopny;
	pixels = stbi_load(filename,&loopnx,&loopny,&bpp,0);

	double expX = log(static_cast<double>(loopnx)) / log(2.0);
	double expY = log(static_cast<double>(loopny)) / log(2.0);

	w_nx = nearPow2(pow(2.0, expX));
	w_ny = nearPow2(pow(2.0, expY));

	w_data.reset(new complex<double>[w_nx * w_ny]);

	Clear();
	
	unsigned int offsetX = (w_nx - loopnx) / 2;
	unsigned int offsetY = (w_ny - loopny) / 2;
	
	unsigned int index;
#pragma omp parallel for private(i, j) num_threads(omp_get_max_threads())
	for (i = 0; i < loopnx; i++)
	{
		for (j = 0; j < loopny; j++)
		{
			index = (loopny - j - 1) * loopnx + i;
			unsigned char r = pixels[index * bpp + 0];
			unsigned char g = pixels[index * bpp + 1];
			unsigned char b = pixels[index * bpp + 2];
			double gray = (double)((r + g + b) / 3.0f / 256.0f);
			SetPixel(i + offsetX, j + offsetY, complex<double>(sqrt(gray), 0));
		}
	}
	stbi_image_free(pixels);

	return *this;
}