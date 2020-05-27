#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#endif
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#endif
#include"WaveFront.h"

// function related to save as image and load real part from image
void WaveFront::SaveBmp(const char* name, Outputformat type)
{
	double val = 0;
	Image* image = new Image(_nx, _ny);
	int i, j;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < _nx; i++)
		for (j = 0; j < _ny; j++)
		{
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
			image->Write(i, (_ny - j - 1), val, val, val, false);// data must be inserted from bottom
		}
	stbi_write_bmp(name, _nx, _ny, sizeof(Image::RGB), image->Pixels());
}
void WaveFront::Normalize()
{
	double max_amp = 0;
	for (int i = 0; i < _nx; i++)
		for (int j = 0; j < _ny; j++)
			max_amp = max<double>(max_amp, GetAmplitude(i, j));// seek max value

	for (int i = 0; i < _nx; i++)
		for (int j = 0; j < _ny; j++)
			SetPixel(i, j, GetPixel(i, j) / max_amp);// normalization
}
WaveFront& WaveFront::LoadBmp(const char* filename)
{
	int i, j;
	int bpp;//byte per pixel
	unsigned char* pixels = nullptr;
	int loopnx, loopny;
	pixels = stbi_load(filename,&loopnx,&loopny,&bpp,0);

	double expX = log((double)loopnx) / log(2.0);
	double expY = log((double)loopny) / log(2.0);

	_nx = nearPow2(pow(2.0, expX));
	_ny = nearPow2(pow(2.0, expY));

	_data.reset(new complex<double>[_nx * _ny]);

	Clear();

	unsigned char r, g, b;
	double gray;
	
	unsigned int offsetX = (_nx - loopnx) / 2;
	unsigned int offsetY = (_ny - loopny) / 2;
	
	unsigned int index;
#pragma omp parallel for private(i, j) num_threads(omp_get_num_threads())
	for (i = 0; i < loopnx; i++)
		for (j = 0; j < loopny; j++)
		{
			index = (loopny - j - 1) * loopnx + i;
			r = pixels[index * bpp + 0];
			g = pixels[index * bpp + 1];
			b = pixels[index * bpp + 2];
			gray = (double)((r + g + b) / 3.0f / 256.0f);
			SetPixel(i + offsetX, j + offsetY, complex<double>(sqrt(gray),0));
		}
	stbi_image_free(pixels);

	return *this;
}