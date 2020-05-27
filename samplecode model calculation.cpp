#include"WaveFront.h"
#include"Model.h"

double lambda = 630e-9;//”g’·[m]
int WIDTH = 1024;//‰æ‘œ‰¡•
int HEIGHT = 1024;//‰æ‘œc•
double Dx = 1e-6;//•W–{ŠÔŠu
double Dy = 1e-6;//•W–{ŠÔŠu
double rw = WIDTH * Dx;
double rh = HEIGHT * Dy;

int main()
{
	WaveFront mfb(WIDTH, HEIGHT, Dx, Dy, lambda);

	vec3 backpos = vec3{ 0.0,0.0,(float)(-1.5 * rw) };

	mfb.SetOrigin(backpos);

	mfb.AllSet(0.5);
	mfb.ModRandomphase();

	BoundingBox bb(rw, rh, 5e-4, vec3(0, 0, -rw));

	char name[200] = "Bunny.mqo";

	MODEL model(name, vec3{ -1,-1,-1 }, SMOOTH, bb, DWIDTH, true);

	//model.SetShieldMethod(SILHOUETTE);
	mat3 id = mat3::identity();
	model.AddObjectField(mfb, 1, id, true, true);

	mfb.Normalize();
	mfb.SaveBmp("•¨‘ÌŒõ”g.bmp", INTENSITY);
}