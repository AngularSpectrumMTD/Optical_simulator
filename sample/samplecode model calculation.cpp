#include"../include/ImagingWaveFront.h"
#include"../include/Model.h"//Model‚ªWaveFront‚ğŠÜ‚ñ‚Å‚¢‚é

double lambda = 630e-9;//”g’·[m]
double Dx = 1e-6;//•W–{ŠÔŠu
double Dy = 1e-6;//•W–{ŠÔŠu

float scale = 16;

int WIDTH = 1024*64/scale;//‰æ‘œ‰¡•
int HEIGHT = 1024*64/scale;//‰æ‘œc•

//•¨‘ÌƒTƒCƒY
double modelwidth = 10e-3/scale;
double modelheight = 10e-3 / scale;
double modeldepth = 10e-3 / scale;

double modelcenterdepth = 50e-3 / scale;//•¨‘Ì’†S‰œs‚«
double backposdepth = 100e-3 / scale;//”wŒi‰œs‚«
double viewpointdepth = 50e-3 / scale;//‹“_ˆÊ’u
//int d = -3;
int d = -5;
double viewpointgap = d * 1e-3 / scale;

int main()//2048
{
	WaveFront mfb(WIDTH, HEIGHT, Dx, Dy, lambda);

	char filenamefield[200];
	char filenameimagee[200];
	char filenameimage[200];
	char filenameimage2[200];
	sprintf(filenamefield, "•¨‘ÌŒõ”g_%d.txt", WIDTH);
	sprintf(filenameimagee, "•¨‘ÌŒõ”g_%d.bmp", WIDTH);
	sprintf(filenameimage, "•¨‘ÌŒõ”gŒ‹‘œ_%d.bmp", WIDTH);
	sprintf(filenameimage2, "•¨‘ÌŒõ”gŒ‹‘œ%d_%d.bmp", d, WIDTH);

	vec3 backpos = vec3{ 0.0,0.0,static_cast<float>(backposdepth) };

	mfb.SetOrigin(backpos);

	BoundingBox bb(modelwidth, modelheight, modeldepth, vec3(0, 0, -modelcenterdepth));

	//char name[200] = "Bunny.mqo";
	char name[200] = "STF DRAGON DOWN.mqo";
	//char name[200] = "kinoko.mqo";

	Model model(name, vec3{ -1,-1,-1 }, SMOOTH, bb, DWIDTH, true);

	model.SetShieldMethod(SILHOUETTE);
	mat3 id = mat3::identity();

	//•¨‘ÌŒõ”gŒvZ•”•ª
	model.AddObjectField(mfb, 1, id, false, false);
	mfb.SaveAsWaveFront(filenamefield);
	mfb.Normalize();
	mfb.SaveBmp(filenameimagee, INTENSITY);
	//

	mfb.LoadAsWaveFront(filenamefield);

	//mfb.DispParam();

	vec3 modelcenter = vec3(0, 0, -modelcenterdepth);
	vec3 view_point = vec3(0, 0, viewpointdepth);
	vec3 view_point2 = vec3(viewpointgap, 0, viewpointdepth);

	ImagingWaveFront eye;
	eye.SetEyeParam();

	printf("\ncenter STRAT");
	eye.SetOrigin(view_point);
	eye.View(mfb, modelcenter);
	eye.TransformforBrainImage();
	eye.Normalize();
	eye.SaveBmp(filenameimage, INTENSITY);
	printf("\ncenter OK");

	printf("\nright START");
	eye.SetOrigin(view_point2);
	eye.View(mfb, modelcenter);
	eye.TransformforBrainImage();
	eye.Normalize();
	eye.SaveBmp(filenameimage2, INTENSITY);
	printf("\nright OK");
}