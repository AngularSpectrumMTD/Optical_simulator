#include"WaveFront.h"
#include"Model.h"

double lambda = 630e-9;//îgí∑[m]
//double lambda = 550e-9;//îgí∑[m]
int SCALE = 2;//max 4!!!!!!!!
int WIDTH = 1024 * SCALE;//âÊëúâ°ïù
int HEIGHT = 1024 * SCALE;//âÊëúècïù
double Dx = 1e-6;
double Dy = 1e-6;
double rw = WIDTH * Dx;
double rh = HEIGHT * Dy;
//double r = 1.0e-4;//äJå˚ïù10[mm]
double r = 1.0e-4;//äJå˚ïù10[mm]

//int main()
//{
//	WaveFront mfb(WIDTH, HEIGHT, Dx, Dy, lambda);
//
//	vec3 backpos = vec3{ 0.0,0.0,(float)(-1.5 * rw) };
//
//	mfb.SetOrigin(backpos);
//
//	mfb.AllSet(0.5);
//	mfb.ModRandomphase();
//
//	BoundingBox bb0(rw / 4, rh / 4, 5e-4, vec3(0, 0, -rw));
//	BoundingBox bb1(rw / 2, rh / 2, 5e-4, vec3(0, 0, -rw));
//	BoundingBox bbfull(rw, rh, 5e-4, vec3(0, 0, -rw));
//	mat3 id = mat3::identity();
//	mat3 rot5 = mat3::rotationY(5 * DEG);
//	mat3 rot10 = mat3::rotationY(10 * DEG);
//	mat3 rot20 = mat3::rotationY(20 * DEG);
//	mat3 rot40 = mat3::rotationY(40 * DEG);
//	mat3 rot60 = mat3::rotationY(60 * DEG);
//	mat3 rot80 = mat3::rotationY(80 * DEG);
//	mat3 rot2020 = mat3::rotationX(20 * DEG) * mat3::rotationY(20 * DEG);
//	mat3 rot2040 = mat3::rotationX(20 * DEG) * mat3::rotationY(40 * DEG);
//	mat3 rot2060 = mat3::rotationX(20 * DEG) * mat3::rotationY(60 * DEG);
//	mat3 rot2080 = mat3::rotationX(20 * DEG) * mat3::rotationY(80 * DEG);
//
//	char name[200] = "drag bunny.mqo";
//	strcpy(name, "sphere.mqo");
//	strcpy(name, "3sphere.mqo");
//	//strcpy(name, "error.mqo");
//	//strcpy(name, "test.mqo");
//	//strcpy(name, "pien.mqo");
//	//strcpy(name, "bunny.mqo");
//	//strcpy(name, "test2.mqo");
//	//strcpy(name, "pera.mqo");
//
//	//strcpy(name, "3tri.mqo");
//
//	//MODEL model(name, vec3{ -1,-1,-1 }, FLAT, bb1, DWIDTH, false);
//	//MODEL model(name, vec3{ -1,-1,-1 }, FLAT, bb1, DWIDTH, true);//É|ÉäÉSÉìÇ‡
//	MODEL model(name, vec3{ -1,-1,-1 }, SMOOTH, bbfull, DWIDTH, true);//É|ÉäÉSÉìÇ‡
//	//MODEL model(name, vec3{ -1,-1,-1 }, SMOOTH, bbfull, DWIDTH, false);//É}ÉXÉNÇÃÇ›
//	//MODEL model(name, vec3{ 0,-1,-1 }, FLAT, bbfull, DWIDTH, false);//É}ÉXÉNÇÃÇ›
//
//	//model.SetShieldMethod(SILHOUETTE);
//
//	//model.AddObjectField(mfb, 1, id, true, false);
//	model.AddObjectField(mfb, 1, id, true, true);
//	
//	//WaveFront mfbcopy(WIDTH, HEIGHT, Dx, Dy, lambda);
//	//mfbcopy = mfb;
//
//	//mfbcopy.ExactAsmProp(-4*model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[1]");
//	//mfb.SaveBmp("ç≈èIåıîg0-1.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[2]");
//	//mfb.SaveBmp("ç≈èIåıîg0-2.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[3]");
//	//mfb.SaveBmp("ç≈èIåıîg0-3.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[4]");
//	//mfb.SaveBmp("ç≈èIåıîg0-4.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[5]");
//	//mfb.SaveBmp("ç≈èIåıîg0-5.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[6]");
//	//mfb.SaveBmp("ç≈èIåıîg0-6.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[7]");
//	//mfb.SaveBmp("ç≈èIåıîg0-7.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[8]");
//	//mfb.SaveBmp("ç≈èIåıîg0-8.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[9]");
//	//mfb.SaveBmp("ç≈èIåıîg0-9.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[10]");
//	//mfb.SaveBmp("ç≈èIåıîg0-10.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[11]");
//	//mfb.SaveBmp("ç≈èIåıîg0-11.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[12]");
//	//mfb.SaveBmp("ç≈èIåıîg0-12.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[13]");
//	//mfb.SaveBmp("ç≈èIåıîg0-13.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[14]");
//	//mfb.SaveBmp("ç≈èIåıîg0-14.bmp", INTENSITY);
//
//	//mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	//mfb = mfbcopy;
//	//mfb.Normalize();
//	//printf("[15]");
//	//mfb.SaveBmp("ç≈èIåıîg0-15.bmp", INTENSITY);
//
//	WaveFront mfbcopy(WIDTH, HEIGHT, Dx, Dy, lambda);
//	mfbcopy = mfb;
//
//	mfb.Normalize();
//	mfb.SaveBmp("ç≈èIåıîg0.bmp", INTENSITY);
//
//	/*mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	mfb = mfbcopy;
//	mfb.Normalize();
//	mfb.SaveBmp("ç≈èIåıîg1.bmp", INTENSITY);
//
//	mfbcopy.ExactAsmProp(model._bbox.Depth() / 8);
//	mfb = mfbcopy;
//	mfb.Normalize();
//	mfb.SaveBmp("ç≈èIåıîg2.bmp", INTENSITY);*/
//}

//int main()
//{
//	WaveFront input(128,128,2e-6,2e-6,lambda);
//	input.LoadBmp("Ç“Ç¶ÇÒ128.bmp");
//	//input.ModRandomphase();
//	//input.SaveAsWaveFront("buriburi.txt");
//
//	
//	input.SaveBmp("input.bmp", INTENSITY);
//	input.Embed();
//	input.Embed();
//	input.Embed();
//
//	double d = 20e-3;
//
//	double f = 10e-3;
//	WaveFront lens(input);
//	lens.Clear();
//	lens.SetGaussian(1e-3, 100);
//	lens.SetQuadraticPhase(f);
//	lens.SaveBmp("lens.bmp", INTENSITY);
//	lens.SaveBmp("lensp.bmp", PHASE);
//
//	input.ExactAsmProp(d);
//
//	WaveFront save(input);
//	save = input;
//	save.Normalize();
//	save.Extract();
//	save.Extract();
//	save.Extract();
//	save.SaveBmp("proped.bmp", INTENSITY);
//
//	input *= lens;
//
//	input.ExactAsmProp(d);
//
//	save.Embed();
//	save.Embed();
//	save.Embed();
//	save = input;
//	save.Extract();
//	save.Extract();
//	save.Extract();
//	save.Normalize();
//	save.SaveBmp("image.bmp",INTENSITY);
//}

//int main()
//{
//	WaveFront input(256, 256, 2e-6, 2e-6, lambda);
//	input.SetGaussian(256*2e-6/8,2);
//	input.ModRandomphase();
//	input.SaveBmp("input.bmp",INTENSITY);
//	input.SaveBmp("inputp.bmp", PHASE);
//	input.SaveAsWaveFront("test.txt");
//	input.LoadAsWaveFront("test.txt");
//	input.SaveBmp("input2.bmp", INTENSITY);
//	input.SaveBmp("inputp2.bmp", PHASE);
//	input.SaveAsWaveFront("buriburi2.txt");
//}