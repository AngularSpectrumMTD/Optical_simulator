#include"ImagingWaveFront.h"

ImagingWaveFront::ImagingWaveFront()
	:WaveFront(4096, 4096, 1e-6, 1e-6, 630e-9), d_eye(20e-3), d_pupil(4e-3)
{
	vec3 origin{ 0.0, 0.0, 20e-3 };
	SetOrigin(origin);
	Clear();
}

WaveFront& ImagingWaveFront::SetOrigin(vec3 p)
{
	WaveFront::SetOrigin(p);
	this->viewpoint = p;
	return *this;
}

void ImagingWaveFront::SetEyeParam()
{
	this->SetNx(4096);
	this->SetNy(4096);
	this->SetPx(1e-6);
	this->SetPy(1e-6);
	this->d_eye = 24e-3;
	this->d_pupil = 6e-3;
	this->Init();
}
void ImagingWaveFront::SetEye(WaveFront& eye, const vec3 & p)
{
	double gd = length(this->GetOrigin() - p);	//注視点までの距離
	//printf("<imaging distance %lf>", gd);
	double f = gd * d_eye / (gd + d_eye);		//レンズの焦点距離

	eye.AllSet(1.0);
	eye.SetQuadraticPhase(f);

	int i, j;							//円形瞳関数の描画
	for (j = 0; j < eye.GetNy(); j++)
	{
		for (i = 0; i < eye.GetNx(); i++)
		{
			double x = eye.itox(i);
			double y = eye.jtoy(j);

			if (x * x + y * y > d_pupil* d_pupil / 4.0)
			{
				eye.SetPixel(i, j, std::complex<double>(0.0, 0.0));
			}
		}
	}
}
void ImagingWaveFront::View(const WaveFront& wf, const vec3 &p)		//光波wfのp点を注視する
{
	if (wf.GetNormal().getX() != 0 && wf.GetNormal().getY() != 0 && wf.GetNormal().getZ() != 1)
	{
		printf(">>ERROR: Normal vector of the WaveField object must be (0, 0, 1)\n");
		printf(">>Process is terminated forcibly...\n");
		exit(0);
	}
	SetLambda(wf.GetLambda());
	SetNormal(vec3(0, 0, 1.0));	//スクリーンを光軸方向に向ける
	Clear();
	ShiftedAsmPropEx(wf);

	//WaveFront test(*this);
	////test.fft2D(-1);
	//test.Normalize();
	//////test.SaveBmp("before imaging center.bmp",INTENSITY);
	//////test.SaveBmp("before imaging right.bmp", INTENSITY);
	//////test.SaveBmp("before imaging center sp.bmp",INTENSITY);
	//test.SaveBmp("before imaging before rotate.bmp", INTENSITY);

	//SaveAsWaveFront("SAVED.txt");
	//system("pause");

	//LoadAsWaveFront("SAVED.txt");
	printf("<imaging START>");
	Imaging(p);
}

void ImagingWaveFront::Imaging(vec3 p)
{
	vec3 reserveViewPoint = this->GetOrigin();
	WaveFront source = *this;
	this->SetNormal(this->GetOrigin() - p);
	vec3 freq;
	//source.TiltedAsmProp(*this, BICUBIC, &freq);
	/*{
		WaveFront save = source;
		save.Normalize();
		save.SaveBmp("origin.bmp", INTENSITY);
	}*/
	TiltedAsmProp(source, BICUBIC, &freq);
	/*{
		if (GetOrigin().getX() == 0 && GetOrigin().getY() == 0)
		{
			WaveFront save = *this;
			save.Normalize();
			save.SaveBmp("rotated center.bmp", INTENSITY);
		}
		else
		{
			WaveFront save = *this;
			save.Normalize();
			save.SaveBmp("rotated right.bmp", INTENSITY);
		}
	}
	vec3 origin = GetOrigin();
	printf("\nEYE origin =");
	DispVec(origin);
	printf("\nCARRIER FREQ =");
	DispVec(freq);*/
	this->MultiplyPlaneWave(freq.getX() * this->GetLambda(), freq.getY() * this->GetLambda());//koko

	source.pitchtrans();
	WaveFront &lens = source;
	SetEye(lens, p);
	*this *= lens;
	/*WaveFront save = *this;
	save.Normalize();
	save.SaveBmp("multipled.bmp",INTENSITY);*/

	this->ExactAsmProp(d_eye);		//スクリーン位置まで伝搬
	this->SetOrigin(reserveViewPoint);
}