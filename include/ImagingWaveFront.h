#include"WaveFront.h"

class ImagingWaveFront : public WaveFront{

protected:
	double w_d_eye;//diameter of eye ball
	double w_d_pupil;//diameter of pupil 

	vec3 viewpoint;//position of viewpoint

public:
	ImagingWaveFront(void);//constructor to imaging by use of default setting
	ImagingWaveFront(double deye, double dpupil, vec3 vp, const WaveFront image);//constructor to imaging by use of default setting

	//override WaveFront::SetOrigin()
	WaveFront& SetOrigin(vec3 p);

	//setter
	void SetImagingDistance(const double dd) { w_d_eye = dd; }
	void SetPupilDiameter(const double dd) { w_d_pupil = dd; }
	
	//getter
	double GetImagingDistance(void) { return w_d_eye; }
	double GetPupilDiameter(void) { return w_d_pupil; }

	void SetEyeParam();

	//imaging
	void Imaging(vec3 p);
	void View(const WaveFront& wf, const vec3 &p);
	void SetEye(WaveFront & eye, const vec3 &p);
};
