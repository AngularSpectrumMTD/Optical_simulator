#include"WaveFront.h"

class ImagingWaveFront : public WaveFront{

protected:
	double d_eye;//diameter of eye ball
	double d_pupil;//diameter of pupil 

	vec3 viewpoint;//position of viewpoint

public:
	ImagingWaveFront(void);//constructor to imaging by use of default setting

	//override WaveFront::SetOrigin()
	WaveFront& SetOrigin(vec3 p);

	//setter
	void SetImagingDistance(double dd) { d_eye = dd; }
	void SetPupilDiameter(double dd) { d_pupil = dd; }
	
	//getter
	double GetImagingDistance(void) const { return d_eye; }
	double GetPupilDiameter(void) const { return d_pupil; }

	void SetEyeParam();

	//imaging
	void Imaging(vec3 p);
	void View(const WaveFront& wf, const vec3 &p);
	void SetEye(WaveFront & eye, const vec3 &p);
};
