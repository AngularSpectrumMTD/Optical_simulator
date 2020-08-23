#include"..\include\WaveFront.h"


int main()
{
    WaveFront input(512,512,1e-6,1e-6,630e-9);
    double radius = input.GetWidth()/8;

    input.GenerateCirc(radius);
    input.SaveBmp("aperture AMPLITUDE.bmp", AMPLITUDE);

    double distance = 0.5e-3;
    distance *= 2;
    //distance *= 3;
    //distance *= 4;

    input.AsmProp(distance);
    WaveFront save(input);
    save = input;
    save.Normalize();
    save.SaveBmp("propagated distribution AMPLITUDE.bmp",AMPLITUDE);

    mat3 rot = mat3::rotationY(60 * DEG);//60“xyŽ²Žü‚è‚É‰ñ“]//

    WaveFront output(input);
    output.SetNormal(rot * input.GetNormal());
    output.TiltedAsmProp(input, BICUBIC);
    //input.TiltedAsmProp(output, BICUBIC);
    output.Normalize();
    output.SaveBmp("rotated distribution AMPLITUDE.bmp", AMPLITUDE);

    output.dispTotalTime();
}