struct VSInput
{
	float4 Position : POSITION;
	float4 UV : TEXCOORD0;
};

#include "param.hlsl"
#include "complex.hlsl"
#include "PS_VS.hlsl"

struct ComputeParameters
{
	float intervalx;
	float intervaly;
	float lambdaR;
	float lambdaG;

	float lambdaB;
	float distance0;
	float gausssigma;
	float gausssigma2;

	float glareintensity;
	float glarelambdasamplenum;
	float threshold;
	float raylength;

	float posx;
	float posy;
	float range;
	float lensradius;

	float focus0;
	int N;
	float r;
	int N1;

	float3 baseColor;
	float r1;

	//float ghostScale;
	float rotAngle;
	float elapsedTime;
	float screenWidth;
	float screenHeight;

	float minColOfDustTex;
	float ghostScale;
};

struct RandomTbl {
	float4 data[RANDNUM/4];
};

struct GhostRotateMatrix {
	float4 elem;//x 00 y 01 z 10 w 11
};

ConstantBuffer<ComputeParameters> computeConstants : register(b0);
ConstantBuffer<RandomTbl> randomTbl : register(b1);
ConstantBuffer<GhostRotateMatrix> ghostGotateMat : register(b2);

Texture2D<float4> sourceImageR : register(t0);
Texture2D<float4> sourceImageI : register(t1);
RWTexture2D<float4> destinationImageR : register(u0);
RWTexture2D<float4> destinationImageI : register(u1);
RWTexture2D<float4> destinationImageR1 : register(u2);
RWTexture2D<float4> destinationImageI1 : register(u3);
RWByteAddressBuffer randomTblIndex : register(u4);

StructuredBuffer<float4> scaleShiftTbl : register(t2);
StructuredBuffer<float4> colorTbl : register(t3);

RWStructuredBuffer<float4> RWscaleShiftTbl : register(u5);
RWStructuredBuffer<float4> RWcolorTbl : register(u6);

SamplerState CSimageSamplerBILINEAR_WRAP : register(s0);
SamplerState CSimageSamplerBILINEAR_CLAMP : register(s1);

float2 clampF2(float2 value)
{
	return float2(clamp(value.x, 0, WIDTH), clamp(value.y, 0, HEIGHT));
}

#include "CMF.hlsl"
#include "FFT.hlsl"
#include "normalize.hlsl"
#include "asm.hlsl"
#include "wavefront.hlsl"
#include "standard.hlsl"
#include "glare.hlsl"
#include "radialblur.hlsl"
#include "quadratic.hlsl"
#include "lensFlare.hlsl"

