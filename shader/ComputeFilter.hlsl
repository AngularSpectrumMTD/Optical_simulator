struct VSInput
{
	float4 Position : POSITION;
	float4 UV : TEXCOORD0;
};

#include "param.hlsl"
#include "complex.hlsl"
#include "PS_VS.hlsl"

//struct ComputeParameters
//{
//	float intervalx;
//	float intervaly;
//	float lambdaR;
//	float lambdaG;
//	float lambdaB;
//	float distance0;
//	float gausssigma;
//	float gausssigma2;
//	float glareintensity;
//	float glarelambdasamplenum;
//	float threshold;
//	float raylength;
//	float posx;
//	float posy;
//	float range;
//	float lensradius;
//	float focus0;
//	int N;
//	float r;
//	int N1;
//	float r1;
//	float ghostScale;
//	float rotAngle;
//	float3 baseColor;
//};

//struct ComputeParameters
//{
//	float intervalx : packoffset(c0.x);
//	float intervaly : packoffset(c0.y);
//	float lambdaR : packoffset(c0.z);
//	float lambdaG : packoffset(c0.w);
//	float lambdaB : packoffset(c1.x);
//	float distance0 : packoffset(c1.y);
//	float gausssigma : packoffset(c1.z);
//	float gausssigma2 : packoffset(c1.w);
//	float glareintensity : packoffset(c2.x);
//	float glarelambdasamplenum : packoffset(c2.y);
//	float threshold : packoffset(c2.z);
//	float raylength : packoffset(c2.w);
//	float posx : packoffset(c3.x);
//	float posy : packoffset(c3.y);
//	float range : packoffset(c3.z);
//	float lensradius : packoffset(c3.w);
//	float focus0 : packoffset(c4.x);
//	int N : packoffset(c4.y);
//	float r : packoffset(c4.z);
//	int N1 : packoffset(c4.w);
//	float r1 : packoffset(c5.x);
//	float ghostScale : packoffset(c5.y);
//	float rotAngle : packoffset(c5.z);
//	float3 baseColor:packoffset(c5.w);
//};

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

	float ghostScale;
	float rotAngle;
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

float2 clampF2(float2 value)
{
	return float2(clamp(value.x, 0, WIDTH), clamp(value.y, 0, HEIGHT));
}

#include "FFT.hlsl"
#include "normalize.hlsl"
#include "asm.hlsl"
#include "wavefront.hlsl"
#include "standard.hlsl"
#include "glare.hlsl"
#include "radialblur.hlsl"
#include "quadratic.hlsl"

