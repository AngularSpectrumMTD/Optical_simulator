struct PSInput
{
	float4 Position : SV_POSITION;
	float4 UV : TEXCOORD0;
};

struct SceneParameters
{
	float4x4 proj;
};

ConstantBuffer<SceneParameters> sceneConstants : register(b0);
Texture2D imageTex : register(t0);
Texture2D ghostImage : register(t1);
Texture2D burstImage : register(t2);
Texture2D backImage : register(t3);
StructuredBuffer<float4> GraphicsScaleShiftTbl : register(t4);
StructuredBuffer<float4> GraphicsScolorTbl : register(t5);
SamplerState imageSampler : register(s0);

PSInput mainVS(VSInput In, uint instanceID : SV_InstanceID)
{
	PSInput result = (PSInput)0;

	result.Position = mul(In.Position, sceneConstants.proj);
	result.UV = In.UV;
	return result;
}

float tonemap(float x)
{
	return x / (1 + x);
}

float3 tonemap(float3 v)
{
	return float3(tonemap(v.x), tonemap(v.y), tonemap(v.z));
}

float4 mainPS(PSInput In) : SV_TARGET
{
  return imageTex.Sample(imageSampler, In.UV.xy);
}

float4 mainPSLensFlare(PSInput In) : SV_TARGET
{
	float4 col = 0.xxxx;

	float2 currentUV = In.UV.xy;

	for (int i = 0; i <= GHOSTCOUNT; i++)
	{
		float4 scaleShift = GraphicsScaleShiftTbl[i];

		float4 colWeight = GraphicsScolorTbl[i];

		float2 scale = scaleShift.xy;
		float2 shift = scaleShift.zw;

		float2 uv = currentUV;

		uv -= shift;
		uv /= scale;

		uv = (uv + float2(1, 1)) * 0.5;

		if (uv.x >= 0 && uv.x <= 1 && uv.y >= 0 && uv.y <= 1)
		{
			col += colWeight * ((i == GHOSTCOUNT) ? burstImage.Sample(imageSampler, uv) : ghostImage.Sample(imageSampler, uv));
		}
	}

	col /= 1.0f * (GHOSTCOUNT + 1);

	return col;
}

float4 mainPSLensFlareAdd(PSInput In) : SV_TARGET
{
	float4 col = 0.xxxx;

	float2 currentUV = In.UV.xy;

	for (int i = 0; i <= GHOSTCOUNT; i++)
	{
		float4 scaleShift = GraphicsScaleShiftTbl[i];

		float4 colWeight = GraphicsScolorTbl[i];

		float2 scale = scaleShift.xy;
		float2 shift = scaleShift.zw;

		float2 uv = currentUV;

		uv -= shift;
		uv /= scale;

		uv = (uv + float2(1, 1)) * 0.5;

		if (uv.x >= 0 && uv.x <= 1 && uv.y >= 0 && uv.y <= 1)
		{
			col += colWeight * ((i == GHOSTCOUNT) ? burstImage.Sample(imageSampler, uv) : ghostImage.Sample(imageSampler, uv));
		}
	}

	col /= 1.0f * (GHOSTCOUNT + 1);

	return col + float4(backImage.Sample(imageSampler, currentUV).rgb, 0);
}