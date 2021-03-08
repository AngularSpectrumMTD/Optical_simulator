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
SamplerState imageSampler : register(s0);

PSInput mainVS(VSInput In, uint instanceID : SV_InstanceID)
{
	PSInput result = (PSInput)0;

	result.Position = mul(In.Position, sceneConstants.proj);
	result.UV = In.UV;
	return result;
}

float4 mainPS(PSInput In) : SV_TARGET
{
  return imageTex.Sample(imageSampler, In.UV.xy);
}