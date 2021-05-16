struct PSInput
{
	float4 Position : SV_POSITION;
	float4 UV : TEXCOORD0;
};

struct SceneParameters
{
	float4x4 proj;
	float r;
	float screenWidth;
	float screenHeight;
	float kerare;
	float disk;
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

float normalDistribution(float x , float u, float sig2)
{
	float X = -(x - u) * (x - u) / (2 * sig2);
	return exp(X);
}

float normalDistribution2D(float x, float u_x, float sig2_x, float y, float u_y, float sig2_y)
{
	return normalDistribution(x, u_x, sig2_x) * normalDistribution(y, u_y, sig2_y);
}

float fade_aperture_edge(float radius, float fade, float signed_distance) {
	float l = radius;
	float u = radius + fade;
	float s = u - l;
	float c = 1.f - saturate(saturate(signed_distance - l) / s);
	return smoothstep(0, 1, c);
}

//float4 mainPSLensFlare(PSInput In) : SV_TARGET
//{
//	float4 col = 0.xxxx;
//
//	float2 currentUV = In.UV.xy;
//
//	for (int i = 0; i <= GHOSTCOUNT; i++)
//	{
//		float4 scaleShift = GraphicsScaleShiftTbl[i];
//
//		float4 colWeight = GraphicsScolorTbl[i];
//
//		float2 scale = scaleShift.xy;
//		float2 shift = scaleShift.zw;
//
//		float2 uv = currentUV;
//
//		uv -= shift;
//		uv /= scale;
//
//		uv = (uv + float2(1, 1)) * 0.5;
//
//		if (uv.x >= 0 && uv.x <= 1 && uv.y >= 0 && uv.y <= 1)
//		{
//			float2 ghostUV = GraphicsScaleShiftTbl[GHOSTCOUNT].zw - 0.5.xx + uv - 0.5.xx;
//
//			const float lengthUV = length(ghostUV);
//
//			float fade = 0.2;
//			float lens_distance = length(ghostUV * (5 + sceneConstants.r));
//			float sun_disk = 1 - saturate((lens_distance - 1.f + fade) / fade);
//			sun_disk = smoothstep(0, 1, sun_disk);
//			sun_disk *= lerp(0.5, 1, saturate(lens_distance));
//			sun_disk /= length(scale);
//
//			float kerarePerGhost = 1;
//
//			col += colWeight * ((i == GHOSTCOUNT) ? burstImage.Sample(imageSampler, uv) :
//				(ghostImage.Sample(imageSampler, uv)
//					* kerarePerGhost));
//		}
//	}
//
//	col /= 1.0f * (GHOSTCOUNT + 1);
//
//	return col;
//}

float kerareMask(float fade, float2 ghostUV, float offset, float2 scale)
{
	float lens_distance = length(ghostUV * (sceneConstants.disk + offset +  1 - sceneConstants.r));
	float sun_disk = 1 - saturate((lens_distance - 1.f + fade) / fade);
	sun_disk = smoothstep(0, 1, sun_disk);
	sun_disk *= lerp(0.5, 1, saturate(lens_distance));
	//sun_disk /= length(scale);
	return sun_disk;
}

float4 compute(float2 currentUV)
{
	float4 col = 0.xxxx;

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
		float2 uvSave = uv;

		float2 direction = GraphicsScaleShiftTbl[GHOSTCOUNT].zw - 0.5.xx;

		float aspect = sceneConstants.screenHeight / sceneConstants.screenWidth;

		direction *= 2 * aspect;

		direction = normalize(direction);
		float2 vertDirection = float2(-direction.y, direction.x);
		if (i < GHOSTCOUNT)
		{
			float2 tmp = shift - 0.5.xx;
			tmp.x *= aspect;
			float Len = length(tmp);
			float scaleFactor = pow(1 + Len, 4);
			uv -= 0.5.xx;
			float rotatedU = dot(uv, direction);
			float rotatedV = dot(uv, vertDirection) * scaleFactor;
			uv = direction * rotatedU + vertDirection * rotatedV;
			uv += 0.5.xx;

			if (uv.x < 0 || uv.x > 1 || uv.y < 0 || uv.y > 1)
			{
				continue;
			}
		}

		if (uv.x >= 0 && uv.x <= 1 && uv.y >= 0 && uv.y <= 1)
		{
			float2 ghostUV = GraphicsScaleShiftTbl[GHOSTCOUNT].zw - 0.5.xx + uvSave - 0.5.xx;

			const float lengthUV = length(ghostUV);

			float fade = 0.2;

			float sun_disk = kerareMask(fade, ghostUV, 0, scale);

			float sun_disk2 = kerareMask(fade, ghostUV, 0.05, scale);

			float sun_disk3 = abs(sun_disk - sun_disk2) + 1;

			//ƒPƒ‰ƒŒ‚³‚¹‚é‚©
			float kerarePerGhost = sun_disk;
			kerarePerGhost = sun_disk3 * sun_disk;
			kerarePerGhost = (sceneConstants.kerare > 0) ? kerarePerGhost : 1;

			col += colWeight * ((i == GHOSTCOUNT) ? burstImage.Sample(imageSampler, uv) :
				(ghostImage.Sample(imageSampler, uv)
					* kerarePerGhost));

			//col += colWeight * ((i == GHOSTCOUNT) ? burstImage.Sample(imageSampler, uv) :
			//	(float4(kerarePerGhost.xxx, 1)));
		}
	}

	col /= 1.0f * (GHOSTCOUNT + 1);

	return col;
}

float4 mainPSLensFlare(PSInput In) : SV_TARGET
{
	float2 currentUV = In.UV.xy;

	float4 col = compute(currentUV);

	return col;
}

float4 mainPSLensFlareAdd(PSInput In) : SV_TARGET
{
	float2 currentUV = In.UV.xy;

	float4 col = compute(currentUV);

	return col + float4(backImage.Sample(imageSampler, currentUV).rgb, 0);
}