static float4 dt[1000];


[numthreads(WIDTH, 1, 1)]
void mainRaiseBottomRealImage(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 inputR = sourceImageR[index].rgb;
	float3 inputI = sourceImageI[index].rgb;

	float3 amplitude = float3(sqrt(inputR.r * inputR.r + inputI.r * inputI.r)
		, sqrt(inputR.g * inputR.g + inputI.g * inputI.g)
		, sqrt(inputR.b * inputR.b + inputI.b * inputI.b));

	float val = (amplitude.r + amplitude.g + amplitude.b) / 3.0f;

	//if (val < 1)// && val > 0.0001)
	//{
		/*inputR = inputR * 10.0f;
		inputI = inputI * 10.0f;*/
		inputR = inputR * computeConstants.glareintensity;
		inputI = inputI * computeConstants.glareintensity;
	//}
	//else if (val < 0.2)
	//{
	//	destinationImageR[index] = float4(0,0,0, 1.0);
	//	destinationImageI[index] = float4(0,0,0, 1.0);
	//}

//	if ((inputR.r + inputR.g + inputR.b) / 3.0f < 1.0)
//{
//	inputR = float4(0, 0, 0, 1.0);
//	inputI = float4(0, 0, 0, 1.0);
//}

	//if (val < 0.00002)// && val > 0.0001)
	//{
	//	/*inputR = inputR * 10.0f;
	//	inputI = inputI * 10.0f;*/
	//	
	//	destinationImageR[index] = float4(0, 0, 0, 1.0);
	//	destinationImageI[index] = float4(0, 0, 0, 1.0);
	//}
	//else if (val < 1)
	//{
	//	inputR = inputR * computeConstants.glareintensity;
	//	inputI = inputI * computeConstants.glareintensity;
	//}

	//if ((inputR.r + inputR.g + inputR.b) / 3.0f >= 1.0)
	//{
	//	inputR = float3(1.0, 1.0, 1.0);
	//}

	//if ((inputI.r + inputI.g + inputI.b) / 3.0f >= 1.0)
	//{
	//	inputI = float3(1.0, 1.0, 1.0);
	//}

	destinationImageR[index] = float4(inputR, 1.0);
	destinationImageI[index] = float4(inputI, 1.0);
}

///////////////////////////////////////////////////////////////////////////////////グレア新
//float lambdafunc(float Imax, float lambdamin,float lambda)//普通
//{
//	return Imax * lambdamin / lambda;
//}

float gauss(float x, float x0, float sigma)
{
	return exp(-(x - x0) * (x - x0) / sigma/ sigma);
}

float lambdafunc(float lambdamin, float lambdamax, float lambda)
{
	float V = (lambda / (lambdamax - lambdamin));
	return  V * V;// *V;// *V* V* V* V* V;
}

float2 indexfunc(float2 IndexStandard, float lambdamax, float lambda)
{
	float2 uvstandard = IndexStandard - float2(0.5 * WIDTH, 0.5 * HEIGHT);
	float2 uv = uvstandard * lambda / lambdamax;

	return uv + float2(0.5 * WIDTH, 0.5 * HEIGHT);
}

float2 indexfunc2(float2 IndexStandard, float scale)
{
	float2 uvstandard = IndexStandard - float2(0.5 * WIDTH, 0.5 * HEIGHT);
	float2 uv = uvstandard * scale;//係数が大きいほどちいさく 波長は大きいほど大きく

	return uv + float2(0.5 * WIDTH, 0.5 * HEIGHT);
}

bool isInRange(float2 index)
{
	if (index.x < 0 || WIDTH < index.x || index.y < 0 || HEIGHT < index.y)
	{
		return false;
	}
	else
	{
		return true;
	}
}

//(https://ja.wikipedia.org/wiki/CIE_1931_%E8%89%B2%E7%A9%BA%E9%96%93)
float3 convertRGBtoXYZ(float3 RGBColor)
{
	float3 colXYZ = 1 / 0.17697 *
		float3(
			0.49 * RGBColor.r + 0.31 * RGBColor.g + 0.2 * RGBColor.b,
			0.17697 * RGBColor.r + 0.81240 * RGBColor.g + 0.010630 * RGBColor.b,
			0 * RGBColor.r + 0.01 * RGBColor.g + 0.99 * RGBColor.b
			);
	return colXYZ;
}

float3 convertXYZtoRGB(float3 XYZColor)
{
	float3 colRGB = 
		float3(
			1.9106 * XYZColor.r + -0.5326 * XYZColor.g -0.2883 * XYZColor.b,
			-0.9843 * XYZColor.r + 1.9984 * XYZColor.g -0.0283 * XYZColor.b,
			0.0584 * XYZColor.r -0.1185* XYZColor.g + 0.8985 * XYZColor.b
			);
	return colRGB;
}

float3 lambdalRGB(float lambda) {
	float3 colRGB;

	if (lambda < 350.0)
		colRGB = float3(0.5, 0.0, 1.0);
	else if ((lambda >= 350.0) && (lambda < 440.0))
		colRGB = float3((440.0 - lambda) / 90.0, 0.0, 1.0);
	else if ((lambda >= 440.0) && (lambda <= 490.0))
		colRGB = float3(0.0, (lambda - 440.0) / 50.0, 1.0);
	else if ((lambda >= 490.0) && (lambda < 510.0))
		colRGB = float3(0.0, 1.0, (-(lambda - 510.0)) / 20.0);
	else if ((lambda >= 510.0) && (lambda < 580.0))
		colRGB = float3((lambda - 510.0) / 70.0, 1.0, 0.0);
	else if ((lambda >= 580.0) && (lambda < 645.0))
		colRGB = float3(1.0, (-(lambda - 645.0)) / 65.0, 0.0);
	else
		colRGB = float3(1.0, 0.0, 0.0);

	if (lambda < 350.0)
		colRGB *= 0.3;
	else if ((lambda >= 350.0) && (lambda < 420.0))
		colRGB *= 0.3 + (0.7 * ((lambda - 350.0) / 70.0));
	else if ((lambda >= 420.0) && (lambda <= 700.0))
		colRGB *= 1.0;
	else if ((lambda > 700.0) && (lambda <= 780.0))
		colRGB *= 0.3 + (0.7 * ((780.0 - lambda) / 80.0));
	else
		colRGB *= 0.3;

	return colRGB;
}

bool Clamped(float2 floatIndex) {
	return floatIndex.x < -0 || floatIndex.x > 1 || floatIndex.y < -0 || floatIndex.y > 1;
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainSpectrumScaling(uint3 dispatchID : SV_DispatchThreadID)
{
	const float starburst_resolution = 1.0f;
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 pos = index / size;
	float2 uv = pos.xy / starburst_resolution - 0.5;
	float d = length(uv) * 2;

	// -ve violet, +v reds
	const float scale = 0.50f;

	const float lambdaRed = 700;
	const float lambdaVio = 380;

	float3 result = 0.f;
	int num_steps = computeConstants.glarelambdasamplenum * 3;
	for (int i = 0; i <= num_steps; ++i) {
		float n = (float)i / (float)num_steps;

		float2 scaled_uv = uv * lerp(1.f + scale, 1.f, n) + 0.5;

		bool clamped = Clamped(scaled_uv);

		float r1 = sourceImageR.SampleLevel(CSimageSamplerBILINEAR_CLAMP, scaled_uv, 0).r * !clamped;

		float starburst = r1 * lerp(0.0f, 25.f, d);

		float lambda = lerp(lambdaVio, lambdaRed, n);

		//float3 rgb = convertToMCF(lambda);
		float3 rgb = lambdalRGB(lambda);

		float cr = 1 / lambda;
		cr *= cr;//対象は強度のため2乗(トーンマッピングを無視するならここは cr = (lambdaRed / lambda)^2)

		rgb = lerp(1.f, rgb, 0.75f);

		result += (cr * r1 * rgb);
	}

	result /= (float)num_steps;

	//result = convertXYZtoRGB(result);

	destinationImageR[index] = float4(result, 1.0);
	destinationImageI[index] = float4(result, 1.0);
}

float2 Rotate(float2 p, float a) {
	float x = p.x;
	float y = p.y;

	float cosa = cos(a);
	float sina = sin(a);

	float x1 = x * cosa - y * sina;
	float y1 = y * cosa + x * sina;

	return float2(x1, y1);
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainSpectrumFilterling(uint3 dispatchID : SV_DispatchThreadID)
{
	const float starburst_resolution = 1.0f;
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 pos = index / size;
	float2 uv = pos.xy / starburst_resolution - 0.5;

	float3 result = 0.f;

	const float lambdaRed = 700;
	const float lambdaVio = 380;

	int num_steps = computeConstants.glarelambdasamplenum * 3;;
	for (int i = 0; i <= num_steps; ++i) {
		float n = (float)i / (float)num_steps;
		float phi = n * 2 * PI * 2.f;

	    //float2 spin = float2(cos(phi), sin(phi)) * n * 0.002f;
		float2 spin = float2(cos(phi), sin(phi)) * n * 0.01f;//maybe include blur effect
		//0中心に回転して0.5加算で01空間へ
		float2 rotated_uv = Rotate(uv + spin, n * 0.05f) + 0.5f;

		bool clamped = Clamped(rotated_uv);

		float3 starburst = sourceImageR.SampleLevel(CSimageSamplerBILINEAR_CLAMP, rotated_uv, 0).rgb * !clamped;

		float lambda = lerp(lambdaVio, lambdaRed, (i % 80) / 80.f);
		float3 rgb = lambdalRGB(lambda);
		rgb = lerp(rgb, 1.f, 0.5f);

		result += starburst * rgb;
	}

	result /= (float)num_steps;

	destinationImageR[index] = float4(result, 1.0);
	//destinationImageR[index] = sourceImageR[index];
}

[numthreads(WIDTH, 1, 1)]
void mainSpectrumAmplitudeAdjustment(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float ampr = sourceImageR[index].r;
	float ampg = sourceImageR[index].g;
	float ampb = sourceImageR[index].b;

	float ratioBG = computeConstants.lambdaB / computeConstants.lambdaG;
	float ratioBR = computeConstants.lambdaB / computeConstants.lambdaR;

	//ratioBG = ratioBG * ratioBG;
	//ratioBR = ratioBR * ratioBR;

	ampr = ampr * ratioBR;
	ampg = ampg * ratioBG;

	destinationImageR[index] = float4(ampr, ampg, ampb, 1.0);
}

float reinhard(float x)
{
	return x / (1.0f + x);
}

float3 ACESFilm(float3 x) {
	float a = 2.51f;
	float b = 0.03f;
	float c = 2.43f;
	float d = 0.59f;
	float e = 0.14f;
	return saturate((x * (a * x + b)) / (x * (c * x + d) + e));
}

[numthreads(WIDTH, 1, 1)]
void mainToneMapping(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	destinationImageR[index] = float4(ACESFilm(sourceImageR[index].rgb), 1.0);
}