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

	if (val < 1)
	{
		/*inputR = inputR * 10.0f;
		inputI = inputI * 10.0f;*/
		inputR = inputR * computeConstants.glareintensity;
		inputI = inputI * computeConstants.glareintensity;
	}
	else if (val < 0.1)
	{
		destinationImageR[index] = float4(0,0,0, 1.0);
		destinationImageI[index] = float4(0,0,0, 1.0);
	}

	if ((inputR.r + inputR.g + inputR.b) / 3.0f >= 1.0)
	{
		inputR = float3(1.0, 1.0, 1.0);
	}

	if ((inputI.r + inputI.g + inputI.b) / 3.0f >= 1.0)
	{
		inputI = float3(1.0, 1.0, 1.0);
	}

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
	return lambda / (lambdamax -lambdamin);
	return 1;

	float red = 800;
	float green = 600;
	float blue = 400;
	float sigma = 16 * 16;
	//等比混合
	return 0.33 * (exp(-(lambda - red) * (lambda - red)/ sigma) + exp(-(lambda - green) * (lambda - green)/ sigma) + exp(-(lambda - blue) * (lambda - blue)/ sigma));
	//赤優勢
	return 0.66 * exp(-(lambda - red) * (lambda - red)/ sigma) + 0.16 * exp(-(lambda - green) * (lambda - green)/ sigma) + 0.16 * exp(-(lambda - blue) * (lambda - blue)/ sigma);
	//緑優勢
	return 0.16 * exp(-(lambda - red) * (lambda - red)/ sigma) + 0.66  * exp(-(lambda - green) * (lambda - green) / sigma) + 0.16 * exp(-(lambda - blue) * (lambda - blue) / sigma);
	//青優勢
	return 0.16 * exp(-(lambda - red) * (lambda - red)/ sigma) + 0.16 * exp(-(lambda - green) * (lambda - green)/ sigma) + 0.66 * exp(-(lambda - blue) * (lambda - blue)/ sigma);
}

float2 indexfunc(float2 IndexStandard, float lambdamax, float lambda)
{
	float2 uvstandard = IndexStandard - float2(0.5 * WIDTH, 0.5 * HEIGHT);
	float2 uv = uvstandard * lambda / lambdamax;

	return uv + float2(0.5 * WIDTH, 0.5 * HEIGHT);
}

[numthreads(WIDTH, 1, 1)]
void mainSpectrumScaling(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 indexR = dispatchID.xy;

	float maxlambda = 800;
	float minlambda = 360;

	//RGB毎に使用するサンプル数
	float samplenum = computeConstants.glarelambdasamplenum;

	float lambdarange = maxlambda - minlambda;

	float lambdaDelta = lambdarange / samplenum / 3;

	float maxr = sourceImageR[indexR].r;

	float3 result = float3(0.0, 0.0, 0.0);

	//スケーリングした先の対応画素値を参照して加算していく
	for (int i = 0; i < samplenum; i++)
	{
		float lamred = maxlambda - i * lambdaDelta;
		float lamgreen = maxlambda - (i + samplenum) * lambdaDelta;//より小さい波長からスタート
		float lamblue = maxlambda - (i + 2 * samplenum) * lambdaDelta;//同上

		float cr = lambdafunc(minlambda, maxlambda, lamred);
		float cg = lambdafunc(minlambda, maxlambda, lamgreen);
		float cb = lambdafunc(minlambda, maxlambda, lamblue);

		float2 indR = indexfunc(indexR, maxlambda, lamred);
		float2 indG = indexfunc(indexR, maxlambda, lamgreen);
		float2 indB = indexfunc(indexR, maxlambda, lamblue);

		cr *= sourceImageR[indR].r;
		cg *= sourceImageR[indG].g;
		cb *= sourceImageR[indB].b;

		result = result + float3(cr, cg, cb);
	}

	//足しただけ割る
	result /= (WIDTH*HEIGHT * samplenum);

	destinationImageR[indexR] = float4(result, 1.0);
	destinationImageI[indexR] = float4(result, 1.0);
}

///////////////////////////////////////////////////////////////////////////////////グレア旧
//[numthreads(WIDTH, 1, 1)]
//void mainSpectrumScaling(uint3 dispatchID : SV_DispatchThreadID)
//{
//	float2 indexR = dispatchID.xy;
//
//	float ratioRG = computeConstants.lambdaG / computeConstants.lambdaR;
//	float ratioRB = computeConstants.lambdaB / computeConstants.lambdaR;
//
//	//float ratioBG = computeConstants.lambdaB / computeConstants.lambdaG;
//	//float ratioBR = computeConstants.lambdaB / computeConstants.lambdaR;
//
//	//ratioBG = ratioBG * ratioBG;
//	//ratioBR = ratioBR * ratioBR;
//
//	float ratioBG = 1;
//	float ratioBR = 1;
//
//	float2 uvR = indexR - float2(0.5 * WIDTH, 0.5 * HEIGHT);
//	float2 uvG = uvR * ratioRG;
//	float2 uvB = uvR * ratioRB;
//
//	float2 indexG = uvG + float2(0.5 * WIDTH, 0.5 * HEIGHT);
//	float2 indexB = uvB + float2(0.5 * WIDTH, 0.5 * HEIGHT);
//
//	float r = sourceImageR[indexR].r * ratioBR;
//	float g = sourceImageR[indexG].g * ratioBG;
//	float b = sourceImageR[indexB].b;
//
//	destinationImageR[indexR] = float4(r, g, b, 1.0);
//
//	r = sourceImageI[indexR].r * ratioBR;
//	g = sourceImageI[indexG].g * ratioBG;
//	b = sourceImageI[indexB].b;
//
//	destinationImageI[indexR] = float4(r, g, b, 1.0);
//}

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

[numthreads(WIDTH, 1, 1)]
void mainToneMapping(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float ampr = sourceImageR[index].r;
	float ampg = sourceImageR[index].g;
	float ampb = sourceImageR[index].b;

	float toneR = reinhard(ampr);
	float toneG = reinhard(ampg);
	float toneB = reinhard(ampb);

	destinationImageR[index] = float4(toneR, toneG, toneB, 1.0);
}