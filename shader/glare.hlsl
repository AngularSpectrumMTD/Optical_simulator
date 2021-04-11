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

	if (val < 1)// && val > 0.0001)
	{
		/*inputR = inputR * 10.0f;
		inputI = inputI * 10.0f;*/
		inputR = inputR * computeConstants.glareintensity;
		inputI = inputI * computeConstants.glareintensity;
	}
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

///////////////////////////////////////////////////////////////////////////////////�O���A�V
//float lambdafunc(float Imax, float lambdamin,float lambda)//����
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

float2 indexfunc2(float2 IndexStandard, float lambdamax, float lambda)
{
	float2 uvstandard = IndexStandard - float2(0.5 * WIDTH, 0.5 * HEIGHT);
	float2 uv = uvstandard * lambdamax / lambda;

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

[numthreads(WIDTH, 1, 1)]
void mainSpectrumScaling(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 indexR = dispatchID.xy;

	float maxlambda = 830;
	float minlambda = 390;

	//RGB���Ɏg�p����T���v����
	float samplenum = computeConstants.glarelambdasamplenum * 4;

	float lambdarange = maxlambda - minlambda;

	float3 result = float3(0.0, 0.0, 0.0);

	float2 size = float2(WIDTH, HEIGHT);
	//�X�P�[�����O������̑Ή���f�l���Q�Ƃ��ĉ��Z���Ă���

	float lambdaDelta = lambdarange / samplenum;

	const float representativeLambda = maxlambda;

	for (int i = 0; i < samplenum; i++)
	{
		float lambda = maxlambda - i * lambdaDelta;
		float cr = representativeLambda / lambda;
		cr *= cr;//�Ώۂ͋��x�̂���2��

		//float2 indR = indexfunc(indexR, maxlambda, representativeLambda + minlambda - lambda);
		float2 indR = indexfunc2(indexR, maxlambda, lambda);// lambda�ɂЂꂢ

		float3 S = (isInRange(indR) == true) ? sourceImageR.SampleLevel(CSimageSampler, indR / size, 0).rgb : float3(0,0,0);

		S *= cr;
		float3 xyzfunc = convertToMCF(lambda);

		S *= xyzfunc;
		S = convertXYZtoRGB(S);
		result = result + S;
	}
	//��������������
	result /= (WIDTH*HEIGHT * samplenum);

	destinationImageR[indexR] = float4(result, 1.0);
	destinationImageI[indexR] = float4(result, 1.0);
}

float3 wl2rgbTannenbaum(float w) {
	float3 r;

	if (w < 350.0)
		r = float3(0.5, 0.0, 1.0);
	else if ((w >= 350.0) && (w < 440.0))
		r = float3((440.0 - w) / 90.0, 0.0, 1.0);
	else if ((w >= 440.0) && (w <= 490.0))
		r = float3(0.0, (w - 440.0) / 50.0, 1.0);
	else if ((w >= 490.0) && (w < 510.0))
		r = float3(0.0, 1.0, (-(w - 510.0)) / 20.0);
	else if ((w >= 510.0) && (w < 580.0))
		r = float3((w - 510.0) / 70.0, 1.0, 0.0);
	else if ((w >= 580.0) && (w < 645.0))
		r = float3(1.0, (-(w - 645.0)) / 65.0, 0.0);
	else
		r = float3(1.0, 0.0, 0.0);

	if (w < 350.0)
		r *= 0.3;
	else if ((w >= 350.0) && (w < 420.0))
		r *= 0.3 + (0.7 * ((w - 350.0) / 70.0));
	else if ((w >= 420.0) && (w <= 700.0))
		r *= 1.0;
	else if ((w > 700.0) && (w <= 780.0))
		r *= 0.3 + (0.7 * ((780.0 - w) / 80.0));
	else
		r *= 0.3;

	return r;
}

//[numthreads(WIDTH, 1, 1)]
//void mainSpectrumScaling(uint3 dispatchID : SV_DispatchThreadID)
//{
//	float2 indexR = dispatchID.xy;
//	float2 size = float2(WIDTH, HEIGHT);
//	float2 uv = indexR / size - 0.5;
//	float d = length(uv) * 2;
//
//	float scale1 = 0.50f;
//	float scale2 = -0.75f;
//	float fft_scale = 0.00001f;
//
//	float3 result = 0.f;
//	int num_steps = 256;
//
//	//�X�P�[�����O������̑Ή���f�l���Q�Ƃ��ĉ��Z���Ă���
//	for (int i = 0; i < num_steps; ++i)
//	{
//		float n = (float)i / (float)num_steps;
//
//		float2 scaled_uv1 = uv * lerp(1.f + scale1, 1.f, n);
//		float2 scaled_uv2 = uv * lerp(1.f + scale2, 1.f, n);
//
//		float r1 = sourceImageR[scaled_uv1 * size].r;
//		float i1 = sourceImageI[scaled_uv1 * size].r;
//
//		float2 p1 = float2(r1, i1);
//
//		float starburst = pow(length(p1), 2.f) * fft_scale * lerp(0.0f, 25.f, d);
//
//		float lambda = lerp(380.f, 700.f, n);
//		float3 rgb = wl2rgbTannenbaum(lambda);
//		rgb = lerp(1.f, rgb, 0.75f);
//
//		result += (starburst + rgb * 0.25f);
//	}
//
//	result /= (float)num_steps;
//	destinationImageR[indexR] = float4(result, 1);
//}

///////////////////////////////////////////////////////////////////////////////////�O���A��
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