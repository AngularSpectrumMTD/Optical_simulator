[numthreads(WIDTH, 1, 1)]
void mainREAL(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 inputR = sourceImageR[index].rgb;

	destinationImageR[index] = float4(inputR, 1.0f);
	destinationImageI[index] = float4(inputR, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainIMAGE(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 inputI = sourceImageI[index].rgb;

	destinationImageR[index] = float4(inputI, 1.0f);
	destinationImageI[index] = float4(inputI, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainAMPLITUDE(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 inputR = sourceImageR[index].rgb;
	float3 inputI = sourceImageI[index].rgb;

	float r = inputR.r * inputR.r + inputI.r * inputI.r;
	float g = inputR.g * inputR.g + inputI.g * inputI.g;
	float b = inputR.b * inputR.b + inputI.b * inputI.b;

	float3 col = float3(sqrt(r), sqrt(g), sqrt(b));

	destinationImageR[index] = float4(col, 1.0f);
	destinationImageI[index] = float4(col, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainINTENSITY(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 inputR = sourceImageR[index].rgb;
	float3 inputI = sourceImageI[index].rgb;

	float r = inputR.r * inputR.r + inputI.r * inputI.r;
	float g = inputR.g * inputR.g + inputI.g * inputI.g;
	float b = inputR.b * inputR.b + inputI.b * inputI.b;

	float3 col = float3(r, g, b);

	destinationImageR[index] = float4(col, 1.0f);
	destinationImageI[index] = float4(col, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainPHASE(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 inputR = sourceImageR[index].rgb;
	float3 inputI = sourceImageI[index].rgb;

	float r = atan2(inputI.r, inputR.r);
	float g = atan2(inputI.g, inputR.g);
	float b = atan2(inputI.b, inputR.b);

	float3 col = float3(r, g, b);

	destinationImageR[index] = float4(col, 1.0f);
	destinationImageI[index] = float4(col, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainREAL_IMAGE4Disp(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input = sourceImageR[index].rgb;

	float r = (input.r + 1.0f) / 2.0f;
	float g = (input.g + 1.0f) / 2.0f;
	float b = (input.b + 1.0f) / 2.0f;

	float3 col = float3(r, g, b);

	destinationImageR[index] = float4(col, 1.0f);
	destinationImageI[index] = float4(col, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainAMPLITUDE_INTENSITY4Disp(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input = sourceImageR[index].rgb;

	destinationImageR[index] = float4(input, 1.0f);
	destinationImageI[index] = float4(input, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainPHASE4Disp(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input = sourceImageR[index].rgb;

	float r = input.r / 2.0f / PI + 1.0f / 2.0f;
	float g = input.g / 2.0f / PI + 1.0f / 2.0f;
	float b = input.b / 2.0f / PI + 1.0f / 2.0f;

	float3 col = float3(r, g, b);

	destinationImageR[index] = float4(col, 1.0f);
	destinationImageI[index] = float4(col, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainDivByMaxAMP(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float2 zero = float2(0.0, 0.0);
	//sourceImageRに最大値格納したテクスチャをセットしておく
	float3 color_max = sourceImageR[zero].rgb;

	//destinatinImageR/Iに正規化したいテクスチャをセットしておく
	float3 colorR = destinationImageR[index].rgb;
	float3 colorI = destinationImageI[index].rgb;

	colorR = colorR / color_max;
	colorI = colorI / color_max;

	float col_max_per = color_max.r;
	col_max_per = max(col_max_per, color_max.g);
	col_max_per = max(col_max_per, color_max.b);

	float3 ratio = color_max / col_max_per;

	colorR = colorR * ratio;
	colorI = colorI * ratio;

	destinationImageR1[index] = float4(colorR, 1.0f);
	destinationImageI1[index] = float4(colorI, 1.0f);
}