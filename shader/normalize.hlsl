[numthreads(1, 1, 1)]
void mainMAXMIN(uint3 dispatchID : SV_DispatchThreadID)
{
	uint2 index = dispatchID.xy;
	float3 color_max = sourceImageR[float2(0, 0)].xyz;
	float3 color_min = color_max;

	int i, j;
	for (i = 0; i < WIDTH; i++)
	{
		for (j = 0; j < HEIGHT; j++)
		{
			uint2 indexx = uint2(i, j);
			float3 color = sourceImageR[indexx].xyz;

			color_max = max(color_max, color);
			color_min = min(color_min, color);
		}
	}

	destinationImageR[index] = float4(color_max, 1.0f);
	destinationImageI[index] = float4(color_min, 1.0f);
}

[numthreads(1, 1, 1)]
void mainMAXMINfirst(uint3 dispatchID : SV_DispatchThreadID)
{
	uint2 index = dispatchID.xy;
	float3 color_max = sourceImageR[float2(0, 0)].xyz;
	float3 color_min = color_max;

	int i, j;
	for (i = 0; i < HEIGHT; i++)
	{
		uint2 indexx = uint2(index.x, i);
		float3 color = sourceImageR[indexx].xyz;

		color_max = max(color_max, color);
		color_min = min(color_min, color);
	}

	destinationImageR[float2(index.x, 0)] = float4(color_max, 1.0f);
	destinationImageR[float2(index.x, 1)] = float4(color_min, 1.0f);
}

[numthreads(1, 1, 1)]
void mainMAXMINsecond(uint3 dispatchID : SV_DispatchThreadID)
{
	uint2 index = dispatchID.xy;
	float3 color_max = sourceImageR[float2(0, 0)].xyz;
	float3 color_min = color_max;

	int i, j;
	for (i = 0; i < HEIGHT; i++)
	{
		uint2 index1 = uint2(i, 0);
		uint2 index2 = uint2(i, 1);
		float3 colormax = sourceImageR[index1].xyz;
		float3 colormin = sourceImageR[index2].xyz;

		color_max = max(color_max, colormax);
		color_min = min(color_min, colormin);
	}

	destinationImageR[float2(0, 0)] = float4(color_max, 1.0f);
	destinationImageI[float2(0, 0)] = float4(color_min, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainNORMALIZE(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float2 zero = float2(0.0, 0.0);
	//sourceImageRに最大値格納したテクスチャをセットしておく
	float3 color_max = sourceImageR[zero].rgb;
	float3 color_min = sourceImageI[zero].rgb;

	//destinatinImageIに正規化したいテクスチャをセットしておく
	float3 color = destinationImageI[index].rgb;

	color = (color - color_min) / (color_max - color_min);

	/*float col_max_per = color_max.r;
	col_max_per = max(col_max_per, color_max.g);
	col_max_per = max(col_max_per, color_max.b);

	float3 ratio = color_max / col_max_per;

	color = color * ratio;*/

	//destinationImageRに正規化したテクスチャを出力
	destinationImageR[index] = float4(color, 1.0f);
}