[numthreads(WIDTH, 1, 1)]
void swap(uint3 dtid : SV_DispatchThreadID)
{
	int i, j;
	if (WIDTH / 2 <= dtid.x && dtid.x < WIDTH)
	{
		i = dtid.x - WIDTH / 2;
	}
	else
	{
		i = dtid.x + WIDTH / 2;
	}

	if (HEIGHT / 2 <= dtid.y && dtid.y < HEIGHT)
	{
		j = dtid.y - HEIGHT / 2;
	}
	else
	{
		j = dtid.y + HEIGHT / 2;
	}

	float2 index = float2(i, j);

	float3 colorR = sourceImageR[index].xyz;
	destinationImageR[dtid.xy] = float4(colorR, 1);
	float3 colorI = sourceImageI[index].xyz;
	destinationImageI[dtid.xy] = float4(colorI, 1);
}

void ComputeWight(uint passIndex, uint x, out uint2 indices, out float2 weights)
{
	uint regionWidth = 2 << passIndex;
	uint halfregionWidth = regionWidth / 2;

	uint regionStartOffset = x & ~(regionWidth - 1);
	uint halfregionOffset = x & (halfregionWidth - 1);
	uint regionOffset = x & (regionWidth - 1);

	float a = 2.0 * PI * float(regionOffset) / float(regionWidth);
	weights.y = -sin(a);
	weights.x = cos(a);

	indices.x = regionStartOffset + halfregionOffset;
	indices.y = regionStartOffset + halfregionOffset + halfregionWidth;

	if (passIndex == 0)
	{
		indices = reversebits(indices) >> (32 - BUTTERFLY_COUNT) & (LENGTH - 1);
	}
}

groupshared float3 ButterflyArray[4][LENGTH];

void ButterflyWeightPass(uint passIndex, uint x, uint t0, uint t1, out float3 resultR, out float3 resultI)
{
	uint2 Indices;
	float2 Weights;
	ComputeWight(passIndex, x, Indices, Weights);

	float3 inputR1 = ButterflyArray[t0][Indices.x];
	float3 inputI1 = ButterflyArray[t1][Indices.x];

	float3 inputR2 = ButterflyArray[t0][Indices.y];
	float3 inputI2 = ButterflyArray[t1][Indices.y];

#if INVERSE
	resultR = (inputR1 + Weights.x * inputR2 + Weights.y * inputI2) * 0.5;
	resultI = (inputI1 - Weights.y * inputR2 + Weights.x * inputI2) * 0.5;
#else
	resultR = inputR1 + Weights.x * inputR2 - Weights.y * inputI2;
	resultI = inputI1 + Weights.y * inputR2 + Weights.x * inputI2;
#endif
}

[numthreads(LENGTH, 1, 1)]
void mainFFT(uint3 dispatchID : SV_DispatchThreadID)
{
	uint2 position = dispatchID.xy;

#if ROW
	uint2 processPos = position.xy;
	uint2 destPos = processPos;
#else
	uint2 processPos = position.yx;
	uint2 destPos = processPos;
#endif

	float4 inputR = sourceImageR[processPos];
	ButterflyArray[0][position.x].xyz = inputR.xyz;

#if ROW && !INVERSE
	ButterflyArray[1][position.x].xyz = (0.0).xxx;
# else
	ButterflyArray[1][position.x].xyz = sourceImageI[processPos].xyz;
#endif

	uint4 processIndices = uint4(0, 1, 2, 3);

	for (int i = 0; i < BUTTERFLY_COUNT - 1; i++)
	{
		GroupMemoryBarrierWithGroupSync();
		ButterflyWeightPass(i, position.x, processIndices.x, processIndices.y, ButterflyArray[processIndices.z][position.x].xyz, ButterflyArray[processIndices.w][position.x].xyz);
		processIndices.xyzw = processIndices.zwxy;
	}

	GroupMemoryBarrierWithGroupSync();

	float3 outputR, outputI;
	ButterflyWeightPass(BUTTERFLY_COUNT - 1, position.x, processIndices.x, processIndices.y, outputR, outputI);

	destinationImageR[destPos] = float4(outputR.rgb, inputR.a);
	destinationImageI[destPos] = float4(outputI.rgb, inputR.a);
}