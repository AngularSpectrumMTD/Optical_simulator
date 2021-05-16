[numthreads(GHOSTCOUNT + 1, 1, 1)]
void mainComputeScaleShiftColor(uint3 dispatchID : SV_DispatchThreadID)
{
	int index = dispatchID.x;

	float2 pos = float2(computeConstants.posx, computeConstants.posy);
	float2 center = 0.5.xx;

	if (index == GHOSTCOUNT)
	{//バースト
		float2 scalingWeight = (1.5 - length(pos - center));// *(2 - computeConstants.r);
		scalingWeight.x *= 2 * computeConstants.screenHeight / computeConstants.screenWidth;

		RWscaleShiftTbl[index].xy = scalingWeight * 0.4f;
		RWscaleShiftTbl[index].zw = pos;
		RWcolorTbl[index] = 1.xxxx;
		//RWcolorTbl[index] = float4(computeConstants.baseColor, 1.0f);
		return;
	}

	{//ゴースト
		float randomValue = randomTbl.data[index / 4][index % 4];// -1 光源と点対象 1 光源の位置
		//randomValue = -1;
		float2 dir = pos - center;

		//スケール
		float weight = (abs(randomValue) * 100) % 10 / 10;
		float2 scalingWeight = 0.1 * (computeConstants.r + 0.5).xx * (length(randomValue * dir + float2(1,1))) * (weight + 0.1);
		scalingWeight.x *= 2 * computeConstants.screenHeight / computeConstants.screenWidth;
		scalingWeight *= randomValue > 0 ? 1 : -1;
		//scalingWeight *= computeConstants.ghostScale;

		//色
		float2 rr = float2(randomValue, randomValue + 1);
		float perlinNoiseR = perlinNoise(float2(randomValue, randomValue + 100));
		float2 gg = rr + 1.xx;
		float2 bb = rr + 2.xx;
		//float amplitudescale = 5 * 1 / (length(dir) + 0.1f);
		float amplitudescale = 0.5 - length(dir);
		amplitudescale = 1 / length(scalingWeight);

		float standardIntensity = computeConstants.glareintensity * 0.001f;

		float3 colorWeight = amplitudescale* computeConstants.baseColor* standardIntensity * float3(perlinNoiseR, perlinNoise(gg), perlinNoise(bb));

		//位置
		float2 move = float2(perlinNoise(rr), perlinNoise(gg));
		float2 dir2 = dir + computeConstants.jitter * weight * float2(dir.y, -dir.x) * ((10 * weight > 5) ? 1 : -1);
		float2 texSize;
		destinationImageR.GetDimensions(texSize.x, texSize.y);
		dir2.y *= texSize.x / texSize.y;
		float2 ghostPos = center + randomValue * dir2;

		RWscaleShiftTbl[index].xy = scalingWeight * computeConstants.ghostScale;
		RWscaleShiftTbl[index].zw = ghostPos;
		RWcolorTbl[index] = float4(colorWeight, 1);
	}
}