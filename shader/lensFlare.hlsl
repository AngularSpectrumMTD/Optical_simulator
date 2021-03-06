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
		return;
	}

	{//ゴースト
		float randomValue = randomTbl.data[index / 4][index % 4];// -1 光源と点対象 1 光源の位置
		float2 dir = pos - center;

		//スケール
		const float weight = (abs(randomValue) * 100) % 10 / 10;
		//const float scaleRadius = computeConstants.r + 0.5;
		const float scaleRadius =1.0;
		const float scalePos = length(float2(abs(dir.x) + 1, abs(dir.y) + 1));
		const float scaleJitter = weight + 0.1;
		const float scaleAdjust = 0.1 * (abs(randomValue) + 0.1);
		float2 scalingWeight = scaleRadius * scalePos * scaleJitter * scaleAdjust;
		scalingWeight.x *= 2 * computeConstants.screenHeight / computeConstants.screenWidth;
		scalingWeight *= randomValue > 0 ? 1 : -1;

		//色
		float2 rr = float2(randomValue, randomValue + 1);
		float perlinNoiseR = perlinNoise(float2(randomValue, randomValue + 100));
		float2 gg = rr + 1.xx;
		float2 bb = rr + 2.xx;
		float amplitudeScale = 0.5 - length(dir);
		amplitudeScale = 1 / length(scalingWeight);
		const float ratioBurstAndGhost = 0.01f;
		float standardIntensity = computeConstants.glareintensity * ratioBurstAndGhost;
		float3 colorWeight = amplitudeScale * computeConstants.baseColor* standardIntensity * float3(perlinNoiseR, perlinNoise(gg), perlinNoise(bb));

		//位置
		float2 dirJitter = dir + computeConstants.jitter * weight * float2(dir.y, -dir.x) * ((10 * weight > 5) ? 1 : -1);
		float2 texSize;
		destinationImageR.GetDimensions(texSize.x, texSize.y);
		dirJitter.y *= texSize.x / texSize.y;
		float2 ghostPos = center + randomValue * dirJitter;

		RWscaleShiftTbl[index].xy = scalingWeight * computeConstants.ghostScale;
		RWscaleShiftTbl[index].zw = ghostPos;
		RWcolorTbl[index] = float4(colorWeight, 1);
	}
}