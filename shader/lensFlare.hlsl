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
		float weight = 1 + (abs(randomValue) * 100) % 10 / 10;
		float2 scalingWeight = (10 * (computeConstants.r + 0.5).xx + 1.xx) * (1.5 - length(randomValue * dir)) * weight;

		float theta = atan2(dir.y, dir.x);
		//scalingWeight = Rotate(scalingWeight, 180);
		//scalingWeight = Rotate(scalingWeight, 90);
		//scalingWeight = Rotate(scalingWeight, 45);
		//scalingWeight = Rotate(scalingWeight, 0);

		//scalingWeight = Rotate(scalingWeight, theta * 180.0 / 3.14);

		scalingWeight.x *= 2 * computeConstants.screenHeight / computeConstants.screenWidth;
		scalingWeight *= 0.005f;

		//scalingWeight =( (theta * 180 / 3.14) % 360) / 360;

		//scalingWeight = float2(0.5, 0.5);
		//float2 dilation = float2(abs(dir.x), abs(dir.y));
		//scalingWeight *= (1 + 100 * length(dilation) * dilation * abs(randomValue));

		//色
		float2 rr = float2(randomValue, randomValue + 1);
		float perlinNoiseR = perlinNoise(float2(randomValue, randomValue + 100));
		float2 gg = float2(perlinNoiseR + randomValue + 1, randomValue + 2);
		float2 bb = float2(randomValue * perlinNoiseR + 3 + 2 * perlinNoiseR, randomValue + 4 + perlinNoiseR);
		float amplitudescale = 5 * 1 / (length(dir) + 0.1f);
		float3 colorWeight = amplitudescale* computeConstants.baseColor* float3(perlinNoiseR, perlinNoise(gg), perlinNoise(bb)) / (0.1 + length(scalingWeight));

		//位置
		float2 move = float2(perlinNoise(rr), perlinNoise(gg));
		float2 dir2 = dir;// +move;
		float2 texSize;
		float  level;
		destinationImageR.GetDimensions(texSize.x, texSize.y);
		dir2.y *= texSize.x / texSize.y;
		float2 ghostPos = center + randomValue * dir2;

		RWscaleShiftTbl[index].xy = scalingWeight;
		RWscaleShiftTbl[index].zw = ghostPos;
		RWcolorTbl[index] = float4(colorWeight, 1);
	}
}