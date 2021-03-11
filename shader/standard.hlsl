
float2 random2(float2 st) {
	st = float2(dot(st, float2(127.1, 311.7)),
		dot(st, float2(269.5, 183.3)));
	return -1.0 + 2.0 * frac(sin(st) * 43758.5453123);
}

float perlinNoise(float2 st)
{
	float2 p = floor(st);
	float2 f = frac(st);
	float2 u = f * f * (3.0 - 2.0 * f);

	float v00 = random2(p + float2(0, 0)).x;
	float v10 = random2(p + float2(1, 0)).x;
	float v01 = random2(p + float2(0, 1)).x;
	float v11 = random2(p + float2(1, 1)).x;

	return lerp(lerp(dot(v00, f - float2(0, 0)), dot(v10, f - float2(1, 0)), u.x),
		lerp(dot(v01, f - float2(0, 1)), dot(v11, f - float2(1, 1)), u.x),
		u.y) + 0.5f;
}

[numthreads(16, 16, 1)]
void mainSepia(uint3 dtid : SV_DispatchThreadID)
{
	if (dtid.x < WIDTH && dtid.y < HEIGHT)
	{
		float3x3 toSepia = float3x3(
			0.393, 0.349, 0.272,
			0.769, 0.686, 0.534,
			0.189, 0.168, 0.131);
		float3 color = mul(sourceImageR[dtid.xy].xyz, toSepia);
		destinationImageR[dtid.xy] = float4(color, 1);
	}
}

[numthreads(1, 1, 1)]
void mainSobel(uint3 dtid : SV_DispatchThreadID)
{
	if (dtid.x < WIDTH && dtid.y < HEIGHT)
	{
		int k = 0;
		float3 pixels[9];
		for (int y = -1; y <= 1; ++y)
		{
			for (int x = -1; x <= 1; ++x)
			{
				float2 index = dtid.xy;
				index += float2(x, y);

				pixels[k] = sourceImageR[index].xyz;
				k++;
			}
		}

		float3 sobelH, sobelV;
		sobelH = pixels[0] * -1 + pixels[2] * 1
			+ pixels[3] * -2 + pixels[5] * 2
			+ pixels[6] * -1 + pixels[8] * 1;
		sobelV = pixels[0] * -1 + pixels[1] * -2 + pixels[2] * -1
			+ pixels[6] * 1 + pixels[7] * 2 + pixels[8] * 1;

		float4 color = float4(sqrt(sobelV * sobelV + sobelH * sobelH), 1);
		destinationImageR[dtid.xy] = color;
	}
}

[numthreads(WIDTH, 1, 1)]
void mainMULTIPLY(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float3 color_real1 = sourceImageR[index].rgb;//元1 実部
	float3 color_image1 = sourceImageI[index].rgb;//元1 虚部

	float3 color_real2 = destinationImageR[index].rgb;//元2 実部
	float3 color_image2 = destinationImageI[index].rgb;//元2 虚部

	//各色からの複素数生成
	float2 compR1 = float2(color_real1.r, color_image1.r);
	float2 compG1 = float2(color_real1.g, color_image1.g);
	float2 compB1 = float2(color_real1.b, color_image1.b);

	float2 compR2 = float2(color_real2.r, color_image2.r);
	float2 compG2 = float2(color_real2.g, color_image2.g);
	float2 compB2 = float2(color_real2.b, color_image2.b);

	//各色毎に乗算
	float2 mulcompR = complex_mul(compR1, compR2);
	float2 mulcompG = complex_mul(compG1, compG2);
	float2 mulcompB = complex_mul(compB1, compB2);

	//代入用の複素数作成
	float3 colorREAL = float3(mulcompR.x, mulcompG.x, mulcompB.x);
	float3 colorIMAGE = float3(mulcompR.y, mulcompG.y, mulcompB.y);

	destinationImageR1[index] = float4(colorREAL, 1.0f);
	destinationImageI1[index] = float4(colorIMAGE, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainColorScaling(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float3 color = sourceImageR[index].rgb;

	color = color / WIDTH / HEIGHT;

	destinationImageR[index] = float4(color, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainBinaryThreshold(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input = sourceImageR[index].rgb;

	float3 col = float3(0.0, 0.0, 0.0);

	if ((input.r + input.g + input.b) / 3.0 >= computeConstants.threshold)
	{
		col = float3(1.0, 1.0, 1.0);
	}

	destinationImageR[index] = float4(col, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainAdd(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input1 = sourceImageR[index].rgb;
	float3 input2 = sourceImageI[index].rgb;

	float3 col = input1 + input2;

	destinationImageR[index] = float4(col, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainDrawGaussian(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float x = index.x - WIDTH / 2.0f;
	float y = index.y - HEIGHT / 2.0f;

	float s = computeConstants.gausssigma;

	float val = exp(-(x * x + y * y) / 2 / s / s) / sqrt(2 * PI) / s;//必ずサイズで割ること

	float3 color_real = float3(val, val, val);

	destinationImageR[index] = float4(color_real, 1.0f);
	destinationImageI[index] = float4(0.0, 0.0, 0.0, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainDrawGaussianNoNormalize(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float x = index.x - WIDTH / 2.0f;
	float y = index.y - HEIGHT / 2.0f;

	float s = computeConstants.gausssigma;

	float val = exp(-(x * x + y * y) / 2 / s / s);//必ずサイズで割ること

	float3 color_real = float3(val, val, val);

	destinationImageR[index] = float4(color_real, 1.0f);
	destinationImageI[index] = float4(0.0, 0.0, 0.0, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainDrawGaussianNoNormalize2(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float x = index.x - WIDTH / 2.0f;
	float y = index.y - HEIGHT / 2.0f;

	

	uint table_index = randomTblIndex.Load(0);

	float randomValue = randomTbl.data[table_index / 4][table_index % 4];
	//randomValue = perlinNoise(float2(randomValue, randomValue + 100));
	randomValue = abs(frac(randomValue));
	randomValue = frac(randomValue);
	randomValue *= randomValue;

	float s = computeConstants.gausssigma2;
	s *= randomValue;

	float length = randomValue * randomValue * x * x + randomValue * randomValue * y * y;
	length = x * x + y * y;

	float val = exp(-length / 2 / s / s);//必ずサイズで割ること

	float3 color_real = float3(val, val, val);

	destinationImageR[index] = float4(color_real, 1.0f);
	destinationImageI[index] = float4(0.0, 0.0, 0.0, 1.0f);
}

[numthreads(WIDTH, 1, 1)]
void mainCopy(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	destinationImageR[index] = sourceImageR[index];
}

[numthreads(WIDTH, 1, 1)]
void mainClear(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	destinationImageR[index] = float4(0.0, 0.0, 0.0, 1.0);
}

[numthreads(WIDTH, 1, 1)]
void mainRaiseBottom(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input = sourceImageR[index].rgb;

	if ((input.r + input.g + input.b) / 3.0f < 0.9)
	{
		input = input * 40.0f;
	}

	if ((input.r + input.g + input.b) / 3.0f >= 1.0)
	{
		input = float3(1.0, 1.0, 1.0);
	}

	destinationImageR[index] = float4(input, 1.0);
}

[numthreads(WIDTH, 1, 1)]//computeConstants.N角形
void mainDrawPolygon(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 uv = index / size - float2(0.5, 0.5);

	//判定したい点の位置
	float pos = length(uv);

	//判定したい点のなす角
	float rad = atan2(uv.x, uv.y) + 2.0 * PI;
	rad = rad % (2.0 * PI / computeConstants.N);

	float r_circ = 0.5 * computeConstants.r;

	//半径r_circの円に内接する正多角形の辺の位置
	float r_polygon = cos(PI / computeConstants.N) / cos(PI / computeConstants.N - rad);
	r_polygon *= r_circ;

	//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0～1)
	//lerp(x,y,s) = x + s(y - x)
	float s = computeConstants.r;
	float r_aperture = lerp(r_polygon, r_circ, s * s * s);

	//判定される点が描きたい多角形の内外かを判定
	float col = step(pos, r_aperture);

	destinationImageR[index] = float4(col, col, col, 1.0);
}

[numthreads(WIDTH, 1, 1)]//computeConstants.N角形
void mainDrawPolygonFixScale(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 uv = index / size - float2(0.5, 0.5);

	//判定したい点の位置
	float pos = length(uv);

	//判定したい点のなす角
	float rad = atan2(uv.x, uv.y) + 2.0 * PI + computeConstants.rotAngle * PI / 180.0f;//ここで角度足したら回る
	rad = rad % (2.0 * PI / computeConstants.N);

	//float r_circ = 0.2;
	float r_circ = 0.4;//max

	//半径r_circの円に内接する正多角形の辺の位置
	float r_polygon = cos(PI / computeConstants.N) / cos(PI / computeConstants.N - rad);
	r_polygon *= r_circ;

	//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0～1)
	//lerp(x,y,s) = x + s(y - x)
	float s = computeConstants.r;
	//float r_aperture = lerp(r_polygon, r_circ, s * s * s);

	float ratio = (1 +  cos(rad)) / 2.0;
	//ratio *= ratio;

	float r_aperture = lerp(r_polygon, r_circ, computeConstants.r * ratio);

	//判定される点が描きたい多角形の内外かを判定
	float col = step(pos, r_aperture);
	//float col = step(pos, r_aperture * (1 + 2 * cos(rad/4.0f)) * 3.0              );
	

	destinationImageR[index] = float4(col, col, col, 1.0);
}

[numthreads(WIDTH, 1, 1)]//computeConstants.N角形
void mainDrawMovingPolygon(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 uv = index / size - float2(computeConstants.posx, computeConstants.posy);

	//判定したい点の位置
	float pos = length(uv);

	//判定したい点のなす角
	float rad = atan2(uv.x, uv.y) + 2.0 * PI;
	rad = rad % (2.0 * PI / computeConstants.N1);

	float r_circ = 0.5 * computeConstants.r1;

	//半径r_circの円に内接する正多角形の辺の位置
	float r_polygon = cos(PI / computeConstants.N1) / cos(PI / computeConstants.N1 - rad);
	r_polygon *= r_circ;

	//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0～1)
	//lerp(x,y,s) = x + s(y - x)
	float s = computeConstants.r;
	float r_aperture = lerp(r_polygon, r_circ, s * s * s);

	//判定される点が描きたい多角形の内外かを判定
	float col = step(pos, r_aperture);

	destinationImageR[index] = float4(col, col, col, 1.0);
}

[numthreads(WIDTH, 1, 1)]
void mainDrawMovingDot(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	uint2 targetIndex =
		uint2
		(
			computeConstants.posx * size.x,
			computeConstants.posy * size.y
			);

	if (targetIndex.x != index.x || targetIndex.y != index.y)
	{
		destinationImageR[index] = float4(0.0, 0.0, 0.0, 1.0);
	}
	else
	{
		destinationImageR[index] = float4(1.0, 1.0, 1.0, 1.0);
	}
}

[numthreads(WIDTH, 1, 1)]
void mainDrawSmallCirc(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 center = float2(WIDTH / 2, HEIGHT / 2);
	float2 vec = index - center;
	
	float length = vec.x * vec.x + vec.y * vec.y;
	
	float range = 3;
	if (length < range)
	{
		destinationImageR[index] = float4(1.0, 1.0, 1.0, 1.0) / range / range;
	}
	else
	{
		destinationImageR[index] = float4(0.0, 0.0, 0.0, 1.0);
	}
}

[numthreads(WIDTH, 1, 1)]
void mainDrawDotByRandomTbl(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	float2 targetPos =
		float2
		(
			computeConstants.posx * size.x,
			computeConstants.posy * size.y
			);

	float2 center =
		float2
		(
			0.5 * size.x,
			0.5 * size.y
			);

	uint table_index = randomTblIndex.Load(0);

	//中心からターゲット方面 係数が負なら逆側へのベクトル
	float2 dir = normalize(targetPos - center);

	float randomValue = randomTbl.data[table_index / 4][table_index % 4];

	//インデックスとして計算(中心基準)
	//float2 dotPoint = center + randomTbl.data[table_index] * dir * size;// randomTbl.data : 0.5で注目点付近 randomTbl.data : -0.5
	float2 dotPoint = center + float2(abs(computeConstants.posx - 0.5), abs(computeConstants.posy - 0.5)) * randomValue * dir * size;// randomTbl.data : 0.5で注目点付近 randomTbl.data : -0.5
	//dotPoint = center - 0.5 * dir * size;
	// dotPoint = center -1 * dir * size;

	if (dotPoint.x >= 0 && dotPoint.x <= WIDTH - 1 && dotPoint.y >= 0 && dotPoint.y <= HEIGHT - 1)
	{
		if ((uint)dotPoint.x == index.x && (uint)dotPoint.y == index.y)
		{
			destinationImageR[index] = float4(1.0, 1.0, 1.0, 1.0);
		}
		else
		{
			destinationImageR[index] = float4(0.0, 0.0, 0.0, 1.0);
		}
	}
	else
	{
		destinationImageR[index] = float4(0.0, 0.0, 0.0, 1.0);
	}
}

[numthreads(1, 1, 1)]
void mainIncrementRandomTblIndex(uint3 dispatchID : SV_DispatchThreadID)
{
	uint tableIndex = randomTblIndex.Load(0);

	tableIndex++;

	randomTblIndex.Store(0, tableIndex);
}

[numthreads(1, 1, 1)]
void mainResetRandomTblIndex(uint3 dispatchID : SV_DispatchThreadID)
{
	randomTblIndex.Store(0, 0);
}

float2 indexfunc2(float2 IndexStandard, float2 scaling)
{
	float2 uvstandard = IndexStandard - float2(0.5 * WIDTH, 0.5 * HEIGHT);
	float2 uv = uvstandard * scaling;

	return uv + float2(0.5 * WIDTH, 0.5 * HEIGHT);
}

float4 sourceImageRValue(uint2 index)
{
	if (index.x >= 0 && index.x < WIDTH && index.y >= 0 && index.y < HEIGHT)
	{
		return sourceImageR[index];
	}
	else
	{
		return float4(0, 0, 0, 1.0);
	}
}

[numthreads(WIDTH, 1, 1)]
void mainRotateByRandomTbl(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	//中心からのベクトル
	float X = (index.x - WIDTH / 2) * ghostGotateMat.elem.x + (index.y - HEIGHT / 2) * ghostGotateMat.elem.z + WIDTH / 2;
	float Y = (index.x - WIDTH / 2) * ghostGotateMat.elem.y + (index.y - HEIGHT / 2) * ghostGotateMat.elem.w + HEIGHT / 2;

	if (X >= 0 && X < WIDTH && Y >= 0 && Y < HEIGHT)
	{
		destinationImageR[index] = sourceImageR[float2(X, Y)];
	}
	else
	{
		destinationImageR[index] = float4(0,0,0,1);
	}
}

[numthreads(WIDTH, 1, 1)]
void mainInverseRotateByRandomTbl(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	//中心からのベクトル
	float X = (index.x - WIDTH / 2) * ghostGotateMat.elem.x + (index.y - HEIGHT / 2) * ghostGotateMat.elem.y + WIDTH / 2;
	float Y = (index.x - WIDTH / 2) * ghostGotateMat.elem.z + (index.y - HEIGHT / 2) * ghostGotateMat.elem.w + HEIGHT / 2;

	if (X >= 0 && X < WIDTH && Y >= 0 && Y < HEIGHT)
	{
		destinationImageR[index] = sourceImageR[float2(X, Y)];
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}
}



float lambdafuncFF(float lambdamax,float lambda)
{
	return ((lambda / lambdamax))/(2 + lambda / lambdamax);
}

[numthreads(WIDTH, 1, 1)]
void mainScalingSizeByRandomTbl(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	float2 targetPos =
		float2
		(
			computeConstants.posx,
			computeConstants.posy
			);

	float2 center =
		float2
		(
			0.5,
			0.5
			);

	//中心からターゲット方面 係数が負なら逆側へのベクトル
	float2 dir = targetPos - center;

	uint table_index = randomTblIndex.Load(0);

	float randomValue = randomTbl.data[table_index / 4][table_index % 4];

	//中心からのベクトル
	float2 vec = randomValue * dir;
	float lengthCenter = length(vec);
	lengthCenter *= lengthCenter;
	float weight = abs(randomValue - uint(randomValue));
	weight *= 100;
	weight = weight % 10;

	float value = (10 * perlinNoise(float2(randomValue, randomValue + 100)) / (abs(randomValue) + 0.1)  / (lengthCenter + weight) / (computeConstants.r + 0.2));//このままでは負の時に小さな値になる
	//float value = 10 * randomValue / (computeConstants.r + 0.2) /( lengthCenter + weight);

	value = abs(value);

	float coef = value / (1 + value);//valueが負の数なら大きくなる
	coef *= (-(lengthCenter - 0.5) * (lengthCenter - 0.5) + 1);

	//0から1の値n をmin maxに変換    y = n(max - min) + min
	float changevalue = lengthCenter / 0.5 / sqrt(2) + 1;
	changevalue *= changevalue;
	changevalue *= changevalue;

	float2 scalingParam = clamp((21 - computeConstants.ghostScale),1, 20) * (1/computeConstants.r) * //ここを大きくすればゴーストは全体的に縮小傾向
		(1 /  (1 + computeConstants.r) / (1 + 0.5 * computeConstants.r)) * value * 
		((((perlinNoise(float2(weight, weight + 3)) * 1000) % 10) > 1) ? float2(coef + 1, coef + changevalue) : float2(coef + changevalue, coef + 1))
		;

	scalingParam += float2(1, 1);//1より小さいとはみ出る 大きいと小さい

	scalingParam *= ((randomValue >= 0) ? 2 : 1);

		float2 texSize;
	float  level;
	destinationImageR.GetDimensions(texSize.x, texSize.y);

	scalingParam.x *= texSize.x / texSize.y;

	scalingParam *= (randomValue > 0) ? -1 : 1;

	float maxlambda = 700;
	float minlambda = 600;

	int samplenumPerRGB = 1;

	float lambdarange = maxlambda - minlambda;

	float lambdaDelta = lambdarange / samplenumPerRGB / 3;

	float3 result = 0.xxx;

	float gapG = 50 * lengthCenter;
	float gapB = 70 * lengthCenter;

	int sampleX = 1;
	int sampleY = 1;

	//for (int x = 0; x < sampleX; x++)
	//{
	//	for (int y = 0; y < sampleY; y++)
	//	{
	//		float2 innnnn = index + float2(x - sampleX/2, y - sampleY/2);

	//		if (innnnn.x >= 0 && innnnn.x < WIDTH && innnnn.y >= 0 && innnnn.y < HEIGHT)
	//		{
	//			for (int i = 0; i < samplenumPerRGB; i++)
	//			{
	//				float lamred = maxlambda - i * lambdaDelta;
	//				float lamgreen = maxlambda - (i + samplenumPerRGB) * lambdaDelta;//より小さい波長からスタート
	//				float lamblue = maxlambda - (i + 2 * samplenumPerRGB) * lambdaDelta;//同上

	//				//波長が大きいほどサイズは大きくする=波長が大きいほどスケーリング係数を小さくする
	//				float2 targetIndexR = indexfunc2(innnnn, scalingParam * (maxlambda / lamred)) + weight.xx;
	//				float2 targetIndexG = indexfunc2(innnnn, scalingParam * (maxlambda / lamgreen)) + weight.xx;
	//				float2 targetIndexB = indexfunc2(innnnn, scalingParam * (maxlambda / lamblue)) + weight.xx;

	//				targetIndexG. x -= gapG;
	//				targetIndexB. x -= gapB;

	//				result += float3(lambdafuncFF(maxlambda, lamred) * sourceImageRValue(targetIndexR).r,
	//					lambdafuncFF(maxlambda, lamgreen) * sourceImageRValue(targetIndexG).g,
	//					lambdafuncFF(maxlambda, lamblue) * sourceImageRValue(targetIndexB).b);
	//			}
	//		}

	//		
	//	}
	//}

	{
		float lamred = maxlambda - lambdaDelta;
		float lamgreen = maxlambda - (samplenumPerRGB)*lambdaDelta;//より小さい波長からスタート
		float lamblue = maxlambda - (2 * samplenumPerRGB) * lambdaDelta;//同上

		//波長が大きいほどサイズは大きくする=波長が大きいほどスケーリング係数を小さくする
		float2 targetIndexR = indexfunc2(index, scalingParam * (maxlambda / lamred)) + weight.xx;
		float2 targetIndexG = indexfunc2(index, scalingParam * (maxlambda / lamgreen)) + weight.xx;
		float2 targetIndexB = indexfunc2(index, scalingParam * (maxlambda / lamblue)) + weight.xx;

		targetIndexG.x -= gapG;
		targetIndexB.x -= gapB;

		result += float3(lambdafuncFF(maxlambda, lamred) * sourceImageRValue(targetIndexR).r,
			lambdafuncFF(maxlambda, lamgreen) * sourceImageRValue(targetIndexG).g,
			lambdafuncFF(maxlambda, lamblue) * sourceImageRValue(targetIndexB).b);
	}



	result /= samplenumPerRGB* sampleX* sampleY;

	float2 rr = float2(randomValue, randomValue + 1);

	float rrr = perlinNoise(rr);

	float2 gg = float2(rrr + randomValue + 1, randomValue + 2);
	float2 bb = float2(randomValue * rrr + 3 + 2 * rrr, randomValue + 4 + rrr);

	float amplitudescale = frac(scalingParam.x) * frac(scalingParam.y) + 2.0;
	amplitudescale *= amplitudescale;
	destinationImageR[index] = float4(amplitudescale * result * float3(rrr, perlinNoise(gg), perlinNoise(bb))//なんかfrac(radis)はだめみたい 1のとき
		, 1.0);
}

[numthreads(WIDTH, 1, 1)]
void mainShiftImageByRandomTbl(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	float2 targetPos =
		float2
		(
			computeConstants.posx * size.x,
			computeConstants.posy * size.y
			);

	float2 center =
		float2
		(
			0.5 * size.x,
			0.5 * size.y
			);

	uint table_index = randomTblIndex.Load(0);

	float randomValue = randomTbl.data[table_index / 4][table_index % 4];

	//中心からターゲット方面 係数が負なら逆側へのベクトル
	//float2 dir = (targetPos - center) / sqrt(size.x * size.x + size.y * size.y);

	float2 rr = float2(randomValue, randomValue + 1);
	float2 gg = float2(randomValue + 1, randomValue + 2);

	//散らばり
	float2 move = float2(perlinNoise(rr), perlinNoise(gg));

	float2 dir = (targetPos - center + move) / sqrt(size.x * size.x + size.y * size.y);

	float2 texSize;
	float  level;
	destinationImageR.GetDimensions(texSize.x, texSize.y);

	dir.y *= texSize.x / texSize.y;

	//インデックスとして計算(中心基準)
	float2 sourcePoint = index - randomValue * dir * size;// randomTbl.data : 0.5で注目点付近 randomTbl.data : -0.5

	if (sourcePoint.x >= 0 && sourcePoint.x < WIDTH && sourcePoint.y >= 0 && sourcePoint.y < HEIGHT)
	{
		destinationImageR[index] = sourceImageR[sourcePoint];
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}
	
}

[numthreads(WIDTH, 1, 1)]
void mainShiftImageByTargetPos(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	float2 targetPos =
		float2
		(
			computeConstants.posx * size.x,
			computeConstants.posy * size.y
			);

	float2 center =
		float2
		(
			0.5 * size.x,
			0.5 * size.y
			);

	//インデックスとして計算(中心基準)
	float2 sourcePoint = index - targetPos + center;

	if (sourcePoint.x >= 0 && sourcePoint.x < WIDTH && sourcePoint.y >= 0 && sourcePoint.y < HEIGHT)
	{
		destinationImageR[index] = sourceImageR[sourcePoint];
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}

}

[numthreads(WIDTH, 1, 1)]
void mainApplyVignettingByRandomTbl(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	float2 targetPos =
		float2
		(
			computeConstants.posx,
			computeConstants.posy
			);

	float2 center =
		float2
		(
			0.5,
			0.5
			);

	uint table_index = randomTblIndex.Load(0);

	float randomValue = randomTbl.data[table_index / 4][table_index % 4];

	float2 rr = float2(randomValue, randomValue + 100);

	float pn = perlinNoise(rr);
	float pnsave = pn;
	pn *= 10;
	pn %= 10;
	pnsave = (pn >= 5) ? -1 : 1;

	float2 distVec = targetPos - center;
	float length = distVec.x * distVec.x + distVec.y * distVec.y;
	//length = sqrt(length);
	length *= sqrt(size.x * size.x + size.y * size.y) * length / 2;

	//インデックスとして計算(中心基準)
	float sourcePointx = index.x - randomValue * pnsave *pn *  length;

	//if (sourcePointx >= 0 && sourcePointx <= HEIGHT - 1) // この条件が原因で切れているように見える
	if ((sourcePointx - WIDTH / 2) * (sourcePointx - WIDTH / 2) + (index.y - HEIGHT / 2) * (index.y - HEIGHT / 2) < (computeConstants.ghostScale) * (WIDTH * WIDTH + HEIGHT * HEIGHT)) // この条件が原因で切れているように見える
	{
		//float weight = smoothstep(0.8, 0.9, 2 * sourcePointx / WIDTH  * length/(1 + length));
		float weight = smoothstep(0.8, 0.9, sourcePointx / (1 + sourcePointx) * abs(randomValue * length) / (abs(randomValue * length) + 1));

		if (pn > 5)
		{
			weight = 0;
		}
		destinationImageR[index] = (1.2 - computeConstants.r) * (1 - weight) * sourceImageI[index];
	}
	else
	{
		destinationImageR[index] = float4(0,0,0,1);
	}

}