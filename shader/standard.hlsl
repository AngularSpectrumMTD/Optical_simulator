
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

[numthreads(THREADNUM, THREADNUM, 1)]
void mainSepia(uint3 dtid : SV_DispatchThreadID)
{
	if (dtid.x < WIDTH && dtid.y < HEIGHT)
	{
		float3x3 toSepia = float3x3(
			0.393, 0.349, 0.272,
			0.769, 0.686, 0.534,
			0.189, 0.16, 0.131);
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

[numthreads(THREADNUM, THREADNUM, 1)]
void mainAdd(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input1 = sourceImageR[index].rgb;
	float3 input2 = sourceImageI[index].rgb;

	float3 col = input1 + input2;

	destinationImageR[index] = float4(col, 1.0f);
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainMul(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 input1 = sourceImageR[index].rgb;
	float3 input2 = sourceImageI[index].rgb;

	float3 col = input1 * input2;

	destinationImageR[index] = float4(col, 1.0f);
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainWhitening(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float3 col = sourceImageR[index].rgb;

	float raiseCol = max((col.r + col.g + col.b) / 3.0f, computeConstants.minColOfDustTex);

	destinationImageR[index] = float4(raiseCol.rrr, 1.0f);
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
	randomValue = abs(randomValue + 0.1);
	randomValue = 1 / randomValue;
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

//[numthreads(WIDTH, 1, 1)]
//void mainCopy(uint3 dispatchID : SV_DispatchThreadID)
//{
//	float2 index = dispatchID.xy;
//
//	destinationImageR[index] = sourceImageR[index];
//}

[numthreads(THREADNUM, THREADNUM, 1)]
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

	//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0〜1)
	//lerp(x,y,s) = x + s(y - x)
	float s = computeConstants.r;
	float r_aperture = lerp(r_polygon, r_circ, s * s * s);

	//判定される点が描きたい多角形の内外かを判定
	float col = step(pos, r_aperture);

	destinationImageR[index] = float4(col, col, col, 1.0);
}

//float fade_aperture_edge(float radius, float fade, float signed_distance) {
//	float l = radius;
//	float u = radius + fade;
//	float s = u - l;
//	float c = 1.f - saturate(saturate(signed_distance - l) / s);
//	return smoothstep(0, 1, c);
//}

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

	//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0〜1)
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

float2 indexfunc2(float2 IndexStandard, float2 scaling, inout bool isinside)
{
	float2 uvstandard = IndexStandard - float2(0.5 * WIDTH, 0.5 * HEIGHT);
	float2 uv = uvstandard * scaling;
	float2 ret = uv + float2(0.5 * WIDTH, 0.5 * HEIGHT);

	isinside = true;

	if (ret.x < 0 || ret.x > WIDTH - 1 || ret.y < 0 || ret.y > HEIGHT - 1)
	{
		isinside = false;
	}

	return ret;
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

float4 sourceImageRValueBilinearClamp(float2 index)
{
	uint i1 = (uint)floor(index.x);
	uint j1 = (uint)floor(index.y);

	float2 texSize;
	sourceImageR.GetDimensions(texSize.x, texSize.y);

	if (i1 >= 0 && i1 < texSize.x && j1 >= 0 && j1 < texSize.y)
	{
		uint i2, j2;

		if (index.x == texSize.x - 1)
		{
			i2 = i1;
		}
		else
		{
			i2 = i1 + 1;
		}

		if (index.y == texSize.y - 1)
		{
			j2 = j1;
		}
		else
		{
			j2 = j1 + 1;
		}

		float t = (index.x - i1);
		float u = (index.y - j1);
		float t1 = 1.0f - t;
		float u1 = 1.0f - u;

		return sourceImageR[float2(i2, j2)] * (t * u)
			+ sourceImageR[float2(i1, j2)] * (t1 * u)
			+ sourceImageR[float2(i1, j1)] * (t1 * u1)
			+ sourceImageR[float2(i2, j1)] * (t * u1);

		//return sourceImageR[index];
	}
	else
	{
		return float4(0, 0, 0, 1.0);
	}
}

void Cubic4(inout float bufx[4], float x1)
{
	float t1 = x1 - 1.0;
	float t2 = x1 * t1;
	float t3 = t2 - 1.0;

	float p[4] = { -t1, t1, -x1, x1 };
	float q[4] = { t2,  t3, t3,  t2 };
	int i;

	for (i = 0; i < 4; i++) {
		bufx[i] = p[i] * q[i];
	}
	return;
}

float4 sourceImageRValueBicubicClamp(float2 index)
{
	float tx = index.x;
	float ty = index.y;

	int i, j;
	int tx_int = (int)(tx);
	int ty_int = (int)(ty);
	float tx_dec = tx - tx_int;
	float ty_dec = ty - ty_int;

	int il = 0; int jl = 0;
	int ir = 4; int jr = 4;
	int i_bufpos = tx_int - 1;
	int j_bufpos = ty_int - 1;

	float2 texSize;
	sourceImageR.GetDimensions(texSize.x, texSize.y);

	if (i_bufpos < 0) {
		il = -i_bufpos;
	}
	if (j_bufpos < 0) {
		jl = -j_bufpos;
	}
	if (ir > texSize.x - i_bufpos) {
		ir = texSize.x - i_bufpos;
	}
	if (jr > texSize.y - j_bufpos) {
		jr = texSize.y - j_bufpos;
	}
	if (il >= ir || jl >= jr) {
		return float4(0, 0, 0, 1.0);
	}

	float bufx[4];
	float bufy[4];
	float dbufR[4][4] = {
		{ 0, 0, 0, 0} ,
		{ 0, 0, 0, 0} ,
		{ 0, 0, 0, 0} ,
		{ 0, 0, 0, 0}
	};

	float dbufG[4][4] = {
	{ 0, 0, 0, 0} ,
	{ 0, 0, 0, 0} ,
	{ 0, 0, 0, 0} ,
	{ 0, 0, 0, 0}
	};

	float dbufB[4][4] = {
	{ 0, 0, 0, 0} ,
	{ 0, 0, 0, 0} ,
	{ 0, 0, 0, 0} ,
	{ 0, 0, 0, 0}
	};

	for (j = jl; j < jr; j++) {
		for (i = il; i < ir; i++) {
			dbufR[j][i] = sourceImageR[float2(i + i_bufpos, j + j_bufpos + jl)].r;
			dbufG[j][i] = sourceImageR[float2(i + i_bufpos, j + j_bufpos + jl)].g;
			dbufB[j][i] = sourceImageR[float2(i + i_bufpos, j + j_bufpos + jl)].b;
		}
	}

	Cubic4(bufx, tx_dec);
	Cubic4(bufy, ty_dec);

	double c1rR;
	double c2rR = 0;

	double c1rG;
	double c2rG = 0;

	double c1rB;
	double c2rB = 0;

	for (j = 0; j < 4; j++) {
		c1rR = 0;
		c1rG = 0;
		c1rB = 0;
		for (i = 0; i < 4; i++) {
			c1rR += dbufR[j][i] * bufx[i];
			c1rG += dbufG[j][i] * bufx[i];
			c1rB += dbufB[j][i] * bufx[i];
		}
		c2rR += c1rR * bufy[j];
		c2rG += c1rG * bufy[j];
		c2rB += c1rB * bufy[j];
	}

	return float4(c2rR, c2rG, c2rB, 1);
}

//[numthreads(WIDTH, 1, 1)]
[numthreads(THREADNUM, THREADNUM, 1)]
void mainRotateByRandomTbl(uint3 dispatchID : SV_DispatchThreadID)//縦より横が長い時 回転するとはみ出てしまうのでやり方考えないといけない
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	//中心からのベクトル
	float X = (index.x - WIDTH / 2) * ghostGotateMat.elem.x + (index.y - HEIGHT / 2) * ghostGotateMat.elem.z + WIDTH / 2;
	float Y = (index.x - WIDTH / 2) * ghostGotateMat.elem.y + (index.y - HEIGHT / 2) * ghostGotateMat.elem.w + HEIGHT / 2;

	if (X >= 0 && X < WIDTH && Y >= 0 && Y < HEIGHT)
	{
		//NN
		destinationImageR[index] = sourceImageR[float2(X, Y)];
		//BL(GOOD)
		//destinationImageR[index] = sourceImageRValueBilinearClamp(float2(X, Y));
		//Bicubic
		//destinationImageR[index] = sourceImageRValueBicubicClamp(float2(X, Y));
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}
}

[numthreads(THREADNUM, THREADNUM, 1)]
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



float lambdafuncFF(float lambdamax, float lambda)
{
	return ((lambda / lambdamax)) / (2 + lambda / lambdamax);
}

[numthreads(WIDTH, 1, 1)]//computeConstants.N角形
void mainDrawPolygonFixScaleForBurst(uint3 dispatchID : SV_DispatchThreadID)
{
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
		float r_circ = 0.48;//max
		//r_circ = 0.1;//max

		//半径r_circの円に内接する正多角形の辺の位置
		float r_polygon = cos(PI / computeConstants.N) / cos(PI / computeConstants.N - rad);
		r_polygon *= r_circ;

		//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0〜1)
		//lerp(x,y,s) = x + s(y - x)
		float s = computeConstants.r;

		uint term = 2;//この周期分はいる

		float ratio = (1 + cos(term * rad)) / 2.0;

		float r_aperture = lerp(r_polygon, r_circ, computeConstants.r * ratio);

		//判定される点が描きたい多角形の内外かを判定
		float col = step(pos, r_aperture);

		destinationImageR[index] = float4(col, col, col, 1.0);
	}
}

[numthreads(WIDTH, 1, 1)]//computeConstants.N角形
void mainDrawPolygonFixScaleForGhost(uint3 dispatchID : SV_DispatchThreadID)
{
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
		float r_circ = 0.48;//max

		//半径r_circの円に内接する正多角形の辺の位置
		float r_polygon = cos(PI / computeConstants.N) / cos(PI / computeConstants.N - rad);
		r_polygon *= r_circ;

		//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0〜1)
		//lerp(x,y,s) = x + s(y - x)
		float s = computeConstants.r;

		uint term = 2;//この周期分はいる

		float ratio = (1 + cos(term * rad)) / 2.0;

		float r_aperture = lerp(r_polygon, r_circ, computeConstants.r * ratio);

		float minRadius = lerp(0.9, 0.8, computeConstants.r * computeConstants.r) * r_aperture;
		//minRadius = 0.9 * r_aperture;

		//判定される点が描きたい多角形の内外かを判定
		float col = smoothstep(1.0 * r_aperture, minRadius, pos);

		destinationImageR[index] = float4(col, col, col, 1.0);
	}
}

//[numthreads(WIDTH, 1, 1)]
[numthreads(THREADNUM, THREADNUM, 1)]
void mainCaustic(uint3 dispatchID : SV_DispatchThreadID)
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
	float r_circ = 0.48;//max

	//半径r_circの円に内接する正多角形の辺の位置
	float r_polygon = cos(PI / computeConstants.N) / cos(PI / computeConstants.N - rad);
	r_polygon *= r_circ;

	//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0〜1)
	//lerp(x,y,s) = x + s(y - x)
	float s = computeConstants.r;
	//float r_aperture = lerp(r_polygon, r_circ, s * s * s);

	float ratio = (1 + cos(rad * 2)) / 2.0;
	//ratio *= ratio;

	float r_aperture = lerp(r_polygon, r_circ, computeConstants.r * ratio);

	//float caustic = smoothstep(0.1, 0.0, abs(r_aperture - pos));//左の値に近いほど0 右の値に近いほど1
	float caustic = 1 + computeConstants.N * computeConstants.r * smoothstep(0.02, 0.01, abs(r_aperture - pos));//左の値に近いほど0 右の値に近いほど1

	float col = step(pos, r_aperture);
	float colsave = col;

	col *= 1.01 - computeConstants.r;
	//caustic += 1;

	if ((sourceImageR[index].x + sourceImageR[index].y + sourceImageR[index].z) / 3.0f > 0.1 && colsave == 1)
	{
		destinationImageR[index] = sourceImageR[index];
	}
	else
	{
		destinationImageR[index] = float4(0,0,0,1);
	}

	//destinationImageR[index] = caustic * col * sourceImageR[index];
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainCutOff(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 uv = index / size - float2(0.5, 0.5);

	//判定したい点の位置
	float pos = length(uv);

	//判定したい点のなす角
	float rad = atan2(uv.x, uv.y) + 2.0 * PI + computeConstants.rotAngle * PI / 180.0f;//ここで角度足したら回る
	rad = rad % (2.0 * PI / computeConstants.N);

	float r_circ = 0.48;

	//半径r_circの円に内接する正多角形の辺の位置
	float r_polygon = cos(PI / computeConstants.N) / cos(PI / computeConstants.N - rad);
	r_polygon *= r_circ;

	//円と多角形の中間(半径を大きくするほど=絞らないほど円形に近づく computeConstants.rは0〜1)
	//lerp(x,y,s) = x + s(y - x)
	float s = computeConstants.r;

	float ratio = (1 + cos(rad * 2)) / 2.0;

	float a = (atan2(uv.x, uv.y) + 0) / 2.0f / PI + 3.f / 4.f;
	float o = frac(a * computeConstants.N + 0.5);
	float w2 = lerp(0.025, 0.001f, saturate((computeConstants.N - 4) / 10.f));
	float s0 = sin(o * 2 * PI);
	float s2 = s0 * w2;

	float r_aperture = lerp(r_polygon, r_circ, computeConstants.r * ratio);

	//float col = step(pos, r_aperture);
	float col = smoothstep(1.0 * r_aperture, 0.9 * r_aperture, pos);
	//float col = smoothstep(1.0 * r_aperture, 0.8 * r_aperture, pos);
	//float col = fade_aperture_edge(1.0 * r_aperture, 0.4, 0.9 * pos);

	destinationImageR[index] = 1 * sourceImageR[index];
}

uint Xorshift(uint seed)
{
	uint x = seed;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return x;
}

float XorshiftZeroToOne(uint seed)
{
	uint x = seed;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;

	float ret = x / (1 + x);

	return ret;
}

//[numthreads(WIDTH, 1, 1)]
[numthreads(THREADNUM, THREADNUM, 1)]
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
	float lengthCenter = length(randomValue * dir);
	float lengthSave = lengthCenter;
	lengthCenter *= lengthCenter;
	float weight = abs(randomValue - uint(randomValue));
	weight *= 100;
	weight = weight % 10;

	float perlinNoiseR = perlinNoise(float2(randomValue, randomValue + 100));

	//float2 scalingParam = float2(51 - computeConstants.ghostScale, 51 - computeConstants.ghostScale) * weight;

	////if (weight < 9)
	////{
	////	scalingParam *= weight;
	////}

	//scalingParam.x = clamp(scalingParam.x, 1, 50);
	//scalingParam.y = clamp(scalingParam.y, 1, 50);

	float2 scalingParam = float2(2 - computeConstants.r, 2 - computeConstants.r) * weight * (lengthSave + 1);
	scalingParam += float2(1, 1);

	float2 texSize;
	destinationImageR.GetDimensions(texSize.x, texSize.y);

	//scalingParam.x *= texSize.x / texSize.y;
	scalingParam.x *= 0.5 * computeConstants.screenWidth / computeConstants.screenHeight;

	scalingParam *= (randomValue > 0) ? -1 : 1;

	float maxlambda = 800;
	float minlambda = 700;

	int samplenumPerRGB = 1;

	float lambdarange = maxlambda - minlambda;

	float lambdaDelta = lambdarange / samplenumPerRGB / 3;

	float3 result = 0.xxx;

	float gapG = 50 * lengthCenter;
	float gapB = 70 * lengthCenter;

	int sampleX = 1;
	int sampleY = 1;

	{
		float lamred = maxlambda - lambdaDelta;
		float lamgreen = maxlambda - (samplenumPerRGB)*lambdaDelta;//より小さい波長からスタート
		float lamblue = maxlambda - (2 * samplenumPerRGB) * lambdaDelta;//同上

		bool judgeR = true;
		bool judgeG = true;
		bool judgeB = true;

		//波長が大きいほどサイズは大きくする=波長が大きいほどスケーリング係数を小さくする         中心からの距離は targetIndex > index
		float2 targetIndexR = indexfunc2(index, scalingParam * (maxlambda / lamred), judgeR) + weight.xx;
		float2 targetIndexG = indexfunc2(index, scalingParam * (maxlambda / lamgreen), judgeG) + weight.xx;
		float2 targetIndexB = indexfunc2(index, scalingParam * (maxlambda / lamblue), judgeB) + weight.xx;

	/*	if (!isBig)
		{*/
			targetIndexG.x -= gapG;
			targetIndexB.x -= gapB;
		//}
		

		if (judgeR && judgeG && judgeB)
		{
			//NN
		/*	result = computeConstants.baseColor * float3(lambdafuncFF(maxlambda, lamred) * sourceImageRValueBilinearClamp(targetIndexR).r,
				lambdafuncFF(maxlambda, lamgreen) * sourceImageRValueBilinearClamp(targetIndexG).g,
				lambdafuncFF(maxlambda, lamblue) * sourceImageRValueBilinearClamp(targetIndexB).b);*/

			

			result = computeConstants.baseColor * float3(lambdafuncFF(maxlambda, lamred) * sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, targetIndexR / texSize, 0).r,
				lambdafuncFF(maxlambda, lamgreen) * sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, targetIndexG / texSize, 0).r,
				lambdafuncFF(maxlambda, lamblue) * sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, targetIndexB / texSize, 0).r);
		}
		else
		{
			result += float3(0, 0, 0);
		}
	}

	result /= samplenumPerRGB * sampleX * sampleY;

	float2 rr = float2(randomValue, randomValue + 1);

	float2 gg = float2(perlinNoiseR + randomValue + 1, randomValue + 2);
	float2 bb = float2(randomValue * perlinNoiseR + 3 + 2 * perlinNoiseR, randomValue + 4 + perlinNoiseR);

	float amplitudescale = length(scalingParam);// / abs(randomValue + 1) / 200 / sqrt(2);
	//amplitudescale = amplitudescale / (1 + amplitudescale);
	destinationImageR[index] = float4(amplitudescale * result * float3(perlinNoiseR, perlinNoise(gg), perlinNoise(bb))//なんかfrac(radis)はだめみたい 1のとき
		, 1.0);
}

//[numthreads(WIDTH, 1, 1)]
[numthreads(THREADNUM, THREADNUM, 1)]
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
	//dir.x *= computeConstants.screenWidth / computeConstants.screenHeight;

	//インデックスとして計算(中心基準)
	float2 sourcePoint = index - randomValue * dir * size;// randomTbl.data : 0.5で注目点付近 randomTbl.data : -0.5

	if (sourcePoint.x >= 0 && sourcePoint.x < WIDTH && sourcePoint.y >= 0 && sourcePoint.y < HEIGHT)
	{
		//destinationImageR[index] = sourceImageRValueBilinearClamp(sourcePoint);
		destinationImageR[index] = sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, sourcePoint / texSize, 0);
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}

}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainShiftImageByRandomTblandAdd(uint3 dispatchID : SV_DispatchThreadID)
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
	//dir.x *= computeConstants.screenWidth / computeConstants.screenHeight;

	//インデックスとして計算(中心基準)
	float2 sourcePoint = index - randomValue * dir * size;// randomTbl.data : 0.5で注目点付近 randomTbl.data : -0.5

	if (sourcePoint.x >= 0 && sourcePoint.x < WIDTH && sourcePoint.y >= 0 && sourcePoint.y < HEIGHT)
	{
		//destinationImageR[index] += sourceImageRValueBilinearClamp(sourcePoint);
		destinationImageR[index] += sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, sourcePoint / texSize, 0);
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
		//destinationImageR[index] = sourceImageRValueBilinearClamp(sourcePoint);
		destinationImageR[index] = sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, sourcePoint / size, 0);
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}

}

//[numthreads(WIDTH, 1, 1)]
[numthreads(THREADNUM, THREADNUM, 1)]
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
	float len = distVec.x * distVec.x + distVec.y * distVec.y;

	//インデックスとして計算(中心基準)
	float sourcePointx = index.x - len * length(size);
	float aperture_mask = fade_aperture_edge(0.7, 0.1, len);

	//if ((sourcePointx - WIDTH / 2) * (sourcePointx - WIDTH / 2) + (index.y - HEIGHT / 2) * (index.y - HEIGHT / 2) < (WIDTH * WIDTH + HEIGHT * HEIGHT) / 4) // この条件が原因で切れているように見える
	//{
	//	//float weight = smoothstep(0.8, 0.9, 2 * sourcePointx / WIDTH  * length/(1 + length));
	//	float weight = smoothstep(0.8, 0.9, sourcePointx / (1 + sourcePointx) * abs(randomValue * len) / (abs(randomValue * len) + 1));

	//	if (pn > 5)
	//	{
	//		weight = 0;
	//	}
	//	//destinationImageR[index] = (1.2 - computeConstants.r) * (1 - weight) * sourceImageI[index];
	//	destinationImageR[index] = aperture_mask *  sourceImageI[index];
	//}
	//else
	//{
	//	destinationImageR[index] = float4(0, 0, 0, 1);
	//}

	float2 texSize;
	float  level;
	destinationImageR.GetDimensions(texSize.x, texSize.y);

	if ((index.x - WIDTH / 2) * (index.x - WIDTH / 2) + (index.y - HEIGHT / 2) * (index.y - HEIGHT / 2) < (abs(randomValue) / (1 + abs(randomValue))) * (WIDTH * WIDTH + HEIGHT * HEIGHT) / 4) // この条件が原因で切れているように見える
	{
		float aaaaaaaa = (abs(randomValue) / (1 + abs(randomValue))) * (WIDTH * WIDTH + HEIGHT * HEIGHT) / 4 - (index.x - WIDTH / 2) * (index.x - WIDTH / 2) + (index.y - HEIGHT / 2) * (index.y - HEIGHT / 2);

		aaaaaaaa /= (WIDTH * WIDTH + HEIGHT * HEIGHT) / 4;
		//float weight = smoothstep(0.0, 0.9, aaaaaaaa);

		float rrrrr = (index.x - WIDTH / 2) * (index.x - WIDTH / 2) + (index.y - HEIGHT / 2) * (index.y - HEIGHT / 2) * (texSize.x / texSize.y) * (texSize.x / texSize.y);
		rrrrr /= (WIDTH * WIDTH + HEIGHT * HEIGHT) / 4;

		float sigsig = 0.5;

		float weight = exp(-(rrrrr * rrrrr) / sigsig * abs(randomValue));
		destinationImageR[index] = 3 * weight *
			sourceImageI[index];
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainBlur(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float3 col = 0.xxx;

	//const float involve = 4 + computeConstants.r * computeConstants.N;
	//const float involve = clamp(2 +computeConstants.r * computeConstants.N, 2.f, 5.f);//このブラーの巻き込み数の変化がちらつきの元
	const float involve = 20;
	//const float involve = 3;

	for (int i = 0; i < involve; i++)
	{
		for (int j = 0; j < involve; j++)
		{
			float2 moved = index + float2(i - involve/2, j - involve / 2);

			col += sourceImageRValue(moved).rgb;
		}
	}

	col /= involve * involve;

	destinationImageR[index] = float4(col, 1.0);
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainFade(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float2 size = float2(WIDTH, HEIGHT);

	float2 uv = index / size - 0.5;
	float d = length(uv) * 2;
	float fade = lerp(1.0f, 0.f, d);

	//fade = 1;

	destinationImageR[index] = float4(fade * sourceImageR[index].rgb, 1.0f);
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainScaleByPos(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 pos = index / size;

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

	//インデックスとして計算(中心基準)
	float2 uv = pos.xy - 0.5;

	float2 scalingWeight = length(targetPos - center) + computeConstants.r * computeConstants.r + 1;
	scalingWeight.x *= 0.5 * computeConstants.screenWidth / computeConstants.screenHeight;

	float2 scaleduv = uv * scalingWeight + 0.5;

	if (scaleduv.x >= 0 && scaleduv.x < 1 && scaleduv.y >= 0 && scaleduv.y < 1)
	{
		//destinationImageR[index] = sourceImageRValueBilinearClamp(sourcePoint);
		destinationImageR[index] = sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, scaleduv, 0);
	}
	else
	{
		destinationImageR[index] = float4(0, 0, 0, 1);
	}
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainDownIntensityByPos(uint3 dispatchID : SV_DispatchThreadID)
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

	//インデックスとして計算(中心基準)
	//float2 sourcePoint = index - targetPos + center;

	float len = length(targetPos - center);
	//float len = 0;

	destinationImageR[index] = sourceImageR[index] / (len * 10 + 1.0f);//シェーダで分数はやばいのか ちらつく
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainShiftSample(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);

	float2 targetIndex = index + 100 * (computeConstants.elapsedTime).xx;

	destinationImageR[index] = sourceImageR.SampleLevel(CSimageSamplerBILINEAR_WRAP, targetIndex / size, 0);

	//destinationImageR[index] = sourceImageR[targetIndex];
}

[numthreads(THREADNUM, THREADNUM, 1)]
void mainMultipleMoveCircle(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;
	float2 size = float2(WIDTH, HEIGHT);
	float2 uv = index / size - 0.5;

	float2 targetPos =
		float2
		(
			computeConstants.posx - 0.5,
			computeConstants.posy - 0.5
			);

	float R = 0.5 + 0.5 * (1 - computeConstants.r);
	R *= R;
	const float R2 = R;

	const float moveRatio = 2.0;

	if (length(float2(uv.x - moveRatio * targetPos.x, uv.y - moveRatio * targetPos.y)) < R2)
	{
		destinationImageR[index] = sourceImageR[index];
	}
	else
	{
		destinationImageR[index] = float4(0,0,0,1);
	}
}