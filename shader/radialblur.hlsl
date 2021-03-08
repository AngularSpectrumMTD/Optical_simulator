

[numthreads(WIDTH, 1, 1)]
void mainRadialBlur(uint3 dispatchID : SV_DispatchThreadID)
{
	float3 color[10];

	float2 index = dispatchID.xy;

	float2 center = float2(computeConstants.posx * WIDTH, computeConstants.posy * HEIGHT);

	float2 direction = index - center;

	//‚±‚ê‚Å•úËó‚ÉL‚ª‚é‚æ‚¤‚É‚È‚é ‚µ‚È‚¯‚ê‚Îû‘©
	direction = -direction;

	float len = length(direction);

	direction = normalize(direction) * computeConstants.raylength;

	float w = 9.0/11.0;//(‚Ó‚Â‚¤‚Í10/11‚Æ‚·‚×‚«‚Å‚·‚ªCƒKƒEƒVƒAƒ“‚Æ‚ÌŒ“‚Ë‡‚¢‚Å1‚ğ’´‚¦‚È‚¢‚æ‚¤‚É‚µ‚Ä‚¢‚Ü‚·)

	//Œ¸Š‚³‚¹‚È‚ª‚çL‚Î‚·
	color[0] = sourceImageR[clampF2(index                            )].rgb * 0.20;
	color[1] = sourceImageR[clampF2(index + direction * 1.0)].rgb * 0.18;
	color[2] = sourceImageR[clampF2(index + direction * 2.0)].rgb * 0.16;
	color[3] = sourceImageR[clampF2(index + direction * 3.0)].rgb * 0.14;
	color[4] = sourceImageR[clampF2(index + direction * 4.0)].rgb * 0.12;
	color[5] = sourceImageR[clampF2(index + direction * 5.0)].rgb * 0.10;
	color[6] = sourceImageR[clampF2(index + direction * 6.0)].rgb * 0.08;
	color[7] = sourceImageR[clampF2(index + direction * 7.0)].rgb * 0.06;
	color[8] = sourceImageR[clampF2(index + direction * 8.0)].rgb * 0.04;
	color[9] = sourceImageR[clampF2(index + direction * 9.0)].rgb * 0.02;

	destinationImageR[index].rgb = w * (color[0] + color[1] + color[2] + color[3] + color[4]
		+ color[5] + color[6] + color[7] + color[8] + color[9]);
}

[numthreads(WIDTH, 1, 1)]
void mainDrawDecayLine(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float2 pos = float2(computeConstants.posx * WIDTH, computeConstants.posy * HEIGHT);
	float2 center = float2(0.5 * WIDTH, 0.5 * HEIGHT);

	float2 requiredirection = pos - center;
	float2 currentdirection = index - center;

	float requirelength = length(requiredirection);
	float currentlength = length(currentdirection);

	if (requirelength < currentlength)
	{
		destinationImageR[index].rgb = float3(0.0,0.0,0.0);
		return;
	}

	if(computeConstants.raylength < currentlength)
	{
		destinationImageR[index].rgb = float3(0.0, 0.0, 0.0);
		return;
	}

	requiredirection = normalize(requiredirection);
	currentdirection = normalize(currentdirection);

	float innerproduct = dot(requiredirection, currentdirection);

	if (innerproduct < computeConstants.range)
	{
		destinationImageR[index].rgb = float3(0.0, 0.0, 0.0);
		return;
	}

	float len = length(currentdirection);
	float val = 1.0;

	if (len > 1e-5)
	{
		val = 1.0 / len;
	}

	float3 col = float3(0.5, 0.5, 0.5);

	destinationImageR[index].rgb = col;
}