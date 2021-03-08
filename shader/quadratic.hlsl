[numthreads(WIDTH, 1, 1)]
void mainDrawQuadraticPhase(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float k_r = 1.0f / computeConstants.lambdaR;
	float k_g = 1.0f / computeConstants.lambdaG;
	float k_b = 1.0f / computeConstants.lambdaB;

	float x = (index.x - WIDTH / 2.0f) * computeConstants.intervalx;
	float y = (index.y - HEIGHT / 2.0f) * computeConstants.intervaly;

	float3 color_real = float3(0.0, 0.0, 0.0);
	float3 color_image = float3(0.0, 0.0, 0.0);

	float p = -(x * x + y * y) / (2 * computeConstants.focus0);//—¼“Ê
	//float p = -(x * x + y * y) / (2 * computeConstants.focus0);//—¼‰š

	float phase_r = 2 * PI * k_r * p;
	float phase_g = 2 * PI * k_g * p;
	float phase_b = 2 * PI * k_b * p;

	float3 val_real = float3(cos(phase_r), cos(phase_g), cos(phase_b));
	float3 val_image = float3(sin(phase_r), sin(phase_g), sin(phase_b));

	color_real = val_real;
	color_image = val_image;

	destinationImageR[index] = float4(color_real, 1.0f);
	destinationImageI[index] = float4(color_image, 1.0f);
}