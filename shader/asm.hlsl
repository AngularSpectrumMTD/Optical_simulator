[numthreads(WIDTH, 1, 1)]
void mainDrawFRF(uint3 dispatchID : SV_DispatchThreadID)
{
	float2 index = dispatchID.xy;

	float k_r = 1.0f / computeConstants.lambdaR / computeConstants.lambdaR;
	float k_g = 1.0f / computeConstants.lambdaG / computeConstants.lambdaG;
	float k_b = 1.0f / computeConstants.lambdaB / computeConstants.lambdaB;

	float fourier_intervalx = 1.0f / (computeConstants.intervalx * WIDTH);
	float fourier_intervaly = 1.0f / (computeConstants.intervaly * HEIGHT);

	float w_x = 2 * fourier_intervalx * computeConstants.distance0;
	w_x *= w_x;
	float halfband_x
		= 1.0f / (sqrt(w_x + 1.0f) * computeConstants.lambdaB);
	float w_y = 2 * fourier_intervaly * computeConstants.distance0;
	w_y *= w_y;
	float halfband_y
		= 1.0f / (sqrt(w_y + 1.0f) * computeConstants.lambdaB);

	float kx = (index.x - WIDTH / 2.0f) * fourier_intervalx;
	float ky = (index.y - HEIGHT / 2.0f) * fourier_intervaly;

	float w_r = k_r - kx * kx - ky * ky;
	float w_g = k_g - kx * kx - ky * ky;
	float w_b = k_b - kx * kx - ky * ky;

	float3 color_real = float3(0.0, 0.0, 0.0);
	float3 color_image = float3(0.0, 0.0, 0.0);

	if (w_b > 0 && abs(kx) < halfband_x && abs(ky) < halfband_y)//‘Ñˆæ§ŒÀ
	{
		float phase_r = 2 * PI * sqrt(w_r) * computeConstants.distance0;
		float phase_g = 2 * PI * sqrt(w_g) * computeConstants.distance0;
		float phase_b = 2 * PI * sqrt(w_b) * computeConstants.distance0;

		float3 val_real = float3(cos(phase_r), cos(phase_g), cos(phase_b));
		float3 val_image = float3(sin(phase_r), sin(phase_g), sin(phase_b));

		color_real = val_real;
		color_image = val_image;
	}

	destinationImageR[index] = float4(color_real, 1.0f);
	destinationImageI[index] = float4(color_image, 1.0f);
}