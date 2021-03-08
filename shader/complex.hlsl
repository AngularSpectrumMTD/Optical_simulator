//ï°ëfêîóp
float2 complex_conjugate(float2 cmp)
{
	return float2(cmp.x, -cmp.y);
}
float complex_sqr(float2 cmp)
{
	return cmp.x * cmp.x + cmp.y * cmp.y;
}
float complex_norm(float2 cmp)
{
	return sqrt(complex_sqr(cmp));
}
float2 complex_add(float2 cmp1, float2 cmp2)
{
	return float2(cmp1.x + cmp2.x, cmp1.y + cmp2.y);
}
float2 complex_sub(float2 cmp1, float2 cmp2)
{
	return float2(cmp1.x - cmp2.x, cmp1.y - cmp2.y);
}
float2 complex_mul(float2 cmp1, float2 cmp2)
{
	return float2(cmp1.x * cmp2.x - cmp1.y * cmp2.y, cmp1.y * cmp2.x + cmp1.x * cmp2.y);
}
float2 complex_div(float2 cmp1, float2 cmp2)
{
	float2 cmp = complex_mul(cmp1, complex_conjugate(cmp2));
	float sqr = complex_sqr(cmp2);
	return float2(cmp.x / sqr, cmp.y / sqr);
}
float2 complex_polar(float amp, float phase)
{
	return float2(amp * cos(phase), amp * sin(phase));
}
//ï°ëfêîóp