#include"..\include\Model.h"

vec3 MODEL::IntersectPoint(const Ray& ray)
{
	vec3 D = ray.origin();
	vec3 F = ray.direction();

	double Dx = D.getX(), Dy = D.getY(), Dz = D.getZ();
	double Fx = F.getX(), Fy = F.getY(), Fz = F.getZ();

	CurrentPolygon p = w_currentpolygon;

	vec3 AB = p.w_vertex[1] - p.w_vertex[0];
	vec3 AC = p.w_vertex[2] - p.w_vertex[0];

	double a = p.w_surfacenormal.getX(), b = p.w_surfacenormal.getY(), c = p.w_surfacenormal.getZ();

	double d = dot(p.w_surfacenormal, p.w_vertex[0]);

	double Bu = AB.getX(), Bv = AB.getY(), Bw = AB.getZ();
	double Cu = AC.getX(), Cv = AC.getY(), Cw = AC.getZ();

	double t = -(a * Dx + b * Dy + c * Dz + d) / (a * Fx + b * Fy + c * Fz);

	vec3 ret;
	ret = ray.at(t);

	return ret;
}