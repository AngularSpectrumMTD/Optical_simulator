#include"Model.h"

vec3 MODEL::IntersectPoint(Ray& ray)
{
	vec3 D = ray.origin();
	vec3 F = ray.direction();

	double Dx = D.getX(), Dy = D.getY(), Dz = D.getZ();
	double Fx = F.getX(), Fy = F.getY(), Fz = F.getZ();

	CurrentPolygon p = _currentpolygon;

	vec3 AB = p._vertex[1] - p._vertex[0];
	vec3 AC = p._vertex[2] - p._vertex[0];

	double a = p._surfacenormal.getX(), b = p._surfacenormal.getY(), c = p._surfacenormal.getZ();

	double d = dot(p._surfacenormal, p._vertex[0]);

	double Bu = AB.getX(), Bv = AB.getY(), Bw = AB.getZ();
	double Cu = AC.getX(), Cv = AC.getY(), Cw = AC.getZ();

	double t = -(a * Dx + b * Dy + c * Dz + d) / (a * Fx + b * Fy + c * Fz);

	vec3 ret;
	ret = ray.at(t);

	return ret;
}