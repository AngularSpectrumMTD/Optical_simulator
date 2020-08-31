#include"../include/Model.h"

vec3 Model::IntersectPoint(const Ray& ray)
{
	vec3 D = ray.origin();
	vec3 F = ray.direction();

	double Dx = D.w_x, Dy = D.w_y, Dz = D.w_z;
	double Fx = F.w_x, Fy = F.w_y, Fz = F.w_z;

	CurrentPolygon p = w_currentpolygon;

	vec3 AB = p.w_vertex[1] - p.w_vertex[0];
	vec3 AC = p.w_vertex[2] - p.w_vertex[0];

	double a = p.w_surfacenormal.w_x, b = p.w_surfacenormal.w_y, c = p.w_surfacenormal.w_z;

	double d = dot(p.w_surfacenormal, p.w_vertex[0]);

	double Bu = AB.w_x, Bv = AB.w_y, Bw = AB.w_z;
	double Cu = AC.w_x, Cv = AC.w_y, Cw = AC.w_z;

	double t = -(a * Dx + b * Dy + c * Dz + d) / (a * Fx + b * Fy + c * Fz);

	vec3 ret;
	ret = ray.at(t);

	return ret;
}