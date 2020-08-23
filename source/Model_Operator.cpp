#include"../include/Model.h"

Model& Model::operator +=(const vec3& vec)
{
#pragma omp parallel for schedule(dynamic, 1) num_threads(omp_get_max_threads())
	for (int n = 0; n < (*this).w_Object.size(); ++n)
	{
		for (int m = 0; m < (*this).w_Object[n].w_Vertex.size(); ++m)
		{
			vec3 v = w_Object[n].w_Vertex[m].w_Coord;
			v = v + vec;
			w_Object[n].w_Vertex[m].w_Coord = v;
		}
	}
	return *this;
}

Model& Model::operator *=(const mat3& mat)
{
#pragma omp parallel for schedule(dynamic, 1) num_threads(omp_get_max_threads())
	for (int n = 0; n < (*this).w_Object.size(); ++n)
	{
		for (int m = 0; m < (*this).w_Object[n].w_Vertex.size(); ++m)
		{
			vec3 v = w_Object[n].w_Vertex[m].w_Coord;
			v = v - w_center;
			v = mat * v;
			v = v + w_center;
			w_Object[n].w_Vertex[m].w_Coord = v;
		}
	}
	return *this;
}

void Model::mul(std::vector<vec3>& vec, const mat3& mat)
{
	for (auto& v : vec)
	{
		v = mat * v;
	}
}

void Model::sub(std::vector<vec3>& vec, const vec3& vv)
{
	for (auto& v : vec)
	{
		v -= vv;
	}
}

void Model::fouriermul(std::vector<vec3>& vec, const mat3& mat)
{
	double invlambda2 = 1 / w_lambda / w_lambda;
	auto itr = vec.begin();
	while (itr != vec.end())
	{
		vec3 v = mat * (*itr);//access the element and transform
		float w = invlambda2 - v.getX() * v.getX() - v.getY() * v.getY();
		if (w < 0) {
			// return the iterator indicate next element of deleted one
			itr = vec.erase(itr);
		}
		// in case of not delete element, get the next indicator
		else {
			*itr = vec3{ v.getX(), v.getY(), sqrt(w) };
			itr++;
		}
	}
}
void Model::fouriersub(std::vector<vec3>& vec, const vec3& vv)
{
	double invlambda2 = 1 / w_lambda / w_lambda;
	auto itr = vec.begin();
	while (itr != vec.end())
	{
		vec3 v = (*itr) - vv;//access the element and transform
		float w = invlambda2 - v.getX() * v.getX() - v.getY() * v.getY();
		if (w < 0) {
			// return the iterator indicate next element of deleted one
			itr = vec.erase(itr);
		}
		// in case of not delete element, get the next indicator
		else {
			*itr = vec3{ v.getX(), v.getY(), sqrt(w) };
			itr++;
		}
	}
}