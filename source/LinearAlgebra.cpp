#include "..\include\LinearAlgebra.h"

LAVector& LAVector::operator = (const LAVector& v)
{
	this->w_x = v.w_x;
	this->w_y = v.w_y;
	this->w_z = v.w_z;
	return *this;
}

LAVector LAVector::operator+()
{
	return *this;
}

LAVector LAVector::operator-()
{
	return LAVector(-w_x, -w_y, -w_z);
}

LAVector& LAVector::operator += (const LAVector& v)
{
	this->w_x += v.w_x;
	this->w_y += v.w_y;
	this->w_z += v.w_z;
	return *this;
}

LAVector& LAVector::operator -= (const LAVector& v)
{
	this->w_x -= v.w_x;
	this->w_y -= v.w_y;
	this->w_z -= v.w_z;
	return *this;
}

LAVector& LAVector::operator *= (const float k)
{
	this->w_x *= k;
	this->w_y *= k;
	this->w_z *= k;
	return *this;
}

LAVector& LAVector::operator /= (const float k)
{
	this->w_x /= k;
	this->w_y /= k;
	this->w_z /= k;
	return *this;
}

LAVector operator +(const LAVector& u, const LAVector& v)
{
	LAVector w;
	w.setX(u.getX() + v.getX());
	w.setY(u.getY() + v.getY());
	w.setZ(u.getZ() + v.getZ());
	return w;
}

LAVector operator -(const LAVector& u, const LAVector& v)
{
	LAVector w;
	w.setX(u.getX() - v.getX());
	w.setY(u.getY() - v.getY());
	w.setZ(u.getZ() - v.getZ());
	return w;
}

float operator *(const LAVector& u, const LAVector& v)//innerproduct
{
	return u.getX() * v.getX() + u.getY() * v.getY() + u.getZ() * v.getZ();
}

LAVector operator *(const float k, const LAVector& v)
{
	LAVector w(k * v.getX(), k * v.getY(), k * v.getZ());
	return w;
}

LAVector operator *(const LAVector& v, const float k)
{
	LAVector w(k * v.getX(), k * v.getY(), k * v.getZ());
	return w;
}

LAVector operator /(const LAVector& v, const float k)
{
	LAVector w(v.getX() / k, v.getY() / k,  v.getZ() / k);
	return w;
}

LAVector operator &&(const LAVector& u, const LAVector& v)
{
	LAVector w(
		u.getY() * v.getZ() - u.getZ() * v.getY()
		, u.getZ() * v.getX() - u.getX() * v.getZ()
		, u.getX() * v.getY() - u.getY() * v.getX());
	return w;
}

LAMatrix& LAMatrix::operator = (const LAMatrix& mat)
{
	w_col0 = mat.w_col0;
	w_col1 = mat.w_col1;
	w_col2 = mat.w_col2;
	return *this;
}

LAMatrix LAMatrix::operator + ()
{
	return *this;
}

LAMatrix LAMatrix::operator - ()
{
	w_col0 = -w_col0;
	w_col1 = -w_col1;
	w_col2 = -w_col2;
	return *this;
}

LAMatrix& LAMatrix::operator += (const LAMatrix& mat)
{
	*this = *this + mat;
	return *this;
}

LAMatrix& LAMatrix::operator -= (const LAMatrix& mat)
{
	*this = *this - mat;
	return *this;
}

LAMatrix& LAMatrix::operator *= (const float k)
{
	*this = *this * k;
	return *this;
}

LAMatrix& LAMatrix::operator *= (const LAMatrix& mat)
{
	*this = *this * mat;
	return *this;
}

LAMatrix operator +(const LAMatrix& mat1, const LAMatrix& mat2)
{
	return LAMatrix(
		mat1.getCol0() + mat2.getCol0(),
		mat1.getCol1() + mat2.getCol1(),
		mat1.getCol2() + mat2.getCol2()
	);
}

LAMatrix operator -(const LAMatrix& mat1, const LAMatrix& mat2)
{
	return LAMatrix(
		mat1.getCol0() - mat2.getCol0(),
		mat1.getCol1() - mat2.getCol1(),
		mat1.getCol2() - mat2.getCol2()
	);
}

LAMatrix operator *(const LAMatrix& mat1, const LAMatrix& mat2)
{
	return LAMatrix(
		mat1 * mat2.getCol0(),
		mat1 * mat2.getCol1(),
		mat1 * mat2.getCol2()
	);
}

LAMatrix operator *(const float k, const LAMatrix& mat)
{
	return LAMatrix(
		mat.getCol0()*k,
		mat.getCol1() *k,
		mat.getCol2() *k
	);
}

LAMatrix operator *(const LAMatrix& mat, const float k)
{
	return LAMatrix(
		mat.getCol0() * k,
		mat.getCol1() * k,
		mat.getCol2() * k
	);
}

LAVector operator *(const LAMatrix& mat ,const LAVector& vec)
{
	return LAVector(
		(((mat.getCol0().getX() * vec.getX()) + (mat.getCol1().getX() * vec.getY())) + (mat.getCol2().getX() * vec.getZ())),
		(((mat.getCol0().getY() * vec.getX()) + (mat.getCol1().getY() * vec.getY())) + (mat.getCol2().getY() * vec.getZ())),
		(((mat.getCol0().getZ() * vec.getX()) + (mat.getCol1().getZ() * vec.getY())) + (mat.getCol2().getZ() * vec.getZ()))
	);
}