#include "..\include\BLAS.h"

BLASVector& BLASVector::operator = (const BLASVector& v)
{
	this->w_x = v.w_x;
	this->w_y = v.w_y;
	this->w_z = v.w_z;
	return *this;
}

BLASVector BLASVector::operator+()
{
	return *this;
}

BLASVector BLASVector::operator-()
{
	return BLASVector(-w_x, -w_y, -w_z);
}

BLASVector& BLASVector::operator += (const BLASVector& v)
{
	this->w_x += v.w_x;
	this->w_y += v.w_y;
	this->w_z += v.w_z;
	return *this;
}

BLASVector& BLASVector::operator -= (const BLASVector& v)
{
	this->w_x -= v.w_x;
	this->w_y -= v.w_y;
	this->w_z -= v.w_z;
	return *this;
}

BLASVector& BLASVector::operator *= (const float k)
{
	this->w_x *= k;
	this->w_y *= k;
	this->w_z *= k;
	return *this;
}

BLASVector& BLASVector::operator /= (const float k)
{
	this->w_x /= k;
	this->w_y /= k;
	this->w_z /= k;
	return *this;
}

BLASVector operator +(const BLASVector& u, const BLASVector& v)
{
	BLASVector w;
	w.setX(u.getX() + v.getX());
	w.setY(u.getY() + v.getY());
	w.setZ(u.getZ() + v.getZ());
	return w;
}

BLASVector operator -(const BLASVector& u, const BLASVector& v)
{
	BLASVector w;
	w.setX(u.getX() - v.getX());
	w.setY(u.getY() - v.getY());
	w.setZ(u.getZ() - v.getZ());
	return w;
}

float operator *(const BLASVector& u, const BLASVector& v)//innerproduct
{
	return u.getX() * v.getX() + u.getY() * v.getY() + u.getZ() * v.getZ();
}

BLASVector operator *(const float k, const BLASVector& v)
{
	BLASVector w(k * v.getX(), k * v.getY(), k * v.getZ());
	return w;
}

BLASVector operator *(const BLASVector& v, const float k)
{
	BLASVector w(k * v.getX(), k * v.getY(), k * v.getZ());
	return w;
}

BLASVector operator /(const BLASVector& v, const float k)
{
	BLASVector w(v.getX() / k, v.getY() / k,  v.getZ() / k);
	return w;
}

BLASVector operator &&(const BLASVector& u, const BLASVector& v)
{
	BLASVector w(
		u.getY() * v.getZ() - u.getZ() * v.getY()
		, u.getZ() * v.getX() - u.getX() * v.getZ()
		, u.getX() * v.getY() - u.getY() * v.getX());
	return w;
}

BLASMatrix& BLASMatrix::operator = (const BLASMatrix& mat)
{
	w_col0 = mat.w_col0;
	w_col1 = mat.w_col1;
	w_col2 = mat.w_col2;
	return *this;
}

BLASMatrix BLASMatrix::operator + ()
{
	return *this;
}

BLASMatrix BLASMatrix::operator - ()
{
	w_col0 = -w_col0;
	w_col1 = -w_col1;
	w_col2 = -w_col2;
	return *this;
}

BLASMatrix& BLASMatrix::operator += (const BLASMatrix& mat)
{
	*this = *this + mat;
	return *this;
}

BLASMatrix& BLASMatrix::operator -= (const BLASMatrix& mat)
{
	*this = *this - mat;
	return *this;
}

BLASMatrix& BLASMatrix::operator *= (const float k)
{
	*this = *this * k;
	return *this;
}

BLASMatrix& BLASMatrix::operator *= (const BLASMatrix& mat)
{
	*this = *this * mat;
	return *this;
}

BLASMatrix operator +(const BLASMatrix& mat1, const BLASMatrix& mat2)
{
	return BLASMatrix(
		mat1.getCol0() + mat2.getCol0(),
		mat1.getCol1() + mat2.getCol1(),
		mat1.getCol2() + mat2.getCol2()
	);
}

BLASMatrix operator -(const BLASMatrix& mat1, const BLASMatrix& mat2)
{
	return BLASMatrix(
		mat1.getCol0() - mat2.getCol0(),
		mat1.getCol1() - mat2.getCol1(),
		mat1.getCol2() - mat2.getCol2()
	);
}

BLASMatrix operator *(const BLASMatrix& mat1, const BLASMatrix& mat2)
{
	return BLASMatrix(
		mat1 * mat2.getCol0(),
		mat1 * mat2.getCol1(),
		mat1 * mat2.getCol2()
	);
}

BLASMatrix operator *(const float k, const BLASMatrix& mat)
{
	return BLASMatrix(
		mat.getCol0()*k,
		mat.getCol1() *k,
		mat.getCol2() *k
	);
}

BLASMatrix operator *(const BLASMatrix& mat, const float k)
{
	return BLASMatrix(
		mat.getCol0() * k,
		mat.getCol1() * k,
		mat.getCol2() * k
	);
}

BLASVector operator *(const BLASMatrix& mat ,const BLASVector& vec)
{
	return BLASVector(
		(((mat.getCol0().getX() * vec.getX()) + (mat.getCol1().getX() * vec.getY())) + (mat.getCol2().getX() * vec.getZ())),
		(((mat.getCol0().getY() * vec.getX()) + (mat.getCol1().getY() * vec.getY())) + (mat.getCol2().getY() * vec.getZ())),
		(((mat.getCol0().getZ() * vec.getX()) + (mat.getCol1().getZ() * vec.getY())) + (mat.getCol2().getZ() * vec.getZ()))
	);
}