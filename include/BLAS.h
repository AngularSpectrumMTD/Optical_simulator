#pragma once

#include <math.h>
//#include <vector>

class BLASVector;
class BLASMatrix;

class BLASVector
{
private:
	float w_x;
	float w_y;
	float w_z;

public:
	
	//constructor
	BLASVector(): w_x(0), w_y(0), w_z(0){}
	BLASVector(float x): w_x(x), w_y(x), w_z(x){}
	BLASVector(float x, float y, float z): w_x(x), w_y(y), w_z(z){}

	//operator
	BLASVector& operator = (const BLASVector& v);
	BLASVector operator +();
	BLASVector operator - ();
	BLASVector& operator += (const BLASVector& v);
	BLASVector& operator -= (const BLASVector& v);
	BLASVector& operator *= (const float k);
	BLASVector& operator /= (const float k);

	float getX() const { return w_x; }
	float getY() const { return w_y; }
	float getZ() const { return w_z; }

	static inline const BLASVector getEx();
	static inline const BLASVector getEy();
	static inline const BLASVector getEz();

	float getElem(int idx) const { return *(&w_x + idx); }

	void setX(float x) { w_x = x; }
	void setY(float y) { w_y = y; } 
	void setZ(float z) { w_z = z; }

	BLASVector& setElem(int idx, float val)
	{
		*(&w_x + idx) = val;
		return *this;
	}
};
inline const BLASVector BLASVector::getEx() { return BLASVector(1.0, 0.0, 0.0); }
inline const BLASVector BLASVector::getEy() { return BLASVector(0.0, 1.0, 0.0); }
inline const BLASVector BLASVector::getEz() { return BLASVector(0.0, 0.0, 1.0); }

inline float dot(const BLASVector & vec0, const BLASVector & vec1)
{
	float result;
	result = (vec0.getX() * vec1.getX());
	result = (result + (vec0.getY() * vec1.getY()));
	result = (result + (vec0.getZ() * vec1.getZ()));
	return result;
}

inline float lengthSqr(const BLASVector & vec)
{
	float result;
	result = (vec.getX() * vec.getX());
	result = (result + (vec.getY() * vec.getY()));
	result = (result + (vec.getZ() * vec.getZ()));
	return result;
}

inline float length(const BLASVector& vec)
{
	return sqrtf(lengthSqr(vec));
}

inline const BLASVector normalize(const BLASVector& vec)
{
	float lenSqr, lenInv;
	lenSqr = lengthSqr(vec);
	lenInv = (1.0f / sqrtf(lenSqr));
	return BLASVector(
		(vec.getX() * lenInv),
		(vec.getY() * lenInv),
		(vec.getZ() * lenInv)
	);
}

inline const BLASVector cross(const BLASVector& vec0, const BLASVector& vec1)
{
	return BLASVector(
		((vec0.getY() * vec1.getZ()) - (vec0.getZ() * vec1.getY())),
		((vec0.getZ() * vec1.getX()) - (vec0.getX() * vec1.getZ())),
		((vec0.getX() * vec1.getY()) - (vec0.getY() * vec1.getX()))
	);
}

inline const BLASVector minPerElem(const BLASVector& vec0, const BLASVector& vec1)
{
	return BLASVector(
		(vec0.getX() < vec1.getX()) ? vec0.getX() : vec1.getX(),
		(vec0.getY() < vec1.getY()) ? vec0.getY() : vec1.getY(),
		(vec0.getZ() < vec1.getZ()) ? vec0.getZ() : vec1.getZ()
	);
}

inline const BLASVector maxPerElem(const BLASVector& vec0, const BLASVector& vec1)
{
	return BLASVector(
		(vec0.getX() > vec1.getX()) ? vec0.getX() : vec1.getX(),
		(vec0.getY() > vec1.getY()) ? vec0.getY() : vec1.getY(),
		(vec0.getZ() > vec1.getZ()) ? vec0.getZ() : vec1.getZ()
	);
}

BLASVector operator +(const BLASVector& u, const BLASVector& v);
BLASVector operator -(const BLASVector& u, const BLASVector& v);
float operator *(const BLASVector& u, const BLASVector& v);//innerproduct
BLASVector operator *(const float k, const BLASVector& v);
BLASVector operator *(const BLASVector& v, const float k);
BLASVector operator /(const BLASVector& v, const float k);
BLASVector operator &&(const BLASVector& u, const BLASVector& v);//outerproduct

class BLASMatrix {
private:
	BLASVector w_col0;
	BLASVector w_col1;
	BLASVector w_col2;

public:
	//constructor
	BLASMatrix()
	{
		w_col0 = BLASVector();
		w_col1 = BLASVector();
		w_col2 = BLASVector();
	}
	BLASMatrix(const BLASMatrix & mat)
	{
		w_col0 = mat.w_col0;
		w_col1 = mat.w_col1;
		w_col2 = mat.w_col2;
	}
	BLASMatrix(const BLASVector& col0, const BLASVector& col1, const BLASVector& col2)
	{
		w_col0 = col0;
		w_col1 = col1;
		w_col2 = col2;
	}

	BLASVector getCol0() const { return w_col0; }
	BLASVector getCol1() const{ return w_col1; }
	BLASVector getCol2() const{ return w_col2; }
	BLASVector getCol(int col) const { return *(&w_col0 + col); }
	BLASVector getRow(int row) const 
	{
		return BLASVector(w_col0.getElem(row), w_col1.getElem(row), w_col2.getElem(row));
	}
	float getElem(int col, int row) const
	{
		return this->getCol(col).getElem(row);
	}

	BLASMatrix& setCol0(const BLASVector& col0)
	{
		w_col0 = col0;
		return *this;
	}
	BLASMatrix& setCol1(const BLASVector& col1)
	{
		w_col1 = col1;
		return *this;
	}
	BLASMatrix& setCol2(const BLASVector& col2)
	{
		w_col2 = col2;
		return *this;
	}
	BLASMatrix& setCol(int col, const BLASVector& vec)
	{
		*(&w_col0 + col) = vec;
		return *this;
	}
	BLASMatrix& setRow(int row, const BLASVector& vec)
	{
		w_col0.setElem(row, vec.getElem(0));
		w_col1.setElem(row, vec.getElem(1));
		w_col2.setElem(row, vec.getElem(2));
		return *this;
	}
	BLASMatrix& setElem(int col, int row, const float val)
	{
		BLASVector tmp;
		tmp = getCol(col);
		tmp.setElem(row, val);
		setCol(col, tmp);
		return *this;
	}

	static inline const BLASMatrix identity();
	static inline const BLASMatrix rotationX(const float rad);
	static inline const BLASMatrix rotationY(const float rad);
	static inline const BLASMatrix rotationZ(const float rad);
	static inline const BLASMatrix rotation(const float rad, const BLASVector & unit_dir);

	//operator
	BLASMatrix& operator = (const BLASMatrix& mat);
	BLASMatrix operator +();
	BLASMatrix operator - ();
	BLASMatrix& operator += (const BLASMatrix& mat);
	BLASMatrix& operator -= (const BLASMatrix& mat);
	BLASMatrix& operator *= (const float k);
	BLASMatrix& operator *= (const BLASMatrix& mat);
};

BLASMatrix operator +(const BLASMatrix & mat1, const BLASMatrix& mat2) ;
BLASMatrix operator -(const BLASMatrix& mat1, const BLASMatrix& mat2);
BLASMatrix operator *(const BLASMatrix& mat1, const BLASMatrix& mat2);
BLASMatrix operator *(const float k, const BLASMatrix& mat);
BLASMatrix operator *(const BLASMatrix& mat, const float k);
BLASVector operator *(const BLASMatrix& mat, const BLASVector& vec);

inline const BLASMatrix BLASMatrix::identity()
{
	return BLASMatrix(
		BLASVector::getEx(),
		BLASVector::getEy(),
		BLASVector::getEz()
	);
}

inline const BLASMatrix BLASMatrix::rotationX(const float rad)
{
	float sinval, cosval;
	sinval = sinf(rad);
	cosval = cosf(rad);
	return BLASMatrix(
		BLASVector::getEx(),
		BLASVector(0.0f, cosval, sinval),
		BLASVector(0.0f, -sinval, cosval)
	);
}
inline const BLASMatrix BLASMatrix::rotationY(const float rad)
{
	float sinval, cosval;
	sinval = sinf(rad);
	cosval = cosf(rad);
	return BLASMatrix(
		BLASVector(cosval, 0.0f, -sinval),
		BLASVector::getEy(),
		BLASVector(sinval, 0.0f, cosval)
	);
}
inline const BLASMatrix BLASMatrix::rotationZ(const float rad)
{
	float sinval, cosval;
	sinval = sinf(rad);
	cosval = cosf(rad);
	return BLASMatrix(
		BLASVector(cosval, sinval, 0.0f),
		BLASVector(-sinval, cosval, 0.0f),
		BLASVector::getEz()
	);
}
inline const BLASMatrix BLASMatrix::rotation(const float rad, const BLASVector& unit_dir)
{
	float sinval, cosval;
	sinval = sinf(rad);
	cosval = cosf(rad);

	float x, y, z;
	x = unit_dir.getX();
	y = unit_dir.getY();
	z = unit_dir.getZ();
	float xy, yz, zx;
	xy = x * y;
	yz = y * z;
	zx = z * x;
	float cosvalone = 1 - cosval;

	return BLASMatrix(
		BLASVector(((x * x) * cosvalone) + cosval, (xy * cosvalone) + (z * sinval), (zx * cosvalone) - (y * sinval)),
		BLASVector((xy * cosvalone) - (z * sinval), (((y * y) * cosvalone) + cosval), ((yz * cosvalone) + (x * sinval))),
		BLASVector(((zx * cosvalone) + (y * sinval)), ((yz * cosvalone) - (x * sinval)), (((z * z) * cosvalone) + cosval))
	);
}

inline const BLASMatrix transpose(const BLASMatrix& mat)
{
	return BLASMatrix(
		BLASVector(mat.getCol0().getX(), mat.getCol1().getX(), mat.getCol2().getX()),
		BLASVector(mat.getCol0().getY(), mat.getCol1().getY(), mat.getCol2().getY()),
		BLASVector(mat.getCol0().getZ(), mat.getCol1().getZ(), mat.getCol2().getZ())
	);
}

inline const BLASMatrix inverse(const BLASMatrix& mat)
{
	BLASVector tmp0, tmp1, tmp2;
	float detinv;
	tmp0 = cross(mat.getCol1(), mat.getCol2());
	tmp1 = cross(mat.getCol2(), mat.getCol0());
	tmp2 = cross(mat.getCol0(), mat.getCol1());
	detinv = (1.0f / dot(mat.getCol2(), tmp2));
	return BLASMatrix(
		BLASVector((tmp0.getX() * detinv), (tmp1.getX() * detinv), (tmp2.getX() * detinv)),
		BLASVector((tmp0.getY() * detinv), (tmp1.getY() * detinv), (tmp2.getY() * detinv)),
		BLASVector((tmp0.getZ() * detinv), (tmp1.getZ() * detinv), (tmp2.getZ() * detinv))
	);
}