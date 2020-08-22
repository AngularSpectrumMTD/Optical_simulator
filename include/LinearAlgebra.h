#pragma once

#include <math.h>
//#include <vector>

class LAVector;
class LAMatrix;

class LAVector
{
private:
	float w_x;
	float w_y;
	float w_z;

public:
	
	//constructor
	LAVector(): w_x(0), w_y(0), w_z(0){}
	LAVector(float x): w_x(x), w_y(x), w_z(x){}
	LAVector(float x, float y, float z): w_x(x), w_y(y), w_z(z){}

	//operator
	LAVector& operator = (const LAVector& v);
	LAVector operator +();
	LAVector operator - ();
	LAVector& operator += (const LAVector& v);
	LAVector& operator -= (const LAVector& v);
	LAVector& operator *= (const float k);
	LAVector& operator /= (const float k);

	float getX() const { return w_x; }
	float getY() const { return w_y; }
	float getZ() const { return w_z; }

	static inline const LAVector getEx();
	static inline const LAVector getEy();
	static inline const LAVector getEz();

	float getElem(int idx) const { return *(&w_x + idx); }

	void setX(float x) { w_x = x; }
	void setY(float y) { w_y = y; } 
	void setZ(float z) { w_z = z; }

	LAVector& setElem(int idx, float val)
	{
		*(&w_x + idx) = val;
		return *this;
	}
};
inline const LAVector LAVector::getEx() { return LAVector(1.0, 0.0, 0.0); }
inline const LAVector LAVector::getEy() { return LAVector(0.0, 1.0, 0.0); }
inline const LAVector LAVector::getEz() { return LAVector(0.0, 0.0, 1.0); }

inline float dot(const LAVector & vec0, const LAVector & vec1)
{
	float result;
	result = (vec0.getX() * vec1.getX());
	result = (result + (vec0.getY() * vec1.getY()));
	result = (result + (vec0.getZ() * vec1.getZ()));
	return result;
}

inline float lengthSqr(const LAVector & vec)
{
	float result;
	result = (vec.getX() * vec.getX());
	result = (result + (vec.getY() * vec.getY()));
	result = (result + (vec.getZ() * vec.getZ()));
	return result;
}

inline float length(const LAVector& vec)
{
	return sqrtf(lengthSqr(vec));
}

inline const LAVector normalize(const LAVector& vec)
{
	float lenSqr, lenInv;
	lenSqr = lengthSqr(vec);
	lenInv = (1.0f / sqrtf(lenSqr));
	return LAVector(
		(vec.getX() * lenInv),
		(vec.getY() * lenInv),
		(vec.getZ() * lenInv)
	);
}

inline const LAVector cross(const LAVector& vec0, const LAVector& vec1)
{
	return LAVector(
		((vec0.getY() * vec1.getZ()) - (vec0.getZ() * vec1.getY())),
		((vec0.getZ() * vec1.getX()) - (vec0.getX() * vec1.getZ())),
		((vec0.getX() * vec1.getY()) - (vec0.getY() * vec1.getX()))
	);
}

inline const LAVector minPerElem(const LAVector& vec0, const LAVector& vec1)
{
	return LAVector(
		(vec0.getX() < vec1.getX()) ? vec0.getX() : vec1.getX(),
		(vec0.getY() < vec1.getY()) ? vec0.getY() : vec1.getY(),
		(vec0.getZ() < vec1.getZ()) ? vec0.getZ() : vec1.getZ()
	);
}

inline const LAVector maxPerElem(const LAVector& vec0, const LAVector& vec1)
{
	return LAVector(
		(vec0.getX() > vec1.getX()) ? vec0.getX() : vec1.getX(),
		(vec0.getY() > vec1.getY()) ? vec0.getY() : vec1.getY(),
		(vec0.getZ() > vec1.getZ()) ? vec0.getZ() : vec1.getZ()
	);
}

LAVector operator +(const LAVector& u, const LAVector& v);
LAVector operator -(const LAVector& u, const LAVector& v);
float operator *(const LAVector& u, const LAVector& v);//innerproduct
LAVector operator *(const float k, const LAVector& v);
LAVector operator *(const LAVector& v, const float k);
LAVector operator /(const LAVector& v, const float k);
LAVector operator &&(const LAVector& u, const LAVector& v);//outerproduct

class LAMatrix {
private:
	LAVector w_col0;
	LAVector w_col1;
	LAVector w_col2;

public:
	//constructor
	LAMatrix()
	{
		w_col0 = LAVector();
		w_col1 = LAVector();
		w_col2 = LAVector();
	}
	LAMatrix(const LAMatrix & mat)
	{
		w_col0 = mat.w_col0;
		w_col1 = mat.w_col1;
		w_col2 = mat.w_col2;
	}
	LAMatrix(const LAVector& col0, const LAVector& col1, const LAVector& col2)
	{
		w_col0 = col0;
		w_col1 = col1;
		w_col2 = col2;
	}

	LAVector getCol0() const { return w_col0; }
	LAVector getCol1() const{ return w_col1; }
	LAVector getCol2() const{ return w_col2; }
	LAVector getCol(int col) const { return *(&w_col0 + col); }
	LAVector getRow(int row) const 
	{
		return LAVector(w_col0.getElem(row), w_col1.getElem(row), w_col2.getElem(row));
	}
	float getElem(int col, int row) const
	{
		return this->getCol(col).getElem(row);
	}

	LAMatrix& setCol0(const LAVector& col0)
	{
		w_col0 = col0;
		return *this;
	}
	LAMatrix& setCol1(const LAVector& col1)
	{
		w_col1 = col1;
		return *this;
	}
	LAMatrix& setCol2(const LAVector& col2)
	{
		w_col2 = col2;
		return *this;
	}
	LAMatrix& setCol(int col, const LAVector& vec)
	{
		*(&w_col0 + col) = vec;
		return *this;
	}
	LAMatrix& setRow(int row, const LAVector& vec)
	{
		w_col0.setElem(row, vec.getElem(0));
		w_col1.setElem(row, vec.getElem(1));
		w_col2.setElem(row, vec.getElem(2));
		return *this;
	}
	LAMatrix& setElem(int col, int row, const float val)
	{
		LAVector tmp;
		tmp = getCol(col);
		tmp.setElem(row, val);
		setCol(col, tmp);
		return *this;
	}

	static inline const LAMatrix identity();
	static inline const LAMatrix rotationX(const float rad);
	static inline const LAMatrix rotationY(const float rad);
	static inline const LAMatrix rotationZ(const float rad);
	static inline const LAMatrix rotation(const float rad, const LAVector & unit_dir);

	//operator
	LAMatrix& operator = (const LAMatrix& mat);
	LAMatrix operator +();
	LAMatrix operator - ();
	LAMatrix& operator += (const LAMatrix& mat);
	LAMatrix& operator -= (const LAMatrix& mat);
	LAMatrix& operator *= (const float k);
	LAMatrix& operator *= (const LAMatrix& mat);
};

LAMatrix operator +(const LAMatrix & mat1, const LAMatrix& mat2) ;
LAMatrix operator -(const LAMatrix& mat1, const LAMatrix& mat2);
LAMatrix operator *(const LAMatrix& mat1, const LAMatrix& mat2);
LAMatrix operator *(const float k, const LAMatrix& mat);
LAMatrix operator *(const LAMatrix& mat, const float k);
LAVector operator *(const LAMatrix& mat, const LAVector& vec);

inline const LAMatrix LAMatrix::identity()
{
	return LAMatrix(
		LAVector::getEx(),
		LAVector::getEy(),
		LAVector::getEz()
	);
}

inline const LAMatrix LAMatrix::rotationX(const float rad)
{
	float sinval, cosval;
	sinval = sinf(rad);
	cosval = cosf(rad);
	return LAMatrix(
		LAVector::getEx(),
		LAVector(0.0f, cosval, sinval),
		LAVector(0.0f, -sinval, cosval)
	);
}
inline const LAMatrix LAMatrix::rotationY(const float rad)
{
	float sinval, cosval;
	sinval = sinf(rad);
	cosval = cosf(rad);
	return LAMatrix(
		LAVector(cosval, 0.0f, -sinval),
		LAVector::getEy(),
		LAVector(sinval, 0.0f, cosval)
	);
}
inline const LAMatrix LAMatrix::rotationZ(const float rad)
{
	float sinval, cosval;
	sinval = sinf(rad);
	cosval = cosf(rad);
	return LAMatrix(
		LAVector(cosval, sinval, 0.0f),
		LAVector(-sinval, cosval, 0.0f),
		LAVector::getEz()
	);
}
inline const LAMatrix LAMatrix::rotation(const float rad, const LAVector& unit_dir)
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

	return LAMatrix(
		LAVector(((x * x) * cosvalone) + cosval, (xy * cosvalone) + (z * sinval), (zx * cosvalone) - (y * sinval)),
		LAVector((xy * cosvalone) - (z * sinval), (((y * y) * cosvalone) + cosval), ((yz * cosvalone) + (x * sinval))),
		LAVector(((zx * cosvalone) + (y * sinval)), ((yz * cosvalone) - (x * sinval)), (((z * z) * cosvalone) + cosval))
	);
}

inline const LAMatrix transpose(const LAMatrix& mat)
{
	return LAMatrix(
		LAVector(mat.getCol0().getX(), mat.getCol1().getX(), mat.getCol2().getX()),
		LAVector(mat.getCol0().getY(), mat.getCol1().getY(), mat.getCol2().getY()),
		LAVector(mat.getCol0().getZ(), mat.getCol1().getZ(), mat.getCol2().getZ())
	);
}

inline const LAMatrix inverse(const LAMatrix& mat)
{
	LAVector tmp0, tmp1, tmp2;
	float detinv;
	tmp0 = cross(mat.getCol1(), mat.getCol2());
	tmp1 = cross(mat.getCol2(), mat.getCol0());
	tmp2 = cross(mat.getCol0(), mat.getCol1());
	detinv = (1.0f / dot(mat.getCol2(), tmp2));
	return LAMatrix(
		LAVector((tmp0.getX() * detinv), (tmp1.getX() * detinv), (tmp2.getX() * detinv)),
		LAVector((tmp0.getY() * detinv), (tmp1.getY() * detinv), (tmp2.getY() * detinv)),
		LAVector((tmp0.getZ() * detinv), (tmp1.getZ() * detinv), (tmp2.getZ() * detinv))
	);
}