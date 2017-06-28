#ifndef _VECTOR_
#define _VECTOR_
struct Vec3 {
	double dx;
	double dy;
	double dz;
	Vec3();
	Vec3(const Vec3 &v);
	Vec3(const double &x, const double &y, const double &z);
	Vec3(const double &data);
	const Vec3 operator=(const Vec3 &v);
	const Vec3 operator-(const Vec3 &v) const;
	const Vec3 operator+(const Vec3 &v) const;
	const Vec3 operator*(const double &coeff) const;
};


struct Vec2 {
	double dx;
	double dy;
	Vec2();
	Vec2(const Vec2 &v);
	Vec2(const double &x, const double &y);
	const Vec2 operator=(const Vec2 &v);
	const Vec2 operator-(const Vec2 &v) const;
	const Vec2 operator+(const Vec2 &v) const;
	const Vec2 operator*(const double &coeff) const;
};

#endif
