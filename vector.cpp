#include "vector.hpp"

Vec3::Vec3():dx(0), dy(0), dz(0) {}
Vec3::Vec3(const Vec3 &v):dx(v.dx), dy(v.dy), dz(v.dz) {}
Vec3::Vec3(const double &x, const double &y, const double &z):dx(x), dy(y), dz(z) {}
Vec3::Vec3(const double &data):dx(data), dy(data), dz(data) {}

const Vec3 Vec3::operator=(const Vec3 &v) {
	dx = v.dx;
	dy = v.dy;
	dz = v.dz;
	return *this;
}

const Vec3 Vec3::operator-(const Vec3 &v) const {
	Vec3 vtemp;
	vtemp.dx = dx - v.dx;
	vtemp.dy = dy - v.dy;
	vtemp.dz = dz - v.dz;
	return vtemp;
}

const Vec3 Vec3::operator+(const Vec3 &v) const {
	Vec3 vtemp;
	vtemp.dx = dx + v.dx;
	vtemp.dy = dy + v.dy;
	vtemp.dz = dz + v.dz;
	return vtemp;
}

const Vec3 Vec3::operator*(const double &coeff) const {
	Vec3 vtemp;
	vtemp.dx = dx * coeff;
	vtemp.dy = dy * coeff;
	vtemp.dz = dz * coeff;
	return vtemp;
}


Vec2::Vec2():dx(0), dy(0) {}
Vec2::Vec2(const Vec2 &v):dx(v.dx), dy(v.dy) {}
Vec2::Vec2(const double &x, const double &y):dx(x), dy(y) {}

const Vec2 Vec2::operator=(const Vec2 &v) {
	dx = v.dx;
	dy = v.dy;
	return *this;
}

const Vec2 Vec2::operator-(const Vec2 &v) const {
	Vec2 vtemp;
	vtemp.dx = dx - v.dx;
	vtemp.dy = dy - v.dy;
	return vtemp;
}

const Vec2 Vec2::operator+(const Vec2 &v) const {
	Vec2 vtemp;
	vtemp.dx = dx + v.dx;
	vtemp.dy = dy + v.dy;
	return vtemp;
}
const Vec2 Vec2::operator*(const double &coeff) const {
	Vec2 vtemp;
	vtemp.dx = dx * coeff;
	vtemp.dy = dy * coeff;
	return vtemp;
}
