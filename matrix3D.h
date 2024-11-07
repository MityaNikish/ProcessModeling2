#pragma once
#include "vector3D.h"

class Matrix3D
{
	Vector3D data_[3];

public:
	Matrix3D() noexcept  = default;
	Matrix3D(const Vector3D& vec1, const Vector3D& vec2, const Vector3D& vec3) noexcept;

	Matrix3D operator+(const Matrix3D& other) const noexcept;
	Matrix3D operator*(double value) const noexcept;
	Matrix3D operator/(double value) const noexcept;
	Vector3D operator*(const Vector3D& vec) const noexcept;
	Matrix3D operator*(const Matrix3D& other) const noexcept;
	Vector3D& operator[](Num index) noexcept;
	Vector3D operator[](Num index) const noexcept;
	Matrix3D operator()(std::function<double(double)> func) const noexcept;

	Matrix3D abs() const noexcept;

	double getMaxElem() noexcept;
};