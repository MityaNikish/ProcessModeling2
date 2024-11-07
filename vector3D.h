#pragma once
#include <functional>

enum class Num
{
	one, two, three
};

class Vector3D
{
	double data_[3];

public:
	Vector3D() noexcept;
	Vector3D(double val1, double val2, double val3) noexcept;

	Vector3D operator+(const Vector3D& other) const noexcept;
	Vector3D operator-(const Vector3D& other) const noexcept;
	Vector3D operator*(double value) const noexcept;
	Vector3D operator/(double value) const noexcept;
	double operator*(const Vector3D& other) const noexcept;
	double& operator[](Num index) noexcept;
	double operator[](Num index) const noexcept;
	Vector3D operator()(std::function<double(double)> func) const noexcept;

	Vector3D abs() const noexcept;
};