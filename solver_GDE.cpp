#include "solver_GDE.h"
#include <fstream>
//#include <iostream>
#include <cmath>

#define NDEBUG
#include <cassert>


namespace 
{
	using enum Num;

	void write(std::filesystem::path file_path, double* arr, size_t size)
	{
		std::ofstream fout(file_path, std::ios::trunc | std::ios::binary);
		fout.write((char*)arr, sizeof(double) * size);
	}

	double minmod(double x, double y)
	{
		const double sign_x = std::copysign(1.0, x);
		return sign_x * std::max(0.0, std::min(abs(x), y * sign_x));
	}

	Vector3D minmod(Vector3D x, Vector3D y)
	{
		return Vector3D(
			minmod(x[one], y[one]),
			minmod(x[two], y[two]),
			minmod(x[three], y[three])
		);
	}
}


//	Показатель адиабаты
double GasDynamicsEquation::gamma = 1.4;	//	1.4 - Воздух

//	Искусственная вязкость Лапидуса
double GasDynamicsEquation::artificial_viscosity = 2.0;	//	2.0 - по умолчанию


GasDynamicsEquation::GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::vector<double>& S) :
	_expanse_grid(expanse_grid),
	_time_grid(time_grid),
	_start_condition(std::move(start_condition)),
	_borderline_condition(std::move(borderline_condition)),
	_ro(2, _expanse_grid.nodes),
	_u(2, _expanse_grid.nodes),
	_p(2, _expanse_grid.nodes),
	_S(S),
	_alpha(_time_grid.tau / _expanse_grid.h)
{ 
	assert(_S.size() == expanse_grid.nodes && "The cross-section array is not of the right lenght!");
}


GasDynamicsEquation::GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::function<double(double)>& S_func) :
	_expanse_grid(expanse_grid),
	_time_grid(time_grid),
	_start_condition(std::move(start_condition)),
	_borderline_condition(std::move(borderline_condition)),
	_ro(2, _expanse_grid.nodes),
	_u(2, _expanse_grid.nodes),
	_p(2, _expanse_grid.nodes),
	_S(_expanse_grid.nodes),
	_alpha(_time_grid.tau / _expanse_grid.h)
{
	for (size_t i = 0; i < _expanse_grid.nodes; i++)
	{
		(*const_cast<std::vector<double>*>(&_S))[i] = S_func(i * _expanse_grid.h);
		//std::cout << _S[i] << "\n";
	}
}


//	Моделирует УГД
void GasDynamicsEquation::solving() noexcept
{
	initConditions();

	for (size_t n = 0; n < _time_grid.nodes - 1; ++n)
	{
		bool not_satisfy_condition = false;
		for (size_t j = 1; j < _expanse_grid.nodes - 1; ++j)
		{
			Vector3D U_star = U_future(n, j);

			const double ro_new = U_star[one];
			const double u_new = U_star[two] / ro_new;
			const double E_new = U_star[three] / ro_new;
			const double p_new = (E_new - u_new * u_new / 2) * ro_new * (gamma - 1);

			_ro.getElement(n + 1, j) = ro_new;
			_u.getElement(n + 1, j) = u_new;
			_p.getElement(n + 1, j) = p_new;

			if (!not_satisfy_condition && !chekStopConditions(n, j, 1.0e-6))
			{
				not_satisfy_condition = !not_satisfy_condition;
			}
		}
		if (!not_satisfy_condition)
		{
			break;
		}
	}

	postProcessing();
}


//	Записывает массив данных плотности в файл
void GasDynamicsEquation::writeRO(const std::filesystem::path& file_path)  const
{
	write(file_path, _ro[0].getData(), _expanse_grid.nodes);
}

//	Записывает массив данных скорости в файл
void GasDynamicsEquation::writeU(const std::filesystem::path& file_path)  const
{
	write(file_path, _u[0].getData(), _expanse_grid.nodes);
}

//	Записывает массив данных давления в файл
void GasDynamicsEquation::writeP(const std::filesystem::path& file_path)  const
{
	write(file_path, _p[0].getData(), _expanse_grid.nodes);
}


//	Вернуть массив данных плотности
std::vector<double> GasDynamicsEquation::getRO() const
{
	Matrix slice = _ro[0];
	return std::vector<double>(slice.getData(), slice.getData() + _expanse_grid.nodes);
}

//	Вернуть массив данных скорости
std::vector<double> GasDynamicsEquation::getU() const
{
	Matrix slice = _u[0];
	return std::vector<double>(slice.getData(), slice.getData() + _expanse_grid.nodes);
}

//	Вернуть массив данных давления
std::vector<double> GasDynamicsEquation::getP() const
{
	Matrix slice = _p[0];
	return std::vector<double>(slice.getData(), slice.getData() + _expanse_grid.nodes);
}


//	Инициализирует сетки начальными и граничными условиями
void GasDynamicsEquation::initConditions() noexcept
{
	for (size_t j = 0; j < _expanse_grid.nodes; ++j)
	{
		_ro.getElement(static_cast<size_t>(0), j) = _S[j] * _start_condition.start_ro(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_u.getElement(static_cast<size_t>(0), j) = _start_condition.start_u(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_p.getElement(static_cast<size_t>(0), j) = _S[j] * _start_condition.start_p(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
	}

	for (size_t n = 0; n < 2; ++n)
	{
		_ro.getElement(n, static_cast<size_t>(0)) = _S[0] * _borderline_condition.left_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, static_cast<size_t>(0)) = _borderline_condition.left_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, static_cast<size_t>(0)) = _S[0] * _borderline_condition.left_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));

		_ro.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.right_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, _expanse_grid.nodes - 1) = _borderline_condition.right_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.right_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
	}
}

//	Проводит послевычеслительные операции
void GasDynamicsEquation::postProcessing() noexcept
{
	for (size_t n = 0; n < 2; ++n)
	{
		for (size_t j = 0; j < _expanse_grid.nodes; ++j)
		{
			_ro.getElement(n, j) /= _S[j];
			_p.getElement(n, j) /= _S[j];
		}
	}
}


//	Проверка на удовлетворение условий остановки вычислений
bool GasDynamicsEquation::chekStopConditions(size_t n, size_t j, double eps) const noexcept
{
	/*
	//	Условие остановки по t
	const double ro_ = _ro.getElement(n, j);
	const double ro_future = _ro.getElement(n + 1, j);

	const double div = (ro_future - ro_) / _time_grid.tau;

	return abs(div) * _S[j] / ro_ < eps;
	*/

	//	Условие остановки по x
	const double ro_ = _ro.getElement(n, j);
	const double ro_next = _ro.getElement(n, j + 1);

	const double u_ = _ro.getElement(n, j);
	const double u_next = _ro.getElement(n, j + 1);

	const double div = (ro_next * u_next - ro_ * u_) / _expanse_grid.h;

	return abs(div) * _S[j] / ro_ < eps;
	

	//return false;
}


//	Удельная полная энергия
double GasDynamicsEquation::E(size_t n, size_t j) const noexcept
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	return p / ro / (gamma - 1) + u * u / 2;
}

//	Удельная польная энтальпия
double GasDynamicsEquation::H(size_t n, size_t j) const noexcept
{
	const double ro = _ro.getElement(n, j);
	const double p = _p.getElement(n, j);

	return E(n, j) + p / ro;
}

double GasDynamicsEquation::g(double z) const noexcept
{
	const double epsl = 0.01 / 0.2;
	return z >= epsl ? abs(z) : (z * z + epsl * epsl) / epsl / 2;
}


//	Консервативная переменная
Vector3D GasDynamicsEquation::U(size_t n, size_t j) const noexcept
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);

	return Vector3D (ro, ro * u, ro * E(n, j));
}

//	Проводит вычисление консервативной переменной, следующей по временной сетке
Vector3D GasDynamicsEquation::U_future(size_t n, size_t j) const noexcept
{
	const Vector3D U_ = U(n, j);

	const Matrix3D R_pref = R(n, j - 1);
	const Matrix3D _R = R(n, j);
	const Matrix3D R_next = R(n, j + 1);

	const Matrix3D R_half_pref = (_R + R_next) / 2;
	const Matrix3D R_half_next = (R_pref + _R) / 2;

	const Vector3D F_pref = F(n, j - 1);
	const Vector3D F_ = F(n, j);
	const Vector3D F_next = F(n, j + 1);

	const Vector3D G_half_pref = G_half(n, j - 1);
	const Vector3D G_half_next = G_half(n, j);

	const Vector3D F_half_pref = (F_ + F_pref - R_half_pref * G_half_pref) / 2;
	const Vector3D F_half_next = (F_next + F_ - R_half_next * G_half_next) / 2;


	return U_ - (F_half_next - F_half_pref) * _alpha + Q(n, j) * _time_grid.tau;
}


//	Значение потока
Vector3D GasDynamicsEquation::F(size_t n, size_t j) const noexcept
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	double ro_u = ro * u;

	return Vector3D(ro_u, ro_u * u + p, ro_u * E(n, j) + u * p);
}

//	Источник
Vector3D GasDynamicsEquation::Q(size_t n, size_t j) const noexcept
{
	const double p = _p.getElement(n, j);
	double dS;
	if (j == 0)
	{
		dS = (_S[2] - _S[0]) / 2 / _expanse_grid.h;
	}
	else if (j == _expanse_grid.nodes - 1)
	{
		dS = (_S[_expanse_grid.nodes - 1] - _S[_expanse_grid.nodes - 3]) / 2 / _expanse_grid.h;
	}
	else
	{
		dS = (_S[j + 1] - _S[j - 1]) / 2 / _expanse_grid.h;
	}

	return Vector3D(0, p / _S[j] * dS, 0);
}


//	d_ksi(n, j+1/2) = minmod(d_W(n, j+1/2), d_W(n, j-1/2)) + minmod(d_W(n, j+1/2), d_W(n, j+3/2)) - d_W(n, j+1/2)
Vector3D GasDynamicsEquation::delta_ksi_half(size_t n, size_t j) const noexcept
{
	const Vector3D dW_half_pref = delta_W_half(n, j - 1);
	const Vector3D dW_half_ = delta_W_half(n, j);
	const Vector3D dW_half_next = delta_W_half(n, j + 1);

	return minmod(dW_half_, dW_half_pref) + minmod(dW_half_, dW_half_next) - dW_half_;
}

//	d_W(n, j+1/2) = L(n, j+1/2) * (U_(n, j+1) - U(n, j))
Vector3D GasDynamicsEquation::delta_W_half(size_t n, size_t j) const noexcept
{
	const Matrix3D L_ = L(n, j);
	const Matrix3D L_next = L(n, j + 1);
	const Matrix3D L_half = (L_next + L_) / 2;

	const Vector3D U_ = U(n, j);
	const Vector3D U_next = U(n, j + 1);
	const Vector3D U_diff = U_next - U_;

	return L_half * U_diff;
}

//	Ф(n, j+1/2) = tau / h * |D(n, j+1/2)^2| * d_ksi(n, j+1/2) + g[D(n, j+1/2)] * (d_W(n, j+1/2) - d_ksi(n, j+1/2))
Vector3D GasDynamicsEquation::G_half(size_t n, size_t j) const noexcept
{
	const Matrix3D D_ = D(n, j);
	const Matrix3D D_next = D(n, j + 1);
	const Matrix3D D_half_ = (D_next + D_) / 2;
	Matrix3D gD = D_half_;
	gD[one][one] = g(gD[one][one]);
	gD[two][two] = g(gD[two][two]);
	gD[three][three] = g(gD[three][three]);

	const Vector3D dW_half_ = delta_W_half(n, j);
	const Vector3D d_ksi_half_ = delta_ksi_half(n, j);

	return 	(D_half_ * D_half_).abs() * _alpha * d_ksi_half_ + gD * (dW_half_ - d_ksi_half_);
}


//	Матрица Якоби (A = dF / dU) [A = L*D*R]
Matrix3D GasDynamicsEquation::A(size_t n, size_t j) const noexcept
{
	const double u = _u.getElement(n, j);

	double H_ = H(n, j);

	return Matrix3D(
		Vector3D{ 0,									1,							0			},
		Vector3D{ (gamma - 3) * u * u / 2,				(3 - gamma) * u,			gamma - 1	},
		Vector3D{ u * ((gamma - 1) * u * u / 2 - H_),	H_ - (gamma - 1) * u * u,	gamma * u	}
	);
}

//	A = L*D*R
Matrix3D GasDynamicsEquation::R(size_t n, size_t j) const noexcept
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	const double c = pow(gamma * p / ro, 0.5);
	const double part = 1 / (2 * gamma * p / ro);

	const double H_ = H(n, j);

	return Matrix3D( 
		Vector3D{ part,					1,			part				},
		Vector3D{ part * (u - c),		u,			part * (u + c)		},
		Vector3D{ part * (H_ - u * c),	u * u / 2,	part * (H_ + u * c)	}
	);
}

//	A = L*D*R
Matrix3D GasDynamicsEquation::L(size_t n, size_t j) const noexcept
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	const double c = pow(gamma * p / ro, 0.5);

	const double k = (gamma - 1);
	const double part_1 = k * u;
	const double part_2 = part_1 * u / 2;

	return Matrix3D(
		Vector3D{ u * c + part_2,	-c - part_1,	k	},
		Vector3D{ c * c - part_2,	part_1,			-k	} / (c * c),
		Vector3D{ -u * c + part_2,	c - part_1,		k	}
	);
}

//	A = L*D*R
Matrix3D GasDynamicsEquation::D(size_t n, size_t j) const noexcept
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	const double c = pow(gamma * p / ro, 0.5);

	return Matrix3D(
		Vector3D{ u - c,	0,		0		},
		Vector3D{ 0,		u,		0		},
		Vector3D{ 0,		0,		u + c	}
	);
}
