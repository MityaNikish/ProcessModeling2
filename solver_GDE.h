﻿#pragma once
#include <functional>
#include <filesystem>
#include "types.h"

//	Cхема Лакс-Вендроффа c TVD модификацией для моделирования уравнение газовой динамики
class GasDynamicsEquation final
{
private:
	//	Данные пространственной сетки
	const ExpanseGrid _expanse_grid;
	//	Данные времянной сетки
	const TimeGrid _time_grid;

	//	Начальные условия
	const StartCondition _start_condition;
	//	Граничные условия
	const BorderlineCondition _borderline_condition;

	//	Пространственно-временая плоскость плотности
	Matrix _ro;
	//	Пространственно-временая плоскость скорости
	Matrix _u;
	//	Пространственно-временая плоскость давления
	Matrix _p;

	//	Площадь поперечного сечения трубы
	const std::vector<double> _S;

	//	Отношение шага по времяни к шагу по пространству
	const double _alpha;

public:
	//	Показатель адиабаты
	static double gamma;
	//	Искусственная вязкость Лапидуса
	static double artificial_viscosity;

public:
	GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::vector<double>& S);
	GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::function<double(double)>& S_func);
	~GasDynamicsEquation() = default;

	GasDynamicsEquation(const GasDynamicsEquation&) = delete;
	GasDynamicsEquation(GasDynamicsEquation&&) = delete;

	GasDynamicsEquation& operator=(const GasDynamicsEquation&) = delete;
	GasDynamicsEquation& operator=(GasDynamicsEquation&&) = delete;

	//	Моделирует УГД
	void solving() noexcept;

	//	Записывает массив данных плотности в файл
	void writeRO(const std::filesystem::path& file_path) const;
	//	Записывает массив данных скорости в файл
	void writeU(const std::filesystem::path& file_path) const;
	//	Записывает массив данных давления в файл
	void writeP(const std::filesystem::path& file_path) const;

	//	Вернуть массив данных плотности
	std::vector<double> getRO() const;
	//	Вернуть массив данных скорости
	std::vector<double> getU() const;
	//	Вернуть массив данных давления
	std::vector<double> getP() const;

private:
	void initConditions() noexcept;
	void postProcessing() noexcept;

	bool chekStopConditions(size_t n, size_t j, double eps) const noexcept;

	double E(size_t n, size_t j) const noexcept;
	double H(size_t n, size_t j) const noexcept;
	double g(double z) const noexcept;
	
	Vector3D U(size_t n, size_t j) const noexcept;
	Vector3D U_future(size_t n, size_t j)  const noexcept;

	Vector3D F(size_t n, size_t j) const noexcept;
	Vector3D Q(size_t n, size_t j) const noexcept;

	Vector3D delta_ksi_half(size_t n, size_t j) const noexcept;
	Vector3D delta_W_half(size_t n, size_t j) const noexcept;
	Vector3D G_half(size_t n, size_t j) const noexcept;

	Matrix3D A(size_t n, size_t j) const noexcept;
	Matrix3D R(size_t n, size_t j) const noexcept;
	Matrix3D L(size_t n, size_t j) const noexcept;
	Matrix3D D(size_t n, size_t j) const noexcept;
};