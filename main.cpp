﻿#include "test.h"
#include <chrono>
#include <iostream>

void main()
{
	auto begin = std::chrono::steady_clock::now();

	modeling1();

	//modeling2();

	//modeling3();

	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "\nLead time: " << elapsed_ms.count() << "ms" << std::endl;

	return;
}