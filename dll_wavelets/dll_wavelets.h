#include <vector>
#pragma once
//extern "C" __declspec(dllexport) int sum(int, int);
extern "C" __declspec(dllexport) void morlet_wavelet(int, int, int, double*, double*, double*);
