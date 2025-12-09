#include <cstdio>
#include <cstdarg>

#include <chrono>

#include "../include/worldgen/Debug.hpp"

namespace DebugNamespace
{
using TimePoint = std::chrono::steady_clock::time_point;

static TimePoint Now()
{
	return std::chrono::steady_clock::now();
}

thread_local TimePoint lastTime = Now();


void RestartDebugTime()
{
	lastTime = Now();
}

void DebugTime(const char *file, int line, const char *fmt, ...)
{
	TimePoint now = Now();
	auto ns = std::chrono::nanoseconds(now - lastTime).count();
	double ms = ns/1'000'000.0;
	
	printf(" DEBUG(%7.3f ms)   ", ms);
	
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	
	printf("     %s : %i \n", file, line);
	fflush(stdout);
	
	lastTime = Now();
}
}
