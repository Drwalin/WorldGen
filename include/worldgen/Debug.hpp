#pragma once

namespace DebugNamespace
{
void RestartDebugTime();
void DebugTime(const char *file, int line, const char *fmt, ...);
} // namespace DebugNamespace

#define RESTART_DEBUG_TIME() DebugNamespace::RestartDebugTime()
#define DEBUG_TIME(...)                                                        \
	DebugNamespace::DebugTime(__FILE__, __LINE__, __VA_ARGS__)
