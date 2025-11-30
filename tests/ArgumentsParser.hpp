#pragma once

#include <cstdio>
#include <cstring>

#include <string>
#include <set>

struct ArgumentParser {
	int argc;
	char **argv;
	
	ArgumentParser(int argc, char **argv) {
		this->argc = argc;
		this->argv = argv;
	}
	
	bool Bool(std::string name)
	{
		if (args.contains(name) == false) {
			args.insert(name);
		}
		for (int i=argc-1; i>0; --i) {
			if (std::string(argv[i]) == ("-" + name)) {
				return true;
			}
			if (std::string(argv[i]).starts_with("-" + name + "=")) {
				std::string n = argv[i] + strlen(name.c_str()) + 2;
				if (n == "0" || n == "false" || n == "off") {
					return false;
				}
				if (n == "1" || n == "true" || n == "on") {
					return true;
				}
			}
		}
		return false;
	}
	
	std::string String(std::string name)
	{
		if (args.contains(name) == false) {
			args.insert(name);
		}
		for (int i=argc-1; i>0; --i) {
			if (std::string(argv[i]).starts_with("-" + name + "=")) {
				return argv[i] + strlen(name.c_str()) + 2;
			}
		}
		return "";
	}
	
	int Int(std::string name)
	{
		return atoi(String(name).c_str());
	}
	
	int Int(std::string name, int min, int max)
	{
		int v = Int(name);
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	int Int(std::string name, int min, int max, int defaultValue)
	{
		std::string s = String(name);
		int v = atoi(s.c_str());
		if (s == "") {
			v = defaultValue;
		}
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	int Float(std::string name)
	{
		return atof(String(name).c_str());
	}
	
	int Float(std::string name, float min, float max)
	{
		float v = Int(name);
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	float Float(std::string name, float min, float max, float defaultValue)
	{
		std::string s = String(name);
		float v = atof(s.c_str());
		if (s == "") {
			v = defaultValue;
		}
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	void PrintHelp() { PrintHelp(false); }
	
	void PrintHelp(bool force)
	{
		if (Bool("help") == false && force == false) {
			return;
		}
		printf("Available program arguments:\n");
		for (auto c : args) {
			printf("   -%s\n", c.c_str());
		}
		printf("\n");
		fflush(stdout);
	}
	
private:
	std::set<std::string> args;
};
