#pragma once

#include <cstdio>
#include <cstring>

#include <string>
#include <map>

class ArgumentParser {
public:
	int argc;
	char **argv;
	
	ArgumentParser(int argc, char **argv) : argc(argc), argv(argv), args({}) {
	}
	
	bool Bool(std::string name)
	{
		if (args.contains(name) == false) {
			args[name] = "(bool)";
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
	
	std::string String(std::string name, std::string defaultValue)
	{
		if (args.contains(name) == false) {
			args[name] = "(string)";
		}
		for (int i=argc-1; i>0; --i) {
			if (std::string(argv[i]).starts_with("-" + name + "=")) {
				return argv[i] + strlen(name.c_str()) + 2;
			}
		}
		return defaultValue;
	}
	
	std::string String(std::string name)
	{
		return String(name, "");
	}
	
	int Int(std::string name)
	{
		std::string s = String(name);
		int v = atoi(s.c_str());
		if (s == "") {
			v = 0.0f;
		}
		args[name] = "(floatwint)";
		return v;
	}
	
	int Int(std::string name, int defaultValue)
	{
		std::string s = String(name);
		int v = atoi(s.c_str());
		if (s == "") {
			v = defaultValue;
		}
		args[name] = "(int) ~(" + std::to_string(defaultValue) + ")";
		return v;
	}
	
	int Int(std::string name, int min, int max)
	{
		int v = Int(name);
		args[name] = "(int) <"+std::to_string(min)+";"+std::to_string(max)+">";
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	int Int(std::string name, int min, int max, int defaultValue)
	{
		int v = Int(name, defaultValue);
		args[name] = "(int) <"+std::to_string(min)+";"+std::to_string(max)+"> ~(" + std::to_string(defaultValue) + ")";
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	int Float(std::string name)
	{
		std::string s = String(name);
		float v = atof(s.c_str());
		if (s == "") {
			v = 0.0f;
		}
		args[name] = "(float)";
		return v;
	}
	
	float Float(std::string name, float defaultValue)
	{
		std::string s = String(name);
		float v = atof(s.c_str());
		if (s == "") {
			v = defaultValue;
		}
		return v;
	}
	
	int Float(std::string name, float min, float max)
	{
		float v = Float(name);
		args[name] = "(float) <"+std::to_string(min)+";"+std::to_string(max)+">";
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	float Float(std::string name, float min, float max, float defaultValue)
	{
		float v = Float(name, defaultValue);
		args[name] = "(float) <"+std::to_string(min)+";"+std::to_string(max)+"> ~(" + std::to_string(defaultValue) + ")";
		v = v > min ? v : min;
		v = v < max ? v : max;
		return v;
	}
	
	void PrintHelp()
	{
		printf("Available program arguments:\n");
		for (auto c : args) {
			printf("   -%s : %s\n", c.first.c_str(), c.second.c_str());
		}
		printf("\n");
		fflush(stdout);
	}
	
private:
	std::map<std::string, std::string> args;
};
