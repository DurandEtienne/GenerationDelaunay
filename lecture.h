#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include<vector>

std::vector<std::vector<double>> lecture(std::string nomdufichier);
void write_mesh (std::vector<double> triangles , std::vector<double> edges ,std::vector<double> vertices, std::string nomdufichier);
void write_sol (std::vector<double> qualite, std::string nomdufichier);
