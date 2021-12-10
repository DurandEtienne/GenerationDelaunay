#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include<vector>


std::vector<double> rotate (std::vector<double> V , double theta);

double scalar (std::vector<double> v1 , std::vector<double> v2);

int getTriangle(std::vector<double> point,std::vector<double> triangles ,std::vector<double> edges ,std::vector<double> vertices);

std::vector<int> getTriangleNeighbors(std::vector<double> triangleIndex,std::vector<double> triangles);

std::vector<int> getTriangleCavity(std::vector<double> point, std::vector<double> triangles, std::vector<double> edges, std::vector<double> vertices);

bool inCircumscribedCircle(std::vector<double> point, int triangleIndex, std::vector<double> triangles, std::vector<double> edges, std::vector<double> vertices);

bool haveOneCommonEdge(std::vector<double> triangle1, std::vector<double> triangle2);
