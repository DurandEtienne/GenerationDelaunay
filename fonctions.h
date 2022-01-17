#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include<vector>

std::vector<double> rotate(std::vector<double> v, double theta);

std::vector<double> rotate (std::vector<double> V , double theta);

double scalar (std::vector<double> v1 , std::vector<double> v2);

double norm(std::vector<double> v);

double distance(std::vector<double> v, std::vector<double> u);

int getIndexOfElementInArray(std::vector<double> v, double element);

void fixEdgesindexing(std::vector<double> &triangles, double ind);

int checkIfAnEdgeAlreadyExist(std::vector<double> edges, std::vector<double> edge);

int getTriangle(std::vector<double> point,std::vector<double> triangles ,std::vector<double> edges ,std::vector<double> vertices);

std::vector<int> getTriangleNeighbors(std::vector<double> triangleIndex,std::vector<double> triangles);

std::vector<int> getTriangleCavity(std::vector<double> point, std::vector<double> triangles, std::vector<double> edges, std::vector<double> vertices);

bool inCircumscribedCircle(std::vector<double> point, int triangleIndex, std::vector<double> triangles, std::vector<double> edges, std::vector<double> vertices);

bool haveOneCommonEdge(std::vector<double> triangle1, std::vector<double> triangle2);

void deleteEdgesOnCavityAndReconnect(std::vector<double> point, std::vector<double> &triangles, std::vector<double> &edges, std::vector<double> &vertices);
