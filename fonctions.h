#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include<vector>

using namespace std;

vector<double> rotate (vector<double> V , double theta);

double scalar (vector<double> v1 , vector<double> v2);

int getTriangle(vector<double> point,vector<double> triangles ,vector<double> edges ,vector<double> vertices);

vector<int> getTriangleNeighbors(vector<double> triangleIndex,vector<double> triangles);

vector<int> getTriangleCavity(vector<double> point, vector<double> triangles, vector<double> edges, vector<double> vertices);
