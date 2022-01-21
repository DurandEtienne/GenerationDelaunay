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

bool edgeAlreadyHere(std::vector<std::vector<double>> edgs,std::vector<double>edg);

bool alreadyConnected(std::vector<double>Points,double id);

void eraseEdge(std::vector<double> &edgs,std::vector<double> edg);

bool inCircumscribedCircle(std::vector<double> point, int triangleIndex, std::vector<double> triangles, std::vector<double> edges, std::vector<double> vertices);

double haveOneCommonEdge(std::vector<double> triangle1, std::vector<double> triangle2);

void deleteEdgesOnCavityAndReconnect(std::vector<double> point, std::vector<double> &triangles, std::vector<double> &edges, std::vector<double> &vertices);

void buildBoiteEnglobante(std::vector<double> &edges, std::vector<double> &vertices);

void InitializeMeshBoite(std::vector<double> &triangles, std::vector<double> &edges, std::vector<double> &vertices);

void AdvanceMeshBoite(int i, std::vector<double> &triangles, std::vector<double> &edges, std::vector<double> &vertices);

void deleteBoiteEnglobante(std::vector<double> &triangles, std::vector<double> &edges, std::vector<double> &vertices, std::vector<double> vertices_initial);

int verifyTriangleInGeometry(int k,std::vector<double> triangles,std::vector<double> vertices_initial);
