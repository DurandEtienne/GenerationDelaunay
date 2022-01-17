#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>

#include "fonctions.h"
#include <bits/stdc++.h>

using namespace std;

vector<double> rotate(vector<double> v, double theta)
{
    vector<double> res(2);
    // convert to radian
    theta *= M_PI / 180.;
    //rotate
    res[0] = cos(theta) * v[0] - sin(theta) * v[1];
    res[1] = sin(theta) * v[0] + cos(theta) * v[1];
    return res;
}
bool haveOneCommonEdge(vector<double> triangle1, vector<double> triangle2)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (triangle1[i] == triangle2[j])
            {
                return triangle1[i];
            }
        }
    }
    return NULL;
}

double scalar(vector<double> v1, vector<double> v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1];
}

double norm(vector<double> v)
{
    double s = 0;

    for (int i = 0; i < v.size(); i++)
    {
        s += v[i] * v[i];
    }

    return sqrt(s);
}
double distance(vector<double> v, vector<double> u)
{
    return sqrt((v[0] - u[0]) * (v[0] - u[0]) + (v[1] - u[1]) * (v[1] - u[1]));
}
int getIndexOfElementInArray(vector<double> v, double element)
{
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] == element)
        {
            return i;
        }
    }
    return NULL;
}
void fixEdgesindexing(vector<double> &triangles, double ind)
{
    for (int i = 0; i < triangles.size(); i++)
    {

        if (triangles[i] % 4 != 3 && triangles[i] > id)
        {
            triangles[i] -= 1;
        }
    }
}
int checkIfAnEdgeAlreadyExist(vector<double> edges, vector<double> edge)
{
    for (int i = 0; i < edges.size(); i += 3)
    {
        if (edges[i] == edge[0] && edges[i + 1] == edge[1])
        {
            return (i / 3) + 1;
        }
    }
    return NULL;
}
int getTriangle(vector<double> point, vector<double> triangles, vector<double> edges, vector<double> vertices)
{
    for (int i = 0; i < triangles.size(); i += 4)
    {
        vector<double> v1(2), v2(2), v3(2);
        // first vertice
        v1[0] = vertices[(triangles[i] - 1) * 3];
        v1[1] = vertices[(triangles[i] - 1) * 3 + 1];
        // second vertice
        v2[0] = vertices[(triangles[i + 1] - 1) * 3];
        v2[1] = vertices[(triangles[i + 1] - 1) * 3 + 1];
        // third vertice
        v3[0] = vertices[(triangles[i + 2] - 1) * 3];
        v3[1] = vertices[(triangles[i + 2] - 1) * 3 + 1];
        // if ((v3[0] == v1[0] && v3[1] == v1[1]) || (v3[0] == v2[0] && v3[1] == v2[1]))
        // {
        //     v3[0] = vertices[edges[triangles[i + 1] * 3 + 1] * 3];
        //     v3[1] = vertices[edges[triangles[i + 1] * 3 + 1] * 3 + 1];
        // }
        vector<double> v1M(2), v1v2(2), v1v3(2);
        v1M[0] = point[0] - v1[0];
        v1M[1] = point[1] - v1[1];
        v1v2[0] = v2[0] - v1[0];
        v1v2[1] = v2[1] - v1[1];
        v1v3[0] = v3[0] - v1[0];
        v1v3[1] = v3[1] - v1[1];
        double lambda2 = scalar(v1M, rotate(v1v3, 90.)) / scalar(v1v2, rotate(v1v3, 90.));
        double lambda3 = scalar(v1M, rotate(v1v2, 90.)) / scalar(v1v3, rotate(v1v2, 90.));
        // checking condition
        if (lambda2 >= 0. && lambda3 >= 0. && lambda2 + lambda3 <= 1.)
        {
            if (lambda2 == 0. || lambda3 == 0. || lambda2 + lambda3 == 1.)
            {
                cout << "You choosed a point in an edge !!! " << endl;
            }

            return i / 4 + 1;
        }
    }
}
vector<int> getTriangleNeighbors(int triangleIndex, vector<double> triangles)
{
    vector<int> res;
    for (int i = 0; i < triangles.size(); i += 4)
    {
        if (haveOneCommonEdge(vector<double>(triangles.begin() + triangleIndex, triangles.begin() + triangleIndex + 2), vector<double>(triangles.begin() + i, triangles.begin() + i + 2)))
            ; //[triangleIndex:(triangleIndex + 3)], triangles [i:i + 3]))
        {
            res.push_back(i / 4);
        }
    }
    return res;
}
bool inCircumscribedCircle(vector<double> point, int triangleIndex, vector<double> triangles, vector<double> edges, vector<double> vertices)
{
    // Retrieving triangle vertices
    int i = (triangleIndex - 1) * 4;
    vector<double> v1(2), v2(2), v3(2);
    // first vertice
    v1[0] = vertices[(triangles[i] - 1) * 3];
    v1[1] = vertices[(triangles[i] - 1) * 3 + 1];
    // second vertice
    v2[0] = vertices[(triangles[i + 1] - 1) * 3];
    v2[1] = vertices[(triangles[i + 1] - 1) * 3 + 1];
    // third vertice
    v3[0] = vertices[(triangles[i + 2] - 1) * 3];
    v3[1] = vertices[(triangles[i + 2] - 1) * 3 + 1];
    // if ((v3[0] == v1[0] && v3[1] == v1[1]) || (v3[0] == v2[0] && v3[1] == v2[1]))
    // {
    //     v3[0] = vertices[edges[triangles[i + 1] * 3 + 1] * 3];
    //     v3[1] = vertices[edges[triangles[i + 1] * 3 + 1] * 3 + 1];
    // }
    // Creating edges vectors
    vector<double> v1M(2), v1v2(2), v1v3(2);
    v1M[0] = point[0] - v1[0];
    v1M[1] = point[1] - v1[1];
    v1v2[0] = v2[0] - v1[0];
    v1v2[1] = v2[1] - v1[1];
    v1v3[0] = v3[0] - v1[0];
    v1v3[1] = v3[1] - v1[1];

    // Finding weight coefficients of our point on our triangle

    double alpha2 = norm(v1v2) * norm(v1v2);
    double alpha3 = norm(v1v3) * norm(v1v3);
    double alpha23 = scalar(v1v2, v1v3);
    double lambda2 = (alpha3 * (alpha2 - alpha23)) / (2 * (alpha2 * alpha3 - alpha23 * alpha23));
    double lambda3 = (alpha2 * (alpha3 - alpha23)) / (2 * (alpha2 * alpha3 - alpha23 * alpha23));
    // Center vector position
    vector<double> v1C(2);
    v1C[0] = lambda2 * v1v2[0] + lambda3 * v1v3[0];
    v1C[1] = lambda2 * v1v2[1] + lambda3 * v1v3[1];
    // Now the position of M relatively to v1 is known let's use it
    // positon relatively to each other
    vector<double> CM(2);
    CM[0] = v1M[0] - v1C[0];
    CM[1] = v1M[1] - v1C[1];
    // verify condition
    return (norm(CM) <= norm(v1C));
}

vector<int> getTriangleCavity(vector<double> point, vector<double> triangles, vector<double> edges, vector<double> vertices)
{
    int nb_triangles = triangles.size() / 4;
    bool cavity;
    vector<int> cavityIndex;
    for (int i = 1; i <= nb_triangles; i++)
    {
        cavity = inCircumscribedCircle(point, i, triangles, edges, vertices);
        if (cavity == true)
        {
            cavityIndex.push_back(i);
        }
    }
    return cavityIndex;
}
void deleteEdgesOnCavityAndReconnect(vector<double> point, vector<double> &triangles, vector<double> &edges, vector<double> &vertices)
{
    vector<int> trianglesInCavity = getTriangleCavity(point, triangles, edges, vertices);
    vector<double> triangle1(3), triangle2(3);
    vector<bool> bTriangles(triangles.size() / 4, true), bEdges(edges.size() / 3, true), bVertices(vertices.size() / 3, true);
    vector<double> pointsToConnect, boundryEdgesToConnect, edgesToDelete;
    int n = trianglesInCavity.size();
    int nbOfVertices = vertices.size() / 3;
    int nbOfEdges = edges.size() / 3;
    int nbOfTriangles = triangles.size() / 4;
    int ibegin, jbegin;
    // deleting elements from cavity and getting vertices of cavity boundary
    for (int i = 0; i < n; i++)
    {
        // constructing first triangle and add edges extremes which are actually the cavity boundry vertices
        ibegin = (int)(4 * (trianglesInCavity[i] - 1));
        for (int r = 0; r < 3; r++)
        {
            triangle1[r] = triangles[ibegin + r];
            // adding boundry edges
            if (!getIndexOfElementInArray(boundryEdgesToConnect, triangle1[r]))
            {
                boundryEdgesToConnect.push_back(triangle1[r]);
            }
            // getting vertices of cavity boundary
            int ibeg = (int)(3 * triangle1[r]);
            // add first boundry vertice if it does not exist already
            if (!getIndexOfElementInArray(pointsToConnect, edges[ibeg]))
            {
                pointsToConnect.push_back(edges[ibeg]);
            }
            // add second boundry vertice if it does not exist already
            if (!getIndexOfElementInArray(pointsToConnect, edges[ibeg + 1]))
            {
                pointsToConnect.push_back(edges[ibeg + 1]);
            }
        }
        // second loop for upcomming triangles to make combinations and get common edges
        for (int j = i + 1; j < n; j++)
        {
            // constructiion of second triangle to check if they have a common edge
            jbegin = (int)(4 * (trianglesInCavity[j] - 1));
            for (int r = 0; r < 3; r++)
            {
                triangle2[r] = triangles[jbegin + r];
            }
            // checking if they have a common edge and erase it in case
            if (haveOneCommonEdge(triangle1, triangle2))
            {
                // getting the beginnig position of the common edge
                double communBndIndex = haveOneCommonEdge(triangle1, triangle2);
                edgesToDelete.push_back(communBndIndex);
                // deleting commun edges from the vector of cavity boundry edges
                if (getIndexOfElementInArray(boundryEdgesToConnect, communBndIndex))
                {
                    boundryEdgesToConnect.erase(boundryEdgesToConnect.begin() + getIndexOfElementInArray(boundryEdgesToConnect, communBndIndex));
                }
            }
        }
    }
    // erasing triangles after using them
    for (int i = 0; i < n; i++)
    {
        ibegin = 4 * (trianglesInCavity[i] - 1);
        triangles.erase(triangles.begin() + ibegin, triangles.begin() + ibegin + 4);
        nbOfTriangles -= 1;
    }
    // adding the new point
    vertices.push_back(point[0]);
    vertices.push_back(point[1]);
    vertices.push_back(0);
    nbOfVertices++;
    // making new edges and triangles
    int edg1, edg2, edg3;
    for (int q = 1; q < boundryEdgesToConnect.size() + 1; q++)
    {
        edg1 = (int)boundryEdgesToConnect[q];
        if (!checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3], nbOfVertices}))
        {
            edges.push_back(edges[(edg1 - 1) * 3]);
            edges.push_back(nbOfVertices);
            edges.push_back(1);
            nbOfEdges++;
            edg2 = nbOfEdges;
        }
        else
        {
            edg2=checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3], nbOfVertices})
        }
        if (!checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3+1], nbOfVertices}))
        {
            edges.push_back(edges[(edg1 - 1) * 3+1]);
            edges.push_back(nbOfVertices);
            edges.push_back(1);
            nbOfEdges++;
            edg3 = nbOfEdges;
        }
        else
        {
            edg3=checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3], nbOfVertices})
        }
        triangles.push_back(edg1);
        triangles.push_back(edg2);
        triangles.push_back(edg3);
        triangles.push_back(2);
        nbOfTriangles++;
    }
    // finally erase edges
    for (int p = 0; p < edgesToDelete.size(); p++)
    {
        ibegin = (int)(edgesToDelete[p] - 1) * 3;
        edges.erase(edges.begin() + ibegin, edges.begin() + ibegin + 3);
        fixEdgesindexing(triangles, edgesToDelete[p]);
    }
}




void buildBoiteEnglobante(vector<double> &edges, vector<double> &vertices)
{
  // Détermination des minimums et maximums en x et y
  int nbOfVertices = vertices.size() / 3;
  int nbOfEdges = edges.size() / 3;
  for (int p = 0; nbOfVertices; p++)
  {
    
  }
}
