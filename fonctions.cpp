#include "fonctions.h"
#include <math.h>
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
                return true;
            }
        }
    }
    return false;
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

int getTriangle(vector<double> point, vector<double> triangles, vector<double> edges, vector<double> vertices)
{
    for (int i = 0; i < triangles.size(); i += 4)
    {
        vector<double> v1(2), v2(2), v3(2);
        // first vertice
        v1[0] = vertices[(triangles[i]-1) * 3];
        v1[1] = vertices[(triangles[i]-1) * 3 + 1];
        // second vertice
        v2[0] = vertices[(triangles[i+1]-1) * 3];
        v2[1] = vertices[(triangles[i+1]-1) * 3 + 1];
        // third vertice
        v3[0] = vertices[(triangles[i+2]-1) * 3];
        v3[1] = vertices[(triangles[i+2]-1) * 3 + 1];
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

            return i / 4+1;
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
    int i = triangleIndex * 4;
    vector<double> v1(2), v2(2), v3(2);
    // first vertice
    v1[0] = vertices[edges[triangles[i] * 3] * 3];
    v1[1] = vertices[edges[triangles[i] * 3] * 3 + 1];
    // second vertice
    v2[0] = vertices[edges[triangles[i] * 3 + 1] * 3];
    v2[1] = vertices[edges[triangles[i] * 3 + 1] * 3 + 1];
    // third vertice
    v3[0] = vertices[edges[triangles[i + 1] * 3] * 3];
    v3[1] = vertices[edges[triangles[i + 1] * 3] * 3 + 1];
    if ((v3[0] == v1[0] && v3[1] == v1[1]) || (v3[0] == v2[0] && v3[1] == v2[1]))
    {
        v3[0] = vertices[edges[triangles[i + 1] * 3 + 1] * 3];
        v3[1] = vertices[edges[triangles[i + 1] * 3 + 1] * 3 + 1];
    }
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
    CM[0] = v1M[0] - v1C[0];
    // verify condition
    return (norm(CM) <= norm(v1C));
}
