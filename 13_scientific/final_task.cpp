#include <iostream>
#include <vector>
#include <cmath>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

// Constants
const int nx = 41;
const int ny = 41;
const int nt = 500;
const int nit = 50;
const double dx = 2.0 / (nx - 1);
const double dy = 2.0 / (ny - 1);
const double dt = 0.01;
const double rho = 1.0;
const double nu = 0.02;

// Function to create a linearly spaced vector
vector<double> linspace(double start, double end, int num) {
    vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

int main() {
    // Variables
    vector<vector<double>> u(ny, vector<double>(nx, 0.0));
    vector<vector<double>> v(ny, vector<double>(nx, 0.0));
    vector<vector<double>> p(ny, vector<double>(nx, 0.0));
    vector<vector<double>> b(ny, vector<double>(nx, 0.0));
    
    vector<double> x = linspace(0.0, 2.0, nx);
    vector<double> y = linspace(0.0, 2.0, ny);

    // Simulation loop
    for (int n = 0; n < nt; ++n) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                b[j][i] = rho * (1.0 / dt * ((u[j][i + 1] - u[j][i - 1]) / (2 * dx) + (v[j + 1][i] - v[j - 1][i]) / (2 * dy)) -
                               pow((u[j][i + 1] - u[j][i - 1]) / (2 * dx), 2) -
                               2 * ((u[j + 1][i] - u[j - 1][i]) / (2 * dy) * (v[j][i + 1] - v[j][i - 1]) / (2 * dx)) -
                               pow((v[j + 1][i] - v[j - 1][i]) / (2 * dy), 2));
            }
        }

        for (int it = 0; it < nit; ++it) {
            vector<vector<double>> pn = p;
            for (int j = 1; j < ny - 1; ++j) {
                for (int i = 1; i < nx - 1; ++i) {
                    p[j][i] = (dy * dy * (pn[j][i + 1] + pn[j][i - 1]) +
                              dx * dx * (pn[j + 1][i] + pn[j - 1][i]) -
                              b[j][i] * dx * dx * dy * dy) /
                              (2 * (dx * dx + dy * dy));
                }
            }
            for (int i = 0; i < ny; ++i) p[i][nx - 1] = p[i][nx - 2];
            for (int i = 0; i < nx; ++i) p[0][i] = p[1][i];
            for (int i = 0; i < ny; ++i) p[i][0] = p[i][1];
            for (int i = 0; i < nx; ++i) p[ny - 1][i] = 0.0;
        }

        vector<vector<double>> un = u;
        vector<vector<double>> vn = v;

        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1]) -
                          vn[j][i] * dt / dy * (un[j][i] - un[j - 1][i]) -
                          dt / (2 * rho * dx) * (p[j][i + 1] - p[j][i - 1]) +
                          nu * dt / (dx * dx) * (un[j][i + 1] - 2 * un[j][i] + un[j][i - 1]) +
                          nu * dt / (dy * dy) * (un[j + 1][i] - 2 * un[j][i] + un[j - 1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1]) -
                          vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i]) -
                          dt / (2 * rho * dx) * (p[j + 1][i] - p[j - 1][i]) +
                          nu * dt / (dx * dx) * (vn[j][i + 1] - 2 * vn[j][i] + vn[j][i - 1]) +
                          nu * dt / (dy * dy) * (vn[j + 1][i] - 2 * vn[j][i] + vn[j - 1][i]);
            }
        }

        for (int i = 0; i < nx; ++i) u[0][i] = 0.0;
        for (int i = 0; i < ny; ++i) u[i][0] = 0.0;
        for (int i = 0; i < ny; ++i) u[i][nx - 1] = 0.0;
        for (int i = 0; i < nx; ++i) u[ny - 1][i] = 1.0;

        for (int i = 0; i < nx; ++i) v[0][i] = 0.0;
        for (int i = 0; i < nx; ++i) v[ny - 1][i] = 0.0;
        for (int i = 0; i < ny; ++i) v[i][0] = 0.0;
        for (int i = 0; i < ny; ++i) v[i][nx - 1] = 0.0;

        // Visualization using OpenCV
        Mat pressureMap(ny, nx, CV_64F, &p[0][0]);
        Mat uMap(ny, nx, CV_64F, &u[0][0]);
        Mat vMap(ny, nx, CV_64F, &v[0][0]);

        Mat colorMap;
        applyColorMap(pressureMap, colorMap, COLORMAP_JET);

        imshow("Pressure", colorMap);

        // Create a quiver plot
        Mat quiverMap(ny, nx, CV_8UC3, Scalar(255, 255, 255));
        for (int j = 0; j < ny; j += 2) {
            for (int i = 0; i < nx; i += 2) {
                Point pt1(i, j);
                Point pt2(i + 5 * u[j][i], j + 5 * v[j][i]);
                arrowedLine(quiverMap, pt1, pt2, Scalar(0, 0, 0), 1, 8, 0, 0.2);
            }
        }

        imshow("Velocity", quiverMap);
        waitKey(10);
    }

    waitKey(0);
    return 0;
}
