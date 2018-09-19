#ifndef PARAMETER_H_INCLUDED
#define PARAMETER_H_INCLUDED

#include <tgmath.h>

namespace p {
    // Universal Constants
    const float pi = 3.141592653589793;
    const float twoPi = 2 * pi;
    const float sqrtTwo = sqrt(2);
    const float radToDeg = (float)180 / pi;

    // Coefficients
    const int maxColloids = 5000;
    const int latRows = 100; // latRows and winHeight decide tileSize
    const int latCols = latRows;

    static float timeStep = 0.01;
    static float diffCo = 0.001;
    static float decayRate = 0.1;
    static float colVelocity = 1;
    static float production = 0.02;
    static float backProd = 0.15;
    static float frontDest = 0;
    static float colShift = 1;
    static float chemCoup = -0.15;
    static float rotDiff = 0.4;

    static float shiftRand = 0.2;

    static bool walls = 1;
    static bool wallsSink = 0;
    static float wallActDist = 0.5;
    static float wallForce = 10;

    // System Space
    const float boxHeight = 10;
    const float boxWidth = (float)boxHeight * latCols / latRows; // Do not change

    static int maxBinRows = 80;
    static int maxBinCols = (int)round(maxBinRows * latCols / latRows);
    static int binRows = 80;
    static int binCols = (int)round(binRows * latCols / latRows);

    // Clustering
    static int minClusSize = 20;
    static float minClusDeny = 120;

    // Timing
    static int simSkip = 25;
    const Time sleepTime = milliseconds(0);
    const Time buttonDelay = milliseconds(200);

    float rand01(); // Random number between [0,1)
    float randpm1(); // Random number between (-1,1)

    // Cosmetic
    const float winScale = 1;
    const int winWidth = (int)round(1200 * winScale);
    const int winHeight = (int)round(400 * winScale);
    static float shadeScale = 2;
    const float colImageSize = 60; // Size of texture [px]
    static float colDiameter = 12 * winScale; // [px]
    const float clusDefaultRad = 15;

    // Dependent
    const float tileSize = boxHeight / latRows;
    const int halfGridSize = (int)floor((float)latRows / 2);
    static float colScale = colDiameter / colImageSize;
    const float sqrtTimeStep = sqrt(timeStep);

    const float boxToWin = (float)winHeight / boxHeight;
}

#endif // PARAMETER_H_INCLUDED
