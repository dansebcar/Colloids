#include "colloid.h"
#include "cluster.h"
#include "parameter.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <SFML/Graphics.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace sf;

int fix(int x, int N); // Fixes an out-of-bounds lattice site
float generateGaussianNoise(const float& mean, const float &stdDev);
int wallSiteCheck(int y);
float getParam(std::string str);
void setParam(std::string str, float val);
float clusAve(float x, int nvc);
int index(int x, int y);
int bindex(int x, int y);
int chromifySigned(float x, int maxVal = 255, int minVal = 0);
int chromifyUnsigned(float x, int maxVal = 255, int minVal = 0);
float dist1d(float x1, float x2, float box);
float dist2d(float x1, float x2, float y1, float y2);

int main()
{
    /* Data gathering
    **********************************************************************/
    bool dataGathering = 0;

    bool routine = 0;
    int routineNo = 0;
    int maxRoutines = 20;
    std::string routineDescription = "";

    // Output type (choose from: aveSep, maxSep, theta, clus, clusHist, collHist, clusPos)
    std::string outputType = "clus";

    int numParams = 9;
    std::vector<std::string> params(numParams);
    params = { "diffCo", "decayRate", "colVelocity", "production", "backProd", "frontDest", "colShift", "chemCoup", "rotDiff" };
    // Don't forget to check the initial value
    std::string indepVbl = "colVelocity";

    float indepVblEnd = 1.2;
    float indepVblStart = getParam(indepVbl);

    int warmUpSteps = 7.5 * 1000;
    int stepsPerDump = 25;
    int numIncs = 15;
    int dumpsPerRun = 60;
    int runsPerInc = 30;

    if (outputType == "clusHist") { numIncs = 0; }
    if (outputType == "collHist") { numIncs = 0; }
    if (outputType == "clusPos") { numIncs = 0; runsPerInc = 1; stepsPerDump = 10; }

    int outHistBins = 200;
    int clusHistMaxColls = 200;

    // Chemical initialisation (choose from: Gauss, unif, rand, wallTest, column, row, tlc, cbar, vgrad) otherwise null
    std::string chemIni = "rand";

    // Colors
    Color bgcolor = Color(240, 190, 20); // was (128, 255, 128)
    Color wallColor = Color(102, 51, 0, 50);
    Color clusGreen = Color(0, 128, 0);

    // Dependent
    int numOfColloids = p::maxColloids;
    int simStep = 0;
    float indepVblInc = 0;
    if (numIncs != 0) indepVblInc = (float)(indepVblEnd - getParam(indepVbl)) / numIncs;
    int dumpsPerInc = dumpsPerRun * runsPerInc;
    int numRuns = runsPerInc * (numIncs + 1);
    int maxClusters = (int)ceil(p::maxBinCols * p::maxBinRows / 2) + 3;
    int curPoint = 0;

    float clusHistWidth = (float)clusHistMaxColls / outHistBins;

    std::vector<float> indepVblVal(numIncs + 1);
    for (int i = 0; i < numIncs + 1; i++)
    {
        indepVblVal[i] = getParam(indepVbl) + i * indepVblInc;
    }

    // Dummies
    int curRun = 0;
    int curInc = 0;
    std::vector<float> depVbl(dumpsPerInc * (numIncs + 1), 0);
    std::vector<float> depVblAve(numIncs + 1, 0);
    std::vector<float> depVblStdDev(numIncs + 1, 0);
    std::vector<float> outHist(outHistBins, 0);
    std::vector<int> dumps(numIncs + 1, 0);

    /* File initialisation
    **********************************************************************/
    std::string fileName = outputType + "_" + indepVbl;
    std::string fileHeader = "";
    if (p::walls) { fileName += "_w"; if (p::wallsSink) fileName = fileName + "s"; }
    fileName += ".txt";

    std::ofstream outputFile;
    if (dataGathering)
    {
        outputFile.open(fileName);

        std::stringstream ss;
        for (int i = 0; i < numParams; i++)
        {
            ss << "#" << params[i] << ": " << getParam(params[i]) << ", ";
        }

        ss << "\n";
        ss << "# box " << p::boxHeight << ", " << p::boxWidth << ", lC " << p::latCols << ", bC " << p::binCols << ", mCS " << p::minClusSize << ", mCD " << p::minClusDeny
        << ", wUS " << warmUpSteps << ", dPR " << dumpsPerInc << ", rPI " << runsPerInc;
        if (p::walls) ss << " wAD " << p::wallActDist << ", wF " << p::wallForce;
        ss << "\n";

        fileHeader = ss.str();
        outputFile << fileHeader;
    }

    // Fields
    std::vector<float> chem(p::latCols * p::latRows, 0);
    std::vector<float> chemProduction(p::latCols * p::latRows, 0);
    std::vector<float> bins(p::maxBinCols * p::maxBinRows, 0);
    std::vector<int> clustogram(p::maxBinCols * p::maxBinRows, 0);
    std::vector<Colloid> Colloids;
    std::vector<Cluster> Clusters;
    std::vector<Color> randColors;

    /* Field initialisation
    **********************************************************************/

    for (int i = 0; i < p::maxColloids; i++)
    {
        Colloids.emplace_back();
    }

    for (int i = 0; i < maxClusters; i++)
    {
        Clusters.emplace_back();
    }

    for (int i = 0; i < maxClusters; i++)
    {
        int maxHue = 240;
        randColors.emplace_back((int)floor(maxHue * p::rand01()), (int)floor(maxHue * p::rand01()), (int)floor(maxHue * p::rand01()));
        if (p::walls && i <= 1) randColors[i] = clusGreen;
    }

    // Seed the random number generator
    srand((unsigned)time(0));

    // Make the window
    RenderWindow window(VideoMode(p::winWidth, p::winHeight), "Janus Colloids");

    // User controlled bools
    bool paused = false;
    bool closing = false;
    bool drawing = true;
    bool histView = false;
    bool clusterFind = true;
    bool drawClusters = true;
    bool drawWalls = false;
    bool drawColloids = true;
    bool resetting = true;
    bool endRun = false;
    bool wallClusters = true;
    if (!p::walls) wallClusters = false;
    bool warmedUp = false;
    bool discardingErroneousClusters = true;
    bool rainbowDanceParty = false;

    // Debugging
    bool outputToConsole = false;
    int chemWentNeg = 0;

    // Texture for colloid
    Texture janus;
    janus.loadFromFile("j60.png");

    // The font for the HUD
    Font font;
    font.loadFromFile("Anonymous.ttf");

    // Results
    float aveChem = 0;
    float chemStdDev = 0;
    float aveSep = 0;
    float maxSep = 0;
    int numClusters = 0;
    int numValidClusters = 0;
    int numPurgedClusters = 0;
    float aveClusSize = 0;
    float clusSizeStdDev = 0;
    float aveClusDeny = 0;
    float clusDenyStdDev = 0;
    int numColloidsInClusters = 0;
    int firstNormClus = 1;
    if (p::walls && wallClusters) firstNormClus = 3;
    float histTileSize = (float)p::boxHeight / p::binRows;

    /* Simulation loop
    **********************************************************************/
    while (window.isOpen())
    {
        if (1)
        {
            /* Handle user input
            **********************************************************************/
            Event event;
            while (window.pollEvent(event))
            {
                if (event.type == Event::Closed)
                    window.close();
            }

            if (Keyboard::isKeyPressed(sf::Keyboard::D))
            {
                drawing = !drawing;
                sf::sleep(p::buttonDelay);
            }

            if (Keyboard::isKeyPressed(sf::Keyboard::H))
            {
                histView = !histView;
                sf::sleep(p::buttonDelay);
            }

            float shadeFactor = 1.3;

            if (Keyboard::isKeyPressed(sf::Keyboard::W))
            {
                p::shadeScale *= shadeFactor;
                sf::sleep(p::buttonDelay);
            }

            if (Keyboard::isKeyPressed(sf::Keyboard::S))
            {
                p::shadeScale /= shadeFactor;
                sf::sleep(p::buttonDelay);
            }

            if (Keyboard::isKeyPressed(sf::Keyboard::Escape))
            {
                closing = true;
            }

            if (!dataGathering)
            {
                if (Keyboard::isKeyPressed(sf::Keyboard::P))
                {
                    paused = !paused;
                    sf::sleep(p::buttonDelay);

                    if (paused)
                    {
                        outputToConsole = true;
                    }
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::C))
                {
                    drawColloids = !drawColloids;
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::V))
                {
                    clusterFind = !clusterFind;
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::L))
                {
                    discardingErroneousClusters = !discardingErroneousClusters;
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::Q))
                {
                    std::string param;
                    float newVal;

                    sf::sleep(p::buttonDelay);

                    std::cout << "Which parameter would you like to update?\n";
                    std::cin >> param;
                    std::cout << "Enter its new value\n";
                    std::cin >> newVal;
                    setParam(param, newVal);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::Up))
                {
                    p::simSkip += 10;
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::Down))
                {
                    p::simSkip -= 10;
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::Right))
                {
                    p::simSkip += 100;
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::Left))
                {
                    p::simSkip -= 100;
                    sf::sleep(p::buttonDelay);
                }

                // Fix skip counter if necessary
                if (p::simSkip % 100 == 1) p::simSkip -= 1;
                if (p::simSkip < 1) p::simSkip = 1;

                float binFactor = 1.2;

                if (Keyboard::isKeyPressed(sf::Keyboard::Y))
                {
                    p::binRows = (int)round(p::binRows * binFactor);
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::G))
                {
                    p::binRows = (int)round(p::binRows / binFactor);
                    sf::sleep(p::buttonDelay);
                }

                if (Keyboard::isKeyPressed(sf::Keyboard::R))
                {
                    resetting = true;
                    sf::sleep(p::buttonDelay);
                }

                // ROUTINE
                if (routine)
                {
                    if (Keyboard::isKeyPressed(sf::Keyboard::E))
                    {
                        resetting = true;
                        routineNo--;
                        sf::sleep(p::buttonDelay);
                    }
                    if (Keyboard::isKeyPressed(sf::Keyboard::T))
                    {
                        resetting = true;
                        routineNo++;
                        sf::sleep(p::buttonDelay);
                    }
                }
            }

            if (closing)
            {
                if (dataGathering) outputFile.close();
                window.close();
            }
        } // End of if "inputting"

        if (resetting)
        {
            /* ROUTINE
            **********************************************************************/
            if (routine)
            {
                int j = 0;
                p::simSkip = 10;
                p::walls = false;
                paused = true;

                if (routineNo >= j)
                {
                    routineDescription = "Diffusion Demonstration (From Gaussian)";

                    numOfColloids = 0;
                    p::walls = false;
                    chemIni = "Gauss";
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Diffusion Demonstration (Randomized)";

                    chemIni = "rand";
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Some Colloids";

                    chemIni = "rand";
                    numOfColloids = 5;
                    p::colDiameter = 100;
                    p::simSkip = 1;

                    p::decayRate = 0.2;
                    p::colVelocity = 3;
                    p::production = 0.5;
                    p::backProd = 10;
                    p::frontDest = 0;
                    p::colShift = 0.5;
                    p::chemCoup = 0.1;
                    p::rotDiff = 0.1;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Many Brownian Colloids";

                    p::simSkip = 10;
                    numOfColloids = 5000;
                    p::colDiameter = 12;

                    p::decayRate = 0.1;
                    p::colVelocity = 1;
                    p::production = 0.02;
                    p::backProd = 0.15;
                    p::frontDest = 0;
                    p::colShift = 1;
                    p::chemCoup = 0;
                    p::rotDiff = 0.4;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Strong Chemoattractive Clustering";

                    p::colShift = 0.1;
                    p::chemCoup = 0.2;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Weak Chemoattractive Clustering";

                    p::chemCoup = 0.05;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Janus Clusters";

                    p::chemCoup = -0.15;
                    p::colShift = 1;
                }
                j++;
/*
                if (routineNo >= j)
                {
                    routineDescription = "Smaller Janus Clusters";

                    p::colShift = 0.5;
                }
                j++;*/

                if (routineNo >= j)
                {
                    routineDescription = "Large Janus Clusters";

                    p::colShift = 2.5;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Untidy Waves (Hopefully)";

                    p::colShift = 0.1;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Vertically Aligned Wave (Hopefully)";

                    chemIni = "column";
                    p::colShift = 0.1;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Walls Tidying Waves";

                    chemIni = "rand";
                    p::walls = true;
                    wallClusters = true;
                    p::decayRate = 0.75;
                    p::colVelocity = 1.1;
                    p::colShift = 0.2;
                    p::chemCoup = -0.3;
                    p::rotDiff = 0.4;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Walls Admit KS Clusters";
                    wallClusters = false;

                    p::chemCoup = 0.1;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Walls Admit Janus Clusters";

                    wallClusters = true;

                    p::decayRate = 0.1;
                    p::colShift = 0.7;
                    p::chemCoup = -0.15;
                }
                j++;

                if (routineNo >= j)
                {
                    routineDescription = "Transitional State";

                    wallClusters = true;

                    p::colVelocity = 1;
                    p::production = 0.02;
                    p::backProd = 0.15;
                    p::frontDest = 0;
                    p::colShift = 0.35;
                    p::chemCoup = -0.2;
                    p::rotDiff = 0.4;
                }
                j++;

                /*
                if (routineNo >= j)
                {
                    p::decayRate = 0.1;
                    p::colVelocity = 1;
                    p::production = 0.02;
                    p::backProd = 0.15;
                    p::frontDest = 0;
                    p::colShift = 1;
                    p::chemCoup = -0.15;
                    p::rotDiff = 0.4;
                }
                j++;
                */

                if (!p::walls) wallClusters = false;
                if (routineNo >= j) routineNo = j - 1;
            }

            /* Lattice initialisations
            **********************************************************************/
            float chemBuffer = 80;
            float chemShift = 0.2;
            int halfLatCols = (int)floor(p::latCols / 2);
            int halfLatRows = (int)floor(p::latRows / 2);

            for (int i = 0; i < p::latCols; i++) {
                for (int j = 0; j < p::latRows; j++) {
                    chem[index(i, j)] = 0;
                }
            }

            // Gaussian
            if (chemIni == "Gauss")
            {
                float x = 0, y = 0, A = 300, B = 50;

                for (int i = 0; i < p::latCols; i++) {
                    for (int j = 0; j < p::latRows; j++) {
                        x = (float)(i - p::halfGridSize) / p::latCols;
                        y = (float)(j - p::halfGridSize) / p::latRows;
                        chem[index(i, j)] = A * exp(-B * (pow(x, 2) + pow(y, 2)));
                    }
                }
            }

            // unif
            if (chemIni == "unif")
            {
                for (int i = 0; i < p::latCols; i++) {
                    for (int j = 0; j < p::latRows; j++) {
                        chem[index(i, j)] = chemBuffer;
                    }
                }
            }

            // rand
            if (chemIni == "rand")
            {
                for (int i = 0; i < p::latCols; i++) {
                    for (int j = 0; j < p::latRows; j++) {
                        chem[index(i, j)] = chemBuffer + chemBuffer * chemShift * p::randpm1();
                    }
                }
            }

            if (chemIni == "wallTest")
            {
                for (int i = 0; i < p::latCols; i++) {
                    for (int j = 0; j < p::latRows; j++) {
                        if (i < (float)p::latCols / 10 || j < (float)p::latRows / 10) chem[index(i, j)] = 5;
                    }
                }
            }

            if (chemIni == "column")
            {
                for (int i = 0; i < std::min(8, p::latCols); i++) {
                    for (int j = 0; j < p::latRows; j++) {
                        chem[index(i, j)] = 100;
                    }
                }
            }

            if (chemIni == "row")
            {
                for (int i = 0; i < p::latCols; i++) {
                    for (int j = 0; j < std::min(8, p::latRows); j++) {
                        chem[index(i, j)] = 100;
                    }
                }
            }

            if (chemIni == "tlc")
            {
                for (int i = 0; i < std::min(30, p::latCols); i++) {
                    for (int j = 0; j < std::min(30, p::latRows); j++) {
                        chem[index(i, j)] = 100;
                    }
                }
            }

            if (chemIni == "cbar")
            {
                for (int i = 0; i < p::latCols; i++) {
                    for (int j = round(0.4 * p::latRows); j < 0.6 * p::latRows; j++) {
                        chem[index(i, j)] = 100;
                    }
                }
            }

            if (chemIni == "vgrad")
            {
                for (int i = 0; i < p::latCols; i++) {
                    for (int j = 0; j < p::latRows; j++) {
                        chem[index(i, j)] = 250 * abs((float)(j - halfLatRows) / p::latRows);
                    }
                }
            }

            /* Colloid initialisation
            **********************************************************************/

            //if (numOfColloids > p::maxColloids) numOfColloids = p::maxColloids;

            for (int i = 0; i < numOfColloids; i++)
            {
                Colloids[i].position.x = p::rand01() * p::boxWidth;
                Colloids[i].position.y = p::rand01() * p::boxHeight;
                Colloids[i].theta = p::rand01() * p::twoPi;
            }

             /* Analysis (as in the paper "Phoretic Interactions Generically Induce Dynamic Clusters and Wave Patterns in Active Colloids" by Benno Liebchen et al)
            **********************************************************************/
            float rho_0 = numOfColloids / (p::boxHeight * p::boxWidth);

            float Pe = p::colVelocity / (p::colShift * p::rotDiff);
            float B = p::chemCoup / (p::rotDiff * pow(p::colShift, 4));
            float K_0 = (p::production + p::backProd) / p::rotDiff;
            float K_d = p::decayRate / p::rotDiff;
            float D = p::diffCo / (pow(p::colShift, 2) * p::rotDiff);

            float benno1 = - B * K_0 * rho_0 * (16 * D + pow(Pe, 2)) / pow((4 * p::sqrtTwo * D + pow(Pe, 2)), 2);
            float benno2 = - B * rho_0 * K_0 * Pe / (2 * D);
            float keller = B * K_0 * rho_0 / (Pe * K_d);

            std::cout << "Pe: " << Pe << ", B: " << B << ", K_0: " << K_0 << ", K_d: " << K_d << ", D: " << D << "\n"
            << "Benno1: " << benno1 << ", Benno2: " << benno2 << ", Keller: " << keller << "\n";

            if (benno1 > 1) std::cout << "Benno predicts Janus Clusters\n";
            if (benno2 > 1) std::cout << "Benno predicts Waves\n";
            if (keller > 1) std::cout << "Keller predicts Clusters\n";

            simStep = 0;
            resetting = false;
            warmedUp = false;
        } // End of if "resetting"

        /* Analysis of data
        **********************************************************************
        **********************************************************************
        **********************************************************************/

        int lowerBound = 0;
        int upperBound = p::binRows;

        // Don't print every step
        if (simStep % stepsPerDump == 0 || simStep % p::simSkip == 0 || paused)
        {
            /* Chemical concentration
            **********************************************************************/
            aveChem = 0;
            chemStdDev = 0;

            for (int i = 0; i < p::latCols; i++) {
                for (int j = 0; j < p::latRows; j++) {
                    aveChem += chem[index(i, j)];
                }
            }

            aveChem = (float)aveChem / (p::latCols * p::latRows);

            for (int i = 0; i < p::latCols; i++) {
                for (int j = 0; j < p::latRows; j++) {
                    chemStdDev += pow(chem[index(i, j)] - aveChem, 2);
                }
            }

            chemStdDev = (float)chemStdDev / (p::latCols * p::latRows);
            chemStdDev = sqrt(chemStdDev);

            /* Make a colloid density histogram
            **********************************************************************/
            // First, clear the bin array and clustogram
            for (int i = 0; i < p::maxBinCols * p::maxBinRows; i++) {
                bins[i] = 0;
                clustogram[i] = - 1;
            }

            // binRows can be changed by the user so we must check it's not too large and update the dependent parameters
            if (p::binRows > p::maxBinRows) p::binRows = p::maxBinRows;
            histTileSize = p::boxHeight / p::binRows;
            p::binCols = (int)round(p::binRows * p::latCols / p::latRows);

            for (int k = 0; k < numOfColloids; k++) {
                bins[bindex((int)floor(Colloids[k].position.x / p::boxWidth * p::binCols),
                            (int)floor(Colloids[k].position.y / p::boxHeight * p::binRows))] += 1;
            }

            if (clusterFind)
            {
                /* Count the clusters
                **********************************************************************/
                // Reset cluster count
                for (int i = 0; i < maxClusters; i++)
                {
                    Clusters[i].wipe();
                }

                numClusters = 0;
                numValidClusters = 0;
                numPurgedClusters = 0;
                numColloidsInClusters = 0;
                int strip = 0;
                bool stripOverlaps = false;
                bool leftOverlap = false;

                if (wallClusters)
                {
                    float scaling = 1.3;
                    lowerBound = (int)floor(scaling * p::wallActDist / histTileSize);
                    upperBound = p::binRows - lowerBound;

                    for (int i = 0; i < p::binCols; i++)
                    {
                        for (int j = 0; j < lowerBound; j++)
                        {
                            if (bins[bindex(i, j)] != 0)
                            {
                                numClusters = 1;
                                clustogram[bindex(i, j)] = numClusters;
                                Clusters[numClusters].update(bins[bindex(i, j)]);
                            }
                        }

                        for (int j = upperBound; j < p::binRows; j++)
                        {
                            if (bins[bindex(i, j)] != 0)
                            {
                                numClusters = 2;
                                clustogram[bindex(i, j)] = numClusters;
                                Clusters[numClusters].update(bins[bindex(i, j)]);
                            }
                        }
                    }
                }

                // Loop over histogram lattice
                bool done = false;
                int sitesAssigned = 0;
                for (int k = 0; k <= 1; k++) {
                    for (int i = 0; i < p::binCols; i++) {
                        for (int j = lowerBound; j < upperBound; j++) {

                            // If the tile is empty and unassigned, assign it
                            if (bins[bindex(i, j)] == 0 && clustogram[bindex(i, j)] == - 1)
                            {
                                clustogram[bindex(i, j)] = 0;
                                sitesAssigned++;
                            }

                            // If the tile is occupied and unassigned, attempt to assign it
                            if (bins[bindex(i, j)] != 0 && clustogram[bindex(i, j)] == - 1)
                            {
                                // Build the strip (column of touching tiles)
                                strip = 1;

                                while (bins[bindex(i, j + strip)] != 0 && strip < p::binRows && (!p::walls || j + strip < upperBound))
                                {
                                    strip++;
                                }

                                // Check if first and last entry in column are connected
                                // In which case we skip this strip and pick up the tiles later
                                // Unless the strip is the entire column
                                if (j == 0 && bins[bindex(i, j - 1)] != 0 && !p::walls && strip <= p::binRows)
                                {
                                    stripOverlaps = true;
                                }

                                // We now check across the other periodic boundary
                                if (!stripOverlaps)
                                {
                                    for (int k = 0; k < strip; k++)
                                    {
                                        if (bins[bindex(i - 1, j + k)] != 0 && clustogram[bindex(i - 1, j + k)] == -1)
                                        {
                                            leftOverlap = true;
                                        }

                                        if (leftOverlap)
                                        {
                                            int horitzStrip = 1;

                                            for (int m = 1; m < p::binCols; m++)
                                            {
                                                if (bins[bindex(i + m, j + k)] != 0) horitzStrip++;
                                            }

                                            if (horitzStrip == p::binCols) leftOverlap = false;
                                        }
                                    }

                                    // If these conditions fail, we skip for now
                                    if (!leftOverlap)
                                    {
                                        int curClus = 0;
                                        int curCollsLeft = 0;

                                        // Look left along the strip for the largest neighbouring cluster
                                        for (int k = 0; k < strip; k++)
                                        {
                                            if (clustogram[bindex(i - 1, j + k)] > 0)
                                            {
                                                if (Clusters[clustogram[bindex(i - 1, j + k)]].numColloids > curCollsLeft)
                                                {
                                                    curClus = clustogram[bindex(i - 1, j + k)];
                                                    curCollsLeft = Clusters[clustogram[bindex(i - 1, j + k)]].numColloids;
                                                }
                                            }
                                        }

                                        if (curClus == 0)
                                        // We didn't find an existing neighbouring cluster
                                        // So we make a new one
                                        {
                                            numClusters++;
                                            if (numClusters >= maxClusters)
                                            {
                                                numClusters = maxClusters - 1;
                                                std::cout << "Too many clusters!\n";
                                            }
                                            curClus = numClusters;
                                        }

                                        // Add the strip to the current cluster
                                        for (int k = 0; k < strip; k++)
                                        {
                                            clustogram[bindex(i, j + k)] = curClus;
                                            Clusters[curClus].update(bins[bindex(i, j + k)]);
                                            sitesAssigned++;
                                        }
                                    }
                                }

                                j += strip - 1; // -1 because j++
                                stripOverlaps = false;
                                leftOverlap = false;

                            } // End of if non-empty

                        if (sitesAssigned > p::binCols * (upperBound - lowerBound)) { done = true; break; }

                        } // End of row-loop

                        if (done) break;

                    } // End of col-loop

                    if (done) break;

                } // End of while loop (and cluster count)

                // Check for big clusters
                for (int i = 1; i <= numClusters; i++)
                {
                    if (Clusters[i].check(histTileSize, p::minClusSize, p::minClusDeny))
                    {
                        Clusters[i].valid = true;
                        Clusters[i].ID = numValidClusters;
                        numColloidsInClusters += Clusters[i].numColloids;
                        numValidClusters++;

                        // Calculate radius
                        Clusters[i].radius = histTileSize * sqrt((float)Clusters[i].numTiles / p::pi);
                    }
                }

                // Find cluster origins
                if (0) // Origin = most occupied
                {
                    for (int i = firstNormClus; i <= numClusters; i++)
                    {
                        int curMax = 0;
                        int curMaxX = 0;
                        int curMaxY = 0;

                        if (Clusters[i].valid)
                        {
                            for (int j = 0; j < p::binCols; j++) {
                                for (int k = 0; k < p::binRows; k++) {
                                    if (clustogram[bindex(j, k)] == i) {
                                        if (bins[bindex(j, k)] > curMax)
                                        {
                                            curMax = bins[bindex(j, k)];
                                            curMaxX = j;
                                            curMaxY = k;
                                        }
                                    }
                                }
                            }

                            Clusters[i].site = Vector2f(curMaxX, curMaxY);
                            Clusters[i].position = histTileSize * Vector2f(curMaxX + 0.5, curMaxY + 0.5);
                        }
                    }
                }

                else // Origin = centre
                {
                    int halfBinCols = (int)floor(p::binCols / 2);
                    int halfBinRows = (int)floor(p::binRows / 2);

                    if (wallClusters)
                    {
                        for (int i = 1; i <= 2; i++)
                        {
                            Clusters[i].site.x = halfBinCols;
                            Clusters[i].position.x = p::boxWidth / 2;
                        }
                        Clusters[1].site.y = std::max(lowerBound - 1, 0);
                        Clusters[2].site.y = std::min(upperBound, p::binRows - 1);
                        Clusters[1].position.y = p::wallActDist;
                        Clusters[2].position.y = p::boxHeight - p::wallActDist;
                    }

                    for (int i = firstNormClus; i <= numClusters; i++) {
                        if (Clusters[i].valid) {

                            float y_Cols = 0, x_Cols = 0, y_Rows = 0, x_Rows = 0;

                            for (int j = 0; j < p::binCols; j++) {
                                for (int k = 0; k < p::binRows; k++) {
                                    if (clustogram[bindex(j, k)] == i) {
                                        float weighting = Clusters[i].numColloids;
                                        float theta_Cols = (float)p::twoPi * (j - halfBinCols) / p::binCols;
                                        y_Cols += sin(theta_Cols) * weighting;
                                        x_Cols += cos(theta_Cols) * weighting;
                                        float theta_Rows = (float)p::twoPi * (k - halfBinRows) / p::binRows;
                                        y_Rows += sin(theta_Rows) * weighting;
                                        x_Rows += cos(theta_Rows) * weighting;
                                    }
                                }
                            }

                            float aveTheta_Cols = atan2(y_Cols, x_Cols);
                            float aveTheta_Rows = atan2(y_Rows, x_Rows);

                            Clusters[i].site.x = (int)floor(halfBinCols * (aveTheta_Cols / p::pi + 1));
                            Clusters[i].site.y = (int)floor(halfBinRows * (aveTheta_Rows / p::pi + 1));
                            Clusters[i].position.x = p::boxWidth * (aveTheta_Cols / p::pi + 1) / 2;
                            Clusters[i].position.y = p::boxHeight * (aveTheta_Rows / p::pi + 1) / 2;
                        }
                    }
                }

                // Discard erroneous clusters
                if ((dataGathering && simStep >= warmUpSteps) && discardingErroneousClusters) {
                    for (int i = firstNormClus; i <= numClusters; i++) {
                        for (int j = firstNormClus; j <= numClusters; j++) {
                            if (Clusters[i].valid && Clusters[j].valid && i != j) {
                                float sepX = dist1d(Clusters[i].position.x, Clusters[j].position.x, p::boxWidth);
                                float sepY = dist1d(Clusters[i].position.y, Clusters[j].position.y, p::boxHeight);
                                if (!p::walls && ((sepX < Clusters[i].radius && sepY < Clusters[i].radius) || (sepX < Clusters[j].radius && sepY < Clusters[j].radius))                                    )
                                {
                                    // Delete the smaller one (if this condition fails then the reverse will get picked up later in the i,j loop)
                                    if (Clusters[i].numColloids > Clusters[j].numColloids)
                                    {
                                        numValidClusters--;
                                        numPurgedClusters++;

                                        Clusters[i].numColloids += Clusters[j].numColloids;
                                        Clusters[i].numTiles += Clusters[j].numTiles;

                                        Clusters[j].wipe();
                                    }
                                }
                            }
                        }
                    }
                }

                // Calculate average cluster size (colloid content) and std dev
                aveClusSize = 0;
                clusSizeStdDev = 0;
                aveClusDeny = 0;
                clusDenyStdDev = 0;

                for (int i = firstNormClus; i <= numClusters; i++)
                {
                    if (Clusters[i].valid)
                    {
                        aveClusSize += Clusters[i].numColloids;
                        aveClusDeny += Clusters[i].density(histTileSize);
                    }
                }

                aveClusSize = clusAve(aveClusSize, numValidClusters);
                aveClusDeny = clusAve(aveClusDeny, numValidClusters);

                for (int i = firstNormClus; i <= numClusters; i++)
                {
                    if (Clusters[i].valid)
                    {
                        clusSizeStdDev += pow(Clusters[i].numColloids - aveClusSize, 2);
                        clusDenyStdDev += pow(Clusters[i].density(histTileSize) - aveClusDeny, 2);
                    }
                }

                clusSizeStdDev = sqrt(clusAve(clusSizeStdDev, numValidClusters));
                clusDenyStdDev = sqrt(clusAve(clusDenyStdDev, numValidClusters));

            } // End of if clusterfind

            /* Output to file
            **********************************************************************/
            if (dataGathering && simStep >= warmUpSteps && simStep % stepsPerDump == 0)
            {
                //if (curRun < 4) drawing = true;

                if (outputType == "aveSep")
                {
                    aveSep = 0;

                    for (int i = 0; i < numOfColloids; i++) {
                        for (int j = 0; j < numOfColloids; j++) {
                            if (i != j) {
                                aveSep += dist2d(
                                    Colloids[i].position.x, Colloids[j].position.x,
                                    Colloids[i].position.y, Colloids[j].position.y);
                            }
                        }
                    }

                    aveSep = (float)aveSep / (numOfColloids * (numOfColloids - 1));

                } // End of aveSep

                if (outputType == "maxSep")
                {
                    maxSep = 0;

                    float curSep = 0;

                    for (int i = 0; i < numOfColloids; i++) {
                        for (int j = 0; j < numOfColloids; j++) {
                            if (i != j) {
                                curSep = dist2d(
                                    Colloids[i].position.x, Colloids[j].position.x,
                                    Colloids[i].position.y, Colloids[j].position.y);
                                if (curSep > maxSep)
                                {
                                    maxSep = curSep;
                                }
                            }
                        }
                    }

                }// End of maxSep

                if (outputType == "theta")
                {
                    outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc); // This clears the file
                    outputFile.close();
                    outputFile.open(fileName);
                    outputFile << fileHeader;

                    for (int i = 0; i < numOfColloids; i++)
                    {
                        outputFile << cos(Colloids[i].theta) << " " << sin(Colloids[i].theta) << "\n";
                    }

                    outputFile.close();
                }// End of theta

                if (outputType == "clus")
                {
                    depVbl[curPoint * dumpsPerInc + dumps[curPoint]] = aveClusSize;
                    dumps[curPoint]++;
                }

                if (outputType == "clusHist")
                {
                    for (int i = firstNormClus; i <= numClusters; i++)
                    {
                        if (Clusters[i].valid)
                        {
                            int bin = (int)floor((float) Clusters[i].numColloids / clusHistWidth);
                            if (bin >= outHistBins) bin = outHistBins - 1;
                            outHist[bin]++;
                        }
                    }
                    dumps[curPoint]++;
                }

                if (outputType == "collHist")
                {
                    for (int i = 0; i < numOfColloids; i++)
                    {
                        int bin = (int)floor(Colloids[i].position.y * outHistBins / p::boxHeight);
                        outHist[bin]++;
                    }
                    dumps[curPoint]++;
                }

                if (outputType == "clusPos")
                {
                    int j = (int)(simStep - warmUpSteps) / stepsPerDump;

                    for (int i = firstNormClus; i < numClusters; i++)
                    {
                        if (Clusters[i].valid)
                        {
                            outputFile << j << " " << Clusters[i].position.x << " " << Clusters[i].radius << "\n";
                        }
                    }

                    dumps[curPoint]++;
                }

                if (dumps[curPoint] % dumpsPerRun == 0) endRun = true;
            }// End of outputting

            /* Check for a reset and update parameters
            **********************************************************************/
            if (dataGathering)
            {
                if (endRun)
                {
                    curRun++;

                    if (outputType == "clus") // etc
                    {
                        outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc); // This clears the file
                        outputFile.close();
                        outputFile.open(fileName);
                        outputFile << fileHeader;

                        // Take average of depVbl
                        for (int i = 0; i < std::min(curRun, numIncs + 1); i++)
                        {
                            float ave = 0;
                            float dev = 0;

                            for (int j = 0; j < dumps[i]; j++) ave += depVbl[i * dumpsPerInc + j];
                            ave = (float)ave / dumps[i];
                            for (int j = 0; j < dumps[i]; j++) dev += pow(depVbl[i * dumpsPerInc + j] - ave, 2);
                            dev = sqrt((float) dev / dumps[i]);

                            outputFile << indepVblVal[i] << " " << ave << " " << dev << "\n";
                        }

                        outputFile.close();

                        // Update independent variable
                        curPoint = curRun % (numIncs + 1);
                        setParam(indepVbl, indepVblVal[curPoint]);
                    }

                    if (outputType == "clusHist")
                    {
                        outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc); // This clears the file
                        outputFile.close();
                        outputFile.open(fileName);
                        outputFile << fileHeader;

                        for (int i = 0; i < outHistBins; i++)
                        {
                            outputFile << i * clusHistWidth << " " << outHist[i] << "\n";
                        }

                        outputFile.close();
                    }

                    if (outputType == "collHist")
                    {
                        outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc); // This clears the file
                        outputFile.close();
                        outputFile.open(fileName);
                        outputFile << fileHeader;

                        for (int i = 0; i < outHistBins; i++)
                        {
                            outputFile << (float)i / outHistBins * p::boxHeight << " " << outHist[i] << "\n";
                        }

                        outputFile.close();
                    }

                    //if (curRun > numIncs) drawing = false;

                    // Decide whether or not to end
                    if (curRun == numRuns)
                    {
                        closing = true;
                    }
                    else
                    {
                        endRun = false;
                        resetting = true;
                    }
                }
            } // End of dataGathering
        }

        if ((float)numColloidsInClusters / numOfColloids > 0.9 && !warmedUp)
        {
            std::cout << "Warmed up at: " << simStep << "\n";
            warmedUp = true;
        }

        /* Draw the frame
        **********************************************************************
        **********************************************************************
        **********************************************************************/
        if (simStep % p::simSkip == 0 || paused)
        {
            // Clear the frame first!
            window.clear(bgcolor);

            // Display
            if (drawing)
            {
                Uint8 shade;
                RectangleShape tile;

                if (!histView)
                {
                    tile.setSize({ (float)p::tileSize * p::boxToWin, (float)p::tileSize * p::boxToWin });

                    // Draw the chemical lattice
                    for (int i = 0; i < p::latCols; i++) {
                        for (int j = 0; j < p::latRows; j++) {
                            shade = chromifySigned((aveChem - chem[index(i, j)]) * p::shadeScale); // (Uint8)floor(255 * - 0.5 * (1 + tanh((chem[index(i, j)] - aveChem) * p::shadeScale)));
                            tile.setPosition(Vector2f(i * p::tileSize * p::boxToWin, j * p::tileSize * p::boxToWin));
                            tile.setFillColor(sf::Color(shade, shade, 255));
                            window.draw(tile);
                        }
                    }

                    // Draw the walls (if necessary)
                    if (p::walls && drawWalls)
                    {
                        RectangleShape wall;
                        wall.setFillColor(wallColor);
                        wall.setSize(Vector2f(p::boxWidth * p::boxToWin, p::wallActDist * p::boxToWin));

                        wall.setPosition(Vector2f(0, 0));
                        window.draw(wall);
                        wall.setPosition(Vector2f(0, (p::boxHeight - p::wallActDist) * p::boxToWin));
                        window.draw(wall);
                    }

                    Sprite coll;

                    // Draw the colloids
                    if (drawColloids)
                    {
                        for (int i = 0; i < numOfColloids; i++)
                        {
                            // Check for colloid size updates
                            p::colScale = p::colDiameter / p::colImageSize;

                            coll = Colloids[i].getShape();

                            coll.setTexture(janus);
                            coll.setScale(Vector2f(p::colScale, p::colScale));

                            coll.setRotation(p::radToDeg * Colloids[i].theta);

                            window.draw(coll);
                        }
                    }

                } // End of if "!histView"

                if (histView)
                {
                    tile.setSize({ (float)histTileSize * p::boxToWin, (float)histTileSize * p::boxToWin });

                    // Draw histogram
                    for (int i = 0; i < p::binCols; i++) {
                        for (int j = 0; j < p::binRows; j++) {

                            shade = 255 - chromifyUnsigned(bins[bindex(i, j)] * 40); //floor(255 * (1 - tanh(((float)bins[bindex(i, j)]) * 0.2)));
                            Color tileColor = Color(255, shade, 255);

                            // If we're looking for clusters, change to green to identify they've been added to a valid one
                            if (clusterFind)
                            {
                                if (!rainbowDanceParty)
                                {
                                    if (clustogram[bindex(i, j)] > 0) {
                                        if (Clusters[clustogram[bindex(i, j)]].valid) {
                                            if (clustogram[bindex(i, j)] > 2 || !wallClusters) {
                                                tileColor = Color(shade, 255, shade);
                                            }
                                            if (clustogram[bindex(i, j)] <= 2 && wallClusters) {
                                                tileColor = Color(255, 255, shade);
                                            }
                                        }
                                    }
                                }
                                if (rainbowDanceParty)
                                {
                                    if (clustogram[bindex(i, j)] > 0) {
                                        if (Clusters[clustogram[bindex(i, j)]].valid) {
                                            if (clustogram[bindex(i, j)] > 0) {
                                                int shade2 = (int)floor(0.8 * shade);
                                                tileColor = randColors[Clusters[clustogram[bindex(i, j)]].ID] + Color(shade2, shade2, shade2);
                                            }
                                        }
                                    }
                                }
                            }

                            tile.setPosition(Vector2f(i * histTileSize * p::boxToWin, j * histTileSize * p::boxToWin));
                            tile.setFillColor(tileColor);

                            window.draw(tile);

                            if (clusterFind)
                            {
                                if (p::binRows <= 32)
                                {
                                    std::stringstream ss;
                                    ss << clustogram[bindex(i, j)];

                                    Text ID;
                                    ID.setFont(font);
                                    ID.setCharacterSize(14);
                                    ID.setFillColor(Color::Black);

                                    ID.setString(ss.str());
                                    ID.setPosition(p::boxToWin * histTileSize * Vector2f(i, j));
                                    window.draw(ID);
                                }
                            }
                        }
                    }

                    if (clusterFind)
                    {
                        // Overwrite tiles at cluster centres
                        for (int i = 0; i <= numClusters && !rainbowDanceParty; i++)
                        {
                            if (Clusters[i].valid)
                            {
                                tile.setPosition(p::boxToWin * histTileSize * Clusters[i].site);
                                tile.setFillColor(clusGreen);
                                window.draw(tile);

                                if (false)
                                {
                                    std::stringstream ss;
                                    ss << Clusters[i].ID;

                                    Text ID;
                                    ID.setFont(font);
                                    ID.setCharacterSize(12);
                                    ID.setFillColor(Color::Black);

                                    ID.setString(ss.str());
                                    ID.setPosition(p::boxToWin * histTileSize * Clusters[i].site);
                                    window.draw(ID);
                                }
                            }
                        }
                    }
                } // End of if "histView"

                if (clusterFind && !histView && drawClusters && drawColloids)
                {
                    // Draw the clusters
                    for (int i = 0; i <= numClusters; i++)
                    {
                        if (Clusters[i].valid)
                        {
                            CircleShape clus;
                            float rad = p::boxToWin * Clusters[i].radius;
                            int alpha = 255 - chromifyUnsigned((float) Clusters[i].numColloids, 75, 185); //(int)floor(60 * (1.3 + tanh(0.01 * (Clusters[i].numColloids - aveClusSize))));
                            if (rainbowDanceParty) clus.setFillColor(randColors[Clusters[i].ID] * Color(255, 255, 255, alpha));
                            if (!rainbowDanceParty) clus.setFillColor(clusGreen * Color(255, 255, 255, alpha));
                            clus.setRadius(rad);
                            clus.setOrigin(rad, rad);
                            clus.setPosition(Clusters[i].position * p::boxToWin);
                            window.draw(clus);

                            if (0)
                            {
                                std::stringstream ss;
                                ss << std::setprecision(3) << Clusters[i].density(histTileSize);

                                Text deny;
                                deny.setFont(font);
                                deny.setCharacterSize(12);
                                deny.setFillColor(Color::White);

                                deny.setString(ss.str());
                                deny.setPosition(Clusters[i].position * p::boxToWin - Vector2f(12, 6));
                                window.draw(deny);
                            }
                        }
                    }
                }

                // The bigbox covers colloids which pass over the edge of the simulation
                RectangleShape bigbox;
                bigbox.setFillColor(bgcolor);
                bigbox.setPosition({ p::boxWidth * p::boxToWin, 0 });
                bigbox.setSize({ (float)p::winWidth - p::boxWidth * p::boxToWin, p::boxHeight * p::boxToWin });

                window.draw(bigbox);
            } // End of if "drawing"

            // Draw the hud
            std::stringstream ss;
            if (routine) ss << routineDescription << "\n\n";

            ss
                << "SimStep : " << simStep << "\n"
                << "AveChem : " << std::setw(4) << aveChem << " pm " << chemStdDev << "\n";

            if (clusterFind)
            {
                ss
                << "Clusters: " << numValidClusters << " (" << numPurgedClusters << ") [" << numClusters << "]\n"
                << "AveSize : " << std::setprecision(4) << aveClusSize << " pm " << clusSizeStdDev << " > " << p::minClusSize << "\n"
                << "AveDeny : " << aveClusDeny << " pm " << clusDenyStdDev << " > " << p::minClusDeny << "\n"
                << "%ofColls: " << (float)numColloidsInClusters / numOfColloids * 100 << "%\n";

                if (wallClusters)
                {
                    ss
                    << "Top, Bot: " << Clusters[2].numColloids << ", " << Clusters[1].numColloids << "\n";
                }
            }

            if (dataGathering)
            {
                ss
                << indepVbl << std::setprecision(3) << ": " << indepVblStart << " -> " << indepVblEnd << " (" << indepVblInc << ")\n"
                << "CurPoint: " << curPoint << " / " << numIncs + 1 << "\n"
                << "CurPass : " << (int)floor(curRun / (numIncs + 1)) << " / " << runsPerInc << "\n"
                << "CurRun  : " << curRun << " / " << numRuns << "\n";
            }

            Text hud;
            hud.setFont(font);
            hud.setCharacterSize((int)round(18 * p::winScale));
            hud.setFillColor(Color::White);

            hud.setString(ss.str());
            hud.setPosition({ p::boxWidth * p::boxToWin + 20 * p::winScale, 20 * p::winScale });
            window.draw(hud);

            ss.str(""); // Clears the string

            ss
                << "NumColloids: " << numOfColloids << "\n"
                << "Lattice    : " << p::latCols << ", " << p::latRows << "\n"
                << "Box        : " << p::boxWidth << ", " << p::boxHeight << "\n"
                << "TimeStep   : " << p::timeStep << "\n"
                << "DiffCo     : " << p::diffCo << "\n"
                << "Decay      : " << p::decayRate << "\n"
                << "ColVelocity: " << p::colVelocity << "\n"
                << "Production : " << p::production << " (" << p::backProd << ", " << p::frontDest << ")\n"
                << "ColShift   : " << p::colShift << " (" << p::shiftRand << ")\n"
                << "ChemCoup   : " << p::chemCoup << "\n"
                << "RotDiff    : " << p::rotDiff << "\n\n"

                << "BinLattice : " << p::binCols << ", " << p::binRows << "\n";

            if (p::walls)
            {
                ss
                << "WallForce  : " << p::wallForce << "\n"
                << "WallActDist: " << p::wallActDist << "\n";
                if (p::wallsSink) ss << "Sink\n";
            }

            if (paused) ss << "Paused\n";
            if (!drawing) ss << "notDrawing\n";
            if (chemWentNeg > 0) ss << "chemWentNeg: " << chemWentNeg << "\n";

            float height1 = hud.getGlobalBounds().height;

            hud.setCharacterSize((int)round(16 * p::winScale));
            hud.setString(ss.str());
            hud.setPosition({ p::boxWidth * p::boxToWin + 20 * p::winScale, height1 + 40 * p::winScale });
            window.draw(hud);

            window.display();
        } // End of if not-skipped

        /* Output to console
        **********************************************************************/
        if (paused && outputToConsole && 1) // DISABLED
        {
            for (int j = 0; j < p::binRows; j++) {
                for (int i = 0; i < p::binCols; i++) {
                    std::cout << std::setw(2) << clustogram[bindex(i, j)] << " ";
                }
                std::cout << "\n";
            }

            std::cout << "\n";

            for (int i = 0; i <= numClusters; i++)
            {
                if (Clusters[i].valid)
                {
                    std::cout << i << ": " << Clusters[i].numColloids << " (" << Clusters[i].numTiles << ")\n";
                }
            }

            std::cout << "\n";

            outputToConsole = false;
        }

        /* Update
        **********************************************************************
        **********************************************************************
        **********************************************************************/

        if (!paused)
        {
            simStep++;

            float curTheta = 0;
            int curX = 0;
            int curY = 0;
            int frontX = 0;
            int frontY = 0;
            int backX = 0;
            int backY = 0;

            for (int i = 0; i < numOfColloids; i++)
            {
                curTheta = Colloids[i].theta;

                // Site containing ith colloid
                float randTheta = p::rand01() * p::twoPi;
                float randRad = p::rand01() * p::shiftRand;
                curX = (int)round((Colloids[i].position.x + randRad * cos(randTheta)) / p::tileSize);
                curY = (int)round((Colloids[i].position.y + randRad * sin(randTheta)) / p::tileSize);

                // Anisotropy
                backX = (int)round((Colloids[i].position.x - p::colShift * cos(curTheta)) / p::tileSize);
                backY = (int)round((Colloids[i].position.y - p::colShift * sin(curTheta)) / p::tileSize);

                frontX = (int)round((Colloids[i].position.x + p::colShift * cos(curTheta)) / p::tileSize);
                frontY = (int)round((Colloids[i].position.y + p::colShift * sin(curTheta)) / p::tileSize);

                if (p::walls)
                {
                    curY = wallSiteCheck(curY);
                    backY = wallSiteCheck(backY);
                    frontY = wallSiteCheck(frontY);
                }

                // Chemical Production
                chemProduction[index(curX, curY)] += p::production;
                chemProduction[index(backX, backY)] += p::backProd;
                chemProduction[index(frontX, frontY)] -= p::frontDest;

                // Rotation in response to chemical (and random wiggle)
                Colloids[i].theta += (float)sqrt(6 * p::rotDiff * p::timeStep) * p::randpm1() //generateGaussianNoise(0, 1)
                    + p::timeStep * p::chemCoup * (
                        cos(curTheta) * (chem[index(curX, curY + 1)] - chem[index(curX, curY - 1)]) -
                        sin(curTheta) * (chem[index(curX + 1, curY)] - chem[index(curX - 1, curY)])
                        )  / (2 * p::tileSize);

                if (Colloids[i].theta < 0) Colloids[i].theta += p::twoPi;
                if (Colloids[i].theta >= p::twoPi) Colloids[i].theta -= p::twoPi;

                // Move colloid in direction its facing
                Colloids[i].update(p::colVelocity, p::walls);
            }

            // Update the chemical lattice
            for (int i = 0; i < p::latCols; i++) {
                for (int j = 0; j < p::latRows; j++) {

                    if (p::walls && (j == 0 || j == p::latRows - 1))
                    {
                        if (p::wallsSink)
                        {
                            chem[index(i, j)] = 0;
                        }
                        else
                        {
                            if (j == 0) chem[index(i, j)] = chem[index(i, j + 1)];
                            if (j == p::latRows - 1) chem[index(i, j)] = chem[index(i, j - 1)];
                        }
                    }
                    else
                    {
                        chem[index(i, j)] += (float)
                            p::timeStep * ((- 4 * p::diffCo * chem[index(i, j)] +
                            p::diffCo * (
                            chem[index(i + 1, j)] + chem[index(i - 1, j)] +
                            chem[index(i, j + 1)] + chem[index(i, j - 1)]) +
                            chemProduction[index(i, j)]) / pow(p::tileSize, 2) - p::decayRate * chem[index(i, j)]);
                    }

                    if (chem[index(i, j)] < 0)
                    {
                        chemWentNeg += 1;
                    }

                    chemProduction[index(i, j)] = 0;
                }
            }

        } // End of if "!paused"

    }// End of the simulation "while" loop

    return 0;
}

int fix(int x, int N)
{
    int a = x;

    if (x < 0) a = x + N;
    if (x >= N) a = x - N;

    return a;
}

float p::rand01()
{
    return (float)(rand() - 1) / RAND_MAX;
}

float p::randpm1()
{
    float a = (float)0.5 - rand01();
    return 2 * a;
}

float clusAve(float x, int nvc)
{
    if (nvc == 0) return 0;
    else return (float)x / nvc;
}

int wallSiteCheck(int y)
{
    int a = y;

    if (y > p::latRows - 1) a = p::latRows - 1;
    if (y < 0) a = 0;

    return a;
}

float getParam(std::string str)
{
    if (str == "diffCo")
    {
        return p::diffCo;
    }
    if (str == "decayRate")
    {
        return p::decayRate;
    }
    if (str == "colVelocity")
    {
        return p::colVelocity;
    }
    if (str == "production")
    {
        return p::production;
    }
    if (str == "backProd")
    {
        return p::backProd;
    }
    if (str == "frontDest")
    {
        return p::frontDest;
    }
    if (str == "colShift")
    {
        return p::colShift;
    }
    if (str == "chemCoup")
    {
        return p::chemCoup;
    }
    if (str == "rotDiff")
    {
        return p::rotDiff;
    }
    return 0;
}

void setParam(std::string str, float val)
{
    if (str == "diffCo")
    {
        p::diffCo = val;
    }
    if (str == "decayRate")
    {
        p::decayRate = val;
    }
    if (str == "colVelocity")
    {
        p::colVelocity = val;
    }
    if (str == "production")
    {
        p::production = val;
    }
    if (str == "backProd")
    {
        p::backProd = val;
    }
    if (str == "frontDest")
    {
        p::frontDest = val;
    }
    if (str == "colShift")
    {
        p::colShift = val;
    }
    if (str == "chemCoup")
    {
        p::chemCoup = val;
    }
    if (str == "rotDiff")
    {
        p::rotDiff = val;
    }
}

// Taken from Wikipedia
float generateGaussianNoise(const float& mean, const float &stdDev)
{
    static bool hasSpare = false;
    static float spare;

    if (hasSpare) {
        hasSpare = false;
        return mean + stdDev * spare;
    }

    hasSpare = true;
    float u, v, s;
    do {
        u = p::rand01() * 2.0 - 1.0;
        v = p::rand01() * 2.0 - 1.0;
        s = u * u + v * v;
    } while ((s >= 1.0) || (s == 0.0));

    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return mean + stdDev * u * s;
}

int index(int x, int y)
{
    return fix(x, p::latCols) + fix(y, p::latRows) * p::latCols;
}

int bindex(int x, int y)
{
    return fix(x, p::binCols) + fix(y, p::binRows) * p::binCols;
}

int chromifySigned(float x, int maxVal, int minVal)
{
    int val = (int)128 + round(x);

    if (val < minVal) return minVal;
    if (val > maxVal) return maxVal;
    return val;
}

int chromifyUnsigned(float x, int maxVal, int minVal)
{
    if (x < 0) return 0;

    int val = (int)round(x);

    if (val < minVal) return minVal;
    if (val > maxVal) return maxVal;
    return val;
}

float dist1d(float x1, float x2, float box)
{
    float hi = std::max(x1, x2), lo = std::min(x1, x2);
    float inner = hi - lo, outer = box - hi + lo;
    return std::min(inner, outer);
}

float dist2d(float x1, float x2, float y1, float y2)
{
    float smallestDistSq = 4 * p::boxWidth * p::boxHeight; // Answer must be smaller than this
    float curDistSq = 0;

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            curDistSq =
                pow(x2 - x1 + (i + j) * p::boxWidth, 2) +
                pow(y2 - y1 + (i + j) * p::boxHeight, 2);

            if (curDistSq < smallestDistSq)
            {
                smallestDistSq = curDistSq;
            }
        }
    }

    return sqrt(smallestDistSq);
}
