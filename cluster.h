#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

#pragma once
#include <SFML/Graphics.hpp>

using namespace sf;

class Cluster
{
public:
    Cluster();

    int ID;
    int numTiles;
    int numColloids;

    bool valid;

    Vector2f site;
    Vector2f position;

    float radius;

    float density(float ts);
    bool check(float ts, int mcs, float md);
    void wipe();
    void update(int x);
};

#endif // CLUSTER_H_INCLUDED
