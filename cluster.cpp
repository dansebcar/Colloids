#include "cluster.h"
#include "parameter.h"

// Constructor
Cluster::Cluster()
{
    numTiles = 0;
    numColloids = 0;

    valid = false;
}

float Cluster::density(float ts)
{
    if (numTiles != 0) return (float)numColloids / (numTiles * pow(ts, 2));
    else return 0;
}

bool Cluster::check(float ts, int mcs, float md)
{
    if (numColloids > mcs && density(ts) > md) return true;
    else return false;
}

void Cluster::wipe()
{
    ID = 0;
    numTiles = 0;
    numColloids = 0;

    valid = false;

    site = Vector2f(0, 0);
    position = Vector2f(0, 0);

    radius = 0;
}

void Cluster::update(int x)
{
    numTiles += 1;
    numColloids += x;
}
