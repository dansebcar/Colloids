#include "colloid.h"
#include "parameter.h"

// Constructor
Colloid::Colloid()
{
    position = { 0, 0 };
    theta = 0;

    colloidShape.setOrigin((p::colImageSize) / 2, (p::colImageSize) / 2); // Sets colloid "origin" as centre
}

Sprite Colloid::getShape()
{
    colloidShape.setPosition(position * p::boxToWin);
    return colloidShape;
}

void Colloid::update(float vel, bool walls)
{
    position.x += p::timeStep * vel * cos(theta);

    float topActDist = p::boxHeight - p::wallActDist;

    if (walls && (position.y < p::wallActDist || position.y > topActDist))
    {
        if (position.y < p::wallActDist) position.y += p::timeStep * (vel * sin(theta) - p::wallForce * (position.y - p::wallActDist));
        if (position.y > topActDist) position.y += p::timeStep * (vel * sin(theta) - p::wallForce * (position.y - topActDist));
    }
    else
    {
        position.y += p::timeStep * vel * sin(theta);
    }

    // Enforce periodic boundary
    if (position.y >= p::boxHeight) position.y -= p::boxHeight;
    if (position.y < 0) position.y += p::boxHeight;
    if (position.x >= p::boxWidth) position.x -= p::boxWidth;
    if (position.x < 0) position.x += p::boxWidth;
}

float rand01()
{
    return (float)(rand() - 1) / RAND_MAX;
}

float randpm1()
{
    float a = (float)0.5 - rand01();
    return 2 * a;
}
