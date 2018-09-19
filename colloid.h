#ifndef COLLOID_H_INCLUDED
#define COLLOID_H_INCLUDED

#pragma once
#include <SFML/Graphics.hpp>

using namespace sf;

class Colloid
{
private:
    Sprite colloidShape;

public:
    Colloid();

    Vector2f position;
    float theta;

    Sprite getShape();
    void update(float vel, bool walls);
};

#endif // COLLOID_H_INCLUDED
