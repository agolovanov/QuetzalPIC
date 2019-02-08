#pragma once

class space2d {
public:
    space2d(double lx, double ly, double dx, double dy);
private:
    int nx, ny;
    double lx, ly;
    double dx, dy;
};

