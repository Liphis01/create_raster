#ifndef __HSL__H_
#define __HSL__H_

static double Min(double a, double b);

static double Max(double a, double b);

void RGBToHSL(int R, int G, int B, int &H, int &S, double &L);

double HueToRGB(double v1, double v2, double vH);

void HSLToRGB(int H, double S, double L, int &r, int &g, int &b);

#endif