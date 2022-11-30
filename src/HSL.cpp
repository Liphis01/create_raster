#include <iostream>

static float Min(float a, float b)
{
	return a <= b ? a : b;
}

static float Max(float a, float b)
{
	return a >= b ? a : b;
}

static void RGBToHSL(int R, int G, int B, int &H, int &S, float &L)
{
	float r = (R / 255.0f);
	float g = (G / 255.0f);
	float b = (B / 255.0f);

	float min = Min(Min(r, g), b);
	float max = Max(Max(r, g), b);
	float delta = max - min;

	L = (max + min) / 2;

	if (delta == 0)
	{
		H = 0;
		S = 0.0f;
	}
	else
	{
		S = (L <= 0.5) ? (delta / (max + min)) : (delta / (2 - max - min));

		float hue;

		if (r == max)
		{
			hue = ((g - b) / 6) / delta;
		}
		else if (g == max)
		{
			hue = (1.0f / 3) + ((b - r) / 6) / delta;
		}
		else
		{
			hue = (2.0f / 3) + ((r - g) / 6) / delta;
		}

		if (hue < 0)
			hue += 1;
		if (hue > 1)
			hue -= 1;

		H = (int)(hue * 360);
	}

	return;
}