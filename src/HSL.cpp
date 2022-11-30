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

static float HueToRGB(float v1, float v2, float vH)
{
	if (vH < 0)
		vH += 1;

	if (vH > 1)
		vH -= 1;

	if ((6 * vH) < 1)
		return (v1 + (v2 - v1) * 6 * vH);

	if ((2 * vH) < 1)
		return v2;

	if ((3 * vH) < 2)
		return (v1 + (v2 - v1) * ((2.0f / 3) - vH) * 6);

	return v1;
}

static void HSLToRGB(int H, float S, float L, int &r, int &g, int &b)
{
	if (S == 0)
	{
		r = g = b = (int)(L * 255);
	}
	else
	{
		float v1, v2;
		float hue = (float)H / 360;

		v2 = (L < 0.5) ? (L * (1 + S)) : ((L + S) - (L * S));
		v1 = 2 * L - v2;

		r = (int)(255 * HueToRGB(v1, v2, hue + (1.0f / 3)));
		g = (int)(255 * HueToRGB(v1, v2, hue));
		b = (int)(255 * HueToRGB(v1, v2, hue - (1.0f / 3)));
	}

	return;
}