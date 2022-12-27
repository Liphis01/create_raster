#ifndef __HSL__H_
#define __HSL__H_

/**
 * @brief Returns the smaller of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The smaller of the two values.
 */
static double Min(double a, double b);

/**
 * @brief Returns the larger of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The larger of the two values.
 */
static double Max(double a, double b);

/**
 * @brief Converts an RGB color to HSL.
 *
 * @param R Red component of the RGB color (0-255).
 * @param G Green component of the RGB color (0-255).
 * @param B Blue component of the RGB color (0-255).
 * @param H Hue component of the HSL color (0-360).
 * @param S Saturation component of the HSL color (0-100).
 * @param L Lightness component of the HSL color (0-1).
 */
void RGBToHSL(int R, int G, int B, int &H, int &S, double &L);

/**
 * @brief Helper function for converting an RGB color to HSL.
 *
 * @param v1 First value.
 * @param v2 Second value.
 * @param vH Third value.
 * @return The result of the conversion.
 */
double hueToRGB(double v1, double v2, double vH);

/**
 * @brief Converts an HSL color to RGB.
 *
 * @param H Hue component of the HSL color (0-360).
 * @param S Saturation component of the HSL color (0-100).
 * @param L Lightness component of the HSL color (0-1).
 * @param r Red component of the RGB color (0-255).
 * @param g Green component of the RGB color (0-255).
 * @param b Blue component of the RGB color (0-255).
 */
void HSLToRGB(int H, double S, double L, int &r, int &g, int &b);

#endif
