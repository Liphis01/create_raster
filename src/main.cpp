#include <iostream>
#include <fstream>
#include <math.h>
#include <proj.h>
#include <delaunator.hpp>
#include <map>
#include "HSL.cpp"

#define M_PI 3.14159265358979323846 /* pi */

using namespace std;

// string filename = "rade.txt";
string filename = "guerledan.txt";

void print_triangle(delaunator::Delaunator &d, int i)
{
    printf(
        "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
        d.coords[2 * d.triangles[i]],         // tx0
        d.coords[2 * d.triangles[i] + 1],     // ty0
        d.coords[2 * d.triangles[i + 1]],     // tx1
        d.coords[2 * d.triangles[i + 1] + 1], // ty1
        d.coords[2 * d.triangles[i + 2]],     // tx2
        d.coords[2 * d.triangles[i + 2] + 1]  // ty2
    );
}

void add_coords(ifstream &stream, PJ *P, vector<double> &coords, map<pair<double, double>, double> &altitudes)
{
    PJ_COORD geo_coord, cartesian_coord;
    double alt;

    stream >> geo_coord.lpzt.phi >> geo_coord.lpzt.lam >> alt;
    geo_coord.lpzt.z = 0.;

    cartesian_coord = proj_trans(P, PJ_FWD, geo_coord);
    coords.push_back(cartesian_coord.xy.x);
    coords.push_back(cartesian_coord.xy.y);
    pair<double, double> xy = make_pair(cartesian_coord.xy.x, cartesian_coord.xy.y);
    altitudes.insert({xy, alt});
}

bool point_in_triangle(double px, double py, double tx0, double ty0, double tx1, double ty1, double tx2, double ty2)
{
    double angle0, angle1, angle2;
    angle0 = (tx0 - px) * (ty1 - py) - (ty0 - py) * (tx1 - px);
    angle1 = (tx1 - px) * (ty2 - py) - (ty1 - py) * (tx2 - px);
    angle2 = (tx2 - px) * (ty0 - py) - (ty2 - py) * (tx0 - px);
    return (angle0 * angle1) >= 0 && (angle1 * angle2) >= 0;
}

int find_triangle(double px, double py, delaunator::Delaunator &d)
{
    for (int i = 0; i < d.triangles.size(); i += 3)
    {
        if (point_in_triangle(
                px, py,
                d.coords[2 * d.triangles[i]],         // tx0
                d.coords[2 * d.triangles[i] + 1],     // ty0
                d.coords[2 * d.triangles[i + 1]],     // tx1
                d.coords[2 * d.triangles[i + 1] + 1], // ty1
                d.coords[2 * d.triangles[i + 2]],     // tx2
                d.coords[2 * d.triangles[i + 2] + 1]  // ty2
                ))
        {
            // print_triangle(d, i);
            return i;
        }
    }
    return -1;
}

double compute_alti(double px, double py, double tx0, double ty0, double tz0, double tx1, double ty1, double tz1, double tx2, double ty2, double tz2)
{
    double cx = (ty1 - ty0) * (tz2 - tz1) - (tz1 - tz0) * (ty2 - ty1);
    double cy = (tz1 - tz0) * (tx2 - tx1) - (tx1 - tx0) * (tz2 - tz1);
    double cz = (tx1 - tx0) * (ty2 - ty1) - (ty1 - ty0) * (tx2 - tx1);
    return tz0 + ((tx0 - px) * cx + (ty0 - py) * cy) / cz;
}

void get_xy_boundaries(const vector<double> &coords, double &min_x, double &min_y, double &max_x, double &max_y)
{
    min_x = INFINITY, min_y = INFINITY, max_x = 0, max_y = 0;
    for (int i = 0; i < coords.size(); i += 2)
    {
        if (coords[i] <= min_x)
            min_x = coords[i];

        if (coords[i] >= max_x)
            max_x = coords[i];
    }

    for (int i = 1; i < coords.size(); i += 2)
    {
        if (coords[i] <= min_y)
            min_y = coords[i];

        if (coords[i] >= max_y)
            max_y = coords[i];
    }
}

float shadowing(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, float azimut_deg = 315., float altitude_deg = 20.)
{
    float zenith_deg = 90 - altitude_deg;
    float zenith_rad = zenith_deg * M_PI / 180;
    // cout << zenith_rad << endl;

    float azimuth_math = fmod(360 - azimut_deg + 90, 360);
    float azimuth_rad = azimuth_math * M_PI / 180;
    // cout << azimuth_rad << endl;

    double t0, t1;
    if (y1 == y2)
        t0 = 0;
    else
        t0 = (y0 - y2) / (y1 - y2);
    if (x1 == x2)
        t1 = 0;
    else
        t1 = (x0 - x2) / (x1 - x2);
    // cout << t0 << " " << t1 << endl;

    double dz_dx = (z2 + t0 * (z1 - z2) - z0) / (x2 + t0 * (x1 - x2) - x0);
    double dz_dy = (z2 + t1 * (z1 - z2) - z0) / (y2 + t1 * (y1 - y2) - y0);
    // cout << dz_dx << " " << dz_dy << endl;

    float slope_rad = atan(1 * sqrt(dz_dx * dz_dx + dz_dy * dz_dy));
    float aspect_rad = fmod(atan2(dz_dy, -dz_dx), 2 * M_PI);
    // cout << slope_rad << " " << aspect_rad << endl;
    // cout << cos(zenith_rad) * cos(slope_rad) << " " << sin(zenith_rad) * sin(slope_rad) * cos(azimuth_rad - aspect_rad) << endl << endl;
    return cos(zenith_rad) * cos(slope_rad) + sin(zenith_rad) * sin(slope_rad) * cos(azimuth_rad - aspect_rad);
}

void draw_raster(delaunator::Delaunator &d, map<pair<double, double>, double> &altitudes, int w, int h)
{
    string raster = filename.substr(8, filename.size() - 12);
    string filename = "../data/generated/raster_" + raster + ".ppm";
    ofstream f(filename);

    int backgroundColor = 0;

    if (!f.is_open())
        cout << "Erreur d'ouverture de " << filename << endl;

    else
    {
        // File characteristics
        f << "P2" << endl
          << w << " " << h << endl
          << 255 << endl;

        double min_x, min_y, max_x, max_y;
        get_xy_boundaries(d.coords, min_x, min_y, max_x, max_y);
        double min_z = min_element(altitudes.begin(), altitudes.end(), [](const auto &x, const auto &y)
                                   { return x.second < y.second; })
                           ->second;
        double max_z = max_element(altitudes.begin(), altitudes.end(), [](const auto &x, const auto &y)
                                   { return x.second < y.second; })
                           ->second;

        // f << min_x << " " << min_y << " " << max_x << " " << max_y << endl;
        double x_step = (max_x - min_x) / w, y_step = (max_y - min_y) / h;

        double x;
        double y = max_y;
        for (int i = 0; i < w; i++)
        {
            x = min_x;
            for (int j = 0; j < h; j++)
            {
                int idx_triangle = find_triangle(x, y, d);
                if (idx_triangle == -1)
                {
                    f << backgroundColor << " ";

                    x += x_step;
                    continue;
                }

                double tx0 = d.coords[2 * d.triangles[idx_triangle]];
                double ty0 = d.coords[2 * d.triangles[idx_triangle] + 1];
                double tz0 = altitudes.find(make_pair(tx0, ty0))->second;
                double tx1 = d.coords[2 * d.triangles[idx_triangle + 1]];
                double ty1 = d.coords[2 * d.triangles[idx_triangle + 1] + 1];
                double tz1 = altitudes.find(make_pair(tx1, ty1))->second;
                double tx2 = d.coords[2 * d.triangles[idx_triangle + 2]];
                double ty2 = d.coords[2 * d.triangles[idx_triangle + 2] + 1];
                double tz2 = altitudes.find(make_pair(tx2, ty2))->second;
                double z = compute_alti(x, y,
                                        tx0, ty0, tz0,
                                        tx1, ty1, tz1,
                                        tx2, ty2, tz2);
                int hueValue = (z - min_z) * 360 / (max_z - min_z);

                int r, g, b;
                float lum = shadowing(tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2);
                // cout << lum << endl;
                HSLToRGB(hueValue, .5f, lum, r, g, b);
                f << max((int)(lum*255), 0) << " ";

                x += x_step;
            }
            f << endl;
            y -= y_step;
        }
    }

    f.close();
}

int main()
{
    // Initialisation des référentiels de coordonnées :
    PJ *P = proj_create_crs_to_crs(
        PJ_DEFAULT_CTX,
        "+proj=longlat +datum=WGS84",
        "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        NULL);
    // ^ peut prendre du temps, à ne faire qu'une seule fois

    // Deux coordonnées à exprimer dans des référentiels différents
    PJ_COORD geo_coord, cartesian_coord;

    // Position géographique en latitude/longitude de l'ENSTA
    geo_coord.lpzt.lam = -4.4720707; // longitude
    geo_coord.lpzt.phi = 48.41908;   // latitude
    geo_coord.lpzt.z = 0.;           // le z dans le référentiel ellipsoidale n'est pas nécessaire pour la projection

    /*
    // Projection géographique
    cartesian_coord = proj_trans(P, PJ_FWD, geo_coord);
    cout << "(" << geo_coord.lpzt.lam << "," << geo_coord.lpzt.phi << ")"
         << " -> "
         << "(" << cartesian_coord.xy.x << "," << cartesian_coord.xy.y << ")" << endl;
    */

    filename = "../data/" + filename;
    ifstream f_data(filename);

    if (!f_data.is_open())
        cout << "Erreur d'ouverture de " << filename << endl;

    else
    {
        int w = 40, h = 40; // Taille de l'image
        vector<double> coords;
        map<pair<double, double>, double> altitudes;
        while (!f_data.eof())
        {
            add_coords(f_data, P, coords, altitudes);
        }

        // Triangulation
        delaunator::Delaunator d(coords);
        /* print_triangle(d, 0);
        double x = 7444000, y = 2219000;
        int i = find_triangle(x, y, d);
        printf("%d\n", i); */
        // print_triangle(d, i);
        draw_raster(d, altitudes, w, h);
    }

    f_data.close();

    return EXIT_SUCCESS;
}