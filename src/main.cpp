#include <iostream>
#include <fstream>
#include <math.h>
#include <proj.h>
#include <map>
#include "delaunator.hpp"
#include "RTree.h"
#include "HSL.h"

#define M_PI 3.14159265358979323846 /* pi */

using namespace std;

void print_triangle(delaunator::Delaunator &d, int i)
{
    printf(
        "Triangle t%d(%f, %f, %f, %f, %f, %f);\n",
        i,
        d.coords[2 * d.triangles[i]],         // tx0
        d.coords[2 * d.triangles[i] + 1],     // ty0
        d.coords[2 * d.triangles[i + 1]],     // tx1
        d.coords[2 * d.triangles[i + 1] + 1], // ty1
        d.coords[2 * d.triangles[i + 2]],     // tx2
        d.coords[2 * d.triangles[i + 2] + 1]  // ty2
    );
    printf("tree.Insert(t%d);\n", i);
}

void ProgressBar(double progress)
{
    int barWidth = 40;

    cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            cout << "=";
        else if (i == pos)
            cout << ">";
        else
            cout << " ";
    }
    cout << "] " << ceil(progress * 100.0) << " %\r";
    // cout.flush();
    fflush(stdout);

    if (progress >= 1)
        cout << endl;
}

void add_coords(ifstream &stream, PJ *P, vector<double> &coords, map<pair<double, double>, double> &altitudes)
{
    PJ_COORD geo_coord, cartesian_coord;
    double alt;

    // phi : latitude, lam : longitude
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

double compute_alti(pair<double, double> point, Triangle triangle, double z0, double z1, double z2)
{
    double px = point.first, py = point.second;
    double x0 = triangle.vertex0.first;
    double y0 = triangle.vertex0.second;
    double x1 = triangle.vertex1.first;
    double y1 = triangle.vertex1.second;
    double x2 = triangle.vertex2.first;
    double y2 = triangle.vertex2.second;

    double cx = (y1 - y0) * (z2 - z1) - (z1 - z0) * (y2 - y1);
    double cy = (z1 - z0) * (x2 - x1) - (x1 - x0) * (z2 - z1);
    double cz = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1);
    return z0 + ((x0 - px) * cx + (y0 - py) * cy) / cz;
}

void compute_derivatives(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double &dz_dx, double &dz_dy)
{
    double cx = (y1 - y0) * (z2 - z1) - (z1 - z0) * (y2 - y1);
    double cy = (z1 - z0) * (x2 - x1) - (x1 - x0) * (z2 - z1);
    double cz = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1);
    dz_dx = - cx / cz;
    dz_dy = - cy / cz;
    // double denom =   (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1);
    // dz_dx =   (z1 - z0) * (y2 - y1) - (y1 - y0) * (z2 - z1) / denom; 
    // dz_dy = - (z1 - z0) * (x2 - x1) - (x1 - x0) * (z2 - z1) / denom;
}

BR get_xyz_boundaries(const vector<double> &coords, std::map<std::pair<double, double>, double> &altitudes)
{
    double min_x = INFINITY, min_y = INFINITY, max_x = 0, max_y = 0;
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

    double min_z = min_element(altitudes.begin(), altitudes.end(), [](const auto &x, const auto &y)
                        { return x.second < y.second; })
                ->second;
    double max_z = max_element(altitudes.begin(), altitudes.end(), [](const auto &x, const auto &y)
                        { return x.second < y.second; })
                ->second;
    
    return {make_pair(min_x, max_x), make_pair(min_y, max_y), make_pair(min_z, max_z)};
}

void InsertTriangles(delaunator::Delaunator &d, RTree &tree)
{
    double area_filter = INFINITY;

    cout << "Taux de triangles insérés dans l'arbre :" << endl;
    for (int i = 0; i < d.triangles.size(); i += 3)
    // for (int i = 0; i < 102; i += 3)
    {
        double x0 = d.coords[2 * d.triangles[i]],      // x0
            y0 = d.coords[2 * d.triangles[i] + 1],     // y0
            x1 = d.coords[2 * d.triangles[i + 1]],     // x1
            y1 = d.coords[2 * d.triangles[i + 1] + 1], // y1
            x2 = d.coords[2 * d.triangles[i + 2]],     // x2
            y2 = d.coords[2 * d.triangles[i + 2] + 1]; // y2

        Triangle triangle(x0, y0, x1, y1, x2, y2);
        double area = triangle.Area();
        if (triangle.Area() <= area_filter)
            tree.Insert(triangle);

        ProgressBar((double)(i + 3) / d.triangles.size());
        // ProgressBar((double)(i + 3) / 102);
    }
}

double hill_shading(Triangle triangle, double z0, double z1, double z2, double azimut_deg = 315., double altitude_deg = 45.)
{
    double zenith_deg = 90 - altitude_deg;
    double zenith_rad = zenith_deg * M_PI / 180;

    double azimuth_math = fmod(360 - azimut_deg + 90, 360);
    double azimuth_rad = azimuth_math * M_PI / 180;

    double dz_dx, dz_dy;
    double neighs[9]; // 9 neighbor cells starting with vertex (x0, y0, z0)
    for (int i = 0; i < 9; i++)
    {
        neighs[i] = compute_alti(make_pair(triangle.vertex0.first + i%3, triangle.vertex0.second + i/3), triangle, z0, z1, z2);
    }
    // compute_derivatives(x0, y0, z0, x1, y1, z1, x2, y2, z2, dz_dx, dz_dy);
    dz_dx = ((neighs[2] + 2*neighs[5] + neighs[8]) - (neighs[0] + 2*neighs[3] + neighs[6])) / 8;
    dz_dy = ((neighs[6] + 2*neighs[7] + neighs[8]) - (neighs[0] + 2*neighs[1] + neighs[2])) / 8;

    double slope_rad = atan(1 * sqrt(dz_dx * dz_dx + dz_dy * dz_dy));
    double aspect_rad = fmod(atan2(dz_dy, -dz_dx) + 2 * M_PI, 2 * M_PI);

    // printf("%f, %f, %f, %f, %f\n", cos(zenith_rad), cos(slope_rad), sin(zenith_rad), sin(slope_rad), cos(azimuth_rad - aspect_rad));
    return cos(zenith_rad) * cos(slope_rad) + sin(zenith_rad) * sin(slope_rad) * cos(azimuth_rad - aspect_rad);
}

void draw_raster(RTree &tree, map<pair<double, double>, double> &altitudes, BR rasterBounds, int w, int h, string filename)
{
    string raster = filename.substr(8, filename.size() - 12);
    cout << "Génération de raster_" << raster << ".ppm" << endl;
    string outputFilename = "../data/generated/raster_" + raster + ".ppm";
    ofstream f(outputFilename);

    int backgroundColor = 0;

    if (!f.is_open())
        cout << "Erreur d'ouverture de " << outputFilename << endl;

    else
    {
        // File characteristics
        f << "P6" << endl
          << w << " " << h << endl
          << 255 << endl;

        double min_x, max_x, min_y, max_y, min_z, max_z;
        min_x = rasterBounds[0].first;
        max_x = rasterBounds[0].second;
        min_y = rasterBounds[1].first;
        max_y = rasterBounds[1].second;
        min_z = rasterBounds[2].first;
        max_z = rasterBounds[2].second;
        double x_step = (max_x - min_x) / w;
        double y_step = (max_y - min_y) / h;

        double x, y = max_y;
        cout << "Avancement de la génération de l'image :" << endl;
        for (int i = 0; i < w; i++)
        {
            x = min_x;
            for (int j = 0; j < h; j++)
            {
                vector<Triangle> resultSearch = tree.Search(make_pair(x, y));
                if (resultSearch.empty())
                {
                    f << (char)backgroundColor
                      << (char)backgroundColor
                      << (char)backgroundColor;
                    x += x_step;
                    continue;
                }

                Triangle triangle = resultSearch[0];

                double z0 = altitudes.find(triangle.vertex0)->second;
                double z1 = altitudes.find(triangle.vertex1)->second;
                double z2 = altitudes.find(triangle.vertex2)->second;

                double z = compute_alti(make_pair(x, y), triangle, z0, z1, z2);
                int hueValue = (z - min_z) * 360 / (max_z - min_z);

                int r, g, b;
                double lum = hill_shading(triangle, z0, z1, z2, 315, 25);
                HSLToRGB(hueValue, .5f, lum, r, g, b);
                if (r < 0 || g < 0 || b < 0)
                {
                    // printf("rgb : %d, %d, %d\nhue, lum : %d, %f\n", r, g, b, hueValue, lum);
                    // printf("compute_derivatives(%f, %f, %f, %f, %f, %f, %f, %f, %f, dz_dx, dz_dy);\n\n", tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2);
                }
                f << (char)max(r, 0)
                  << (char)max(g, 0)
                  << (char)max(b, 0);

                x += x_step;
            }
            y -= y_step;
            ProgressBar(double(i+1)/w);
        }
    }

    f.close();
}

int main(int argc, char *argv[])
{
    // Initialisation des référentiels de coordonnées :
    PJ *P = proj_create_crs_to_crs(
        PJ_DEFAULT_CTX,
        "+proj=longlat +datum=WGS84",
        "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        NULL);
    // ^ peut prendre du temps, à ne faire qu'une seule fois

    if (0 == P) {
        fprintf(stderr, "Failed to create transformation object.\n");

        return 1;
    }

    if (argc < 3)
    {
        printf("Not enough arguments.\n");
        // return -1;
        argv[1] = (char *)"rade.txt";
        argv[2] = (char *)"100";
    }

    int w = stoi(argv[2]), h = stoi(argv[2]); // Image size
    string filename = "../data/" + (string)argv[1];
    ifstream f_data(filename);

    vector<double> coords;
    map<pair<double, double>, double> altitudes;

    if (!f_data.is_open())
        cout << "Erreur d'ouverture de " << filename << endl;

    else
    {
        // Storing file's data in a vector and a map
        while (!f_data.eof())
            add_coords(f_data, P, coords, altitudes);
        
        proj_destroy(P);
        f_data.close();
    }

    // Triangulation
    delaunator::Delaunator d(coords);
    
    // Insertion of the triangles in a R-Tree
    RTree tree(2, 5);
    InsertTriangles(d, tree);

    // Draw raster
    double min_x, max_x, min_y, max_y, min_z, max_z;
    BR boundaries = get_xyz_boundaries(d.coords, altitudes);
    min_x = boundaries[0].first;
    max_x = boundaries[0].second;
    min_y = boundaries[1].first;
    max_y = boundaries[1].second;
    min_z = boundaries[2].first;
    max_z = boundaries[2].second;
    cout << min_x << " " << min_y << " " << max_x << " " << max_y << endl;
    BR rasterBounds = {make_pair(min_x, max_x), make_pair(min_y, max_y), make_pair(min_z, max_z)};
    draw_raster(tree, altitudes, rasterBounds, w, h, filename);

    return EXIT_SUCCESS;
}


