#include <iostream>
#include <fstream>
#include <math.h>
#include <proj.h>
#include <map>
#include "delaunator.hpp"
#include "RTree.h"
#include "HSL.h"

using namespace std;

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

void ReadFile(ifstream &stream, vector<double> &coords, map<pair<double, double>, double> &altitudes)
{
    // Initialisation des référentiels de coordonnées :
    PJ *P = proj_create_crs_to_crs(
        PJ_DEFAULT_CTX,
        "+proj=longlat +datum=WGS84",
        "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        NULL);
    // ^ peut prendre du temps, à ne faire qu'une seule fois

    if (0 == P)
    {
        fprintf(stderr, "Failed to create transformation object.\n");
    }
    // Storing file's data in a vector and a map
    while (!stream.eof())
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
        altitudes.insert({xy, abs(alt)});
    }

    proj_destroy(P);
}

double ComputeAltitudes(pair<double, double> point, Triangle triangle, double z0, double z1, double z2)
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

BR GetBR(const vector<double> &coords, std::map<std::pair<double, double>, double> &altitudes)
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

void InsertTriangles(delaunator::Delaunator &d, RTree &tree, double areaFilter = INFINITY)
{
    int nbSkipped = 0;

    cout << "Taux de triangles insérés dans l'arbre :" << endl;
    for (int i = 0; i < d.triangles.size(); i += 3)
    // for (int i = 0; i < 102; i += 3)
    {
        double x0 = d.coords[2 * d.triangles[i]],
               y0 = d.coords[2 * d.triangles[i] + 1],
               x1 = d.coords[2 * d.triangles[i + 1]],
               y1 = d.coords[2 * d.triangles[i + 1] + 1],
               x2 = d.coords[2 * d.triangles[i + 2]],
               y2 = d.coords[2 * d.triangles[i + 2] + 1];

        Triangle triangle(x0, y0, x1, y1, x2, y2);
        double area = triangle.Area();
        if (0 < area && area <= areaFilter)
            tree.Insert(triangle);

        else
            nbSkipped++;

        if (i % (d.triangles.size() / 300) == (d.triangles.size() - 3) % (d.triangles.size() / 300))
            ProgressBar((double)(i + 3) / d.triangles.size());
    }
    cout << "Nombre de triangles supprimés car trop larges : " << nbSkipped << "/" << d.triangles.size() / 3 << endl;
}

double HillShading(Triangle triangle, double z0, double z1, double z2, double azimut_deg = 315., double altitude_deg = 45.)
{
    double zenith_deg = 90 - altitude_deg;
    double zenith_rad = zenith_deg * M_PI / 180;

    double azimuth_math = fmod(360 - azimut_deg + 90, 360);
    double azimuth_rad = azimuth_math * M_PI / 180;

    double dz_dx, dz_dy;
    double neighs[9]; // 9 neighbor cells starting with vertex (x0, y0, z0)
    for (int i = 0; i < 9; i++)
    {
        neighs[i] = ComputeAltitudes(make_pair(triangle.vertex0.first + i % 3, triangle.vertex0.second + i / 3), triangle, z0, z1, z2);
    }
    dz_dx = ((neighs[2] + 2 * neighs[5] + neighs[8]) - (neighs[0] + 2 * neighs[3] + neighs[6])) / 8;
    dz_dy = ((neighs[6] + 2 * neighs[7] + neighs[8]) - (neighs[0] + 2 * neighs[1] + neighs[2])) / 8;

    double slope_rad = atan(1 * sqrt(dz_dx * dz_dx + dz_dy * dz_dy));
    double aspect_rad = fmod(atan2(dz_dy, -dz_dx) + 2 * M_PI, 2 * M_PI);

    return cos(zenith_rad) * cos(slope_rad) + sin(zenith_rad) * sin(slope_rad) * cos(azimuth_rad - aspect_rad);
}

void CreateRaster(ofstream &f, RTree &tree, map<pair<double, double>, double> &altitudes, BR rasterBounds, int w)
{
    const int backgroundColor = 0;
    const float azimut_deg = 315, altitude_deg = 35;

    const int h = w * (rasterBounds[1].second - rasterBounds[1].first) /
                  (rasterBounds[0].second - rasterBounds[0].first);
    printf("Génération de l'image en taille : (%d, %d)\n", w, h);
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
    const double x_step = (max_x - min_x) / w;
    const double y_step = (max_y - min_y) / h;

    double x, y = max_y;
    for (int i = 0; i < h; i++)
    {
        x = min_x;
        for (int j = 0; j < w; j++)
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

            double z = ComputeAltitudes(make_pair(x, y), triangle, z0, z1, z2);
            int hueValue = (max_z - z) * 360 / (max_z - min_z);

            int r, g, b;
            double lum = HillShading(triangle, z0, z1, z2, azimut_deg, altitude_deg);
            HSLToRGB(hueValue, .5f, lum, r, g, b);

            f << (char)max(r, 0)
              << (char)max(g, 0)
              << (char)max(b, 0);

            x += x_step;
        }
        y -= y_step;
        if (i % (h / 100) == (h - 1) % (h / 100))
            ProgressBar(double(i + 1) / h);
    }
    cout << endl;
}

int main(int argc, char const *argv[])
{
    if (argc < 3)
    {
        printf("Not enough arguments.\n");
        return -1;
    }

    int imageWidth = stoi(argv[2]); // Image size
    string filename = "../data/" + (string)argv[1];
    ifstream f_data(filename);

    vector<double> coords;
    map<pair<double, double>, double> altitudes;

    if (!f_data.is_open())
        cout << "Erreur d'ouverture de " << filename << endl;
    else
    {
        ReadFile(f_data, coords, altitudes);
        f_data.close();
    }

    // Triangulation
    delaunator::Delaunator d(coords);

    // Insertion of the triangles in a R-Tree
    RTree tree(2, 7);
    if (filename == "../data/rade_1m_IM.txt")
        InsertTriangles(d, tree, .6);
    else if (filename == "../data/Guerledan_Feb19_50cm_wgs84.txt")
        InsertTriangles(d, tree, .2);
    else
        InsertTriangles(d, tree);
    cout << tree << endl;

    // Create raster
    string raster = filename.substr(8, filename.size() - 12);
    string outputFilename = "../data/generated/" + raster + ".ppm";
    ofstream f(outputFilename);
    if (!f.is_open())
        cout << "Erreur d'ouverture de " << outputFilename << endl;
    else
    {
        BR rasterBounds = GetBR(d.coords, altitudes);
        CreateRaster(f, tree, altitudes, rasterBounds, imageWidth);
    }
    f.close();

    return EXIT_SUCCESS;
}
