#include "viewer.hpp"

namespace V = COL781::Viewer;
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//#include <bits/stdc++.h>
#include <random>
#include<map>
#include<string>

#define DEG_TO_RAD (M_PI / 180.0)

using namespace glm;
using namespace std;

bool add_noise = false;

glm::vec3 addRandomNoise(glm::vec3 input, float noiseStrength) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(-noiseStrength, noiseStrength);

    return input + glm::vec3(dist(gen), dist(gen), dist(gen));
}

void writeOBJ(const std::string& filename, const std::vector<vec3>& vertices, const std::vector<vector<int>>& faces) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    for (const auto& v : vertices) {
        vec3 temp = v;
        if(add_noise){
            temp = addRandomNoise(v, 0.02);
        }
        file << "v " << temp.x << " " << temp.y << " " << temp.z << "\n";
    }
    for (const auto& f : faces) {
        file << "f ";
        for(int i = 0; i < f.size(); i++){
            file << f[i] << " ";
        }file << "\n";
    }
    file.close();
    std::cout << "OBJ file saved: " << filename << std::endl;
}

void create_square(int m, int n){
    string filename = "./meshes/square_";
    // filename += to_string(m);
    // filename += "_";
    // filename += to_string(n);
    vector<vec3> vertices;
    vector<vector<int>> faces;
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            vertices.push_back(vec3(j*float(1.0/float(n)), i*float(1.0/float(m)), 0.0));
            if(i != 0 && j != 0){
                vector<int> temp = {i*n+j+1, i*n+j, (i-1)*n+j, (i-1)*n+j+1};
                faces.push_back(temp);
            }
        }
    }
    filename += ".obj";
    writeOBJ(filename, vertices, faces);
}

vec3 spherical_coord(float lat, float longi){
    lat *= DEG_TO_RAD;
    longi *= DEG_TO_RAD;
    float x = 0.50*sin(lat)*cos(longi);
    float y = 0.50*sin(lat)*sin(longi);
    float z = 0.50*cos(lat);
    return vec3(x, y, z);
}

void create_sphere(int m, int n){
    string filename = "./meshes/sphere_";
    // filename += to_string(m);
    // filename += "_";
    // filename += to_string(n);
    vector<vec3> vertices;
    vector<vector<int>> faces;
    for(int i = 0; i <= n; i++){
        float lat = (180.0/n)*float(i);
        float longi = 0.0;
        if(i == 0 || i == n){
            // if(i == 0){lat = 0.0;}
            // else lat = 180.0;
            vertices.push_back(spherical_coord(lat, longi));
            // cout << spherical_coord(lat, longi).z << endl;
            // cout << cos(180.0) << endl;
        }
        if(i == 0){continue;}
        if(i == n){
            int start = 1+(n-2)*m;
            for(int j = 0; j < m; j++){
                int k = (j+1)%m;
                vector<int> temp = {m*(n-1)+2, start+j+1, start+k+1};
                faces.push_back(temp);
            }
            continue;
        }
        for(int j = 0; j < m; j++){
            int k = (j-1+m)%m;
            longi = (-180.0+(360.0/m)*float(j));
            vertices.push_back(spherical_coord(lat, longi));
            if(i == 1){
                vector<int> temp = {1, j+2, k+2};
                faces.push_back(temp);
            }
            else{
                int start = (i-1)*m;
                vector<int> temp = {start+j+2, start+k+2, start+k+2-m, start+j+2-m};
                faces.push_back(temp);
            }
        }
    }
    filename += ".obj";
    writeOBJ(filename, vertices, faces);
}

void create_pot(int m, int n){
    string filename = "./meshes/pot";
    // filename += to_string(m);
    // filename += "_";
    // filename += to_string(n);
    vector<vec3> vertices;
    vector<vector<int>> faces;
    for(int i = 5; i <= n; i++){
        float lat = (180.0/n)*float(i);
        float longi = 0.0;
        if(i == 0 || i == n){
            // if(i == 0){lat = 0.0;}
            // else lat = 180.0;
            vertices.push_back(spherical_coord(lat, longi));
            // cout << spherical_coord(lat, longi).z << endl;
            // cout << cos(180.0) << endl;
        }
        if(i == 0){continue;}
        if(i == n){
            int start = 1+(n-2)*m;
            for(int j = 0; j < m; j++){
                int k = (j+1)%m;
                vector<int> temp = {m*(n-1)+2 - 81, start+j+1 - 81, start+k+1 - 81};
                faces.push_back(temp);
            }
            continue;
        }
        for(int j = 0; j < m; j++){
            int k = (j-1+m)%m;
            longi = (-180.0+(360.0/m)*float(j));
            vertices.push_back(spherical_coord(lat, longi));
            if(i == 1){
                vector<int> temp = {1, j+2, k+2};
                faces.push_back(temp);
            }
            else{
                int start = (i-1)*m;
                if(std::min(start+k+2-m - 81, start+j+2-m - 81) >= 1) {
                    vector<int> temp = {start+j+2 - 81, start+k+2 - 81, start+k+2-m - 81, start+j+2-m - 81};
                    faces.push_back(temp);
                }
            }
        }
    }
    filename += ".obj";
    writeOBJ(filename, vertices, faces);
}

struct Vec3Compare {
    bool operator()(const glm::vec3& a, const glm::vec3& b) const {
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.z < b.z;
    }
};

vec3 cube_coord(int i, int j, int k, int m, int n, int o){
    return vec3(-0.5+float(i)*(1/float(m)), -0.5+float(j)*(1/float(n)), -0.5+float(k)*(1/float(o)));
}
void create_cube(int m, int n, int o){
    string filename = "./meshes/cube_";
    // filename += to_string(m);
    // filename += "_";
    // filename += to_string(n);
    // filename += "_";
    // filename += to_string(o);
    vector<vec3> vertices;
    vector<vector<int>> faces;
    map<vec3, int, Vec3Compare> index;
    for(int i = 0; i <= m; i++){
        for(int j = 0; j <= n; j++){
            for(int k = 0; k <= o; k++){
                if(i == 0 || i == m || j == 0 || j == n || k == 0 || k == o){
                    if(index[vec3(-0.5+float(i)*(1/float(m)), -0.5+float(j)*(1/float(n)), -0.5+float(k)*(1/float(o)))]){continue;}
                    vertices.push_back(vec3(-0.5+float(i)*(1/float(m)), -0.5+float(j)*(1/float(n)), -0.5+float(k)*(1/float(o))));
                    index[vec3(-0.5+float(i)*(1/float(m)), -0.5+float(j)*(1/float(n)), -0.5+float(k)*(1/float(o)))] = vertices.size();
                }
            }
        }
    }
    for(int i = 0; i <= m; i++){
        for(int j = 0; j <= n; j++){
            for(int k = 0; k <= o; k++){
                if((i == 0 || i == m) && j != 0 && k != 0){
                    vec3 a = cube_coord(i, j, k, m, n, o);
                    vec3 b = cube_coord(i, j-1, k, m, n, o);
                    vec3 c = cube_coord(i, j-1, k-1, m, n, o);
                    vec3 d = cube_coord(i, j, k-1, m, n, o);
                    vector<int> temp = {index[a], index[b], index[c], index[d]};
                    faces.push_back(temp);
                }
                if((j == 0 || j == n) && i != 0 && k != 0){
                    vec3 a = cube_coord(i, j, k, m, n, o);
                    vec3 b = cube_coord(i-1, j, k, m, n, o);
                    vec3 c = cube_coord(i-1, j, k-1, m, n, o);
                    vec3 d = cube_coord(i, j, k-1, m, n, o);
                    vector<int> temp = {index[a], index[b], index[c], index[d]};
                    faces.push_back(temp);
                }
                if((k == 0 || k == o) && i != 0 && j != 0){
                    vec3 a = cube_coord(i, j, k, m, n, o);
                    vec3 b = cube_coord(i, j-1, k, m, n, o);
                    vec3 c = cube_coord(i-1, j-1, k, m, n, o);
                    vec3 d = cube_coord(i-1, j, k, m, n, o);
                    vector<int> temp = {index[a], index[b], index[c], index[d]};
                    faces.push_back(temp);
                }
            }
        }
    }
    filename += ".obj";
    writeOBJ(filename, vertices, faces);
}

void create_noisy_cube(int m, int n, int o){
    string filename = "./meshes/noisy_cube_";
    // filename += to_string(m);
    // filename += "_";
    // filename += to_string(n);
    // filename += "_";
    // filename += to_string(o);
    vector<vec3> vertices;
    vector<vector<int>> faces;
    map<vec3, int, Vec3Compare> index;
    for(int i = 0; i <= m; i++){
        for(int j = 0; j <= n; j++){
            for(int k = 0; k <= o; k++){
                if(i == 0 || i == m || j == 0 || j == n || k == 0 || k == o){
                    if(index[vec3(-0.5+float(i)*(1/float(m)), -0.5+float(j)*(1/float(n)), -0.5+float(k)*(1/float(o)))]){continue;}
                    vertices.push_back(vec3(-0.5+float(i)*(1/float(m)), -0.5+float(j)*(1/float(n)), -0.5+float(k)*(1/float(o))));
                    index[vec3(-0.5+float(i)*(1/float(m)), -0.5+float(j)*(1/float(n)), -0.5+float(k)*(1/float(o)))] = vertices.size();
                }
            }
        }
    }
    for(int i = 0; i <= m; i++){
        for(int j = 0; j <= n; j++){
            for(int k = 0; k <= o; k++){
                if((i == 0 || i == m) && j != 0 && k != 0){
                    vec3 a = cube_coord(i, j, k, m, n, o);
                    vec3 b = cube_coord(i, j-1, k, m, n, o);
                    vec3 c = cube_coord(i, j-1, k-1, m, n, o);
                    vec3 d = cube_coord(i, j, k-1, m, n, o);
                    vector<int> temp = {index[a], index[b], index[c], index[d]};
                    faces.push_back(temp);
                }
                if((j == 0 || j == n) && i != 0 && k != 0){
                    vec3 a = cube_coord(i, j, k, m, n, o);
                    vec3 b = cube_coord(i-1, j, k, m, n, o);
                    vec3 c = cube_coord(i-1, j, k-1, m, n, o);
                    vec3 d = cube_coord(i, j, k-1, m, n, o);
                    vector<int> temp = {index[a], index[b], index[c], index[d]};
                    faces.push_back(temp);
                }
                if((k == 0 || k == o) && i != 0 && j != 0){
                    vec3 a = cube_coord(i, j, k, m, n, o);
                    vec3 b = cube_coord(i, j-1, k, m, n, o);
                    vec3 c = cube_coord(i-1, j-1, k, m, n, o);
                    vec3 d = cube_coord(i-1, j, k, m, n, o);
                    vector<int> temp = {index[a], index[b], index[c], index[d]};
                    faces.push_back(temp);
                }
            }
        }
    }
    filename += ".obj";
    writeOBJ(filename, vertices, faces);
}

int main(int argc, char* argv[]) {


    string shape_type = argv[1];
    if(shape_type == "square") {
        int l = atoi(argv[2]);
        int b = atoi(argv[3]);
        create_square(l, b);
    }
    else if(shape_type == "cube") {
        int l = atoi(argv[2]);
        int b = atoi(argv[3]);
        int h = atoi(argv[4]);
        //create_cube(l, b, h);
        add_noise = true;
        create_noisy_cube(l, b, h);
    }
    else if(shape_type == "sphere") {
        int l = atoi(argv[2]);
        int b = atoi(argv[3]);
        create_sphere(l, b);
    }
    else if(shape_type == "pot") {
        int l = atoi(argv[2]);
        int b = atoi(argv[3]);
        create_pot(l, b);
    }

    return 0;
}
