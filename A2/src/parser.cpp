#include "viewer.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include<set>
// #include </opt/homebrew/include/c++/14/bits>

namespace V = COL781::Viewer;
using namespace glm;
using namespace std;

struct Polygon {vector<int> vertexIndices, texCoordIndices, normalIndices;};

struct HEdge;
// Face structure
struct Face {
    int index;
    HEdge* h = nullptr;  // One of the half-edges on the face
};

// Vertex structure
struct Vertex {
    int index;
    HEdge* h = nullptr;  // One of the outgoing half-edges
};

// Half-Edge structure
struct HEdge {
    HEdge* pair = nullptr;  // Opposite half-edge
    HEdge* next = nullptr;  // Next half-edge in the face cycle
    Vertex* v = nullptr;    // Associated vertex
    Face* f = nullptr;      // Associated face
};

std::vector<vec3> vertices;
std::vector<Vertex*> Vert_Pointer;
std::vector<vec2> texCoords;
std::vector<vec3> normals;
std::vector<Polygon> faces;
std::vector<Face*> Face_Pointer;
std::vector<vec2> edges;
set<pair<int,int>> actual_edges;
map<pair<int, int>, HEdge*> Init;
map<pair<int, int>, vector<int>> Neigh;

bool smoothening = false;
bool extrusion = false;


void loadOBJ(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if (prefix == "v") {
            vec3 vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            vertices.push_back(vertex);
            Vertex *temp = new Vertex;
            temp->index = Vert_Pointer.size();
            Vert_Pointer.push_back(temp);
        } else if (prefix == "vt") {
            vec2 texCoord;
            iss >> texCoord.x >> texCoord.y;
            texCoords.push_back(texCoord);
        } else if (prefix == "vn") {
            vec3 normal;
            iss >> normal.x >> normal.y >> normal.z;
            normals.push_back(normal);
        } 
        else if (prefix == "f") {
            // continue;
            Polygon face;
            std::string vertex;
            
            while (iss >> vertex) {
                std::istringstream vss(vertex);
                std::string index;
                int vi = -1, ti = -1, ni = -1; // Initialize to -1 to handle missing values

                // Read vertex index
                if (std::getline(vss, index, '/') && !index.empty()) {
                    vi = std::stoi(index) - 1;
                }

                // Read texture index (optional)
                if (std::getline(vss, index, '/') && !index.empty()) {
                    ti = std::stoi(index) - 1;
                }

                // Read normal index (optional)
                if (std::getline(vss, index, '/') && !index.empty()) {
                    ni = std::stoi(index) - 1;
                }

                // Store indices in the face structure
                face.vertexIndices.push_back(vi);
                face.texCoordIndices.push_back(ti);
                face.normalIndices.push_back(ni);
            }
            // faces.push_back(face);
            for(int i = 1; i < (int)(face.vertexIndices.size())-1; i++){
                Polygon temp = Polygon();
                int index = faces.size();
                temp.vertexIndices.push_back(face.vertexIndices[0]);
                temp.vertexIndices.push_back(face.vertexIndices[i]);
                temp.vertexIndices.push_back(face.vertexIndices[i+1]);
                temp.texCoordIndices.push_back(face.texCoordIndices[0]);
                temp.texCoordIndices.push_back(face.texCoordIndices[i]);
                temp.texCoordIndices.push_back(face.texCoordIndices[i+1]);
                temp.normalIndices.push_back(face.normalIndices[0]);
                temp.normalIndices.push_back(face.normalIndices[i]);
                temp.normalIndices.push_back(face.normalIndices[i+1]);
                faces.push_back(temp);
                Neigh[{face.vertexIndices[0], face.vertexIndices[i]}].push_back(index);
                Neigh[{face.vertexIndices[i+1], face.vertexIndices[0]}].push_back(index);
                Neigh[{face.vertexIndices[i+1], face.vertexIndices[i]}].push_back(index);
                if(i != 1){
                    int i1 = face.vertexIndices[0];
                    int i2 = face.vertexIndices[i];
                    HEdge* e1 = new HEdge; 
                    HEdge* e2 = new HEdge;
                    e1 -> pair = e2;
                    e2 -> pair = e1;
                    e1->v = Vert_Pointer[i1];
                    Vert_Pointer[i1]->h = e1;
                    e2->v = Vert_Pointer[i2];
                    Vert_Pointer[i2]->h = e2;
                    Init[{i1, i2}] = e1;
                    Init[{i2, i1}] = e2;
                }
            }
            for(int i = 0;i < face.vertexIndices.size(); i++){
                int j = (i+1)%(int)(face.vertexIndices.size());
                vec2 temp;
                temp.x = face.vertexIndices[i];
                temp.y = face.vertexIndices[j];
                edges.push_back(temp);
                actual_edges.insert({temp.x, temp.y});
                actual_edges.insert({temp.y, temp.x});
                HEdge* e1 = new HEdge; HEdge* e2 = new HEdge;
                e1->pair = e2;
                e2->pair = e1;
                e1->v = Vert_Pointer[temp.x];
                Vert_Pointer[temp.x]->h = e1;
                e2->v = Vert_Pointer[temp.y];
                Vert_Pointer[temp.y]->h = e2;
                Init[{temp.x, temp.y}] = e1;
                Init[{temp.y, temp.x}] = e2;
            }
        }
    }
    file.close();
}

void Process(int index, bool first, vector<bool> &Vis);

void Next(int index, vector<bool> &Vis){
    int i1 = faces[index].vertexIndices[0];
    int i2 = faces[index].vertexIndices[1];
    int i3 = faces[index].vertexIndices[2];
    vector<int> indices = {i1, i2, i3};
    for(int i = 0; i < 3; i++){
        int j = (i+1)%3;
        auto it1 = &Neigh[{indices[i], indices[j]}];
        for(int i = 0; i < (*it1).size(); i++){
            if(!Vis[(*it1)[i]]){
                Process((*it1)[i], false, Vis);
            }
        }
        auto it2 = &Neigh[{indices[j], indices[i]}];
        for(int i = 0; i < (*it2).size(); i++){
            if(!Vis[(*it2)[i]]){
                Process((*it2)[i], false, Vis);
            }
        }
    }
}

void Process(int index, bool first, vector<bool> &Vis){
    Face *face = new Face;
    Vis[index] = true;
    face->index = index;
    int i1 = faces[index].vertexIndices[0];
    int i2 = faces[index].vertexIndices[1];
    int i3 = faces[index].vertexIndices[2];
    Vert_Pointer[i1]->h = Init[{i1, i2}];
    Vert_Pointer[i2]->h = Init[{i2, i3}];
    Vert_Pointer[i3]->h = Init[{i3, i1}];
    vec3 vert1 = vertices[i1];
    vec3 vert3 = vertices[i2];
    vec3 vert2 = vertices[i3];
    if(first){
        vec3 temp1 = vert2-vert1;
        vec3 temp2 = vert3-vert2;
        vec3 check = cross(temp1, temp2);
        if(check.z >= 0){
            face->h = Init[{i1, i2}];
            Init[{i1, i2}]->f = face;
            Init[{i2, i3}]->f = face;
            Init[{i3, i1}]->f = face;
            Init[{i1, i2}]->next = Init[{i3, i1}];
            Init[{i2, i3}]->next = Init[{i1, i2}];
            Init[{i3, i1}]->next = Init[{i2, i3}];
        }
        else{
            face->h = Init[{i2, i1}];
            Init[{i2, i1}]->f = face;
            Init[{i1, i3}]->f = face;
            Init[{i3, i2}]->f = face;
            Init[{i2, i1}]->next = Init[{i3, i2}];
            Init[{i1, i3}]->next = Init[{i2, i1}];
            Init[{i3, i2}]->next = Init[{i1, i3}];
        }
        Next(index, Vis);
        return;
    }

    vector<int> indices = {i1, i2, i3};
    bool clock = false;
    if(Init[{i2, i1}]->f != nullptr && Init[{i1, i2}]->f == nullptr){clock = true;}
    if(Init[{i3, i2}]->f != nullptr && Init[{i2, i3}]->f == nullptr){clock = true;}
    if(Init[{i1, i3}]->f != nullptr && Init[{i3, i1}]->f == nullptr){clock = true;}
    if(clock){
        face->h = Init[{i1, i2}];
        Init[{i1, i2}]->f = face;
        Init[{i2, i3}]->f = face;
        Init[{i3, i1}]->f = face;
        Init[{i1, i2}]->next = Init[{i3, i1}];
        Init[{i2, i3}]->next = Init[{i1, i2}];
        Init[{i3, i1}]->next = Init[{i2, i3}];
    }
    else{
        face->h = Init[{i2, i1}];
        Init[{i2, i1}]->f = face;
        Init[{i1, i3}]->f = face;
        Init[{i3, i2}]->f = face;
        Init[{i2, i1}]->next = Init[{i3, i2}];
        Init[{i1, i3}]->next = Init[{i2, i1}];
        Init[{i3, i2}]->next = Init[{i1, i3}];
    }
    Next(index, Vis);
}

vec3 compute_normal(Vertex *vert){
    vec3 nor = vec3(0);
    vec3 ver = vertices[vert->index];
    assert(vert -> h != nullptr);

    HEdge *edge = vert->h;
    int temp = 0;
    do{
        temp++;
        if(edge->next == nullptr)break;
        if(edge->f == nullptr)break;
        vec3 vert1 = vertices[edge->pair->v->index];
        vec3 vert2 = vertices[edge->next->v->index];
        nor += cross(vert1-ver, vert2-ver)/((length(vert1-ver)*length(vert1-ver))*(length(vert2-ver)*length(vert2-ver)));
        edge = edge->next->pair;
    }while(edge != nullptr && edge != vert->h);

    nor /= length(nor);
    return nor;
}

void Check_Boundary(){
    for(int i = 0; i < vertices.size(); i++){
        int temp1 = 0, temp2 = 0;
        HEdge *h1;
        HEdge *h2;
        for(int j = 0; j < vertices.size(); j++){
            if(Init[{i,j}] != nullptr && Init[{i,j}]->f == nullptr){temp1++;h1 = Init[{i,j}];}
        }
        for(int j = 0; j < vertices.size(); j++){
            if(Init[{j,i}] != nullptr && Init[{j,i}]->f == nullptr){temp2++;h2 = Init[{j,i}];}
        }
        if(temp1 != 0 && temp2 != 0) {
            (h1->pair)->pair = h2->pair;
            Vert_Pointer[i]->h = h2->pair;
        }

        if(temp1 > 1 || temp2 > 1) {
            cout << "This node has weird half-edges : " << vertices[i].x << ' ' 
            << vertices[i].y << ' ' << vertices[i].z << endl;
        }
    }
}

long double dis(int index, vec3 point){
    vec3 vert1 = vertices[faces[index].vertexIndices[0]];
    vec3 vert2 = vertices[faces[index].vertexIndices[1]];
    vec3 vert3 = vertices[faces[index].vertexIndices[2]];
    vec3 normal = cross(vert2-vert1, vert3-vert1);
    normal /= length(normal);
    return abs(dot(point-vert1, normal));
}

vector<int> nearest_faces(vec3 point){
    long double min = 99999.0;
    vector<int> ans;
    map<int, bool> vis;
    for(int i = 0; i < vertices.size(); i++){
        Vertex *vert = Vert_Pointer[i];
        HEdge *h = vert->h;
        do{
            if(h->next == nullptr || h->f == nullptr){break;}
            int index = h->f->index;
            long double distance = dis(index, point);
            if(abs(distance-min) < 1e-6){
                if(!vis[index])ans.push_back(index);
                vis[index] = true;
            }
            else if(distance < min - 1e-6){
                min = distance;
                ans.clear();
                vis.clear();
                ans.push_back(index);
                vis[index] = true;
            }
            h = h->next->pair;
        }
        while(h != nullptr && h != vert->h);
    }
    return ans;
}

void Smoothen(float lambda, int iter){
    smoothening = true;

    vector<vec3> disp;
    while(iter){
        for(int i = 0; i < vertices.size(); i++){
            vec3 temp = vec3(0);
            Vertex *vert = Vert_Pointer[i];
            vec3 curr = vertices[i];
            HEdge *edge = vert->h;
            int num = 0;
            do{
                if(edge->next == nullptr)break;
                if(edge->f == nullptr)break;
                vec3 vert1 = vertices[edge->pair->v->index];
                temp += (vert1-curr);
                edge = edge->next->pair;
                num++;
            }while(edge != vert->h);
            disp.push_back(temp/float(num));
        }
        for(int i = 0; i < vertices.size(); i++){
            vertices[i] += lambda*disp[i];
        }
        iter--;
    }
}

int create_vertex(vec3 position) {
    Vertex* ver = new Vertex;
    ver->index = (int)vertices.size();
    Vert_Pointer.push_back(ver);
    vertices.push_back(position);
    return ver->index;
}

void create_edge(int id1, int id2, bool show) {
    //always both side edges are created
    vec2 ed;
    ed.x = id1, ed.y = id2;

    if(Init.find({id1, id2}) == Init.end()) {
        Init[{id1, id2}] = new HEdge;
        Init[{id2, id1}] = new HEdge;

        Init[{id1, id2}] -> pair = Init[{id2, id1}];
        Init[{id2, id1}] -> pair = Init[{id1, id2}];
    
        Init[{id1, id2}] -> v = Vert_Pointer[id1];
        Init[{id2, id1}] -> v = Vert_Pointer[id2];

        Vert_Pointer[id1] -> h = Init[{id1, id2}];
        Vert_Pointer[id2] -> h = Init[{id2, id1}];
    
        if(show) {
            edges.push_back(ed); 
            actual_edges.insert({ed.x, ed.y});
            actual_edges.insert({ed.y, ed.x});
        } 
    }
}

vector<int> get_vertex_orientation(vector<int> vertex_ids, vec3 point) {
    vec3 v0 = vertices[vertex_ids[0]];
    vec3 v1 = vertices[vertex_ids[1]];
    vec3 v2 = vertices[vertex_ids[2]];

    vec3 normal = cross(v1 - v0, v2 - v0);
    normal /= length(normal);

    if(dot(point - v0, normal) > 0) {
        return {vertex_ids[0], vertex_ids[1], vertex_ids[2]};
    }
    else {
        return {vertex_ids[0], vertex_ids[2], vertex_ids[1]};
    }
}

void create_face(vector<int> vertex_ids) {
    //give ids in the order you want to create the face in
    int id0 = vertex_ids[0];
    int id1 = vertex_ids[1];
    int id2 = vertex_ids[2];

    //id0 ke peeche id1 ke peeche id2
    Init[{id0, id2}] -> next = Init[{id1, id0}];
    Init[{id1, id0}] -> next = Init[{id2, id1}];
    Init[{id2, id1}] -> next = Init[{id0, id2}];

    Face* f = new Face;
    Init[{id0, id2}] -> f = f;
    Init[{id1, id0}] -> f = f;
    Init[{id2, id1}] -> f = f;

    Polygon fac;
    fac.vertexIndices = vertex_ids;
    faces.push_back(fac);

    f -> index = Face_Pointer.size();
    f -> h = Init[{id0, id2}];
    Face_Pointer.push_back(f);
}


vector<HEdge*> get_half_edges_of_face(int face_id) {
    Polygon fc = faces[face_id];
    if((Init[{fc.vertexIndices[0], fc.vertexIndices[1]}] -> f) -> index == face_id) {
        return {Init[{fc.vertexIndices[0], fc.vertexIndices[1]}], Init[{fc.vertexIndices[1], fc.vertexIndices[2]}], Init[{fc.vertexIndices[2], fc.vertexIndices[0]}]};
    }
    else {
        return {Init[{fc.vertexIndices[1], fc.vertexIndices[0]}], Init[{fc.vertexIndices[0], fc.vertexIndices[2]}], Init[{fc.vertexIndices[2], fc.vertexIndices[1]}]};
    }
}

vector<vector<HEdge*>> get_boundary_for_extrude(vector<int> face_ids) {
    set<int> all_face_ids;
    for(int i = 0; i < (int)face_ids.size(); i ++) {
        all_face_ids.insert(face_ids[i]);
    }

    Vertex* start_vertex = nullptr;
    map<Vertex*, HEdge*> all_edges;
    map<Vertex*, bool> visited;
    vector<Vertex*> all_vertices;

    for(int i = 0; i < (int)face_ids.size(); i ++) {
        vector<HEdge*> hs = get_half_edges_of_face(face_ids[i]);
        for(int j = 0; j < 3; j ++) {
            if(((hs[j] -> pair) -> f == nullptr) 
            || (all_face_ids.find(((hs[j] -> pair) -> f) -> index) == all_face_ids.end())) {
                all_edges[hs[j] -> v] = hs[j];
                visited[hs[j] -> v] = false;
                all_vertices.push_back(hs[j] -> v);
                start_vertex = hs[j] -> v;
            }
        }
    }

    vector<vector<HEdge*>> all_boundary_edges;

    while(start_vertex != nullptr) {
        vector<HEdge*> temp;
        temp.push_back(all_edges[start_vertex]);
        visited[start_vertex] = true;

        Vertex* cur_vertex;
        cur_vertex = (all_edges[start_vertex] -> pair) -> v;

        while(cur_vertex != start_vertex) {
            temp.push_back(all_edges[cur_vertex]);
            visited[cur_vertex] = true;
            cur_vertex = (all_edges[cur_vertex] -> pair) -> v;
        }

        all_boundary_edges.push_back(temp);
        start_vertex = nullptr;
        for(int i = 0; i <(int)all_vertices.size(); i ++) {
            if(! visited[all_vertices[i]]) start_vertex = all_vertices[i];
        }
    }

    return all_boundary_edges;
}

void extrude_boundary(vector<Vertex*> original_vertices, vector<Vertex*> new_vertices) {

    //this function just assumes the new vertices are created, nothing else
    int i = 0;
    int s = (int)original_vertices.size();
    assert(s == (int)new_vertices.size());

    while(i < (int)original_vertices.size()) {
        create_edge(new_vertices[i] -> index, original_vertices[i] -> index, true);
        create_edge(original_vertices[i] -> index, new_vertices[(i + 1) % s] -> index, false);
        i ++;
    }

    i = 0;
    while(i < (int)original_vertices.size()) {
        create_face({original_vertices[i] -> index, new_vertices[i] -> index, 
            new_vertices[(i + 1) % s] -> index});
        create_face({original_vertices[i] -> index, new_vertices[(i + 1) % s] -> index, 
            original_vertices[(i + 1) % s] -> index});
        i ++;
    }
}

void extrude_helper(vector<int> face_ids, vec3 displacement, vec3 point) {

    //creating the new vertices first
    map<int,int> extended_vertices;
    for(int i = 0; i < (int)face_ids.size(); i ++) {
        vec3 pt1 = vertices[faces[face_ids[i]].vertexIndices[0]];
        vec3 pt2 = vertices[faces[face_ids[i]].vertexIndices[1]];
        vec3 pt3 = vertices[faces[face_ids[i]].vertexIndices[2]];

        vec3 normal = cross(pt2 - pt1, pt3 - pt1);
        normal = normal / length(normal);
        vec3 disp = displacement;

        if(dot((point - pt1), normal) * dot(disp, normal) <= 0) {
            disp = - disp;
        }

        for(int j = 0; j < 3; j ++) {
            int id = faces[face_ids[i]].vertexIndices[j];
            vec3 vv = vertices[id];
            if(extended_vertices.find(id) == extended_vertices.end()) {
                extended_vertices[id] = create_vertex(vv + disp);
            }
        }
    }

    vector<vector<HEdge*>> boundary_vectors = get_boundary_for_extrude(face_ids);

    for(int idx = 0; idx < (int)boundary_vectors.size(); idx ++) {
        vector<HEdge*> boundary = boundary_vectors[idx];
        for(int i = 0; i < (int)boundary.size(); i ++) {
            int id1 = boundary[i] -> v -> index;
            int id2 = (boundary[(i + 1) % (int)boundary.size()] -> v) -> index;
            create_edge(extended_vertices[id1], extended_vertices[id2], true);
        }
    }


    //extruding the faces now
    for(int i = 0; i < (int)face_ids.size(); i ++) {
        vector<HEdge*> h = get_half_edges_of_face(face_ids[i]);
        assert((int)h.size() == 3);
        vector<int> old_v(3);
        vector<int> new_v(3);
        for(int j = 0; j < 3; j ++) {
            assert(extended_vertices.find((h[j] -> v) -> index) != extended_vertices.end());
            old_v[j] = (h[j] -> v) -> index;
            new_v[j] = extended_vertices[(h[j] -> v) -> index];
        }

        for(int j = 0; j < 3; j ++) {
            if(actual_edges.find({old_v[j], old_v[(j + 1) % 3]}) != actual_edges.end()
            || actual_edges.find({old_v[(j + 1) % 3], old_v[j]}) != actual_edges.end()) {
                create_edge(new_v[j], new_v[(j + 1) % 3], true);
            }
            else create_edge(new_v[j], new_v[(j + 1) % 3], false);
        }

        create_face({new_v[0], new_v[2], new_v[1]}); 
        // this order is important remember that
    }


    //extruding the boundary now
    for(int idx = 0; idx < (int)boundary_vectors.size(); idx ++) {
        vector<HEdge*> boundary = boundary_vectors[idx];
        vector<Vertex*> original_vertices;
        vector<Vertex*> new_vertices;
        for(int i = 0; i < (int)boundary.size(); i ++) {
            original_vertices.push_back(boundary[i] -> v);
        }
        for(int i = 0; i < (int)boundary.size(); i ++) {
            new_vertices.push_back(
                Vert_Pointer[
                    extended_vertices[boundary[i] -> v -> index]
                ]
            );
        }

        extrude_boundary(original_vertices, new_vertices);
    }
}

vec3 average_face_normal(vector<int> face_ids) {
    vec3 normal = vec3(0);
    for(int i = 0; i < (int)face_ids.size(); i ++) {
        vector<HEdge*> hs = get_half_edges_of_face(face_ids[i]);
        vec3 v1 = vertices[(hs[1] -> v) -> index] - vertices[(hs[0] -> v) -> index];
        vec3 v2 = vertices[(hs[2] -> v) -> index] - vertices[(hs[0] -> v) -> index];
        vec3 normal_area = cross(v1, v2);
        normal_area = normal_area / length(normal_area);
        normal += normal_area;
    }

    normal = normal / length(normal);
    return normal;
}

void extrude(vector<int> faces, double mag, vec3 point) {
    extrusion = true;
    vec3 disp = vec3(mag) * average_face_normal(faces);
    extrude_helper(faces, disp, point);
}

vector<int> get_all_face_ids(vector<int> vertex_ids) {
    set<int> vids;
    vector<int> all_face_ids;
    for(int i = 0; i < (int)vertex_ids.size(); i ++) {
        vids.insert(vertex_ids[i]);
    }

    for(int i = 0; i < (int)faces.size(); i ++) {
        if(vids.find(faces[i].vertexIndices[0]) != vids.end() 
        && (vids.find(faces[i].vertexIndices[1]) != vids.end())
        && (vids.find(faces[i].vertexIndices[2]) != vids.end())) {
            all_face_ids.push_back(i);
        }
    }

    return all_face_ids;
}

// void catmull_clark() {
    
// }

void show(){
    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        cout << "FAILURE" << endl;
        return ;
    }
    vec3 vertices_new[(int)Vert_Pointer.size()];
    vec3 normals_new[(int)Vert_Pointer.size()];
    vector<ivec3> triangles;
    ivec2 edges_new[(int)edges.size()];
    for(int i = 0; i < edges.size(); i++){
        edges_new[i] = edges[i];
    }
    for(int i = 0; i < vertices.size(); i++){
        vertices_new[i] = vertices[i];
    }

    if((! smoothening) && (! extrusion) && ((int)normals.size() == (int)vertices.size())) {
        for(int i = 0; i < (int)vertices.size(); i ++) {
            normals_new[i] = normals[i];
        }
    }
    else {
        for(int i = 0; i < (int)vertices.size(); i ++) {
            normals_new[i] = compute_normal(Vert_Pointer[i]);
        }
    }

    for(int i = 0; i < (int)vertices.size(); i ++) {
        HEdge *h = Vert_Pointer[i]->h;
        do {
            ivec3 temp = ivec3(0);
            temp.x = i;
            if(h->next == nullptr){break;}
            temp.y = h->pair->v->index;
            temp.z = h->next->v->index;
            triangles.push_back(temp);
            h = h->next->pair;
        }   while(h != nullptr && h != Vert_Pointer[i]->h);
    }


    ivec3 faces_new[(int)triangles.size()];
    for(int i = 0; i < triangles.size(); i++){
        faces_new[i] = triangles[i];
    }
    v.setMesh(vertices.size(), triangles.size(), edges.size(), vertices_new, faces_new, edges_new, normals_new);
    v.view();
}

int main(int argc, char* argv[]) {

    string filename = argv[1];
    int smooth = atoi(argv[2]);
    int iter = atoi(argv[3]);
    int ext = atoi(argv[4]);

    double x = atof(argv[5]);
    double y = atof(argv[6]);
    double z = atof(argv[7]);
    double scaling = atof(argv[8]);

    loadOBJ(filename);

    // Print loaded data
    std::cout << "Loaded " << vertices.size() << " vertices, " 
              << texCoords.size() << " texture coordinates, " 
              << normals.size() << " normals, " 
              << faces.size() << " faces." << std::endl;

    vector<bool> Face_Visited((int)faces.size());
    Process(0, true, Face_Visited);
    Check_Boundary();

    if(ext == 1) {
        //vector<int> extrude_faces = nearest_faces(vec3(x, y, z));
        // cout << "Number of nearest faces is : " << (int)extrude_faces.size() << endl;
        // extrude(extrude_faces, scaling * dis(extrude_faces[0], vec3(x, y, z)), vec3(x, y, z));


        //get the vector of planes for extrude

        //extruding the center faces of a 3 x 3 cube
        // vector<int> all_face_ids = get_all_face_ids({5, 6, 9, 10});
        // extrude(all_face_ids, 0.3, vec3(-0.6, 0.0, 0.0));
        // all_face_ids = get_all_face_ids({17, 18, 29, 30});
        // extrude(all_face_ids, 0.3, vec3(0.0, -0.6, 0.0));
        // all_face_ids = get_all_face_ids({20, 22, 32, 34});
        // extrude(all_face_ids, 0.3, vec3(0.0, 0.0, -0.6));
        // all_face_ids = get_all_face_ids({21, 23, 33, 35});
        // extrude(all_face_ids, 0.3, vec3(0.0, 0.0, 0.6));
        // all_face_ids = get_all_face_ids({25, 26, 37, 38});
        // extrude(all_face_ids, 0.3, vec3(0.0, 0.6, 0.0));
        // all_face_ids = get_all_face_ids({45, 46, 49, 50});
        // extrude(all_face_ids, 0.3, vec3(0.6, 0.0, 0.0));

        //making a pot
        // vector<int> vc;
        // for(int i = 0; i < 40; i ++) vc.push_back(i);
        // vector<int> all_face_ids = get_all_face_ids(vc);
        // cout << "Number of faces is : " << (int)all_face_ids.size() << endl;
        // extrude(all_face_ids, scaling, vec3(x, y, z));

        //extruding a connected region of a sphere
        // vector<int> vc;
        // for(int i = 81; i <= 120; i ++) vc.push_back(i);
        // vector<int> all_face_ids = get_all_face_ids(vc);
        // cout << "Number of faces is : " << (int)all_face_ids.size() << endl;
        // extrude(all_face_ids, scaling, vec3(x, y, z));
    }

    if(smooth == 1) {
        Smoothen(0.01, iter);
    }

    show(); //show the object

    return 0;
}
