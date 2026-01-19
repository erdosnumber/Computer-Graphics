#include "camera.hpp"

#include <iostream>
#include<vector>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;
using namespace std;

GL::Rasterizer r;
GL::ShaderProgram program;

const int L = 15; // L is the number of vertices

glm::vec3 vel[L * L];
bool fix[L * L];

const float structural_constant = 1000.0;
const float shear_constant = 700.0;
const float bending_constant = 200.0;
const float damping_structural_constant = 100.0;
const float damping_shear_constant = 100.0;
const float damping_bending_constant = 100.0;
const float gravity = 0.98;
const float mass = 1.0;

CameraControl camCtl;

class Object {
	private:
	GL::Object object;
	GL::AttribBuf vertexBuf, normalBuf;
	vector<vec3> vertices;
	vector<vec3> normals;
	vector<ivec3> triangles;
	vector<ivec2> edges;

	public:
	Object(int nv, int ne, int nt) {
		vertices.resize(nv);
		normals.resize(nv);
		edges.resize(ne);
		triangles.resize(nt);
		object = r.createObject();
	};

	void add_vertices(vec3 v) {
		vertices.push_back(v);
	}
	void add_normals(vec3 v) {
		normals.push_back(v);
	}
	void add_edges(ivec2 v) {
		edges.push_back(v);
	}
	void add_triangles(ivec3 v) {
		triangles.push_back(v);
	}

	void set_vertex(int id, vec3 v) {
		vertices[id] = v;
	}
	void set_normal(int id, vec3 v) {
		normals[id] = v;
	}
	void set_edge(int id, ivec2 v) {
		edges[id] = v;
	}
	void set_triangle(int id, ivec3 v) {
		triangles[id] = v;
	}

	vec3 get_vertex(int id) {
		return vertices[id];
	}

	vec3 get_normal(int id) {
		return normals[id];
	}

	ivec2 get_edge(int id) {
		return edges[id];
	}

	ivec3 get_triangle(int id) {
		return triangles[id];
	}

	GL::Object& get_object() {
		return object;
	}

	void build_vertices() {
		vertexBuf = r.createVertexAttribs(object, 0, (int)vertices.size(), &vertices[0]);
	}
	void build_normals() {
		normalBuf = r.createVertexAttribs(object, 1, (int)normals.size(), &normals[0]);
	}
	void build_triangles() {
		r.createTriangleIndices(object, (int)triangles.size(), &triangles[0]);
	}
	void build_edges() {
		r.createEdgeIndices(object, (int)edges.size(), &edges[0]);
	}

	void update_vertices() {
		r.updateVertexAttribs(vertexBuf, (int)vertices.size(), &vertices[0]);
	}

	void update_normals() {
		r.updateVertexAttribs(normalBuf, (int)normals.size(), &normals[0]);
	}
};

Object* cloth;

void initializeScene() {
    //need to create the edges
	cloth = new Object(L * L, 2 * L * (L - 1), 2 * (L - 1) * (L - 1));

    //need to create a L x L sized square
    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            cloth -> set_vertex(L * i + j, vec3((float)(j * 1.0), 0.0, (float)(i * 1.0)));
        }
    }

    //need to create the normals 
    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            cloth -> set_normal(L * i + j, vec3(0.0, -1.0, 0.0));
        }
    }

    //need to create the edges
    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L - 1; j ++) {
            cloth -> add_edges(ivec2(L * i + j, L * i + (j + 1)));
        }
    }
    for(int j = 0; j < L; j ++) {
        for(int i = 0; i < L - 1; i ++) {
            cloth -> add_edges(ivec2(L * i + j, L * (i + 1) + j));
        }
    }

    for(int i = 0; i < L - 1; i ++) {
        for(int j = 0; j < L - 1; j ++) {
            cloth -> add_triangles(ivec3(L * i + j, L * i + (j + 1), L * (i + 1) + j));
        }
    }

    for(int i = 1; i < L; i ++) {
        for(int j = 1; j < L; j ++) {
            cloth -> add_triangles(ivec3(L * i + j, L * i + (j - 1), L * (i - 1) + j));
        }
    }

    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            fix[L * i + j] = false;
        }
    }

    fix[L * (L - 1)] = true;
    fix[L * (L - 1) + (L - 1)] = true;

    cloth -> build_vertices();
    cloth -> build_normals();
    cloth -> build_edges();
    cloth -> build_triangles();
}

vec3 get_force(int src, int dest, float nat_length, float spring_constant, float damping_constant) {
    //rest length is 1 for all springs
    vec3 dir = cloth -> get_vertex(dest) - cloth -> get_vertex(src);
    float lt = length(dir);
    dir /= length(dir);
    vec3 f1 = spring_constant * (lt - nat_length) * dir;


    vec3 vdest = vel[dest];
    vec3 vsrc = vel[src];
    vec3 vrel = dot((vdest - vsrc), dir) * dir;
    vec3 f2 = damping_constant * vrel;

    vec3 f = -f1 - f2;
    return f;
}

vec3 get_all_forces(int i, int j) {
    //get all structural forces
    vec3 f = vec3(0.0, 0.0, 0.0);
    if(i > 0) {
        f += get_force(L * (i - 1) + j, L * i + j, 1.0, structural_constant, damping_structural_constant);
    }
    if(i < L - 1) {
        f += get_force(L * (i + 1) + j, L * i + j, 1.0, structural_constant, damping_structural_constant);
    }
    if(j > 0) {
        f += get_force(L * i + (j - 1), L * i + j, 1.0, structural_constant, damping_structural_constant);
    }
    if(j < L - 1) {
        f += get_force(L * i + (j + 1), L * i + j, 1.0, structural_constant, damping_structural_constant);
    }

    //get all shear forces
    if(i > 0 && j > 0) {
        f += get_force(L * (i - 1) + (j - 1), L * i + j, sqrt(2.0), shear_constant, damping_shear_constant);
    }
    if(i > 0 && j < L - 1) {
        f += get_force(L * (i - 1) + (j + 1), L * i + j, sqrt(2.0), shear_constant, damping_shear_constant);
    }
    if(i < L - 1 && j > 0) {
        f += get_force(L * (i + 1) + (j - 1), L * i + j, sqrt(2.0), shear_constant, damping_shear_constant);
    }
    if(i < L - 1 && j < L - 1) {
        f += get_force(L * (i + 1) + (j + 1), L * i + j, sqrt(2.0), shear_constant, damping_shear_constant);
    }

    //get all bending forces
    if(i > 1 && j > 0) {
        f += get_force(L * (i - 2) + (j - 1), L * i + j, sqrt(5.0), bending_constant, damping_bending_constant);
    }
    if(i > 1 && j < L - 1) {
        f += get_force(L * (i - 2) + (j + 1), L * i + j, sqrt(5.0), bending_constant, damping_bending_constant);
    }
    if(i < L - 2 && j > 0) {
        f += get_force(L * (i + 2) + (j - 1), L * i + j, sqrt(5.0), bending_constant, damping_bending_constant);
    }
    if(i < L - 2 && j < L - 1) {
        f += get_force(L * (i + 2) + (j + 1), L * i + j, sqrt(5.0), bending_constant, damping_bending_constant);
    }

    //finally the gravity force, acting in -Y direction
    f += vec3(0.0, -gravity, 0.0);

    return f;
}

vec3 compute_normal(int i, int j) {
    vector<vec3> edges;
    if(i > 0) {
        edges.push_back(cloth -> get_vertex(L * (i - 1) + j) - cloth -> get_vertex(L * i + j));
    }
    if(j < L - 1) {
        edges.push_back(cloth -> get_vertex(L * i + (j + 1)) - cloth -> get_vertex(L * i + j));
    }
    if(i < L - 1) {
        edges.push_back(cloth -> get_vertex(L * (i + 1) + j) - cloth -> get_vertex(L * i + j));
    }
    if(j > 0) {
        edges.push_back(cloth -> get_vertex(L * i + (j - 1)) - cloth -> get_vertex(L * i + j));
    }

    vec3 norm = vec3(0.0, 0.0, 0.0);

    for(int i = 0; i < (int)edges.size(); i ++) {
        vec3 v1 = edges[i];
        vec3 v2 = edges[(i + 1) % (int)edges.size()];
        norm += 
        cross(v1, v2) / (length(v1) * length(v1) * length(v2) * length(v2));
    }

    return norm;
}

void updateScene(float delt) {
    vector<vec3> new_vertices(L * L);
    vector<vec3> new_vel(L * L);

    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            if(fix[L * i + j] == false) {
                vec3 f = get_all_forces(i, j);
                new_vel[L * i + j] = vel[L * i + j] + (1 / mass) * f * delt;
                new_vertices[L * i + j] = cloth -> get_vertex(L * i + j) + new_vel[L * i + j] * delt;
            }
            else {
                new_vertices[L * i + j] = cloth -> get_vertex(L * i + j);
                new_vel[L * i + j] = vel[L * i + j];
            }
        }
    }

    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            cloth -> set_vertex(L * i + j, new_vertices[L * i + j]);
            vel[L * i + j] = new_vel[L * i + j];
        }
    }

    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            cloth -> set_normal(L * i + j, compute_normal(i, j));
        }
    }

    cloth -> update_vertices();
    cloth -> update_normals();
}

int main(int argc, char *argv[]) {
	int width = 640, height = 480;
	if (!r.initialize("Animation", width, height)) {
		return EXIT_FAILURE;
	}
	camCtl.initialize(width, height);
	camCtl.camera.setCameraView(vec3(L / 2, - 3 * L / 4 , 5 * L / 2), vec3(L / 2, -3 * L / 4, 0.0), vec3(0.0, 1.0, 0.0));
	//camCtl.camera.setCameraView(vec3(0.5, -0.5, 1.5), vec3(0.5, -0.5, 0.0), vec3(0.5, 0.0, 1.5));
	program = r.createShaderProgram(
		r.vsBlinnPhong(),
		r.fsBlinnPhong()
	);

	initializeScene();

	while (!r.shouldQuit()) {
        float t = SDL_GetTicks64()*1e-2;
		updateScene(0.001);

		camCtl.update();
		Camera &camera = camCtl.camera;

		r.clear(vec4(0.4, 0.4, 0.4, 1.0));
		r.enableDepthTest();
		r.useShaderProgram(program);

		r.setUniform(program, "model", glm::mat4(1.0));
		r.setUniform(program, "view", camera.getViewMatrix());
		r.setUniform(program, "projection", camera.getProjectionMatrix());
		r.setUniform(program, "lightPos", camera.position);
		r.setUniform(program, "viewPos", camera.position);
		r.setUniform(program, "lightColor", vec3(1.0f, 1.0f, 1.0f));

		r.setupFilledFaces();
        glm::vec3 orange(1.0f, 0.6f, 0.2f);
        glm::vec3 white(1.0f, 1.0f, 1.0f);
		glm::vec3 blue(0.2f, 0.2f, 1.0f);
        r.setUniform(program, "ambientColor", 0.2f*white);
        r.setUniform(program, "extdiffuseColor", 0.9f*orange);
        r.setUniform(program, "intdiffuseColor", blue);
        r.setUniform(program, "specularColor", 0.6f*white);
        r.setUniform(program, "phongExponent", 20.f);
		r.drawTriangles(cloth -> get_object());

		r.setupWireFrame();
        glm::vec3 black(0.0f, 0.0f, 0.0f);
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(cloth -> get_object());

		r.show();
	}
}

//command to execute
//cmake -B build && cmake --build build && ./build/part2 1000 50 1 50 2.5 0.05
