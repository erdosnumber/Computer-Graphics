#include "camera.hpp"

#include <iostream>
#include<vector>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;
using namespace std;

GL::Rasterizer r;
GL::ShaderProgram program;

const int L = 30; // L is the number of vertices
glm::vec3 vel[L * L];
glm::vec3 force[L * L];
bool fix[L * L];

const float structural_constant = 1000.0;
const float shear_constant = 50.0;
const float bending_constant = 1.0;
const float damping_structural_constant = 50.0;
const float damping_shear_constant = 2.5;
const float damping_bending_constant = 0.05;
const float gravity = 9.8;
const float mass = 1.0;

const int cv = L * L;
const int ce =  2 * L * (L - 1);
const int ct =  2 * (L - 1) * (L - 1);

// const float radius = 7.0;
const int v_div = 30;
const int h_div = 30;
vec3 axis = vec3(0.0, 1.0, 0.0);

const int sv = (h_div - 1) * v_div + 2;
const int se = h_div * v_div;
const int st = 2 * (h_div - 2) * v_div + 2 * v_div;

const float base = -20.0;
const float plane_edge = 25.0;

CameraControl camCtl;
const float DEL = 0.01;
const float EPS = 1.0;
const float MU = 0.0;

class Object {
	private:
	GL::Object object;
	GL::AttribBuf vertexBuf, normalBuf;
	vector<vec3> vertices;
	vector<vec3> normals;
	vector<ivec3> triangles;
	vector<ivec2> edges;
    vec3 center;
    float radius;
    vec3 speed;
    vec3 angular_velocity;

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
    void set_center(vec3 v) {
        center = v;
    }
    void set_radius(float r) {
        radius = r;
    }
    void set_speed(vec3 v) {
        speed = v;
    }
    void set_angular_velocity(vec3 v) {
        angular_velocity = v;
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
    vec3 get_center() {
        return center;
    }
    float get_radius() {
        return radius;
    }
    vec3 get_speed() {
        return speed;
    }
    vec3 get_angular_velocity() {
        return angular_velocity;
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
Object* moving_sphere;
Object* static_sphere;
Object* plane;

void create_sphere(Object* object) {
    object -> set_vertex(0, object -> get_center() + object -> get_radius() * vec3(0.0, 1.0, 0.0));
    for(int i = 1; i < h_div; i ++) {
        //north pole corresponds to i = 0
        //south pole corresponds to i = L
        for(int j = 1; j <= v_div; j ++) {
            float theta = (float)(j * (2 * M_PI / v_div));
            float phi = (float)(i * M_PI / h_div);
            vec3 offset = vec3(object -> get_radius() * sin(phi) * cos(theta), object -> get_radius() * cos(phi), object -> get_radius() * sin(phi) * sin(theta)); //axis is along the Y
            vec3 new_pos = object -> get_center() + offset;
            object -> set_vertex((i - 1) * v_div + j, new_pos);
        }
    }
    object -> set_vertex((h_div - 1) * v_div + 1, object -> get_center() - object -> get_radius() * vec3(0.0, 1.0, 0.0));

    for(int i = 0; i < h_div; i ++) {
        //edges along the same longitude
        if(i == 0) {
            for(int j = 1; j <= v_div; j ++) {
                object -> add_edges(ivec2(0, j));
            }
        }
        else if(i > 0 && i < h_div - 1) {
            for(int j = 1; j <= v_div; j ++) {
                object -> add_edges(ivec2((i - 1) * v_div + j, i * v_div + j));
            }
        }
        else {
            for(int j = 1; j <= v_div; j ++) {
                object -> add_edges(ivec2((i - 1) * v_div + j, (h_div - 1) * v_div + 1));
            }
        }

        if(i > 0) {
            //edges along the same latitude
            for(int j = 1; j <= v_div; j ++) {
                object -> add_edges(ivec2((i - 1) * v_div + j, (i - 1) * v_div + (j % v_div) + 1));
            }
        }
    }

    for(int i = 0; i < sv; i ++) {
        vec3 diff = object -> get_vertex(i) - object -> get_center();
        diff /= length(diff);
        object -> set_normal(i, diff);
    }

    for(int i = 0; i < h_div; i ++) {
        if(i == 0) {
            for(int j = 1; j <= v_div; j ++) {
                object -> add_triangles(ivec3(0, (j % v_div) + 1, j));
            }
        }
        else if(i > 0 && i < h_div - 1) {
            for(int j = 1; j <= v_div; j ++) {
                object -> add_triangles(ivec3((i - 1) * v_div + j, (i - 1) * v_div + (j % v_div) + 1, i * v_div + j));
                object -> add_triangles(ivec3((i - 1) * v_div + (j % v_div) + 1, i * v_div + (j % v_div) + 1, i * v_div + j));
            }
        }
        else {
            for(int j = 1; j <= v_div; j ++) {
                object -> add_triangles(ivec3((i - 1) * v_div + j, (i - 1) * v_div + (j % v_div) + 1, (h_div - 1) * v_div + 1));
            }
        }
    }

    object -> build_vertices();
    object -> build_normals();
    object -> build_edges();
    object -> build_triangles();
}

void initializeScene() {
    //need to create the edges
	cloth = new Object(cv, ce, ct);

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

    // fix[L * (L - 1)] = true;
    // fix[L * (L - 1) + (L - 1)] = true;

    cloth -> build_vertices();
    cloth -> build_normals();
    cloth -> build_edges();
    cloth -> build_triangles();

    //create a sphere
    //(11.0, -11.0, 21.0) contact point of sphere and cloth
    //keep center of sphere at (11.0, -11.0, 21.0 + R)
    //+ Y axis is up

    moving_sphere = new Object(sv, se, st);
    moving_sphere -> set_center(vec3(L / 2, -L, 5 * L / 4));
    moving_sphere -> set_radius(8.0);
    moving_sphere -> set_speed(vec3(0.0, 0.0, -8.0));
    moving_sphere -> set_angular_velocity(vec3(0.0, 20.0, 0.0));

    create_sphere(moving_sphere);

    static_sphere = new Object(sv, se, st);
    static_sphere -> set_center(vec3(L / 2, -L / 2 , 0.0));
    static_sphere -> set_radius(4.0);

    create_sphere(static_sphere);

    plane = new Object(4, 4, 2);
    plane -> set_vertex(0, vec3(1000.0, base, plane_edge));
    plane -> set_vertex(1, vec3(-1000.0, base, plane_edge));
    plane -> set_vertex(2, vec3(-1000.0, base, 1000.0));
    plane -> set_vertex(3, vec3(1000.0, base, 1000.0));

    plane -> set_normal(0, vec3(0.0, 1.0, 0.0));
    plane -> set_normal(1, vec3(0.0, 1.0, 0.0));
    plane -> set_normal(2, vec3(0.0, 1.0, 0.0));
    plane -> set_normal(3, vec3(0.0, 1.0, 0.0));

    plane -> set_edge(0, ivec2(0, 1));
    plane -> set_edge(1, ivec2(1, 2));
    plane -> set_edge(2, ivec2(2, 3));
    plane -> set_edge(3, ivec2(3, 0));

    plane -> set_triangle(0, ivec3(0, 1, 2));
    plane -> set_triangle(1, ivec3(0, 2, 3));
    
    plane -> build_vertices();
    plane -> build_normals();
    plane -> build_edges();
    plane -> build_triangles();
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

    vec3 f = - f1 - f2;
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

vec3 compute_cloth_normal(int i, int j) {
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

vec4 mult_q(vec4 q1, vec4 q2) {
	float a = q1.x;
	vec3 b = vec3(q1.y, q1.z, q1.w);
	float c = q2.x;
	vec3 d = vec3(q2.y, q2.z, q2.w);

	float e = a * c - dot(b, d);
	vec3 f = a * d + c * b + cross(b, d);
	return vec4(e, f.x, f.y, f.z);
}

vec3 compute_sphere_normal(Object* obj,int id) {
    vec3 diff = obj -> get_vertex(id) - obj -> get_center();
    diff /= length(diff);
    return diff;
}

void rotate_sphere(Object* obj, float angle, vec3 axis) {
    for(int i = 0; i < sv; i ++) {
        vec3 v = obj -> get_vertex(i) - obj -> get_center();
        float r = length(v);
        v /= length(v);

        vec4 q1 = vec4(cos(angle / 2), axis.x * sin(angle / 2),
                axis.y * sin(angle / 2), axis.z * sin(angle / 2));
        vec4 q2 = vec4(0.0, v.x, v.y, v.z);
        vec4 q1_inv = vec4(q1.x, -q1.y, -q1.z, -q1.w);

        vec4 resv = mult_q(mult_q(q1, q2), q1_inv);
        vec3 new_v = vec3(resv.y, resv.z, resv.w);
        new_v *= r;
        new_v += obj -> get_center();
        obj -> set_vertex(i, new_v);
    }
}

void collision_sphere_forces(Object* sph, int i, float tol) {
    vec3 c = sph -> get_center();
    vec3 v = cloth -> get_vertex(i);
    vec3 diff = v - c;
    float d = length(diff);
    if(d < sph -> get_radius() + tol) {
        vec3 norm = diff / d;
        vec3 new_v = c + (sph -> get_radius() + tol) * norm;

        //get the velocity of the cloth at this point
        vec3 cloth_vel = vel[i];
        vec3 sphere_vel = sph -> get_speed() + cross(sph -> get_angular_velocity(), new_v - c);

        vec3 rel_vel = cloth_vel - sphere_vel;
        vec3 nvel = dot(rel_vel, norm) * norm;
        vec3 tvel = rel_vel - dot(rel_vel, norm) * norm;

        if(dot(force[i], norm) < 0) {
            vec3 fn = - EPS * (dot(force[i], norm) * norm);
            vec3 ft = vec3(0.0);
            if(length(tvel) > 0) ft = - MU * abs(dot(force[i], norm)) * tvel / length(tvel);

            force[i] = force[i] + fn + ft;
        }
    }
}

void collision_plane_forces(int i, float tol, float base) {
    float d = length((cloth -> get_vertex(i)).y - base);
    if(d < tol && (cloth -> get_vertex(i)).z > plane_edge) {
        vec3 norm = vec3(0.0, 1.0, 0.0);
        vec3 rel_vel = vel[i];
        vec3 nvel = dot(rel_vel, norm) * norm;
        vec3 tvel = rel_vel - dot(rel_vel, norm) * norm;

        if(dot(force[i], norm) < 0) {
            vec3 fn = -EPS * (dot(force[i], norm) * norm);
            vec3 ft = vec3(0.0);
            if(length(tvel) > 0) ft = - MU * abs(dot(force[i], norm)) * tvel / length(tvel);
            force[i] = force[i] + fn + ft;
        }
    }
}

void collision_sphere_velocities(Object* sph, int i, float tol) {
    vec3 c = sph -> get_center();
    vec3 v = cloth -> get_vertex(i);
    vec3 diff = v - c;
    float d = length(diff);
    if(d < sph -> get_radius() + tol) {
        vec3 norm = diff / d;
        vec3 new_v = c + (sph -> get_radius() + tol) * norm;

        //get the velocity of the cloth at this point
        vec3 cloth_vel = vel[i];
        vec3 sphere_vel = sph -> get_speed() + cross(sph -> get_angular_velocity() * axis, new_v - c);

        vec3 rel_vel = cloth_vel - sphere_vel;
        vec3 nvel = dot(rel_vel, norm) * norm;
        vec3 tvel = rel_vel - dot(rel_vel, norm) * norm;

        if(dot(rel_vel, norm) < 0) {
            //veclocity toh is always changing
            vec3 new_nvel = -EPS * nvel;
            vec3 new_tval = vec3(0.0);
            if(length(tvel) > 0) new_tval = tvel - (1 + EPS) * MU * length(nvel) * (tvel / length(tvel));
            vec3 new_rel_vel = new_nvel + new_tval;

            vel[i] = new_rel_vel + sphere_vel;
        }
        //now we know the normal and the relative velocity of the
        //particle with respect to the sphere
    }
}

void collision_plane_velocities(int i, float tol, float base) {
    float d = length((cloth -> get_vertex(i)).y - base);
    if(d < tol && (cloth -> get_vertex(i)).z > plane_edge) {
        vec3 norm = vec3(0.0, 1.0, 0.0);
        vec3 rel_vel = vel[i];
        vec3 nvel = dot(rel_vel, norm) * norm;
        vec3 tvel = rel_vel - dot(rel_vel, norm) * norm;

        if(dot(rel_vel, norm) < 0) {
            vec3 new_nvel = -EPS * nvel;
            vec3 new_tval = vec3(0.0);
            if(length(tvel) > 0) new_tval = tvel - (1 + EPS) * MU * length(nvel) * (tvel / length(tvel));
            vec3 new_rel_vel = new_nvel + new_tval;

            vel[i] = new_rel_vel;
        }
    }
}

void updateScene(float delt) {

    vector<vec3> new_vertices(cv);
    vector<vec3> new_vel(cv);

    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            if(fix[L * i + j] == false) {
                vec3 f = get_all_forces(i, j);
                force[L * i + j] = f;
                collision_sphere_forces(moving_sphere, L * i + j, 0.5);
                collision_sphere_forces(static_sphere, L * i + j, 0.5);
                collision_plane_forces(L * i + j, 0.5, base);
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
            collision_sphere_velocities(moving_sphere, L * i + j, 0.5);
            collision_sphere_velocities(static_sphere, L * i + j, 0.5);
            collision_plane_velocities(L * i + j, 0.5, base);
        }
    }

    for(int i = 0; i < L; i ++) {
        for(int j = 0; j < L; j ++) {
            cloth -> set_normal(L * i + j, compute_cloth_normal(i, j));
        }
    }

    cloth -> update_vertices();
    cloth -> update_normals();

    for(int i = 0; i < sv; i ++) {
        moving_sphere -> set_vertex(i, moving_sphere -> get_vertex(i) + moving_sphere -> get_speed() * delt);
    }
    moving_sphere -> set_center(moving_sphere -> get_center() + moving_sphere -> get_speed() * delt);

    rotate_sphere(moving_sphere, length(moving_sphere -> get_angular_velocity()) * delt, vec3(0.0, 1.0, 0.0));

    for(int i = 0; i < sv; i ++) {
        moving_sphere -> set_normal(i, compute_sphere_normal(moving_sphere, i));
    }

    moving_sphere -> update_vertices();
    moving_sphere -> update_normals();
}

int main(int argc, char *argv[]) {
	int width = 640, height = 480;
	if (!r.initialize("Animation", width, height)) {
		return EXIT_FAILURE;
	}
	camCtl.initialize(width, height);
	camCtl.camera.setCameraView(vec3(2 * L, - L / 2, L / 2), vec3(0.0, - L / 2, L / 2), vec3(0.0, 1.0, 0.0));
	//camCtl.camera.setCameraView(vec3(0.5, -0.5, 1.5), vec3(0.5, -0.5, 0.0), vec3(0.5, 0.0, 1.5));
	program = r.createShaderProgram(
		r.vsBlinnPhong(),
		r.fsBlinnPhong()
	);

	initializeScene();

	while (!r.shouldQuit()) {
        float t = SDL_GetTicks64()*1e-2;
		updateScene(0.0005);

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
        glm::vec3 green(0.2f, 1.0f, 0.2f);
        glm::vec3 red(1.0f, 0.2f, 0.2f);
        glm::vec3 brown(0.6f, 0.4f, 0.2f);
        glm::vec3 black(0.0f, 0.0f, 0.0f);


        r.setUniform(program, "ambientColor", 0.2f*white);
        r.setUniform(program, "extdiffuseColor", 0.9f*orange);
        r.setUniform(program, "intdiffuseColor", blue);
        r.setUniform(program, "specularColor", 0.6f*white);
        r.setUniform(program, "phongExponent", 20.f);
		r.drawTriangles(cloth -> get_object());

		r.setupWireFrame();
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(cloth -> get_object());


        r.setUniform(program, "ambientColor", 0.2f*white);
        r.setUniform(program, "extdiffuseColor", 0.9f*red);
        r.setUniform(program, "intdiffuseColor", blue);
        r.setUniform(program, "specularColor", 0.6f*white);
        r.setUniform(program, "phongExponent", 20.f);
		r.drawTriangles(moving_sphere -> get_object());

		r.setupWireFrame();
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(moving_sphere -> get_object());

        r.setUniform(program, "ambientColor", 0.2f*white);
        r.setUniform(program, "extdiffuseColor", green);
        r.setUniform(program, "intdiffuseColor", 0.8f * blue);
        r.setUniform(program, "specularColor", 0.6f*white);
        r.setUniform(program, "phongExponent", 20.f);
		r.drawTriangles(static_sphere -> get_object());

		r.setupWireFrame();
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(static_sphere -> get_object());


        r.setUniform(program, "ambientColor", 0.2f*white);
        r.setUniform(program, "extdiffuseColor", brown);
        r.setUniform(program, "intdiffuseColor", 0.8f * blue);
        r.setUniform(program, "specularColor", 0.6f*white);
        r.setUniform(program, "phongExponent", 20.f);
		r.drawTriangles(plane -> get_object());

		r.setupWireFrame();
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(plane -> get_object());

        //display the object
		r.show();
	}
}

//command to execute
//cmake -B build && cmake --build build && ./build/part2 1000 50 1 50 2.5 0.05
