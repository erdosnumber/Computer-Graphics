#include "camera.hpp"

#include <iostream>
#include<vector>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;
using namespace std;

GL::Rasterizer r;
GL::ShaderProgram program;
const float thickness = 0.25;

const int nv = 56, nt = 28, ne = 56;

vec3 vertices[nv]
{
vec3(0.0, 2.0, 0.0), vec3(0.7, 2.7, 0.0), vec3(0.0, 3.4, 0.0), vec3(-0.7, 2.7, 0.0), //face vertices
vec3(1.0,2.0,0.0), vec3(-1.0, 2.0, 0.0), vec3(-1.0, -2.0, 0.0), vec3(1.0, -2.0, 0.0), //torso vertices
vec3(1.0, -2.0, 0.0), vec3(0.5, -2.0, 0.0), vec3(1.0, -3.0, 0.0), vec3(1.5, -3.0, 0.0), //upper left leg vertices
vec3(1.5, -3.0, 0.0), vec3(1.0, -3.0, 0.0), vec3(1.0, -5.0, 0.0), vec3(1.5, -5.0, 0.0), //lower left leg vertices
vec3(-0.5, -2.0, 0.0), vec3(-1.0, -2.0, 0.0), vec3(-1.5, -3.0, 0.0), vec3(-1, -3.0, 0.0), //upper right leg vertices
vec3(-1.0, -3.0, 0.0), vec3(-1.5, -3.0, 0.0), vec3(-1.5, -5.0, 0.0), vec3(-1.0, -5.0, 0.0), //lower right leg vertices
vec3(1.0, 2.0, 0.0), vec3(1.0, 1.5, 0.0), vec3(2.0, 1.0, 0.0), vec3(2.0, 1.5, 0.0), //upper left arm vertices
vec3(2.0, 1.5, 0.0), vec3(2.0, 1.0, 0.0), vec3(4.0, 1.0, 0.0), vec3(4.0, 1.5, 0.0), //lower left arm vertices
vec3(-1.0, 2.0, 0.0), vec3(-2.0, 1.5, 0.0), vec3(-2.0, 1.0, 0.0), vec3(-1.0, 1.5, 0.0), //upper right arm vertices
vec3(-2.0, 1.5, 0.0), vec3(-4.0, 1.5, 0.0), vec3(-4.0, 1.0, 0.0), vec3(-2.0, 1.0, 0.0), //lower right arm vertices
vec3(4.0, 1.6, 0.0), vec3(4.0, 0.9, 0.0), vec3(4.5, 0.9, 0.0), vec3(4.5, 1.6, 0.0), //left hand vertices
vec3(-4.0, 1.6, 0.0), vec3(-4.5, 1.6, 0.0), vec3(-4.5, 0.9, 0.0), vec3(-4.0, 0.9, 0.0), //right hand vertices
vec3(1.75, -5.0, 0.0), vec3(0.75, -5.0, 0.0), vec3(0.75, -5.5, 0.0), vec3(1.75, -5.5, 0.0), //left foot vertices
vec3(-0.75, -5.0, 0.0), vec3(-1.75, -5.0, 0.0), vec3(-1.75, -5.5, 0.0), vec3(-0.75, -5.5, 0.0) //right foot vertices
};

vec3 normals[nv]
{
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //face normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //torso normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //upper left leg normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //lower left leg normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //upper right leg normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //lower right leg normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //upper left arm normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //lower left arm normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //upper right arm normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //lower right arm normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //left hand normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //right hand normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), //left foot normals
vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0) //right foot normals
};

ivec3 triangles[nt]
{
ivec3(0, 1, 2), ivec3(0, 2, 3), //face triangles
ivec3(4, 5, 6), ivec3(4, 6, 7), //torso triangles
ivec3(8, 9, 10), ivec3(8, 10, 11), //upper left leg triangles
ivec3(12, 13, 14), ivec3(12, 14, 15), //lower left leg triangles
ivec3(16, 17, 18), ivec3(16, 18, 19), //upper right leg triangles
ivec3(20, 21, 22), ivec3(20, 22, 23), //lower right leg triangles
ivec3(24, 25, 26), ivec3(24, 26, 27), //upper left arm triangles
ivec3(28, 29, 30), ivec3(28, 30, 31), //lower left arm triangles
ivec3(32, 33, 34), ivec3(32, 34, 35), //upper right arm triangles
ivec3(36, 37, 38), ivec3(36, 38, 39), //lower right arm triangles
 ivec3(40, 41, 42), ivec3(40, 42, 43), //left hand triangles
ivec3(44, 45, 46), ivec3(44, 46, 47), //right hand triangles
ivec3(48, 49, 50), ivec3(48, 50, 51), //left foot triangles
ivec3(52, 53, 54), ivec3(52, 54, 55) //right foot triangles
};

ivec2 edges[nv]
{
ivec2(0, 1), ivec2(1, 2), ivec2(2, 3), ivec2(3, 0), //edges of face
ivec2(4, 5), ivec2(5, 6), ivec2(6, 7), ivec2(7, 4), //edges of torso
ivec2(8, 9), ivec2(9, 10), ivec2(10, 11), ivec2(11, 8), //edges of upper left leg
ivec2(12, 13), ivec2(13, 14), ivec2(14, 15), ivec2(15, 12), //edges of lower left leg
ivec2(16, 17), ivec2(17, 18), ivec2(18, 19), ivec2(19, 16), //edges of upper right leg
ivec2(20, 21), ivec2(21, 22), ivec2(22, 23), ivec2(23, 20), //edges of lower right leg
ivec2(24, 25), ivec2(25, 26), ivec2(26, 27), ivec2(27, 24), //edges of upper left arm
ivec2(28, 29), ivec2(29, 30), ivec2(30, 31), ivec2(31, 28), //edges of lower left arm
ivec2(32, 33), ivec2(33, 34), ivec2(34, 35), ivec2(35, 32), //edges of upper right arm
ivec2(36, 37), ivec2(37, 38), ivec2(38, 39), ivec2(39, 36), //edges of lower right arm
ivec2(40, 41), ivec2(41, 42), ivec2(42, 43), ivec2(43, 40), //edges of left hand
ivec2(44, 45), ivec2(45, 46), ivec2(46, 47), ivec2(47, 44), //edges of right hand
ivec2(48, 49), ivec2(49, 50), ivec2(50, 51), ivec2(51, 48), //edges of left foot
ivec2(52, 53), ivec2(53, 54), ivec2(54, 55), ivec2(55, 52) //edges of right foot
};

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

CameraControl camCtl;

struct bone {
	bone *parent = nullptr;
	vec3 origin; //keeps on changing
	vec3 parent_coord_pos; //co-ordinates wrt parent where it is attached, always going to be fixed
	vec3 rotation_axis; //keeps on changing
	std::vector<bone*> children = {};
	std::vector<int> vertex_indices = {};
	std::vector<vec3> vertex_pos = {}; //always going to be fixed
	vec3 initial_origin; //always going to be fixed
	vec3 initial_rotation_axis; //always going to be fixed
	string name;

	float theta = 0.0; // the value that would be updated at each iteration

	bone(bone* parent,vec3 par_pos, vec3 org, vec3 rot_axis, string name) {
		(this->parent) = parent;
		(this->origin) = org;
		(this -> initial_origin) = org;
		(this->parent_coord_pos) = par_pos;
		(this->rotation_axis) = rot_axis;
		(this->initial_rotation_axis) = rot_axis;
		(this->name) = name; 
	}

	void set_vertices(vector<int> vertex_ids) {
		(this -> vertex_indices) = vertex_ids;
		for(int i = 0; i < vertex_ids.size(); i ++) {
			(this -> vertex_pos).push_back(vertices[vertex_ids[i]] - (this -> initial_origin));
		}
	}

	//the below 2 functions would be used while updating at each timestamp
	void set_actual_vertices(Object *object) {
		for(int i = 0; i < vertex_indices.size(); i ++) {
			int id = vertex_indices[i];
			object -> set_vertex(id, vertices[id] + thickness * normals[id]);
			object -> set_vertex(id + nv, vertices[id] - thickness * normals[id]);
		}
	}

	void set_actual_normals(Object *object) {
		for(int i = 0; i < vertex_indices.size(); i ++) {	
			//compute normals for the first half		
			int id0 = vertex_indices[i];
			int id1 = vertex_indices[(i + 1) % (int)vertex_indices.size()];
			int id2 = vertex_indices[(i + (int)vertex_indices.size() - 1) % (int)vertex_indices.size()];
			int id3 = vertex_indices[i] + nv;
			vec3 v0 = object -> get_vertex(id0);
			vec3 v1 = object -> get_vertex(id1);
			vec3 v2 = object -> get_vertex(id2);
			vec3 v3 = object -> get_vertex(id3);

			v3 -= v0;
			v2 -= v0;
			v1 -= v0;

			vec3 norm = cross(v1, v2) + cross(v2, v3) + cross(v3, v1);
			// vec3 norm = cross(v1, v2);

			norm /= length(norm);
			object -> set_normal(id0, norm);

			//compute normals for the second half
			id0 = vertex_indices[i] + nv;
			id1 = vertex_indices[(i + (int)vertex_indices.size() - 1) % (int)vertex_indices.size()] + nv;
			id2 = vertex_indices[(i + 1) % (int)vertex_indices.size()] + nv;
			id3 = vertex_indices[i];
			v0 = object -> get_vertex(id0);
			v1 = object -> get_vertex(id1);
			v2 = object -> get_vertex(id2);
			v3 = object -> get_vertex(id3);

			v1 -= v0;
			v2 -= v0;
			v3 -= v0;

			norm = cross(v1, v2) + cross(v2, v3) + cross(v3, v1);
			// norm = cross(v1, v2);

			norm /= length(norm);
			object -> set_normal(id0, norm);
		}
	}

	void get_lateral_triangles(Object *object) {
		for(int i = 0; i < vertex_indices.size(); i ++) {
			int id0 = vertex_indices[i];
			int id1 = vertex_indices[(i + 1) % (int)vertex_indices.size()];
			int id2 = id0 + nv;
			int id3 = id1 + nv;

			object -> add_triangles(ivec3(id0, id3, id1));
			object -> add_triangles(ivec3(id0, id2, id3));
		}
	}
};

bone *torso,
*head,
*lu_leg,
*ru_leg,
*lu_arm,
*ru_arm,
*ll_leg,
*rl_leg,
*ll_arm,
*rl_arm,
*l_foot,
*r_foot,
*l_hand,
*r_hand;

Object* skeleton;
Object* bottom_plane;

//general structure of an update is go to parent, find its origin,
//then find the position where the child is attached and find its position there as it's
//already updated

vec4 mult_q(vec4 q1, vec4 q2) {
	float a = q1.x;
	vec3 b = vec3(q1.y, q1.z, q1.w);
	float c = q2.x;
	vec3 d = vec3(q2.y, q2.z, q2.w);

	float e = a * c - dot(b, d);
	vec3 f = a * d + c * b + cross(b, d);
	return vec4(e, f.x, f.y, f.z);
}

vec3 rotate(vec3 src, vec3 dest, vec4 rot_q) {
	//we use the quaternion formula here to find the new vector
	//here rot_q is a unit quarternion
	vec4 src_q(0.0, dest.x - src.x, dest.y - src.y, dest.z - src.z);
	vec4 rot_q_inv(rot_q.x, -rot_q.y, -rot_q.z, -rot_q.w);
	vec4 dest_q  = mult_q(mult_q(rot_q, src_q), rot_q_inv);
	return vec3(dest_q.y + src.x, dest_q.z + src.y, dest_q.w + src.z);
}

void update_bone(bone* bone, vec3 disp, vec4 rot_q) {
	//disp is the displacement of the attached point ->
	bone -> origin = disp;
	bone -> rotation_axis = rotate(vec3(0.0, 0.0, 0.0), bone -> initial_rotation_axis, rot_q);

	float ang = bone -> theta;
	vec4 new_rot_q = mult_q(vec4(cos(ang / 2), sin(ang / 2) * (bone -> rotation_axis).x, 
	sin(ang / 2) * (bone -> rotation_axis).y, sin(ang / 2) * (bone -> rotation_axis).z), rot_q);
	new_rot_q = normalize(new_rot_q); //this step is important

	//displacing the vertices
	for(int i = 0; i < (int)(bone->vertex_indices.size()); i ++) {
		int e = bone -> vertex_indices[i];
		vertices[e] = rotate(bone->origin, bone -> origin + bone -> vertex_pos[i],
			new_rot_q);
	}

	for(auto &bp : bone->children) {
		vec3 tmp = bp->parent_coord_pos;
		tmp = rotate(bone->origin, bp -> parent_coord_pos + bone -> origin, 
			new_rot_q);
		update_bone(bp, tmp, new_rot_q);
	}
}

void update_bone_normals(bone* bone) {
	for(int i = 0; i < (int)((bone -> vertex_indices).size()); i ++) {
		if(i == 0) {
			vec3 norm = cross(vertices[(bone -> vertex_indices)[1]] - 
			vertices[(bone -> vertex_indices)[0]], vertices[(bone -> vertex_indices)[2]]
			- vertices[(bone -> vertex_indices)[0]]);
			norm /= length(norm);
			normals[(bone -> vertex_indices)[i]] = norm;
		}
		else if(i == 1) {
			vec3 norm = cross(vertices[(bone -> vertex_indices)[2]] - 
			vertices[(bone -> vertex_indices)[1]], vertices[(bone -> vertex_indices)[0]]
			- vertices[(bone -> vertex_indices)[1]]);
			norm /= length(norm);
			normals[(bone -> vertex_indices)[i]] = norm;
		}
		else {
			vec3 norm = cross(vertices[(bone -> vertex_indices)[0]] - 
			vertices[(bone -> vertex_indices)[i]], vertices[(bone -> vertex_indices)[1]]
			- vertices[(bone -> vertex_indices)[i]]);
			norm /= length(norm);
			normals[(bone -> vertex_indices)[i]] = norm;
		}
	}

	for(int i = 0; i < (int)((bone -> children).size()); i ++) {
		update_bone_normals((bone -> children)[i]);
	}
}

void initialize_actual_bone(bone* bone, Object* object) {
	bone -> set_actual_vertices(object);
	bone -> set_actual_normals(object);
	bone -> get_lateral_triangles(object);
	for(int i = 0; i < (int)(bone -> children).size(); i ++) {
		initialize_actual_bone((bone -> children)[i], object);
	}
}
void update_actual_bone(bone* bone, Object* object) {
	//only vertices and normals need to be considered here
	bone -> set_actual_vertices(object);
	bone -> set_actual_normals(object); 
	for(int i = 0; i < (int)((bone -> children).size()); i ++) {
		update_actual_bone((bone -> children)[i], object);
	}
}

float sine_func(float t, float amp, float w, float phase) {
	return amp * sin(w * t + phase);
}

void initializeScene() {

	//initialize the bones
	torso = new bone(nullptr, vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), vec3(0.0,0.0, 1.0), "torso");
	head = new bone(torso, vec3(0.0, 2.0, 0.0) - torso->origin ,vec3(0.0,2.0,0.0), vec3(0.0, 1.0, 0.0), "head");
	lu_leg = new bone(torso, vec3(0.75, -2.0, 0.0) - torso->origin , vec3(0.75, -2.0, 0.0), vec3(1.0, 0.0, 0.0), "lu_leg");
	ru_leg = new bone(torso, vec3(-0.75, -2.0, 0.0) - torso->origin , vec3(-0.75, -2.0, 0.0), vec3(1.0, 0.0, 0.0), "ru_leg");
	lu_arm = new bone(torso, vec3(1.0, 1.75, 0.0) - torso->origin , vec3(1.0, 1.75, 0.0), vec3(0.0, 1.0, .0), "lu_arm");
	ru_arm = new bone(torso, vec3(-1.0, 1.75, 0.0) - torso->origin , vec3(-1.0, 1.75, 0.0), vec3(0.0, 1.0, 0.0), "ru_arm");
	ll_leg = new bone(lu_leg, vec3(1.25, -3.0, 0.0) - lu_leg->origin , vec3(1.25, -3.0, 0.0), vec3(1.0, 0.0, 0.0), "ll_leg");
	rl_leg = new bone(ru_leg, vec3(-1.25, -3.0, 0.0) - ru_leg->origin , vec3(-1.25, -3.0, 0.0), vec3(1.0, 0.0, 0.0), "rl_leg");
	ll_arm = new bone(lu_arm, vec3(2.0, 1.25, 0.0) - lu_arm->origin , vec3(2.0, 1.25, 0.0), vec3(0.0, 1.0, 0.0), "ll_arm");
	rl_arm = new bone(ru_arm, vec3(-2.0, 1.25, 0.0) - ru_arm->origin , vec3(-2.0, 1.25, 0.0), vec3(0.0, 1.0, 0.0), "rl_arm");
	l_foot = new bone(ll_leg, vec3(1.25, -5.0, 0.0) - ll_leg->origin , vec3(1.25, -5.0, 0.0), vec3(0.0, 0.0, 1.0), "l_foot");
	r_foot = new bone(rl_leg, vec3(-1.25, -5.0, 0.0) - rl_leg->origin , vec3(-1.25, -5.0, 0.0), vec3(0.0, 0.0, 1.0), "r_foot");
	l_hand = new bone(ll_arm, vec3(4.0, 1.25, 0.0) - ll_arm->origin , vec3(4.0, 1.25, 0.0), vec3(0.0, 0.0, 1.0), "l_hand");
	r_hand = new bone(rl_arm, vec3(-4.0, 1.25, 0.0) - rl_arm->origin , vec3(-4.0, 1.25, 0.0), vec3(0.0, 0.0, 1.0), "r_hand");

	//set the children
	torso->children = {head, lu_arm, ru_arm, lu_leg, ru_leg};
	lu_arm->children = {ll_arm};
	ru_arm->children = {rl_arm};
	lu_leg->children = {ll_leg};
	ru_leg->children = {rl_leg};
	ll_arm->children = {l_hand};
	rl_arm->children = {r_hand};
	ll_leg->children = {l_foot};
	rl_leg->children = {r_foot};

	//set the vertex indices, as these would be needed for updating their positions
	head->set_vertices({0, 1, 2, 3});
	torso->set_vertices({4, 5, 6, 7});
	lu_leg->set_vertices({8, 9, 10, 11});
	ll_leg->set_vertices({12, 13, 14, 15});
	ru_leg->set_vertices({16, 17, 18, 19});
	rl_leg->set_vertices({20, 21, 22, 23});
	lu_arm->set_vertices({24, 25, 26, 27});
	ll_arm->set_vertices({28, 29, 30, 31});
	ru_arm->set_vertices({32, 33, 34, 35});
	rl_arm->set_vertices({36, 37, 38, 39});
	l_hand->set_vertices({40, 41, 42, 43});
	r_hand->set_vertices({44, 45, 46, 47});
	l_foot->set_vertices({48, 49, 50, 51});
	r_foot->set_vertices({52, 53, 54, 55});

	skeleton = new Object(2 * nv, 3 * nv, 2 * nt);

	//create actual_edges -> 3nv edges
	for(int i = 0; i < nv; i ++) {
		skeleton->add_edges(edges[i]); //front edges
		skeleton->add_edges(ivec2(edges[i].x + nv, edges[i].y + nv)); //back edges
		skeleton->add_edges(ivec2(i, i + nv)); //cross edges
	}

	//create actual_triangles
	for(int i = 0; i < nt; i ++) {
		skeleton->add_triangles(triangles[i]);
		skeleton->add_triangles(ivec3(triangles[i].x + nv, triangles[i].z + nv, triangles[i].y + nv));
	}

	initialize_actual_bone(torso, skeleton);

	skeleton -> build_vertices();
	skeleton -> build_normals();
	skeleton -> build_edges();
	skeleton -> build_triangles();


	bottom_plane = new Object(4, 4, 2);
	bottom_plane -> set_vertex(0, vec3(5000.0, -5.501, 500.0));
	bottom_plane -> set_vertex(1, vec3(5000.0, -5.501, -500.0));
	bottom_plane -> set_vertex(2, vec3(-5000.0, -5.501, -500.0));
	bottom_plane -> set_vertex(3, vec3(-5000.0, -5.501, 500.0));

	bottom_plane -> set_normal(0, vec3(0.0, 1.0, 0.0));
	bottom_plane -> set_normal(1, vec3(0.0, 1.0, 0.0));
	bottom_plane -> set_normal(2, vec3(0.0, 1.0, 0.0));
	bottom_plane -> set_normal(3, vec3(0.0, 1.0, 0.0));

	bottom_plane -> set_edge(0, ivec2(0, 1));
	bottom_plane -> set_edge(1, ivec2(1, 2));
	bottom_plane -> set_edge(2, ivec2(2, 3));
	bottom_plane -> set_edge(3, ivec2(3, 0));

	bottom_plane -> set_triangle(0, ivec3(0, 1, 2));
	bottom_plane -> set_triangle(1, ivec3(0, 2, 3));

	bottom_plane -> build_vertices();
	bottom_plane -> build_normals();
	bottom_plane -> build_edges();
	bottom_plane -> build_triangles();
}

void updateScene(float t) {
	head ->theta = sine_func(t, M_PI / 3, 0.01, 0.00);
	lu_arm -> theta = -abs(sine_func(t, M_PI / 4, 0.01, 0.00));
	ru_arm -> theta = abs(sine_func(t, M_PI / 4, 0.01, 0.00));
	ll_arm->theta = - abs(sine_func(t, M_PI / 4, 0.01, 0.00));
	rl_arm->theta = abs(sine_func(t, M_PI / 4, 0.01, 0.00));
	lu_leg->theta = -sine_func(t, M_PI / 10, 0.01, 0.00);
	ru_leg->theta = sine_func(t, M_PI / 10, 0.01, 0.00);
	ll_leg->theta = abs(sine_func(t, M_PI / 20, 0.01, 0.00));
	rl_leg->theta = abs(sine_func(t, M_PI / 20, 0.01, 0.00));
	vec3 or_disp = vec3(0.0, 0.0, 0.001 * t);

	update_bone(torso, or_disp, vec4(1.0, 0.0, 0.0, 0.0));
	update_bone_normals(torso);
	update_actual_bone(torso, skeleton);

	skeleton -> update_vertices();
	skeleton -> update_normals();

	bottom_plane -> build_vertices();
	bottom_plane -> build_normals();
}

int main() {
	int width = 640, height = 480;
	if (!r.initialize("Animation", width, height)) {
		return EXIT_FAILURE;
	}
	camCtl.initialize(width, height);
	camCtl.camera.setCameraView(vec3(0.0, -1.0, 12.0), vec3(0.0, -1.0, 0.0), vec3(0.0, 1.0, 0.0));
	//camCtl.camera.setCameraView(vec3(0.0, 10.0, 0.0), vec3(0.0, -1.0 , 0.0), vec3(0.0, 0.0, 1.0));
	program = r.createShaderProgram(
		r.vsBlinnPhong(),
		r.fsBlinnPhong()
	);

	initializeScene();

	while (!r.shouldQuit()) {
        float t = SDL_GetTicks64()*1e-1;
		updateScene(t);

		camCtl.update();
		Camera &camera = camCtl.camera;

		glm::vec3 orange(1.0f, 0.6f, 0.2f);
        glm::vec3 white(1.0f, 1.0f, 1.0f);
		glm::vec3 blue(0.2f, 0.2f, 1.0f);
		glm::vec3 yellow(1.0f, 1.0f, 0.2f);
		glm::vec3 red(1.0f, 0.2f, 0.2f);
		glm::vec3 black(0.0f, 0.0f, 0.0f);

		r.clear(vec4(0.4, 0.4, 0.4, 1.0));
		r.enableDepthTest();
		r.useShaderProgram(program);

		r.setUniform(program, "model", glm::mat4(1.0));
		r.setUniform(program, "view", camera.getViewMatrix());
		r.setUniform(program, "projection", camera.getProjectionMatrix());
		r.setUniform(program, "lightPos", camera.position);
		r.setUniform(program, "viewPos", camera.position);
		r.setUniform(program, "lightColor", white);

		r.setupFilledFaces();

		//this is for the color of the skeleton
        r.setUniform(program, "ambientColor", 0.2f * white);
        r.setUniform(program, "extdiffuseColor", 0.9f * orange);
        r.setUniform(program, "intdiffuseColor", blue);
        r.setUniform(program, "specularColor", white);
        r.setUniform(program, "phongExponent", 20.0f);
		r.drawTriangles(skeleton -> get_object());

		r.setupWireFrame();
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(skeleton -> get_object());


		//this is for the color of the bottom_plane
		r.setUniform(program, "ambientColor", 0.2f * white);
        r.setUniform(program, "extdiffuseColor", yellow);
        r.setUniform(program, "intdiffuseColor", blue);
        r.setUniform(program, "specularColor", white);
        r.setUniform(program, "phongExponent", 20.0f);
		r.drawTriangles(bottom_plane -> get_object());

		r.setupWireFrame();
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(bottom_plane -> get_object());


		//show the image
		r.show();

	}
}
