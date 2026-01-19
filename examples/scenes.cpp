#include "../src/a1.hpp"
#include "/opt/homebrew/include/glm/gtc/matrix_transform.hpp"
#include <iostream>

namespace R = COL781::Software;
// namespace R = COL781::Hardware;
using namespace std;
using namespace glm;

// x-coordinate of plane at a given time
float time_x(float time){
    return (time)-1.0f;
}

// y-coordinate of plane at a given time
float time_y(float time){
    return sqrt((1-time_x(time))/8.0f);
}

// slope of plane-path at a given time
float slope_t(float time){
    return -1*atan(sqrt(1/(2*(1-time_x(time))))/4);
}

int main() {
	R::Rasterizer r;
	int width = 640, height = 480;
    if (!r.initialize("Scene", width, height))
        return EXIT_FAILURE;

    R::ShaderProgram program = r.createShaderProgram(
        r.vsColorTransform(),
        r.fsIdentity()
    );

    // vertices of plane's triangles
    float vertices1[] = {
        -0.2, 0.0, 0.005, 1.0, 
        0.3, 0.0, 0.005, 1.0,
        -0.3, -0.1, 0.2, 1.0,

        -0.2, 0.0, -0.005, 1.0, 
        0.3, 0.0, -0.005, 1.0,
        -0.3, -0.1, -0.2, 1.0,

        -0.2, -0.005, 0.0, 1.0,
        0.3, -0.005, 0.0, 1.0,
        -0.2, -0.1, 0.0, 1.0
    };

    // vertices of square pyramid
    float vertices2[] = {
        0.5, 0.0, -0.5, 1.0,
        0.5, 0.0, 0.5, 1.0,
        -0.5, 0.0, -0.5, 1.0,
        -0.5, 0.0, 0.5, 1.0,
        0.0, 1.0, 0.0, 1.0
    };

    // colors of plane's triangles
    float colors1[] = {
        0.05, 0.9, 0.2, 1.0,
        0.05, 0.9, 0.2, 1.0,
        0.05, 0.9, 0.6, 1.0,

        0.05, 0.9, 0.2, 1.0,
        0.05, 0.9, 0.2, 1.0,
        0.4, 0.9, 0.2, 1.0,

        0.73, 0.1, 0.2, 1.0,
        0.73, 0.1, 0.2, 1.0,
        0.73, 0.1, 0.2, 1.0
    };

    // triangles for plane
    int triangles1[] = {
        0, 1, 2,
        3, 4, 5, 
        6, 7, 8
    };

    // triangles for pyramid
    int triangles2[] = {
        0, 1, 2, 
        1, 2, 3, 
        0, 1, 4,
        1, 3, 4,
        3, 2, 4,
        2, 0, 4
    };

    // colors of pyramid
    float colors2[] = {
        0.55, 0.98, 0.02, 1.0,
        0.02, 0.87, 0.98, 1.0,
        0.02, 0.87, 0.98, 1.0,
        0.55, 0.98, 0.02, 1.0,
        0.98, 0.75, 0.02, 1.0
    };

	R::Object shape = r.createObject();
	{
        r.setVertexAttribs(shape, 0, 9, 4, vertices1);
        r.setVertexAttribs(shape, 1, 9, 4, colors1);
        r.setTriangleIndices(shape, 3, triangles1);
        r.enableDepthTest();
    }

    R::Object pyramid = r.createObject();
    {
        r.setVertexAttribs(pyramid, 0, 5, 4, vertices2);
        r.setVertexAttribs(pyramid, 1, 5, 4, colors2);
        r.setTriangleIndices(pyramid, 6, triangles2);
    }

    mat4 model = mat4(1.0f);
    // Camera view
	mat4 view = translate(mat4(1.0f), vec3(0.0f, 0.0f, -2.3f)); 
    view =  rotate(view, radians(15.0f), vec3(1.0, 0.0, 0.0));
    view =  rotate(view, radians(45.0f), vec3(0.0, 1.0, 0.0));
    // Projection matrix
    mat4 projection = perspective(radians(60.0f), (float)width/(float)height, 0.1f, 100.0f);
    while (!r.shouldQuit()) {
        float time = SDL_GetTicks64()*1e-3;
        time = fmod(time/10, 1.5f);
        r.clear(vec4(0.96, 0.9, 0.75, 1.0));
        r.useShaderProgram(program);
        model = rotate(mat4(1.0f), slope_t(time), vec3(0.0, 0.0, 1.0));
        // Plane moves on a specified path
        mat4 translater = translate(mat4(1.0f), vec3(time_x(time)-0.3, time_y(time), 0.0f));
        mat4 rotater = rotate(mat4(1.0f), radians(0.0f), vec3(0.0, 1.0, 0.0));
        r.setUniform(program, "transform", projection * view * translater * rotater * model);
		r.drawObject(shape);

        time = SDL_GetTicks64()*1e-3;
        translater = translate(mat4(1.0f), vec3(1.0, -1.0, 0.0));
        rotater = rotate(mat4(1.0f), radians(15*time), vec3(0.0, 1.0, 0.0));
        // upscaler magnifies the object
        mat4 upscaler = mat4(1.5f);
        upscaler[3][3] = 1;
        r.setUniform(program, "transform", projection * view * translater * rotater * upscaler * model);
        r.drawObject(pyramid);
        r.show();
    }
    r.deleteShaderProgram(program);
    return EXIT_SUCCESS;
}
