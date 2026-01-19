#include "../src/a1.hpp"
#include <iostream>
#include "/opt/homebrew/include/glm/gtc/matrix_transform.hpp"
// #include <chrono>
#include <ctime>
namespace R = COL781::Software;
// namespace R = COL781::Hardware;
using namespace glm;
using namespace std;

int main() {
	R::Rasterizer r;
	int width = 640, height = 480;
    if (!r.initialize("Clock", width, height))
        return EXIT_FAILURE;

    R::ShaderProgram program = r.createShaderProgram(
        r.vsColorTransform(),
        r.fsIdentity()
    );

    // vertices for clock hands
    float vertices[] = {
        -0.025, -0.1, 0.0, 1.0,
        0.025, -0.1, 0.0, 1.0,
        0.0, 0.5, 0.0, 1.0,

        -0.015, -0.1, 0.0, 1.0,
        0.015, -0.1, 0.0, 1.0,
        0.0, 0.5, 0.0, 1.0
    };

    // vertices for markings on the clock
    float new_vertices[] = {
        0.025, 0.05, 0.0, 1.0,
        0.025, -0.05, 0.0, 1.0,
        -0.025, 0.05, 0.0, 1.0,
        -0.025, -0.05, 0.0, 1.0
    };

    // colors of clocks hands
    float colors[] = {
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0,

        1.0, 0.0, 0.0, 1.0,
        1.0, 0.0, 0.0, 1.0,
        1.0, 0.0, 0.0, 1.0
    };

    // Triangles for clock hands
    int triangles[] = {
        0, 1, 2
    };

    // Triangles for clock markings
    int new_triangles[] = {
        0, 1, 2,
        1, 2, 3
    };

    // Color for clock markings
    float new_colors[] = {
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0
    };

    R::Object shape = r.createObject();
    r.enableDepthTest();
    r.useShaderProgram(program);

    while (!r.shouldQuit()) {
        r.clear(vec4(0.2, 0.89, 0.99, 1.0));
        std::time_t t = std::time(0);   
        // Fethching current time
        std::tm* now = std::localtime(&t);
        float hours = fmod(now->tm_hour+5.0f, 12.0f);
        if(now->tm_min >= 30.0f){hours += 1.0f;hours = fmod(hours, 12.0f);}
        float minutes = fmod(now->tm_min+30, 60.0f);
        float seconds = now->tm_sec;
        // Computing angles
        float hour_angle = 30*hours+minutes/2;
        float minute_angle = 6*minutes;
        float second_angle = 6*seconds;
        mat4 rotater = rotate(mat4(1.0f), -1*radians(hour_angle), vec3(0.0f,0.0f,1.0f));
        // upscale magnifies the object
        mat4 upscale = mat4(0.8f);
        upscale[3][3] = 1;
        mat4 translater = translate(mat4(1.0f), vec3(0.0f, 0.9f, 0.0f)); 
        translater = mat4(1.0f);
        
	    r.setVertexAttribs(shape, 0, 6, 4, vertices);
        r.setVertexAttribs(shape, 1, 6, 4, colors);
        r.setTriangleIndices(shape, 1, triangles);
        r.setUniform(program, "transform", rotater * translater * upscale);
		r.drawObject(shape);

        rotater = rotate(mat4(1.0f), -1*radians(minute_angle), vec3(0.0f,0.0f,1.0f));
        // upscale magnifies the object
        upscale = mat4(1.3f);
        upscale[3][3] = 1;
        translater = translate(mat4(1.0f), vec3(0.0f, 0.0f, -0.1f)); 
        shape = r.createObject();
	    r.setVertexAttribs(shape, 0, 6, 4, vertices);
        r.setVertexAttribs(shape, 1, 6, 4, colors);
        r.setTriangleIndices(shape, 1, triangles);
        r.setUniform(program, "transform", rotater * translater * upscale);
		r.drawObject(shape);

        int new_indices[] = {3, 4, 5};
        rotater = rotate(mat4(1.0f), -1*radians(second_angle), vec3(0.0f,0.0f,1.0f));
        // upscale magnifies the object
        upscale = mat4(1.75f);
        upscale[3][3] = 1;
        translater = translate(mat4(1.0f), vec3(0.0f, 0.0f, -0.2f)); 
        shape = r.createObject();
	    r.setVertexAttribs(shape, 0, 6, 4, vertices);
        r.setVertexAttribs(shape, 1, 6, 4, colors);
        r.setTriangleIndices(shape, 1, new_indices);
        r.setUniform(program, "transform", rotater * translater * upscale);
		r.drawObject(shape);

        for(int i = 0; i < 12; i++){
            rotater = rotate(mat4(1.0f), -1*radians(i*30.0f), vec3(0.0f,0.0f,1.0f));
            // upscale magnifies the object
            upscale = mat4(1.0);
            upscale[3][3] = 1;
            translater = translate(mat4(1.0f), vec3(0.0f, 0.9f, 0.1f)); 
            shape = r.createObject();
            r.setVertexAttribs(shape, 0, 4, 4, new_vertices);
            r.setVertexAttribs(shape, 1, 4, 4, new_colors);
            r.setTriangleIndices(shape, 2, new_triangles);
            r.setUniform(program, "transform", rotater * translater * upscale);
            r.drawObject(shape);
        }
        r.show();
    }

    r.deleteShaderProgram(program);
    return EXIT_SUCCESS;
}