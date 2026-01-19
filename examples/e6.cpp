#include "../src/a1.hpp"
#include <glm/gtc/matrix_transform.hpp>
// Program with perspective correct interpolation of vertex attributes.

namespace R = COL781::Software;
// namespace R = COL781::Hardware;
using namespace glm;
int main() {
	R::Rasterizer r;
	int width = 640, height = 480;
    if (!r.initialize("Example 6", width, height))
        return EXIT_FAILURE;

    R::ShaderProgram program = r.createShaderProgram(
        r.vsNormalTransform(),
        r.fsDiffuseLighting()
        // r.fsSpecularLighting()
    );

    float vertices[] = {
        -0.5,  -0.5, -0.5, 1.0,
         0.5,  -0.5, -0.5, 1.0,
        -0.5,   0.5, -0.5, 1.0,
         0.5,   0.5, -0.5, 1.0,

        -0.5,  -0.5,  0.5, 1.0,
         0.5,  -0.5,  0.5, 1.0,
        -0.5,   0.5,  0.5, 1.0,
         0.5,   0.5,  0.5, 1.0,

        -0.5,  -0.5, -0.5, 1.0,
        -0.5,  -0.5,  0.5, 1.0,
        -0.5,   0.5, -0.5, 1.0,
        -0.5,   0.5,  0.5, 1.0,

         0.5,  -0.5, -0.5, 1.0,
         0.5,  -0.5,  0.5, 1.0,
         0.5,   0.5, -0.5, 1.0,
         0.5,   0.5,  0.5, 1.0,

        -0.5,  -0.5, -0.5, 1.0,
         0.5,  -0.5, -0.5, 1.0,
        -0.5,  -0.5,  0.5, 1.0,
         0.5,  -0.5,  0.5, 1.0,

        -0.5,   0.5, -0.5, 1.0,
         0.5,   0.5, -0.5, 1.0,
        -0.5,   0.5,  0.5, 1.0,
         0.5,   0.5,  0.5, 1.0
    };

    float normals[] = {
        0.0,  0.0, -1.0,
        0.0,  0.0, -1.0,
        0.0,  0.0, -1.0,
        0.0,  0.0, -1.0,

        0.0,  0.0,  1.0,
        0.0,  0.0,  1.0,
        0.0,  0.0,  1.0,
        0.0,  0.0,  1.0,

        -1.0,  0.0,  0.0,
        -1.0,  0.0,  0.0,
        -1.0,  0.0,  0.0,
        -1.0,  0.0,  0.0,

        1.0,  0.0,  0.0,
        1.0,  0.0,  0.0,
        1.0,  0.0,  0.0,
        1.0,  0.0,  0.0,

        0.0, -1.0,  0.0,
        0.0, -1.0,  0.0,
        0.0, -1.0,  0.0,
        0.0, -1.0,  0.0,

        0.0,  1.0,  0.0,
        0.0,  1.0,  0.0,
        0.0,  1.0,  0.0,
        0.0,  1.0,  0.0
    };

    int triangles[] = {
        0, 1, 2,
        1, 2, 3,

        4, 5, 6,
        5, 6, 7,

        8, 9, 10,
        9, 10, 11,

        12, 13, 14,
        13, 14, 15,

        16, 17, 18,
        17, 18, 19,

        20, 21, 22,
        21, 22, 23
    };
	R::Object shape = r.createObject();
	r.setVertexAttribs(shape, 0, 24, 4, vertices);
	r.setVertexAttribs(shape, 1, 24, 3, normals);
	r.setTriangleIndices(shape, 12, triangles);
    r.enableDepthTest();
    vec3 objectColor(0.8f, 0.4f, 0.248f);
    vec3 ambientColor(0.2f, 0.1f, 0.062f);
    vec3 specularColor(0.8f, 0.8f, 0.8f);

    // The transformation matrix.
    mat4 model = mat4(1.0f);
	mat4 view = translate(mat4(1.0f), vec3(0.0f, 0.0f, -3.0f));
    view = rotate(view, radians(45.0f), vec3(1.0f, 0.0f, 0.0f));
    mat4 projection = perspective(radians(60.0f), (float)width/(float)height, 0.1f, 100.0f);

    // Lighting parameters.
    vec3 lightColor(1.0f, 1.0f, 1.0f);
    vec3 lightDir(0.0f, 1.0f, -1.0f);

    while (!r.shouldQuit()) {
        float time = SDL_GetTicks64()*1e-3;
        r.clear(vec4(0.1, 0.1, 0.1, 1.0));
        r.useShaderProgram(program);
        model = rotate(mat4(1.0f), radians(45.0f), vec3(0.0f,1.0f,0.0f));
        lightDir = vec3(cos(time), 0.5f, sin(time));
        r.setUniform(program, "transform", projection * view * model);
        r.setUniform(program, "wsTransform", model);
        r.setUniform(program, "lightColor", lightColor);
        r.setUniform(program, "lightDir", lightDir);
        r.setUniform(program, "objectColor", objectColor);
        r.setUniform(program, "ambientColor", ambientColor);
        // for specular lighting
        r.setUniform(program, "viewPos", vec3(inverse(view) * vec4(0.0f, 0.0f, 0.0f, 1.0f)));
        r.setUniform(program, "specularColor", specularColor);
        r.setUniform(program, "blinnpow", 32);
		r.drawObject(shape);
        r.show();
    }
    r.deleteShaderProgram(program);
    return EXIT_SUCCESS;
}
