#include "../src/a1.hpp"

namespace R = COL781::Software;
// namespace R = COL781::Hardware;
using namespace glm;

int main() {
	R::Rasterizer r;
    if (!r.initialize("Example 2", 640, 480))
        return EXIT_FAILURE;
    R::ShaderProgram program = r.createShaderProgram(
        r.vsColor(),
        r.fsIdentity()
    );

    float vertices[] = {
        0.0, -0.8, 0.0, 1.0,
        0.8,  0.0, 0.0, 1.0,
        0.0,  0.8, 0.0, 1.0,
        -0.8,  0.0, 0.0, 1.0
    };

    float colors[] = {
        0.0, 0.0, 0.0, 1.0,
        1.0, 0.0, 0.0, 1.0,
        1.0, 1.0, 1.0, 1.0,
        0.0, 0.0, 1.0, 1.0
    };

    int triangles[] = {
        0, 1, 2,
        0, 2, 3
    };

	R::Object shape = r.createObject();
	r.setVertexAttribs(shape, 0, 4, 4, vertices);
	r.setVertexAttribs(shape, 1, 4, 4, colors);
	r.setTriangleIndices(shape, 2, triangles);
    while (!r.shouldQuit()) {
        r.clear(vec4(1.0, 1.0, 1.0, 1.0));
        r.useShaderProgram(program);
		r.drawObject(shape);
        r.show();
    }
    r.deleteShaderProgram(program);
    return EXIT_SUCCESS;
}
