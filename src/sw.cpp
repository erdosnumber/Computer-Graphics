#include "sw.hpp"

#include <iostream>
#include <vector>

using namespace glm;
using namespace std;
namespace COL781 {
	namespace Software {

		// Forward declarations

		template <> float Attribs::get(int index) const;
		template <> glm::vec2 Attribs::get(int index) const;
		template <> glm::vec3 Attribs::get(int index) const;
		template <> glm::vec4 Attribs::get(int index) const;

		template <> void Attribs::set(int index, float value);
		template <> void Attribs::set(int index, glm::vec2 value);
		template <> void Attribs::set(int index, glm::vec3 value);
		template <> void Attribs::set(int index, glm::vec4 value);

		// Built-in shaders

		VertexShader Rasterizer::vsIdentity() {
			return [](const Uniforms &uniforms, const Attribs &in, Attribs &out) {
				glm::vec4 vertex = in.get<glm::vec4>(0);
				return vertex;
			};
		}

		VertexShader Rasterizer::vsTransform() {
			return [](const Uniforms &uniforms, const Attribs &in, Attribs &out) {
				glm::vec4 vertex = in.get<glm::vec4>(0);
				glm::mat4 transform = uniforms.get<glm::mat4>("transform");
				return transform * vertex;
			};
		}

		VertexShader Rasterizer::vsColor() {
			return [](const Uniforms &uniforms, const Attribs &in, Attribs &out) {
				glm::vec4 vertex = in.get<glm::vec4>(0);
				glm::vec4 color = in.get<glm::vec4>(1);
				out.set<glm::vec4>(0, color);
				return vertex;
			};
		}

		FragmentShader Rasterizer::fsConstant() {
			return [](const Uniforms &uniforms, const Attribs &in) {
				// cout << 'A' << endl;
				glm::vec4 color = uniforms.get<glm::vec4>("color");
				return color;
			};
		}

		FragmentShader Rasterizer::fsIdentity() {
			return [](const Uniforms &uniforms, const Attribs &in) {
				glm::vec4 color = in.get<glm::vec4>(0);
				return color;
			};
		}

		// Implementation of Attribs and Uniforms classes

		void checkDimension(int index, int actual, int requested) {
			if (actual != requested) {
				std::cout << "Warning: attribute " << index << " has dimension " << actual << " but accessed as dimension " << requested << std::endl;
			}
		}

		template <> float Attribs::get(int index) const {
			checkDimension(index, dims[index], 1);
			return values[index].x;
		}

		template <> glm::vec2 Attribs::get(int index) const {
			checkDimension(index, dims[index], 2);
			return glm::vec2(values[index].x, values[index].y);
		}

		template <> glm::vec3 Attribs::get(int index) const {
			checkDimension(index, dims[index], 3);
			return glm::vec3(values[index].x, values[index].y, values[index].z);
		}

		template <> glm::vec4 Attribs::get(int index) const {
			checkDimension(index, dims[index], 4);
			return values[index];
		}

		void expand(std::vector<int> &dims, std::vector<glm::vec4> &values, int index) {
			if (dims.size() < index+1)
				dims.resize(index+1);
			if (values.size() < index+1)
				values.resize(index+1);
		}

		template <> void Attribs::set(int index, float value) {
			expand(dims, values, index);
			dims[index] = 1;
			values[index].x = value;
		}

		template <> void Attribs::set(int index, glm::vec2 value) {
			expand(dims, values, index);
			dims[index] = 2;
			values[index].x = value.x;
			values[index].y = value.y;
		}

		template <> void Attribs::set(int index, glm::vec3 value) {
			expand(dims, values, index);
			dims[index] = 3;
			values[index].x = value.x;
			values[index].y = value.y;
			values[index].z = value.z;
		}

		template <> void Attribs::set(int index, glm::vec4 value) {
			expand(dims, values, index);
			dims[index] = 4;
			values[index] = value;
		}

		template <typename T> T Uniforms::get(const std::string &name) const {
			return *(T*)values.at(name);
		}

		template <typename T> void Uniforms::set(const std::string &name, T value) {
			auto it = values.find(name);
			if (it != values.end()) {
				delete it->second;
			}
			values[name] = (void*)(new T(value));
		}

		vec4 Convert_Color(vec4 color){
			Uint8 r = static_cast<Uint8>(color.r * 255);
			Uint8 g = static_cast<Uint8>(color.g * 255);
			Uint8 b = static_cast<Uint8>(color.b * 255);
			Uint8 a = static_cast<Uint8>(color.a * 255);
			return vec4(r, g, b, a);
		}

		void clearScreen(SDL_Renderer* renderer, const glm::vec4& color) {
			vec4 new_color = Convert_Color(color);
			SDL_SetRenderDrawColor(renderer, new_color.r, new_color.g, new_color.b, new_color.a);
			SDL_RenderClear(renderer);
		}

		bool Rasterizer::initialize(const std::string &title, int width, int height, int spp){
			this->width = width;
			this->height = height;
			if (SDL_Init(SDL_INIT_VIDEO) < 0) {
				SDL_Log("SDL could not initialize! SDL_Error: %s", SDL_GetError());
				this->quit = true;
				return 0;
    		}
			const char *title_char = (&title)->c_str();
    		SDL_Window* window = SDL_CreateWindow(title_char,
                                          SDL_WINDOWPOS_CENTERED,
                                          SDL_WINDOWPOS_CENTERED,
                                          width, height,
                                          SDL_WINDOW_SHOWN);
			if (!window) {
				SDL_Log("Window could not be created! SDL_Error: %s", SDL_GetError());
				SDL_Quit();
				this->quit = true;
				return 0;
			}
			this->window = window;
			this->renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
			if (!renderer) {
				SDL_Log("Renderer could not be created! SDL_Error: %s", SDL_GetError());
				SDL_DestroyWindow(window);
				SDL_Quit();
				this->quit = true;
				return 0;
			}
			return true;
		}

		void Rasterizer::clear(glm::vec4 color){
			this->used.clear();
			this->z_buffer.clear();
			clearScreen(this->renderer, color);
		}
		
		bool Rasterizer::shouldQuit(){
			SDL_Event event;
			while (SDL_PollEvent(&event)) {
				if (event.type == SDL_QUIT) {
					this->quit = true;  
					break;
				}
			}
			return this->quit;
		}

		Object Rasterizer::createObject(){
			Object temp = Object();
			return temp;
		}

		void Rasterizer::setVertexAttribs(Object &object, int attribIndex, int n, int d, const float* data){
			int iter = 0;
			int N = 0;
			object.vertexAttributes.resize(n);
			while(N < n){
				vec4 value = vec4(1.0, 1.0, 1.0, 1.0);
				if(d >= 1){value.x = data[iter];iter++;}
				if(d >= 2){value.y = data[iter];iter++;}
				if(d >= 3){value.z = data[iter];iter++;}
				if(d >= 4){value.w = data[iter];iter++;}
				object.vertexAttributes[N].set<vec4>(attribIndex, value);
				N++;
			}
		}

		void Rasterizer::setTriangleIndices(Object &object, int n, int* indices){
			int iter = 0;
			for(int i = 0; i < n; i++){
				ivec3 triangle;
				triangle.x = indices[iter];
				triangle.y = indices[iter+1];
				triangle.z = indices[iter+2];
				iter+=3;
				object.indices.push_back(triangle);
			}
		}

		ShaderProgram Rasterizer::createShaderProgram(const VertexShader &vs, const FragmentShader &fs){
			ShaderProgram shaderProgram  = ShaderProgram();
			shaderProgram.fs = fs;
			shaderProgram.vs = vs;
			return shaderProgram;
		}

		void Rasterizer::useShaderProgram(const ShaderProgram &program){
			this->program = &program;
			this->shader_flag = true;
		}

		void Rasterizer::deleteShaderProgram(ShaderProgram &program){
			this->shader_flag = false;
		}
		
		void Rasterizer::show(){
			SDL_RenderPresent(this->renderer);
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, vec4 value){
			(program.uniforms).set<vec4>(name, value);
		}
		
		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, mat4 value){
			(program.uniforms).set<mat4>(name, value);
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, vec3 value){
			(program.uniforms).set<vec3>(name, value);
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, int value){
			(program.uniforms).set<float>(name, value);
		}

		vec2 glm_to_screen_space(vec4 point, float width, float height){
			float max_x = 1.5;
			float max_y = max_x*(height/width);
			float x = point.x;
			float y = -1*point.y;
			vec2 screen_point;
			screen_point.x = x*width/(2*max_x);
			screen_point.x += width/2;
			screen_point.y = y*height/(2*max_y);
			screen_point.y += height/2; 
			return screen_point;
		}

		float Area(vec4 point1, vec4 point2, vec4 point3){
			float area = abs(length(cross(vec3(point1)-vec3(point2), vec3(point2)-vec3(point3))))/2.0;
			return area;
		}
		
		float find_z(vec4 point1, vec4 point2, vec4 point3, float x, float y){
			vec3 normal = cross(vec3(point1-point2), vec3(point1-point3));
			float dot = normal.x*(x-point1.x)+normal.y*(y-point1.y);
			float z = (-1*dot)/normal.z+point1.z;
			return z;
		}

		void Rasterizer::drawObject(const Object &object){
			float Increment = 0.002;
			float x_start = -1.0, y_start = -1.0;
			float x_end = 1.0, y_end = 1.0;
			for(int i = 0; i < object.indices.size(); i++){
				Attribs temp;
				vec4 point1 = this->program->vs((this->program)->uniforms, object.vertexAttributes[(object.indices)[i].x], temp);
				vec4 color1 = this->program->fs((this->program)->uniforms, temp);
				vec4 point2 = this->program->vs((this->program)->uniforms, object.vertexAttributes[(object.indices)[i].y], temp);
				vec4 color2 = this->program->fs((this->program)->uniforms, temp);
				vec4 point3 = this->program->vs((this->program)->uniforms, object.vertexAttributes[(object.indices)[i].z], temp);
				vec4 color3 = this->program->fs((this->program)->uniforms, temp);
				float w1 = point1.w, w2 = point2.w, w3 = point3.w;
				point1 /= w1;point2 /= w2;point3 /= w3;
				point1.w = point1.w/w1;point2.w = point2.w/w2;point3.w = point3.w/w3;
				color1 /= w1;color2 /= w2;color3 /= w3;
				float xx = x_start;
				while(xx <= x_end){
					float yy = y_start;
					while(yy <= y_end){
						int check = 0;
						if(((xx-point1.x)*(point1.y-point2.y)-(point1.x-point2.x)*(yy-point1.y))*((point3.x-point1.x)*(point1.y-point2.y)-(point1.x-point2.x)*(point3.y-point1.y)) >= 0){
							check++;
						}
						if(((xx-point2.x)*(point3.y-point2.y)-(point3.x-point2.x)*(yy-point2.y))*((point1.x-point2.x)*(point3.y-point2.y)-(point3.x-point2.x)*(point1.y-point2.y)) >= 0){
							check++;
						}
						if(((xx-point1.x)*(point3.y-point1.y)-(point3.x-point1.x)*(yy-point1.y))*((point2.x-point1.x)*(point3.y-point1.y)-(point3.x-point1.x)*(point2.y-point1.y)) >= 0){
							check++;
						}
						if(check == 3){
							float zz = find_z(point1, point2, point3, xx, yy);
							float weight1 = Area(point2, point3, vec4(xx, yy, zz, 0));
							float weight2 = Area(point1, point3, vec4(xx, yy, zz, 0));
							float weight3 = Area(point1, point2, vec4(xx, yy, zz, 0));
							float w_inverse = (weight1*point1.w+weight2*point2.w+weight3*point3.w)/(weight1+weight2+weight3);
							vec4 color = (weight1*color1+weight2*color2+weight3*color3)/(weight1+weight2+weight3);
							color /= w_inverse;
							vec4 new_color = Convert_Color(color);
							SDL_SetRenderDrawColor(renderer, new_color.r, new_color.g, new_color.b, new_color.a);
							vec2 screen_point = glm_to_screen_space(vec4(xx, yy, 0, 0), this->width, this->height);
							if(abs(zz) > 1){yy += Increment;continue;}
							if(!this->depth_test || !this->used[{xx,yy}] || this->z_buffer[{xx, yy}] > zz){SDL_RenderDrawPoint(renderer, screen_point.x, screen_point.y);}
							if(!this->used[{xx, yy}] || this->z_buffer[{xx,yy}] > zz){this->z_buffer[{xx, yy}] = zz;}
							this->used[{xx, yy}] = true;
						}
						yy += Increment;
					}
					xx += Increment;
				}
			}
		}
		
		VertexShader Rasterizer::vsColorTransform() {
			return [](const Uniforms &uniforms, const Attribs &in, Attribs &out) {
				glm::vec4 vertex = in.get<glm::vec4>(0);
				glm::vec4 color = in.get<glm::vec4>(1);
				glm::mat4 transform = uniforms.get<glm::mat4>("transform");
				vec4 temp = transform * vertex;
				out.set<glm::vec4>(0, color);
				return temp;
			};
		}
		
		VertexShader Rasterizer::vsNormalTransform(){
			return [](const Uniforms &uniforms, const Attribs &in, Attribs &out) {
				glm::vec4 vertex = in.get<glm::vec4>(0);
				glm::vec4 color = in.get<glm::vec4>(1);
				glm::mat4 transform = uniforms.get<glm::mat4>("transform");
				mat4 wsTransform = uniforms.get<mat4>("wsTransform");
				vec4 temp = transform * vertex;
				color = vec4(transpose(inverse(mat3(wsTransform))) * vec3(color), 0.0);
				out.set<glm::vec4>(0, color);
				return temp;
			};
		}

		FragmentShader Rasterizer::fsDiffuseLighting(){
			return [](const Uniforms &uniforms, const Attribs &in) {
				vec3 objectColor = uniforms.get<vec3>("objectColor");
				vec3 ambientColor = uniforms.get<vec3>("ambientColor");
				vec3 lightDir = uniforms.get<vec3>("lightDir");
				vec3 lightColor = uniforms.get<vec3>("lightColor");
				mat4 wsTransform = uniforms.get<mat4>("wsTransform");
				vec4 normal = in.get<vec4>(0);
				vec3 light = lightDir;
				vec3 new_normal = normalize(vec3(normal));
				light = normalize(light);
				float factor = 0.0;
				if(factor < dot(new_normal, light)){factor = dot(new_normal, light);}
				vec3 diffuse = lightColor * objectColor;
				diffuse *= factor;
				vec3 ambient = ambientColor * objectColor;
				vec3 final_light = diffuse+ambient;
				vec4 final_color = vec4(final_light.x, final_light.y, final_light.z, 1.0);
				return final_color;
			};
		}

		FragmentShader Rasterizer::fsSpecularLighting(){
			return [](const Uniforms &uniforms, const Attribs &in) {
				vec3 objectColor = uniforms.get<vec3>("objectColor");
				vec3 ambientColor = uniforms.get<vec3>("ambientColor");
				vec3 lightDir = uniforms.get<vec3>("lightDir");
				vec3 lightColor = uniforms.get<vec3>("lightColor");
				vec3 specularColor = uniforms.get<vec3>("specularColor");
				float blinnpow = uniforms.get<float>("blinnpow");
				vec3 viewPos = uniforms.get<vec3>("viewPos");
				viewPos = normalize(viewPos);
				mat4 wsTransform = uniforms.get<mat4>("wsTransform");
				vec4 normal = in.get<vec4>(0);
				vec3 light = lightDir;
				vec3 new_normal = normalize(vec3(normal));
				light = normalize(light);
				vec3 H = (light+viewPos)/length(light+viewPos);
				H = normalize(H);
				float factor = 0.0;
				if(factor < dot(new_normal, light)){factor = dot(new_normal, light);}
				vec3 diffuse = lightColor * objectColor;
				diffuse *= factor;
				vec3 ambient = ambientColor * objectColor;
				float spec = 0.0;
				if(spec < dot(H, new_normal)){spec = dot(H, new_normal);}
				float specular_factor = pow(spec, blinnpow);
				vec3 specular = specularColor * specular_factor;
				vec3 final_light = diffuse+ambient+specular;
				vec4 final_color = vec4(final_light.x, final_light.y, final_light.z, 1.0);
				return final_color;
			};
		}
		
		void Rasterizer::enableDepthTest(){
			this->depth_test = true;
		}
	}
}
