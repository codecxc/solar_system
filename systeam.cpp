#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
struct Vertex {
	float x,y,z;
	float r,g,b;
};
struct Planet {
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
};

void createPlanet(Planet& planet,float radius) {
	int s=400;
	for(int i=0;i<s;++i) {
		float phi=M_PI*i/s;
		for(int j=0;j<=s;++j) {
			float theta=2*M_PI*j/s;
			Vertex v;
			v.x=radius*sin(phi)*cos(theta);
			v.y=radius*cos(phi);
			v.z=radius*sin(phi)*sin(theta);
            float nx = sin(phi) * cos(theta);
            float ny = cos(phi);
            float nz = sin(phi) * sin(theta);
            
            float light = (ny + 1.0f) * 0.5f;
            
            v.r = light * 0.8f;
            v.g = light * 0.5f;
            v.b = light * 0.2f;
			planet.vertices.push_back(v);
		}
	}
	for(int i=0;i<s;++i) {
		for(int j=0;j<s;++j) {
			int first=i*(s+1)+j;
			int second=first+1;
			int third=first+(s+1);
			int fourth=third+1;
			planet.indices.push_back(first);
            		planet.indices.push_back(second);
            		planet.indices.push_back(third);
            
            		planet.indices.push_back(second);
            		planet.indices.push_back(fourth);
            		planet.indices.push_back(third);
		}
	
	}

}


void drawPlanet(const Planet& planet) {
    glBegin(GL_TRIANGLES);
    for(int i = 0; i < planet.indices.size(); i++) {
        int vertexIndex = planet.indices[i];
        const Vertex& v = planet.vertices[vertexIndex];
        glColor3f(v.r, v.g, v.b);
        glVertex3f(v.x, v.y, v.z);
    }
    glEnd();
}

int main() {
if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
	GLFWwindow* window = glfwCreateWindow(800, 600, "3D Planet", NULL, NULL);
	glfwMakeContextCurrent(window);
if (glewInit() != GLEW_OK) {
    std::cerr << "Failed to initialize GLEW" << std::endl;
    return -1;
}	
Planet planet;
    	createPlanet(planet, 1.0f);
	glEnable(GL_DEPTH_TEST);
	float rotationAngle = 0.0f;
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)width / (float)height;
  glFrustum(-0.1 * aspect, 0.1 * aspect, -0.1, 0.1, 0.1, 100.0);      
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        glTranslatef(0.0f, 0.0f, -5.0f);
	glRotatef(30.0f, 1.0f, 0.0f, 0.0f);        
         rotationAngle += 10.0f;
        glRotatef(rotationAngle, 0.0f, 1.0f, 0.0f);
	        drawPlanet(planet);
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;

}
