#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct Cube {
    Vertex vertices[8];
};

void initCube(Cube& cube) {
    float vertices[8][3] = {
        {-0.5f, -0.5f, -0.5f},
        { 0.5f, -0.5f, -0.5f},
        { 0.5f,  0.5f, -0.5f},
        {-0.5f,  0.5f, -0.5f},
        {-0.5f, -0.5f,  0.5f},
        { 0.5f, -0.5f,  0.5f},
        { 0.5f,  0.5f,  0.5f},
        {-0.5f,  0.5f,  0.5f}
    };
    
    float colors[8][3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f},
        {1.0f, 1.0f, 0.0f},
        {1.0f, 0.0f, 1.0f},
        {0.0f, 1.0f, 1.0f},
        {0.5f, 0.5f, 0.5f},
        {1.0f, 0.5f, 0.0f}
    };
    
    for (int i = 0; i < 8; i++) {
        cube.vertices[i].x = vertices[i][0];
        cube.vertices[i].y = vertices[i][1];
        cube.vertices[i].z = vertices[i][2];
        cube.vertices[i].r = colors[i][0];
        cube.vertices[i].g = colors[i][1];
        cube.vertices[i].b = colors[i][2];
    }
}

void drawCube(const Cube& cube) {
    int indices[36] = {
        0, 1, 2, 2, 3, 0,
        4, 5, 6, 6, 7, 4,
        0, 3, 7, 7, 4, 0,
        1, 2, 6, 6, 5, 1,
        0, 1, 5, 5, 4, 0,
        3, 2, 6, 6, 7, 3
    };
    
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < 36; i++) {
        int vertexIndex = indices[i];
        glColor3f(cube.vertices[vertexIndex].r, 
                  cube.vertices[vertexIndex].g, 
                  cube.vertices[vertexIndex].b);
        glVertex3f(cube.vertices[vertexIndex].x, 
                   cube.vertices[vertexIndex].y, 
                   cube.vertices[vertexIndex].z);
    }
    glEnd();
}

struct Cube {
	Mesh verticales;
	std::vector<unsigned int> indexes;
};


void createCone(Mesh& mesh, hight, radius, segments) {
	mesh.verticales.push_back(hight/2);

}
int main() {
    glfwInit();
    GLFWwindow* window = glfwCreateWindow(800, 600, "3D Cube", NULL, NULL);
    glfwMakeContextCurrent(window);
    glewInit();
    
    Cube cube;
    initCube(cube);
    glEnable(GL_DEPTH_TEST);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    float angleX = 0.0f;
    float angleY = 0.0f;    
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)width / (float)height;
        gluPerspective(45.0f, aspect, 0.1f, 100.0f);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -3.0f);

	angleX += 0.5f;
        angleY += 0.3f;
        glRotatef(angleX, 1.0f, 0.0f, 0.0f);
        glRotatef(angleY, 0.0f, 1.0f, 0.0f);        
        drawCube(cube);
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
