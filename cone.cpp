#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <cmath>

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct Cone {
    Vertex* vertices;
    int vertexCount;
};

void initCone(Cone& cone, int segments = 64) {
    cone.vertexCount = segments * 3 + segments * 3;
    cone.vertices = new Vertex[cone.vertexCount];
    
    float height = 1.0f;
    float radius = 0.5f;
    int index = 0;
    
    for (int i = 0; i < segments; i++) {
        float angle1 = 2.0f * M_PI * i / segments;
        float angle2 = 2.0f * M_PI * (i + 1) / segments;
        
        cone.vertices[index].x = 0.0f;
        cone.vertices[index].y = height/2;
        cone.vertices[index].z = 0.0f;
        cone.vertices[index].r = 1.0f;
        cone.vertices[index].g = 0.0f;
        cone.vertices[index].b = 0.0f;
        index++;
        
g++ program.cpp -o program -lGL -lGLU -lglfw -lGLEW        cone.vertices[index].x = radius * cos(angle1);
        cone.vertices[index].y = -height/2;
        cone.vertices[index].z = radius * sin(angle1);
        cone.vertices[index].r = 1.0f;
        cone.vertices[index].g = 0.0f;
        cone.vertices[index].b = 0.0f;
        index++;
        
        cone.vertices[index].x = radius * cos(angle2);
        cone.vertices[index].y = -height/2;
        cone.vertices[index].z = radius * sin(angle2);
        cone.vertices[index].r = 1.0f;
        cone.vertices[index].g = 0.0f;
        cone.vertices[index].b = 0.0f;
        index++;
    }
    
    for (int i = 0; i < segments; i++) {
        float angle1 = 2.0f * M_PI * i / segments;
        float angle2 = 2.0f * M_PI * (i + 1) / segments;
        
        cone.vertices[index].x = 0.0f;
        cone.vertices[index].y = -height/2;
        cone.vertices[index].z = 0.0f;
        cone.vertices[index].r = 1.0f;
        cone.vertices[index].g = 0.0f;
        cone.vertices[index].b = 0.0f;
        index++;
        
        cone.vertices[index].x = radius * cos(angle1);
        cone.vertices[index].y = -height/2;
        cone.vertices[index].z = radius * sin(angle1);
        cone.vertices[index].r = 1.0f;
        cone.vertices[index].g = 0.0f;
        cone.vertices[index].b = 0.0f;
        index++;
        
        cone.vertices[index].x = radius * cos(angle2);
        cone.vertices[index].y = -height/2;
        cone.vertices[index].z = radius * sin(angle2);
        cone.vertices[index].r = 1.0f;
        cone.vertices[index].g = 0.0f;
        cone.vertices[index].b = 0.0f;
        index++;
    }
}

void drawCone(const Cone& cone) {
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < cone.vertexCount; i++) {
        glColor3f(cone.vertices[i].r, cone.vertices[i].g, cone.vertices[i].b);
        glVertex3f(cone.vertices[i].x, cone.vertices[i].y, cone.vertices[i].z);
    }
    glEnd();
}

void cleanupCone(Cone& cone) {
    delete[] cone.vertices;
}

int main() {
    glfwInit();
    GLFWwindow* window = glfwCreateWindow(800, 600, "3D Cone", NULL, NULL);
    glfwMakeContextCurrent(window);
    glewInit();
    
    Cone cone;
    initCone(cone);
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
        
        drawCone(cone);
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    cleanupCone(cone);
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
