#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <algorithm>
struct Triangle {
    float size;
};
void updateTriangleSize(Triangle& triangle, int width, int height) {
    float scale = std::min((float)width / 800.0f, (float)height / 600.0f);
    triangle.size = 0.5f * scale/2;
}
void drawTriangle(const Triangle& triangle) {
    glBegin(GL_TRIANGLES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, triangle.size, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(-triangle.size, -triangle.size, 0.0f);
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(triangle.size, -triangle.size, 0.0f);
    glEnd();
}
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}




int main() {
    glfwInit();
    GLFWwindow* window = glfwCreateWindow(800, 600, "Мое окно", NULL, NULL);
    glfwMakeContextCurrent(window);

	Triangle triangle;
	triangle.size=0.5f;

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    updateTriangleSize(triangle, width, height);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glfwGetFramebufferSize(window, &width, &height);
        updateTriangleSize(triangle, width, height);
	drawTriangle(triangle);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
