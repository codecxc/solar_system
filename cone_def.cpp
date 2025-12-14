#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <cmath>
#include <vector>

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct Cone {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    float height;
    float baseRadius;
    int segments;
};

void initCone(Cone& cone, int segments = 64, float height = 1.0f, float baseRadius = 0.5f) {
    cone.segments = segments;
    cone.height = height;
    cone.baseRadius = baseRadius;
    cone.vertices.clear();
    cone.indices.clear();
    
    // Вершина конуса
    Vertex tip;
    tip.x = 0.0f; tip.y = height/2; tip.z = 0.0f;
    tip.r = 1.0f; tip.g = 0.0f; tip.b = 0.0f;
    cone.vertices.push_back(tip);
    
    // Точки основания
    for (int i = 0; i <= segments; i++) {
        float angle = 2.0f * M_PI * i / segments;
        Vertex v;
        v.x = baseRadius * cos(angle);
        v.y = -height/2;
        v.z = baseRadius * sin(angle);
        v.r = 0.0f; v.g = 1.0f; v.b = 0.0f;
        cone.vertices.push_back(v);
    }
    
    // Центр основания
    Vertex center;
    center.x = 0.0f; center.y = -height/2; center.z = 0.0f;
    center.r = 0.0f; center.g = 0.0f; center.b = 1.0f;
    cone.vertices.push_back(center);
    
    int centerIndex = cone.vertices.size() - 1;
    
    // Индексы для боковой поверхности
    for (int i = 1; i <= segments; i++) {
        cone.indices.push_back(0); // вершина
        cone.indices.push_back(i);
        cone.indices.push_back(i + 1);
    }
    
    // Индексы для основания
    for (int i = 1; i <= segments; i++) {
        cone.indices.push_back(centerIndex);
        cone.indices.push_back(i + 1);
        cone.indices.push_back(i);
    }
}

// Функции для деформации конуса

void addProtrusion(Cone& cone, float angle, float heightPos, float radius, float width) {
    // Добавляем выступ на боковой поверхности
    int baseIndex = 1 + (int)(angle * cone.segments / (2 * M_PI)) % cone.segments;
    
    Vertex protrusion;
    protrusion.x = (cone.baseRadius + radius) * cos(angle);
    protrusion.y = heightPos;
    protrusion.z = (cone.baseRadius + radius) * sin(angle);
    protrusion.r = 1.0f; protrusion.g = 1.0f; protrusion.b = 0.0f;
    
    int protIndex = cone.vertices.size();
    cone.vertices.push_back(protrusion);
    
    // Добавляем треугольники для выступа
    int nextBase = baseIndex + 1;
    if (nextBase > cone.segments) nextBase = 1;
    
    cone.indices.push_back(baseIndex);
    cone.indices.push_back(protIndex);
    cone.indices.push_back(nextBase);
}

void addCutout(Cone& cone, float startAngle, float endAngle) {
    // Вырезаем сегмент из основания
    // Упрощенная реализация - удаляем треугольники в диапазоне углов
    int startSeg = (int)(startAngle * cone.segments / (2 * M_PI));
    int endSeg = (int)(endAngle * cone.segments / (2 * M_PI));
    
    // Здесь нужно перестроить индексы, исключив треугольники в этом диапазоне
    // Для простоты просто изменим радиус в этом секторе
    for (int i = startSeg; i <= endSeg; i++) {
        if (i <= cone.segments) {
            float angle = 2.0f * M_PI * i / cone.segments;
            cone.vertices[i + 1].x = cone.baseRadius * 0.3f * cos(angle);
            cone.vertices[i + 1].z = cone.baseRadius * 0.3f * sin(angle);
            cone.vertices[i + 1].r = 0.5f; cone.vertices[i + 1].g = 0.5f; cone.vertices[i + 1].b = 0.5f;
        }
    }
}

void resizeCone(Cone& cone, float newHeight, float newRadius) {
    cone.height = newHeight;
    cone.baseRadius = newRadius;
    
    // Обновляем вершину
    cone.vertices[0].y = newHeight/2;
    
    // Обновляем точки основания
    for (int i = 1; i <= cone.segments + 1; i++) {
        float angle = 2.0f * M_PI * (i - 1) / cone.segments;
        cone.vertices[i].x = newRadius * cos(angle);
        cone.vertices[i].z = newRadius * sin(angle);
        cone.vertices[i].y = -newHeight/2;
    }
    
    // Обновляем центр основания
    cone.vertices[cone.segments + 2].y = -newHeight/2;
}

void bendCone(Cone& cone, float bendAngle) {
    // Изгибаем конус
    for (int i = 1; i <= cone.segments + 1; i++) {
        float yPos = cone.vertices[i].y;
        float bend = bendAngle * (yPos + cone.height/2) / cone.height;
        cone.vertices[i].x += bend;
    }
}

void twistCone(Cone& cone, float twistAngle) {
    // Скручиваем конус
    for (int i = 1; i <= cone.segments + 1; i++) {
        float yPos = cone.vertices[i].y;
        float twist = twistAngle * (yPos + cone.height/2) / cone.height;
        float currentAngle = atan2(cone.vertices[i].z, cone.vertices[i].x);
        float radius = sqrt(cone.vertices[i].x * cone.vertices[i].x + 
                           cone.vertices[i].z * cone.vertices[i].z);
        cone.vertices[i].x = radius * cos(currentAngle + twist);
        cone.vertices[i].z = radius * sin(currentAngle + twist);
    }
}

void drawCone(const Cone& cone) {
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < cone.indices.size(); i++) {
        const Vertex& v = cone.vertices[cone.indices[i]];
        glColor3f(v.r, v.g, v.b);
        glVertex3f(v.x, v.y, v.z);
    }
    glEnd();
}

int main() {
    glfwInit();
    GLFWwindow* window = glfwCreateWindow(800, 600, "Editable Cone", NULL, NULL);
    glfwMakeContextCurrent(window);
    glewInit();
    
    Cone cone;
    initCone(cone, 64, 1.5f, 0.7f);
    
    // Примеры деформаций (раскомментируй нужные):
//     addProtrusion(cone, M_PI/4, 0.0f, 0.3f, 0.2f); // Выступ под 45°
     addCutout(cone, M_PI/2, 3*M_PI/4);            // Вырез от 90° до 135°
    // resizeCone(cone, 2.0f, 0.3f);                 // Вытянутый тонкий конус
    // bendCone(cone, 0.5f);                         // Изогнутый конус
     twistCone(cone, M_PI/3);                      // Скрученный конус
    
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
    
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
