#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>

struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct Planet {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    float orbitRadius;//радиус
    float orbitSpeed;//орбитальная скорость
    float rotationSpeed;//вокруг оси вращение
    float orbitAngle;     // текущий угол орбиты
    float rotationAngle;  // текущий угол вращения
    float size;           // размер планеты
};

void createPlanet(Planet& planet, float radius, float orbitRadius, float orbitSpeed, float rotationSpeed, float size) {
    planet.orbitRadius = orbitRadius;
    planet.orbitSpeed = orbitSpeed;
    planet.rotationSpeed = rotationSpeed;
    planet.orbitAngle = 0.0f;
    planet.rotationAngle = 0.0f;
    planet.size = size;
    
    int segments = 32; // Уменьшил для производительности
    
    // Создаем вершины сферы
    for(int i = 0; i <= segments; ++i) {
        float phi = M_PI * i / segments;
        for(int j = 0; j <= segments; ++j) {
            float theta = 2 * M_PI * j / segments;
            Vertex v;
            
            // Координаты сферы
            v.x = radius * sin(phi) * cos(theta);
            v.y = radius * cos(phi);
            v.z = radius * sin(phi) * sin(theta);
            
            // Цвета на основе положения
            v.r = 0.2f + 0.8f * (float)i / segments;  // красный
            v.g = 0.1f + 0.4f * (float)j / segments;  // зеленый
            v.b = 0.3f + 0.7f * (1.0f - (float)i / segments); // синий
            
            planet.vertices.push_back(v);
        }
    }
    
    // Создаем индексы для треугольников
    for(int i = 0; i < segments; ++i) {
        for(int j = 0; j < segments; ++j) {
            int first = i * (segments + 1) + j;
            int second = first + 1;
            int third = first + (segments + 1);
            int fourth = third + 1;
            
            // Два треугольника образуют квадрат
            planet.indices.push_back(first);
            planet.indices.push_back(second);
            planet.indices.push_back(third);
            
            planet.indices.push_back(second);
            planet.indices.push_back(fourth);
            planet.indices.push_back(third);
        }
    }
}void updatePlanet(Planet& planet, float deltaTime) {
    planet.orbitAngle += planet.orbitSpeed * deltaTime;
    planet.rotationAngle += planet.rotationSpeed * deltaTime;
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

void updatePlanet(Planet& planet, float deltaTime) {
    planet.orbitAngle += planet.orbitSpeed * deltaTime;
    planet.rotationAngle += planet.rotationSpeed * deltaTime;
}

void renderPlanet(const Planet& planet) {
    glPushMatrix();
    
    // Вращение по орбите
    glRotatef(planet.orbitAngle, 0.0f, 1.0f, 0.0f);
    glTranslatef(planet.orbitRadius, 0.0f, 0.0f);
    
    // Собственное вращение
    glRotatef(planet.rotationAngle, 0.0f, 1.0f, 0.0f);
    
    // Масштаб планеты
    glScalef(planet.size, planet.size, planet.size);
    
    drawPlanet(planet);
    
    glPopMatrix();
}

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    GLFWwindow* window = glfwCreateWindow(800, 600, "Solar System", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }
    
    Planet planets[3];
    
    createPlanet(planets[0], 0.3f, 2.0f, 50.0f, 100.0f, 0.3f);
    planets[0].vertices[0].r = 1.0f; // Красная планета
    
    // Планета 2: средняя, средняя скорость, средний размер
    createPlanet(planets[1], 0.4f, 3.5f, 30.0f, 70.0f, 0.4f);
    for(auto& v : planets[1].vertices) v.g = 1.0f; // Зеленая планета
    
    // Планета 3: дальняя, медленная, большая
    createPlanet(planets[2], 0.5f, 5.0f, 15.0f, 40.0f, 0.5f);
    for(auto& v : planets[2].vertices) v.b = 1.0f; // Синяя планета
    
    // Настройки OpenGL
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    
    float lastTime = glfwGetTime();
    
    // Главный цикл
    while (!glfwWindowShouldClose(window)) {
        // Вычисляем время между кадрами
        float currentTime = glfwGetTime();
        float deltaTime = currentTime - lastTime;
        lastTime = currentTime;
        
        // Очистка экрана
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Получаем размеры окна
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        
        // Настройка проекции
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)width / (float)height;
        gluPerspective(45.0f, aspect, 0.1f, 100.0f);
        
        // Настройка видовой матрицы
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        // Позиция камеры
        glTranslatef(0.0f, 0.0f, -10.0f);
        
        // Наклон для лучшего обзора
        glRotatef(25.0f, 1.0f, 0.0f, 0.0f);
        
        // Обновляем и рисуем все планеты
        for(int i = 0; i < 3; i++) {
            updatePlanet(planets[i], deltaTime);
            renderPlanet(planets[i]);
        }
        
        // Обмен буферов и обработка событий
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    // Очистка
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
glRotatef(25.0f, 1.0f, 0.0f, 0.0f);
glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
