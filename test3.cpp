#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <unistd.h>

struct Vertex {
    float x, y, z;  // OpenGL использует float
    float r, g, b;
};

struct Planet {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    float visual_radius;  // Для визуализации (float)
    double mass;          // Масса (double)
    double x, y, z;       // Позиция (double)
    double vx, vy, vz;    // Скорость (double)
    double ax, ay, az;    // Ускорение (double)
};

void create_sphere(Planet& planet, float visual_radius, float r, float g, float b) {
    planet.visual_radius = visual_radius;
    int s = 30;
    
    for (int i = 0; i <= s; ++i) {
        float phi = M_PI * i / s;
        for (int j = 0; j <= s; ++j) {
            float theta = 2.0f * M_PI * j / s;
            Vertex v;
            v.x = visual_radius * sin(phi) * cos(theta);
            v.y = visual_radius * cos(phi);
            v.z = visual_radius * sin(phi) * sin(theta);
            v.r = r;
            v.g = g;
            v.b = b;
            planet.vertices.push_back(v);
        }
    }
    
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            int first = i * (s + 1) + j;
            int second = first + 1;
            int third = first + (s + 1);
            int fourth = third + 1;
            planet.indices.push_back(first);
            planet.indices.push_back(second);
            planet.indices.push_back(third);
            planet.indices.push_back(second);
            planet.indices.push_back(fourth);
            planet.indices.push_back(third);
        }
    }
}

void draw_planet(const Planet& planet) {
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < planet.indices.size(); ++i) {
        const Vertex& v = planet.vertices[planet.indices[i]];
        glColor3f(v.r, v.g, v.b);
        glVertex3f(v.x, v.y, v.z);
    }
    glEnd();
}

void draw_coordinate() {
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();
}

// Используем double для точных расчетов
void update_physics(Planet& planet, const Planet& sun, double dt) {
    const double G = 1.0;  // Упрощенная гравитационная постоянная
    
    double dx = sun.x - planet.x;
    double dy = sun.y - planet.y;
    double dz = sun.z - planet.z;
    
    double r2 = dx*dx + dy*dy + dz*dz;
    
    // Защита от деления на ноль
    if (r2 < 1e-6) r2 = 1e-6;
    
    double r = sqrt(r2);
    
    // Ускорение от гравитации: a = G * M / r^2
    double a = G * sun.mass / r2;
    
    // Направление ускорения (к Солнцу)
    planet.ax = a * dx / r;
    planet.ay = a * dy / r;
    planet.az = a * dz / r;
    
    // Обновляем скорость: v = v0 + a * dt
    planet.vx += planet.ax * dt;
    planet.vy += planet.ay * dt;
    planet.vz += planet.az * dt;
    
    // Обновляем позицию: x = x0 + v * dt
    planet.x += planet.vx * dt;
    planet.y += planet.vy * dt;
    planet.z += planet.vz * dt;
}

// Преобразуем double в float для OpenGL
void render_planet(const Planet& planet) {
    glPushMatrix();
    // Кастим double в float для OpenGL
    glTranslatef(static_cast<float>(planet.x), 
                 static_cast<float>(planet.y), 
                 static_cast<float>(planet.z));
    draw_planet(planet);
    glPopMatrix();
}

int main() {
    if (!glfwInit()) {
        return -1;
    }
    
    GLFWwindow* window = glfwCreateWindow(800, 600, "Solar System (Double Precision)", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent(window);
    
    if (glewInit() != GLEW_OK) {
        return -1;
    }
    
    Planet sun, earth;
    
    // Солнце
    sun.mass = 100.0;       // double
    sun.x = 0.0;
    sun.y = 0.0;
    sun.z = 0.0;
    sun.vx = 0.0;
    sun.vy = 0.0;
    sun.vz = 0.0;
    sun.ax = 0.0;
    sun.ay = 0.0;
    sun.az = 0.0;
    create_sphere(sun, 0.5f, 1.0f, 1.0f, 0.0f);  // float для визуализации
    
    // Земля - начинаем ближе и с меньшей скоростью
    earth.mass = 1.0;
    earth.x = 2.0;         // Начинаем на расстоянии 2.0
    earth.y = 0.0;
    earth.z = 0.0;
    earth.vx = 0.0;
    earth.vy = 0.5;        // Начальная орбитальная скорость
    earth.vz = 0.0;
    earth.ax = 0.0;
    earth.ay = 0.0;
    earth.az = 0.0;
    create_sphere(earth, 0.2f, 0.0f, 0.5f, 1.0f);
    
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.1f, 1.0f);
    
    double last_time = glfwGetTime();
    double fixed_dt = 0.016;  // Фиксированный шаг времени (double)
    
    while (!glfwWindowShouldClose(window)) {
        double current_time = glfwGetTime();
        double delta_time = current_time - last_time;
        last_time = current_time;
        
        // Используем фиксированный шаг времени для стабильности
        update_physics(earth, sun, fixed_dt);
        
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
        
        // Камера
        gluLookAt(0.0f, 5.0f, 8.0f,
                  0.0f, 0.0f, 0.0f,
                  0.0f, 1.0f, 0.0f);
        
        draw_coordinate();
        
        // Солнце
        render_planet(sun);
        
        // Земля
        render_planet(earth);
        
        // Рисование орбиты
        glColor3f(0.3f, 0.3f, 0.3f);
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < 100; i++) {
            float angle = 2.0f * M_PI * i / 100.0f;
            glVertex3f(2.0f * cos(angle), 0.0f, 2.0f * sin(angle));
        }
        glEnd();
        
        // Вывод информации (double значения)
        std::cout << "Pos: (" << earth.x << ", " << earth.y 
                  << ") | Vel: (" << earth.vx << ", " << earth.vy 
                  << ") | Acc: (" << earth.ax << ", " << earth.ay << ")\n";
        
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            break;
        }
        
        // Уменьшил задержку для более плавного движения
        usleep(16000);  // ~60 FPS (16ms)
    }
    
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
