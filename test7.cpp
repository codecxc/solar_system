#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <algorithm>
#include <string>

// Константы
const double AU = 1.496e11;
const double G = 6.67430e-11;
const double SCALE = 1.0 / AU * 15.0;
const double TIME_SCALE = 86400.0 * 10;

// Массы
const double MASS_SUN = 1.989e30;
const double MASS_MERCURY = 3.285e23;
const double MASS_VENUS = 4.867e24;
const double MASS_EARTH = 5.972e24;
const double MASS_MARS = 6.39e23;
const double MASS_JUPITER = 1.898e27;
const double MASS_SATURN = 5.683e26;
const double MASS_MIPT_1 = 1.0e3; // Масса спутника (1000 кг)

// Орбитальные параметры
const float EX_SUN = 0.0f;
const float EX_MERCURY = 0.206f;
const float EX_VENUS = 0.007f;
const float EX_EARTH = 0.017f;
const float EX_MARS = 0.093f;
const float EX_JUPITER = 0.049f;
const float EX_SATURN = 0.057f;

const float A_SUN = 0.0f;
const float A_MERCURY = 0.387 * AU;
const float A_VENUS = 0.723 * AU;
const float A_EARTH = 1.000 * AU;
const float A_MARS = 1.524 * AU;
const float A_JUPITER = 5.2044 * AU;
const float A_SATURN = 9.5826 * AU;

const float TILT_MERCURY = 7.01f * M_PI / 180.0f;
const float TILT_VENUS = 3.39f * M_PI / 180.0f;
const float TILT_EARTH = 0.0f;
const float TILT_MARS = 1.85f * M_PI / 180.0f;
const float TILT_JUPITER = 1.31f * M_PI / 180.0f;
const float TILT_SATURN = 2.49f * M_PI / 180.0f;

const double PERIHELION_MERCURY = A_MERCURY * (1 - EX_MERCURY);
const double PERIHELION_VENUS = A_VENUS * (1 - EX_VENUS);
const double PERIHELION_EARTH = A_EARTH * (1 - EX_EARTH);
const double PERIHELION_MARS = A_MARS * (1 - EX_MARS);

// Визуальные радиусы
const float RADIUS_SUN = 0.4f;
const float RADIUS_MERCURY = 0.08f;
const float RADIUS_VENUS = 0.12f;
const float RADIUS_EARTH = 0.15f;
const float RADIUS_MARS = 0.10f;
const float RADIUS_JUPITER = 0.35f;
const float RADIUS_SATURN = 0.3f;
const float RADIUS_MIPT_1 = 0.05f;

// Цвета
const float COLOR_SUN[3] = {1.0f, 1.0f, 0.0f};
const float COLOR_MERCURY[3] = {0.7f, 0.7f, 0.7f};
const float COLOR_VENUS[3] = {1.0f, 0.8f, 0.4f};
const float COLOR_EARTH[3] = {0.0f, 0.5f, 1.0f};
const float COLOR_MARS[3] = {1.0f, 0.3f, 0.1f};
const float COLOR_JUPITER[3] = {0.9f, 0.7f, 0.5f};
const float COLOR_SATURN[3] = {0.95f, 0.85f, 0.6f};
const float COLOR_MIPT_1[3] = {1.0f, 0.0f, 0.0f};

// Параметры орбиты спутника вокруг Земли
const double MIPT_1_ORBIT_RADIUS = 6.371e6 + 4.0e5; // 400 км над поверхностью Земли
const double MIPT_1_ORBITAL_VELOCITY = sqrt(G * MASS_EARTH / MIPT_1_ORBIT_RADIUS);

// Переменные для скоростей планет
double VEL_MERCURY, VEL_VENUS, VEL_EARTH, VEL_MARS;

struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct Planet {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    float visual_radius;
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    std::string name;
    float a;
    float orbit_eccentricity;
    float orbit_tilt;
    
    // Для вычисления орбитальных параметров
    double last_angle;
    double period_start_time;
    double sectorial_velocity;
    double angular_momentum;
    double semi_major_axis;
    double period;
    bool period_measured;
    int orbit_count;
    double last_period_check_time;
    
    // Для спутника - родительская планета
    Planet* parent;
};

// Вычисление начальных скоростей планет
void vel_calculate(const std::vector<Planet*>& planets) {
    for (const auto& planet : planets) {
        if (planet->name == "Mercury")
            VEL_MERCURY = sqrt(2 * ((-G * MASS_SUN / (2 * A_MERCURY)) + (G * MASS_SUN / PERIHELION_MERCURY)));
        else if (planet->name == "Venus")
            VEL_VENUS = sqrt(2 * ((-G * MASS_SUN / (2 * A_VENUS)) + (G * MASS_SUN / PERIHELION_VENUS)));
        else if (planet->name == "Earth")
            VEL_EARTH = sqrt(2 * ((-G * MASS_SUN / (2 * A_EARTH)) + (G * MASS_SUN / PERIHELION_EARTH)));
        else if (planet->name == "Mars")
            VEL_MARS = sqrt(2 * ((-G * MASS_SUN / (2 * A_MARS)) + (G * MASS_SUN / PERIHELION_MARS)));
    }
}

// Создание сферы
void create_sphere(Planet& planet, float visual_radius, float r, float g, float b) {
    planet.visual_radius = visual_radius;
    int sectors = 30;
    int stacks = 30;
    
    for (int i = 0; i <= stacks; ++i) {
        float phi = M_PI * i / stacks;
        for (int j = 0; j <= sectors; ++j) {
            float theta = 2.0f * M_PI * j / sectors;
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
    
    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < sectors; ++j) {
            int first = i * (sectors + 1) + j;
            int second = first + 1;
            int third = first + (sectors + 1);
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

// Отрисовка сферы
void draw_planet(const Planet& planet) {
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < planet.indices.size(); ++i) {
        const Vertex& v = planet.vertices[planet.indices[i]];
        glColor3f(v.r, v.g, v.b);
        glVertex3f(v.x, v.y, v.z);
    }
    glEnd();
}

// Отрисовка орбиты
void draw_orbit(float a, float eccentricity, float tilt, float r, float g, float b) {
    glColor3f(r, g, b);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 200; i++) {
        float angle = 2.0f * M_PI * i / 200.0f;
        float r_orbit = a * (1 - eccentricity * eccentricity) / (1 + eccentricity * cos(angle));
        
        r_orbit *= SCALE;
        float x = r_orbit * cos(angle);
        float z = r_orbit * sin(angle);
        float y = 0.0f;
        
        float y_tilted = z * sin(tilt);
        float z_tilted = z * cos(tilt);
        glVertex3f(x, y_tilted, z_tilted);
    }
    glEnd();
}

// Отрисовка орбиты спутника вокруг Земли
void draw_satellite_orbit(const Planet& earth, const Planet& satellite, float r, float g, float b) {
    glColor3f(r, g, b);
    glBegin(GL_LINE_LOOP);
    int segments = 100;
    
    // Рисуем орбиту в системе координат Земли
    for (int i = 0; i < segments; i++) {
        float angle = 2.0f * M_PI * i / segments;
        float orbit_radius = MIPT_1_ORBIT_RADIUS * SCALE;
        
        // Координаты в системе Земли
        float local_x = orbit_radius * cos(angle);
        float local_z = orbit_radius * sin(angle);
        
        // Преобразуем в глобальные координаты
        float global_x = (earth.x + local_x / SCALE) * SCALE;
        float global_z = (earth.z + local_z / SCALE) * SCALE;
        float global_y = earth.y * SCALE;
        
        glVertex3f(global_x, global_y, global_z);
    }
    glEnd();
}

// Обновление физики
void update_physics(Planet* planet, const std::vector<Planet*>& all_planets, double dt, double simulation_time) {
    double total_ax = 0.0, total_ay = 0.0, total_az = 0.0;
    
    for (const auto& other : all_planets) {
        if (planet == other) continue;
        
        double dx = other->x - planet->x;
        double dy = other->y - planet->y;
        double dz = other->z - planet->z;
        
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < 1e-12) continue;
        
        double r = sqrt(r2);
        double a = G * other->mass / r2;
        
        total_ax += a * dx / r;
        total_ay += a * dy / r;
        total_az += a * dz / r;
    }
    
    // Сохраняем старые скорости для интегрирования
    double old_vx = planet->vx;
    double old_vy = planet->vy;
    double old_vz = planet->vz;
    
    // Обновляем скорости
    planet->vx += total_ax * dt;
    planet->vy += total_ay * dt;
    planet->vz += total_az * dt;
    
    // Обновляем позиции (используем среднюю скорость для лучшей точности)
    planet->x += (old_vx + planet->vx) * 0.5 * dt;
    planet->y += (old_vy + planet->vy) * 0.5 * dt;
    planet->z += (old_vz + planet->vz) * 0.5 * dt;
    
    planet->ax = total_ax;
    planet->ay = total_ay;
    planet->az = total_az;
    
    // Вычисление орбитальных параметров
    if (planet->name != "Sun" && planet->parent == nullptr) {
        // Для планет вычисляем относительно Солнца
        double rx = planet->x;
        double ry = planet->y;
        double rz = planet->z;
        
        double vx = planet->vx;
        double vy = planet->vy;
        double vz = planet->vz;
        
        // Угловой момент
        double Lx = ry * vz - rz * vy;
        double Ly = rz * vx - rx * vz;
        double Lz = rx * vy - ry * vx;
        planet->angular_momentum = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
        planet->sectorial_velocity = planet->angular_momentum / (2.0 * planet->mass);
        
        // Большая полуось из энергии
        double r = sqrt(rx * rx + ry * ry + rz * rz);
        double v2 = vx * vx + vy * vy + vz * vz;
        double energy = v2 / 2.0 - G * MASS_SUN / r;
        if (energy < 0) {
            planet->semi_major_axis = -G * MASS_SUN / (2.0 * energy);
        }
        
        // Определение периода
        if (!planet->period_measured) {
            double current_angle = atan2(planet->z, planet->x);
            
            if (planet->orbit_count == 0) {
                planet->last_angle = current_angle;
                planet->period_start_time = simulation_time;
                planet->orbit_count = 1;
            } else {
                // Отслеживаем изменения угла
                double angle_diff = current_angle - planet->last_angle;
                
                // Корректируем переход через 0
                if (angle_diff > M_PI) angle_diff -= 2 * M_PI;
                if (angle_diff < -M_PI) angle_diff += 2 * M_PI;
                
                // Если накопилось примерно 2π, значит прошел период
                static double accumulated_angle = 0;
                accumulated_angle += fabs(angle_diff);
                
                if (accumulated_angle > 6.0) { // Немного меньше 2π для надежности
                    planet->period = simulation_time - planet->period_start_time;
                    planet->period_measured = true;
                    accumulated_angle = 0;
                }
                
                planet->last_angle = current_angle;
            }
        }
    }
}

// Рендер планеты
void render_planet(const Planet* planet) {
    glPushMatrix();
    glTranslatef(static_cast<float>(planet->x * SCALE),
                 static_cast<float>(planet->y * SCALE),
                 static_cast<float>(planet->z * SCALE));
    draw_planet(*planet);
    glPopMatrix();
}

// Отрисовка координатной системы
void draw_coordinate_system() {
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(2.0f, 0.0f, 0.0f);
    
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 2.0f, 0.0f);
    
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 2.0f);
    
    glEnd();
    glLineWidth(1.0f);
}

// Инициализация планет
void initialize_planets(std::vector<Planet*>& planets) {
    for (auto p : planets) delete p;
    planets.clear();
    
    // Солнце
    Planet* sun = new Planet();
    sun->name = "Sun";
    sun->mass = MASS_SUN;
    sun->x = 0.0; sun->y = 0.0; sun->z = 0.0;
    sun->vx = 0.0; sun->vy = 0.0; sun->vz = 0.0;
    sun->ax = 0.0; sun->ay = 0.0; sun->az = 0.0;
    sun->orbit_eccentricity = 0.0f;
    sun->a = 0.0f;
    sun->orbit_tilt = 0.0f;
    sun->sectorial_velocity = 0.0;
    sun->angular_momentum = 0.0;
    sun->semi_major_axis = 0.0;
    sun->period = 0.0;
    sun->period_measured = false;
    sun->orbit_count = 0;
    sun->parent = nullptr;
    create_sphere(*sun, RADIUS_SUN, COLOR_SUN[0], COLOR_SUN[1], COLOR_SUN[2]);
    planets.push_back(sun);
    
    // Меркурий
    Planet* mercury = new Planet();
    mercury->name = "Mercury";
    mercury->mass = MASS_MERCURY;
    mercury->x = PERIHELION_MERCURY; mercury->y = 0.0; mercury->z = 0.0;
    mercury->vx = 0.0; mercury->vy = 0.0; mercury->vz = 0.0; // Вычислится ниже
    mercury->ax = 0.0; mercury->ay = 0.0; mercury->az = 0.0;
    mercury->orbit_eccentricity = EX_MERCURY;
    mercury->a = A_MERCURY;
    mercury->orbit_tilt = TILT_MERCURY;
    mercury->sectorial_velocity = 0.0;
    mercury->angular_momentum = 0.0;
    mercury->semi_major_axis = 0.0;
    mercury->period = 0.0;
    mercury->period_measured = false;
    mercury->orbit_count = 0;
    mercury->parent = nullptr;
    create_sphere(*mercury, RADIUS_MERCURY, COLOR_MERCURY[0], COLOR_MERCURY[1], COLOR_MERCURY[2]);
    planets.push_back(mercury);
    
    // Венера
    Planet* venus = new Planet();
    venus->name = "Venus";
    venus->mass = MASS_VENUS;
    venus->x = PERIHELION_VENUS; venus->y = 0.0; venus->z = 0.0;
    venus->vx = 0.0; venus->vy = 0.0; venus->vz = 0.0; // Вычислится ниже
    venus->ax = 0.0; venus->ay = 0.0; venus->az = 0.0;
    venus->orbit_eccentricity = EX_VENUS;
    venus->a = A_VENUS;
    venus->orbit_tilt = TILT_VENUS;
    venus->sectorial_velocity = 0.0;
    venus->angular_momentum = 0.0;
    venus->semi_major_axis = 0.0;
    venus->period = 0.0;
    venus->period_measured = false;
    venus->orbit_count = 0;
    venus->parent = nullptr;
    create_sphere(*venus, RADIUS_VENUS, COLOR_VENUS[0], COLOR_VENUS[1], COLOR_VENUS[2]);
    planets.push_back(venus);
    
    // Земля
    Planet* earth = new Planet();
    earth->name = "Earth";
    earth->mass = MASS_EARTH;
    earth->x = PERIHELION_EARTH; earth->y = 0.0; earth->z = 0.0;
    earth->vx = 0.0; earth->vy = 0.0; earth->vz = 0.0; // Вычислится ниже
    earth->ax = 0.0; earth->ay = 0.0; earth->az = 0.0;
    earth->orbit_eccentricity = EX_EARTH;
    earth->a = A_EARTH;
    earth->orbit_tilt = TILT_EARTH;
    earth->sectorial_velocity = 0.0;
    earth->angular_momentum = 0.0;
    earth->semi_major_axis = 0.0;
    earth->period = 0.0;
    earth->period_measured = false;
    earth->orbit_count = 0;
    earth->parent = nullptr;
    create_sphere(*earth, RADIUS_EARTH, COLOR_EARTH[0], COLOR_EARTH[1], COLOR_EARTH[2]);
    planets.push_back(earth);
    
    // Марс
    Planet* mars = new Planet();
    mars->name = "Mars";
    mars->mass = MASS_MARS;
    mars->x = PERIHELION_MARS; mars->y = 0.0; mars->z = 0.0;
    mars->vx = 0.0; mars->vy = 0.0; mars->vz = 0.0; // Вычислится ниже
    mars->ax = 0.0; mars->ay = 0.0; mars->az = 0.0;
    mars->orbit_eccentricity = EX_MARS;
    mars->a = A_MARS;
    mars->orbit_tilt = TILT_MARS;
    mars->sectorial_velocity = 0.0;
    mars->angular_momentum = 0.0;
    mars->semi_major_axis = 0.0;
    mars->period = 0.0;
    mars->period_measured = false;
    mars->orbit_count = 0;
    mars->parent = nullptr;
    create_sphere(*mars, RADIUS_MARS, COLOR_MARS[0], COLOR_MARS[1], COLOR_MARS[2]);
    planets.push_back(mars);
    
    // Спутник MIPT_1 (на круговой орбите вокруг Земли)
    Planet* mipt_1 = new Planet();
    mipt_1->name = "MIPT_1";
    mipt_1->mass = MASS_MIPT_1;
    
    // ПРАВИЛЬНЫЕ НАЧАЛЬНЫЕ УСЛОВИЯ ДЛЯ ОРБИТАЛЬНОГО ДВИЖЕНИЯ ВОКРУГ ЗЕМЛИ:
    // 1. Расположим спутник на расстоянии MIPT_1_ORBIT_RADIUS от Земли
    // 2. Зададим ему перпендикулярную скорость для круговой орбиты
    
    // Начальное положение: справа от Земли (вдоль оси X)
    mipt_1->x = PERIHELION_EARTH + MIPT_1_ORBIT_RADIUS;
    mipt_1->y = 0.0;
    mipt_1->z = 0.0;
    
    // Начальная скорость: перпендикулярно радиус-вектору (вдоль оси Z)
    // Это создаст круговую орбиту в плоскости XZ
    mipt_1->vx = 0.0;  // Скорость Земли по X = 0
    mipt_1->vy = 0.0;
    
    // Скорость спутника = скорость Земли + орбитальная скорость вокруг Земли
    // Орбитальная скорость направлена перпендикулярно радиус-вектору
    mipt_1->vz = VEL_EARTH + MIPT_1_ORBITAL_VELOCITY;
    
    mipt_1->ax = 0.0; mipt_1->ay = 0.0; mipt_1->az = 0.0;
    mipt_1->sectorial_velocity = 0.0;
    mipt_1->angular_momentum = 0.0;
    mipt_1->semi_major_axis = 0.0;
    mipt_1->period = 0.0;
    mipt_1->period_measured = false;
    mipt_1->orbit_count = 0;
    mipt_1->parent = earth;  // Указываем родительскую планету
    
    create_sphere(*mipt_1, RADIUS_MIPT_1, COLOR_MIPT_1[0], COLOR_MIPT_1[1], COLOR_MIPT_1[2]);
    planets.push_back(mipt_1);
}

// Вывод информации о планетах
void print_planet_info(const std::vector<Planet*>& planets, double simulation_time) {
    static int counter = 0;
    static double last_print_time = 0;
    
    // Выводим информацию каждые 5 секунд симуляции
    if (simulation_time - last_print_time < 86400.0 * 30) return; // Каждые 30 дней
    
    last_print_time = simulation_time;
    counter++;
    
    std::cout << "\n════════════════════════════════════════════════════════════════════\n";
    std::cout << " ВРЕМЯ СИМУЛЯЦИИ: " << simulation_time / 86400.0 << " дней\n";
    std::cout << "════════════════════════════════════════════════════════════════════\n";
    
    // Находим Землю и спутник
    Planet* earth = nullptr;
    Planet* mipt_1 = nullptr;
    
    for (const auto& planet : planets) {
        if (planet->name == "Earth") earth = planet;
        if (planet->name == "MIPT_1") mipt_1 = planet;
    }
    
    for (const auto& planet : planets) {
        if (planet->name == "Sun") continue; // Пропускаем Солнце для краткости
        
        std::cout << "\n[" << planet->name << "]\n";
        std::cout << "  Позиция: x=" << planet->x << " м, y=" << planet->y << " м, z=" << planet->z << " м\n";
        std::cout << "  Скорость: v=" << sqrt(planet->vx*planet->vx + planet->vy*planet->vy + planet->vz*planet->vz) << " м/с\n";
        
        if (planet->parent == nullptr && planet->name != "Sun") {
            // Для планет вычисляем относительно Солнца
            double distance_from_sun = sqrt(planet->x*planet->x + planet->y*planet->y + planet->z*planet->z);
            std::cout << "  Расстояние от Солнца: " << distance_from_sun << " м\n";
            
            if (planet->period_measured) {
                std::cout << "  Секториальная скорость: " << planet->sectorial_velocity << " м²/с\n";
                std::cout << "  Большая полуось (a): " << planet->semi_major_axis << " м\n";
                std::cout << "  Период (T): " << planet->period / 86400.0 << " дней\n";
                
                if (planet->semi_major_axis > 0) {
                    double t2_a3 = (planet->period * planet->period) / 
                                  (planet->semi_major_axis * planet->semi_major_axis * planet->semi_major_axis);
                    std::cout << "  T²/a³: " << t2_a3 << " с²/м³\n";
                    
                    // Теоретическое значение по третьему закону Кеплера
                    double theoretical = (4.0 * M_PI * M_PI) / (G * MASS_SUN);
                    std::cout << "  Теоретическое T²/a³: " << theoretical << " с²/м³\n";
                }
            }
        }
        
        // Для спутника вычисляем относительно Земли
        if (planet->name == "MIPT_1" && earth) {
            std::cout << "\n  [ОТНОСИТЕЛЬНО ЗЕМЛИ]\n";
            
            // Относительные координаты и скорость
            double dx = planet->x - earth->x;
            double dy = planet->y - earth->y;
            double dz = planet->z - earth->z;
            double distance = sqrt(dx*dx + dy*dy + dz*dz);
            
            double dvx = planet->vx - earth->vx;
            double dvy = planet->vy - earth->vy;
            double dvz = planet->vz - earth->vz;
            double rel_velocity = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            
            std::cout << "  Расстояние до Земли: " << distance << " м\n";
            std::cout << "  Относительная скорость: " << rel_velocity << " м/с\n";
            std::cout << "  Теоретическая орбитальная скорость: " << MIPT_1_ORBITAL_VELOCITY << " м/с\n";
            
            // Угловой момент относительно Земли
            double Lx = dy * dvz - dz * dvy;
            double Ly = dz * dvx - dx * dvz;
            double Lz = dx * dvy - dy * dvx;
            double angular_momentum = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
            double sectorial_velocity = angular_momentum / (2.0 * planet->mass);
            
            std::cout << "  Секториальная скорость (отн. Земли): " << sectorial_velocity << " м²/с\n";
            
            // Орбитальный период вокруг Земли
            double orbital_period = 2.0 * M_PI * distance / rel_velocity;
            std::cout << "  Орбитальный период: " << orbital_period / 3600.0 << " часов\n";
        }
    }
    std::cout << "════════════════════════════════════════════════════════════════════\n\n";
}

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    
    GLFWwindow* window = glfwCreateWindow(1200, 800, "Solar System Simulation with Satellite", NULL, NULL);
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
    
    std::vector<Planet*> planets;
    initialize_planets(planets);
    
    // Вычисляем скорости планет
    vel_calculate(planets);
    
    // Устанавливаем вычисленные скорости для планет
    for (auto planet : planets) {
        if (planet->name == "Mercury") planet->vz = VEL_MERCURY;
        else if (planet->name == "Venus") planet->vz = VEL_VENUS;
        else if (planet->name == "Earth") planet->vz = VEL_EARTH;
        else if (planet->name == "Mars") planet->vz = VEL_MARS;
    }
    
    // Для спутника устанавливаем правильную скорость
    for (auto planet : planets) {
        if (planet->name == "MIPT_1") {
            // Ищем Землю
            Planet* earth = nullptr;
            for (auto p : planets) {
                if (p->name == "Earth") {
                    earth = p;
                    break;
                }
            }
            
            if (earth) {
                // Правильная скорость для орбитального движения вокруг Земли
                planet->vx = earth->vx;  // Совпадает со скоростью Земли по X
                planet->vy = earth->vy;  // Совпадает со скоростью Земли по Y
                planet->vz = earth->vz + MIPT_1_ORBITAL_VELOCITY;  // Плюс орбитальная скорость
            }
        }
    }
    
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.05f, 1.0f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    double last_time = glfwGetTime();
    double simulation_time = 0.0;
    
    while (!glfwWindowShouldClose(window)) {
        double current_time = glfwGetTime();
        double delta_time = current_time - last_time;
        last_time = current_time;
        
        double dt = delta_time * TIME_SCALE;
        simulation_time += dt;
        
        // Обновляем физику
        for (auto& planet : planets) {
            update_physics(planet, planets, dt, simulation_time);
        }
        
        // Очистка экрана
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Настройка видового экрана
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)width / (float)height;
        gluPerspective(60.0f, aspect, 0.1f, 100.0f);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(30.0f, 20.0f, 30.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
        
        // Рисуем координатную систему
        draw_coordinate_system();
        
        // Рисуем орбиты планет
        glDisable(GL_DEPTH_TEST);
        glLineWidth(0.5f);
        for (const auto& planet : planets) {
            if (planet->name == "Mercury") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_MERCURY, 0.5f, 0.5f, 0.5f);
            else if (planet->name == "Venus") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_VENUS, 1.0f, 0.8f, 0.4f);
            else if (planet->name == "Earth") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_EARTH, 0.3f, 0.3f, 0.8f);
            else if (planet->name == "Mars") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_MARS, 0.8f, 0.3f, 0.2f);
        }
        glEnable(GL_DEPTH_TEST);
        glLineWidth(1.0f);
        
        // Рисуем планеты
        for (const auto& planet : planets) {
            render_planet(planet);
        }
        
        // Рисуем орбиту спутника вокруг Земли (если визуально нужно)
        glDisable(GL_DEPTH_TEST);
        glLineWidth(0.3f);
        Planet* earth = nullptr;
        Planet* mipt_1 = nullptr;
        for (const auto& planet : planets) {
            if (planet->name == "Earth") earth = planet;
            if (planet->name == "MIPT_1") mipt_1 = planet;
        }
        if (earth && mipt_1) {
            draw_satellite_orbit(*earth, *mipt_1, 1.0f, 0.0f, 0.0f);
        }
        glEnable(GL_DEPTH_TEST);
        glLineWidth(1.0f);
        
        // Вывод информации
        print_planet_info(planets, simulation_time);
        
        // Обмен буферов и обработка событий
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        usleep(20000);
    }
    
    // Очистка
    for (auto p : planets) delete p;
    glfwDestroyWindow(window);
    glfwTerminate();
    
    std::cout << std::endl << "Симуляция завершена" << std::endl;
    return 0;
}
