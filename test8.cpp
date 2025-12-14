#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <algorithm>
#include <map>

// Константы
const double AU = 1.496e11;
const double G = 6.67430e-11;
const double SCALE = 1.0 / AU * 30;
const double TIME_SCALE = 80000.0;

const double MASS_SUN = 1.989e30;
const double MASS_MERCURY = 3.285e23;
const double MASS_VENUS = 4.867e24;
const double MASS_EARTH = 5.972e24;
const double MASS_MARS = 6.39e23;
const double MASS_JUPITER = 1.898e27;
const double MASS_SATURN = 5.683e26;
const double MASS_MIPT_1 = 440075;

const double MIPT_1_ORBIT_RADIUS = 6.371e6 + 418200;
const double VEL_MIPT_1 = sqrt(G * MASS_EARTH / MIPT_1_ORBIT_RADIUS);

double VEL_MERCURY, VEL_VENUS, VEL_EARTH, VEL_MARS;

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

const float RADIUS_SUN = 0.4f;
const float RADIUS_MERCURY = 0.08f;
const float RADIUS_VENUS = 0.12f;
const float RADIUS_EARTH = 0.15f;
const float RADIUS_MARS = 0.10f;
const float RADIUS_JUPITER = 0.35f;
const float RADIUS_SATURN = 0.3f;
const float RADIUS_MIPT_1 = 0.1;

const float COLOR_SUN[3] = {1.0f, 1.0f, 0.0f};
const float COLOR_MERCURY[3] = {0.7f, 0.7f, 0.7f};
const float COLOR_VENUS[3] = {1.0f, 0.8f, 0.4f};
const float COLOR_EARTH[3] = {0.0f, 0.5f, 1.0f};
const float COLOR_MARS[3] = {1.0f, 0.3f, 0.1f};
const float COLOR_JUPITER[3] = {0.9f, 0.7f, 0.5f};
const float COLOR_SATURN[3] = {0.95f, 0.85f, 0.6f};
const float COLOR_MIPT_1[3] = {0.5f, 0.55f, 0.7f};

const double MASS_ASTEROID = 1.0e19;
const float RADIUS_ASTEROID = 0.05f;
const float COLOR_ASTEROID[3] = {0.8f, 0.4f, 0.1f};

// Структуры
struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct Camera {
    glm::vec3 position;
    glm::vec3 target;
    glm::vec3 up;
    float distance;
    float angleX;
    float angleY;
    
    Camera() : position(30.0f, 20.0f, 30.0f), 
               target(0.0f, 0.0f, 0.0f), 
               up(0.0f, 1.0f, 0.0f),
               distance(50.0f), 
               angleX(0.0f), 
               angleY(0.0f) {}
    
    void update() {
        position.x = target.x + distance * cos(angleY) * sin(angleX);
        position.y = target.y + distance * sin(angleY);
        position.z = target.z + distance * cos(angleY) * cos(angleX);
    }
    
    glm::mat4 getViewMatrix() {
        return glm::lookAt(position, target, up);
    }
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
    float a; // большая полуось
    float orbit_eccentricity; // эксцентриситет
    float orbit_tilt; // угол наклона орбиты
    
    // Для расчетов
    double last_perihelion_time; // время последнего прохождения перигелия
    double orbital_period; // период обращения (сек)
    double area_swept; // секториальная скорость (площадь/время)
    double last_area_calc_time;
    glm::vec3 last_position;
};

// Глобальные переменные
Camera camera;
bool camera_rotating = false;
double last_mouse_x = 0.0, last_mouse_y = 0.0;
double simulation_time = 0.0;
std::map<std::string, double> planet_periods;

// Функции
void vel_calculate(const std::vector<Planet*>& planets) {
    for(const auto& planet : planets) {
        if(planet->name == "Mercury") 
            VEL_MERCURY = sqrt(2 * ((-G * MASS_SUN / (2 * A_MERCURY)) + (G * MASS_SUN / PERIHELION_MERCURY)));
        else if(planet->name == "Venus") 
            VEL_VENUS = sqrt(2 * ((-G * MASS_SUN / (2 * A_VENUS)) + (G * MASS_SUN / PERIHELION_VENUS)));
        else if(planet->name == "Earth") 
            VEL_EARTH = sqrt(2 * ((-G * MASS_SUN / (2 * A_EARTH)) + (G * MASS_SUN / PERIHELION_EARTH)));
        else if(planet->name == "Mars") 
            VEL_MARS = sqrt(2 * ((-G * MASS_SUN / (2 * A_MARS)) + (G * MASS_SUN / PERIHELION_MARS)));
    }
}

void create_sphere(Planet& planet, float visual_radius, float r, float g, float b) {
    planet.visual_radius = visual_radius;
    int sectors = 30;
    int stacks = 30;
    
    for(int i = 0; i <= stacks; ++i) {
        float phi = M_PI * i / stacks;
        for(int j = 0; j <= sectors; ++j) {
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
    
    for(int i = 0; i < stacks; ++i) {
        for(int j = 0; j < sectors; ++j) {
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

void draw_planet(const Planet& planet) {
    glBegin(GL_TRIANGLES);
    for(size_t i = 0; i < planet.indices.size(); ++i) {
        const Vertex& v = planet.vertices[planet.indices[i]];
        glColor3f(v.r, v.g, v.b);
        glVertex3f(v.x, v.y, v.z);
    }
    glEnd();
}

void draw_orbit(float a, float eccentricity, float tilt, float r, float g, float b) {
    glColor3f(r, g, b);
    glBegin(GL_LINE_LOOP);
    for(int i = 0; i < 200; i++) {
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

void update_physics(Planet* planet, const std::vector<Planet*>& all_planets, double dt) {
    planet->vx += 0.5 * planet->ax * dt;
    planet->vy += 0.5 * planet->ay * dt;
    planet->vz += 0.5 * planet->az * dt;
    planet->x += planet->vx * dt;
    planet->y += planet->vy * dt;
    planet->z += planet->vz * dt;
    
    double total_ax = 0.0, total_ay = 0.0, total_az = 0.0;
    for(const auto& other : all_planets) {
        if(planet == other) continue;
        
        double dx = other->x - planet->x;
        double dy = other->y - planet->y;
        double dz = other->z - planet->z;
        
        double r2 = dx * dx + dy * dy + dz * dz;
        if(r2 < 1e-12) continue;
        
        double r = sqrt(r2);
        double a = G * other->mass / r2;
        
        total_ax += a * dx / r;
        total_ay += a * dy / r;
        total_az += a * dz / r;
    }
    
    planet->vx += 0.5 * total_ax * dt;
    planet->vy += 0.5 * total_ay * dt;
    planet->vz += 0.5 * total_az * dt;
    planet->ax = total_ax;
    planet->ay = total_ay;
    planet->az = total_az;
}

void render_planet(const Planet* planet) {
    glPushMatrix();
    glTranslatef(static_cast<float>(planet->x * SCALE),
                 static_cast<float>(planet->y * SCALE),
                 static_cast<float>(planet->z * SCALE));
    draw_planet(*planet);
    glPopMatrix();
}

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

void initialize_planets(std::vector<Planet*>& planets) {
    for(auto p : planets) delete p;
    planets.clear();
    
    Planet* sun = new Planet();
    sun->name = "Sun";
    sun->mass = MASS_SUN;
    sun->x = 0.0; sun->y = 0.0; sun->z = 0.0;
    sun->vx = 0.0; sun->vy = 0.0; sun->vz = 0.0;
    sun->ax = 0.0; sun->ay = 0.0; sun->az = 0.0;
    sun->orbit_eccentricity = 0.0f;
    sun->a = 0.0f;
    sun->orbit_tilt = 0.0f;
    sun->orbital_period = 0.0;
    sun->area_swept = 0.0;
    sun->last_perihelion_time = -1.0;
    sun->last_area_calc_time = 0.0;
    sun->last_position = glm::vec3(0.0f);
    create_sphere(*sun, RADIUS_SUN, COLOR_SUN[0], COLOR_SUN[1], COLOR_SUN[2]);
    planets.push_back(sun);
    
    Planet* mercury = new Planet();
    mercury->name = "Mercury";
    mercury->mass = MASS_MERCURY;
    mercury->x = PERIHELION_MERCURY; mercury->y = 0.0; mercury->z = 0.0;
    mercury->vx = 0.0; mercury->vy = 0.0; mercury->vz = VEL_MERCURY;
    mercury->ax = 0.0; mercury->ay = 0.0; mercury->az = 0.0;
    mercury->orbit_eccentricity = EX_MERCURY;
    mercury->a = A_MERCURY;
    mercury->orbit_tilt = TILT_MERCURY;
    mercury->orbital_period = 0.0;
    mercury->area_swept = 0.0;
    mercury->last_perihelion_time = 0.0;
    mercury->last_area_calc_time = 0.0;
    mercury->last_position = glm::vec3(mercury->x, mercury->y, mercury->z);
    create_sphere(*mercury, RADIUS_MERCURY, COLOR_MERCURY[0], COLOR_MERCURY[1], COLOR_MERCURY[2]);
    planets.push_back(mercury);
    
    Planet* venus = new Planet();
    venus->name = "Venus";
    venus->mass = MASS_VENUS;
    venus->x = PERIHELION_VENUS; venus->y = 0.0; venus->z = 0.0;
    venus->vx = 0.0; venus->vy = 0.0; venus->vz = VEL_VENUS;
    venus->ax = 0.0; venus->ay = 0.0; venus->az = 0.0;
    venus->orbit_eccentricity = EX_VENUS;
    venus->a = A_VENUS;
    venus->orbit_tilt = TILT_VENUS;
    venus->orbital_period = 0.0;
    venus->area_swept = 0.0;
    venus->last_perihelion_time = 0.0;
    venus->last_area_calc_time = 0.0;
    venus->last_position = glm::vec3(venus->x, venus->y, venus->z);
    create_sphere(*venus, RADIUS_VENUS, COLOR_VENUS[0], COLOR_VENUS[1], COLOR_VENUS[2]);
    planets.push_back(venus);
    
    Planet* earth = new Planet();
    earth->name = "Earth";
    earth->mass = MASS_EARTH;
    earth->x = PERIHELION_EARTH; earth->y = 0.0; earth->z = 0.0;
    earth->vx = 0.0; earth->vy = 0.0; earth->vz = VEL_EARTH;
    earth->ax = 0.0; earth->ay = 0.0; earth->az = 0.0;
    earth->orbit_eccentricity = EX_EARTH;
    earth->a = A_EARTH;
    earth->orbit_tilt = TILT_EARTH;
    earth->orbital_period = 0.0;
    earth->area_swept = 0.0;
    earth->last_perihelion_time = 0.0;
    earth->last_area_calc_time = 0.0;
    earth->last_position = glm::vec3(earth->x, earth->y, earth->z);
    create_sphere(*earth, RADIUS_EARTH, COLOR_EARTH[0], COLOR_EARTH[1], COLOR_EARTH[2]);
    planets.push_back(earth);
    
    Planet* mars = new Planet();
    mars->name = "Mars";
    mars->mass = MASS_MARS;
    mars->x = PERIHELION_MARS; mars->y = 0.0; mars->z = 0.0;
    mars->vx = 0.0; mars->vy = 0.0; mars->vz = VEL_MARS;
    mars->ax = 0.0; mars->ay = 0.0; mars->az = 0.0;
    mars->orbit_eccentricity = EX_MARS;
    mars->a = A_MARS;
    mars->orbit_tilt = TILT_MARS;
    mars->orbital_period = 0.0;
    mars->area_swept = 0.0;
    mars->last_perihelion_time = 0.0;
    mars->last_area_calc_time = 0.0;
    mars->last_position = glm::vec3(mars->x, mars->y, mars->z);
    create_sphere(*mars, RADIUS_MARS, COLOR_MARS[0], COLOR_MARS[1], COLOR_MARS[2]);
    planets.push_back(mars);
    
    Planet* mipt_1 = new Planet();
    mipt_1->name = "MIPT_1";
    mipt_1->mass = MASS_MIPT_1;
    mipt_1->x = earth->x + MIPT_1_ORBIT_RADIUS; mipt_1->y = 0.0; mipt_1->z = earth->z;
    mipt_1->vx = earth->vx; mipt_1->vy = VEL_MIPT_1; mipt_1->vz = VEL_MIPT_1 + earth->vz;
    mipt_1->ax = 0.0; mipt_1->ay = 0.0; mipt_1->az = 0.0;
    mipt_1->orbital_period = 0.0;
    mipt_1->area_swept = 0.0;
    mipt_1->last_perihelion_time = 0.0;
    mipt_1->last_area_calc_time = 0.0;
    mipt_1->last_position = glm::vec3(mipt_1->x, mipt_1->y, mipt_1->z);
    create_sphere(*mipt_1, RADIUS_MIPT_1, COLOR_MIPT_1[0], COLOR_MIPT_1[1], COLOR_MIPT_1[2]);
    planets.push_back(mipt_1);
}

void update_orbital_data(Planet* planet, double dt) {
    // Расчет секториальной скорости
    glm::vec3 current_pos(planet->x, planet->y, planet->z);
    glm::vec3 sun_pos(0.0f, 0.0f, 0.0f);
    
    // Вектор от Солнца к планете
    glm::vec3 r_vec = current_pos - sun_pos;
    glm::vec3 v_vec(planet->vx, planet->vy, planet->vz);
    
    // Секториальная скорость (площадь, заметаемая радиус-вектором в единицу времени)
    // dA/dt = 0.5 * |r x v|
    glm::vec3 cross_product = glm::cross(r_vec, v_vec);
    planet->area_swept = 0.5 * glm::length(cross_product);
    
    // Определение прохождения перигелия для расчета периода
    double distance_to_sun = glm::length(r_vec);
    double perihelion_distance = planet->a * (1.0 - planet->orbit_eccentricity);
    
    // Если близко к перигелию и скорость направлена правильно
    if (distance_to_sun <= perihelion_distance * 1.01 && 
        distance_to_sun >= perihelion_distance * 0.99) {
        
        // Проверяем, что радиальная скорость близка к 0 (в перигелии/афелии)
        glm::vec3 radial_dir = glm::normalize(r_vec);
        double radial_velocity = glm::dot(v_vec, radial_dir);
        
        if (fabs(radial_velocity) < 100 && planet->last_perihelion_time > 0) {
            double current_period = simulation_time - planet->last_perihelion_time;
            
            // Усредняем значение периода
            if (planet->orbital_period == 0.0) {
                planet->orbital_period = current_period;
            } else {
                planet->orbital_period = 0.9 * planet->orbital_period + 0.1 * current_period;
            }
            
            planet_periods[planet->name] = planet->orbital_period;
            planet->last_perihelion_time = simulation_time;
        } else if (planet->last_perihelion_time < 0) {
            planet->last_perihelion_time = simulation_time;
        }
    }
    
    planet->last_position = current_pos;
}

void print_planet_info(const std::vector<Planet*>& planets) {
    static int frame_count = 0;
    frame_count++;
    
    // Выводим информацию только каждые 100 кадров
    if (frame_count % 100 != 0) return;
    
    std::cout << "\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "                                  ПАРАМЕТРЫ ПЛАНЕТ (t = " << simulation_time << " с)\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════════════════\n";
    
    double kepler_constant_sum = 0.0;
    int planet_count = 0;
    
    for(const auto& planet : planets) {
        if (planet->name == "Sun" || planet->name == "MIPT_1") continue;
        
        std::cout << "\n[" << planet->name << "]\n";
        std::cout << "  Позиция: x=" << planet->x << " м, y=" << planet->y << " м, z=" << planet->z << " м\n";
        std::cout << "  Скорость: v=" << sqrt(planet->vx*planet->vx + planet->vy*planet->vy + planet->vz*planet->vz) << " м/с\n";
        
        // Расчет расстояния до Солнца
        double distance_to_sun = sqrt(planet->x*planet->x + planet->y*planet->y + planet->z*planet->z);
        std::cout << "  Расстояние до Солнца: " << distance_to_sun << " м (" << distance_to_sun/AU << " AU)\n";
        
        // Период обращения
        if (planet->orbital_period > 0) {
            double period_days = planet->orbital_period / (24.0 * 3600.0);
            std::cout << "  Период обращения: " << planet->orbital_period << " с (" << period_days << " дней)\n";
            
            // Теоретический период по 3-му закону Кеплера
            double theoretical_period = 2.0 * M_PI * sqrt(pow(planet->a, 3) / (G * MASS_SUN));
            double theoretical_days = theoretical_period / (24.0 * 3600.0);
            std::cout << "  Теоретический период (Кеплер): " << theoretical_period << " с (" << theoretical_days << " дней)\n";
            
            // Отношение T²/a³
            double T2_a3 = (planet->orbital_period * planet->orbital_period) / pow(planet->a, 3);
            double expected_T2_a3 = (4.0 * M_PI * M_PI) / (G * MASS_SUN);
            std::cout << "  T²/a³: " << T2_a3 << " (ожидается: " << expected_T2_a3 << ")\n";
            std::cout << "  Отклонение: " << fabs((T2_a3 - expected_T2_a3)/expected_T2_a3 * 100.0) << "%\n";
            
            kepler_constant_sum += T2_a3;
            planet_count++;
        }
        
        // Секториальная скорость
        std::cout << "  Секториальная скорость: " << planet->area_swept << " м²/с\n";
        
        // Орбитальные параметры
        std::cout << "  Орбитальные параметры:\n";
        std::cout << "    Большая полуось: " << planet->a << " м (" << planet->a/AU << " AU)\n";
        std::cout << "    Эксцентриситет: " << planet->orbit_eccentricity << "\n";
        std::cout << "    Наклон орбиты: " << planet->orbit_tilt * 180.0/M_PI << "°\n";
    }
    
    if (planet_count > 0) {
        double avg_kepler_constant = kepler_constant_sum / planet_count;
        double expected_constant = (4.0 * M_PI * M_PI) / (G * MASS_SUN);
        std::cout << "\n════════════════════════════════════════════════════════════════════════════════════════\n";
        std::cout << "  СРЕДНЕЕ T²/a³: " << avg_kepler_constant << "\n";
        std::cout << "  ОЖИДАЕМОЕ: " << expected_constant << "\n";
        std::cout << "  ПОГРЕШНОСТЬ: " << fabs((avg_kepler_constant - expected_constant)/expected_constant * 100.0) << "%\n";
        std::cout << "════════════════════════════════════════════════════════════════════════════════════════\n";
    }
    
    // Информация о камере
    std::cout << "\n[КАМЕРА]\n";
    std::cout << "  Позиция: (" << camera.position.x << ", " << camera.position.y << ", " << camera.position.z << ")\n";
    std::cout << "  Направление: (" << camera.target.x << ", " << camera.target.y << ", " << camera.target.z << ")\n";
    std::cout << "  Управление: ЛКМ + движение мыши - вращение, Прокрутка - приближение/отдаление\n";
}

// Callback функции для камеры
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            camera_rotating = true;
            glfwGetCursorPos(window, &last_mouse_x, &last_mouse_y);
        } else if (action == GLFW_RELEASE) {
            camera_rotating = false;
        }
    }
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    if (camera_rotating) {
        double dx = xpos - last_mouse_x;
        double dy = ypos - last_mouse_y;
        
        camera.angleX += dx * 0.01f;
        camera.angleY -= dy * 0.01f;
        
        // Ограничиваем угол Y чтобы не переворачивать камеру
camera.angleY = std::max<float>(-M_PI/2.0f + 0.1f, std::min<float>(M_PI/2.0f - 0.1f, camera.angleY));
        
        camera.update();
        
        last_mouse_x = xpos;
        last_mouse_y = ypos;
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    camera.distance -= yoffset * 2.0f;
    camera.distance = std::max(10.0f, std::min(200.0f, camera.distance));
    camera.update();
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        float camera_speed = 1.0f;
        
        if (key == GLFW_KEY_W) camera.target.y += camera_speed;
        if (key == GLFW_KEY_S) camera.target.y -= camera_speed;
        if (key == GLFW_KEY_A) camera.target.x -= camera_speed;
        if (key == GLFW_KEY_D) camera.target.x += camera_speed;
        if (key == GLFW_KEY_Q) camera.target.z -= camera_speed;
        if (key == GLFW_KEY_E) camera.target.z += camera_speed;
        
        if (key == GLFW_KEY_R) {
            // Сброс камеры
            camera.position = glm::vec3(30.0f, 20.0f, 30.0f);
            camera.target = glm::vec3(0.0f, 0.0f, 0.0f);
            camera.distance = 50.0f;
            camera.angleX = 0.0f;
            camera.angleY = 0.0f;
            camera.update();
        }
        
        if (key == GLFW_KEY_F) {
            // Переключение фокуса на Солнце
            camera.target = glm::vec3(0.0f, 0.0f, 0.0f);
            camera.update();
        }
        
        camera.update();
    }
}

int main() {
    if(!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    
    GLFWwindow* window = glfwCreateWindow(1200, 800, "Solar System Simulation", NULL, NULL);
    if(!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent(window);
    
    if(glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }
    
    // Устанавливаем callback функции
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);
    
    std::vector<Planet*> planets;
    initialize_planets(planets);
    
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.05f, 1.0f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    vel_calculate(planets);
    for(auto planet : planets) {
        if(planet->name == "Mercury") planet->vz = VEL_MERCURY;
        else if(planet->name == "Venus") planet->vz = VEL_VENUS;
        else if(planet->name == "Earth") planet->vz = VEL_EARTH;
        else if(planet->name == "Mars") planet->vz = VEL_MARS;
    }
    
    double last_time = glfwGetTime();
    camera.update();
    
    while(!glfwWindowShouldClose(window)) {
        double current_time = glfwGetTime();
        double delta_time = current_time - last_time;
        last_time = current_time;
        
        double dt = delta_time * TIME_SCALE;
        simulation_time += dt;
        
        // Обновляем физику
        for(auto& planet : planets) {
            update_physics(planet, planets, dt);
        }
        
        // Обновляем орбитальные данные
        for(auto& planet : planets) {
            if(planet->name != "Sun") {
                update_orbital_data(planet, dt);
            }
        }
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)width / (float)height;
        gluPerspective(60.0f, aspect, 0.1f, 1000.0f);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        // Используем матрицу вида из камеры
        glm::mat4 view = camera.getViewMatrix();
        glLoadMatrixf(glm::value_ptr(view));
        
        draw_coordinate_system();
        
        glDisable(GL_DEPTH_TEST);
        glLineWidth(0.5f);
        
        // Рисуем орбиты
        for(const auto& planet : planets) {
            if(planet->name == "Mercury") 
                draw_orbit(planet->a, planet->orbit_eccentricity, TILT_MERCURY, 0.5f, 0.5f, 0.5f);
            else if(planet->name == "Venus") 
                draw_orbit(planet->a, planet->orbit_eccentricity, TILT_VENUS, 1.0f, 0.8f, 0.4f);
            else if(planet->name == "Earth") 
                draw_orbit(planet->a, planet->orbit_eccentricity, TILT_EARTH, 0.3f, 0.3f, 0.8f);
            else if(planet->name == "Mars") 
                draw_orbit(planet->a, planet->orbit_eccentricity, TILT_MARS, 0.8f, 0.3f, 0.2f);
        }
        
        glEnable(GL_DEPTH_TEST);
        glLineWidth(1.0f);
        
        // Рисуем планеты
        for(const auto& planet : planets) {
            render_planet(planet);
        }
        
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        print_planet_info(planets);
        usleep(20000);
    }
    
    for(auto p : planets) delete p;
    glfwDestroyWindow(window);
    glfwTerminate();
    
    std::cout << std::endl << "Симуляция завершена" << std::endl;
    return 0;
}
