#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>

// ==================== КОНСТАНТЫ ====================
const double AU = 1.496e11;
const double G = 6.67430e-11;
const double SCALE = 2500.0 / AU;  // Увеличил масштаб для лучшей видимости
const double TIME_SCALE = 86400.0;
const double RADIUS_SCALE = 2.0;  // Увеличил масштаб радиусов
const unsigned int TRIAL_STEPS = 50000;

// Массы
const double MASS_SUN = 1.9891e30;
const double MASS_MERCURY = 3.302e23;
const double MASS_VENUS = 4.868e24;
const double MASS_EARTH = 5.9736e24;
const double MASS_MARS = 6.4185e23;
const double MASS_JUPITER = 1.8986e27;
const double MASS_SATURN = 5.68460e26;
const double MASS_URAN = 8.6832e25;
const double MASS_NEPTUN = 1.02430e26;
const double MASS_LUNA = 7.342e22;
const double MASS_IO = 8.9319e22;
const double MASS_EUROPA = 4.7998e22;
const double MASS_GANYMEDE = 1.4819e23;
const double MASS_CALLISTO = 1.0759e23;
const double MASS_TITAN = 1.3452e23;
const double MASS_RHEA = 2.306e21;
const double MASS_TETHYS = 6.174e20;
const double MASS_DIONE = 1.095e21;
const double MASS_TITANIA = 3.527e21;
const double MASS_OBERON = 3.014e21;
const double MASS_TRITON = 2.14e22;
const double MASS_PHOBOS = 1.0659e16;
const double MASS_DEIMOS = 1.4762e15;

// Орбитальные параметры
const float EX_MERCURY = 0.206f;
const float EX_VENUS = 0.007f;
const float EX_EARTH = 0.017f;
const float EX_MARS = 0.093f;
const float EX_JUPITER = 0.049f;
const float EX_SATURN = 0.057f;
const float EX_URAN = 0.046f;
const float EX_NEPTUN = 0.011f;

const float A_MERCURY = 0.387f * AU;
const float A_VENUS = 0.723f * AU;
const float A_EARTH = 1.000f * AU;
const float A_MARS = 1.524f * AU;
const float A_JUPITER = 5.2044f * AU;
const float A_SATURN = 9.5826f * AU;
const float A_URAN = 19.21840f * AU;
const float A_NEPTUN = 30.11000f * AU;

// Наклоны орбит
const float TILT_MERCURY = 7.01f * M_PI / 180.0f;
const float TILT_VENUS = 3.39f * M_PI / 180.0f;
const float TILT_EARTH = 0.0f;
const float TILT_MARS = 1.85f * M_PI / 180.0f;
const float TILT_JUPITER = 1.31f * M_PI / 180.0f;
const float TILT_SATURN = 2.49f * M_PI / 180.0f;
const float TILT_URAN = 0.77f * M_PI / 180.0f;
const float TILT_NEPTUN = 1.77f * M_PI / 180.0f;

// Перигелии
const double PERIHELION_MERCURY = A_MERCURY * (1 - EX_MERCURY);
const double PERIHELION_VENUS = A_VENUS * (1 - EX_VENUS);
const double PERIHELION_EARTH = A_EARTH * (1 - EX_EARTH);
const double PERIHELION_MARS = A_MARS * (1 - EX_MARS);
const double PERIHELION_JUPITER = A_JUPITER * (1 - EX_JUPITER);
const double PERIHELION_SATURN = A_SATURN * (1 - EX_SATURN);
const double PERIHELION_URAN = A_URAN * (1 - EX_URAN);
const double PERIHELION_NEPTUN = A_NEPTUN * (1 - EX_NEPTUN);

// Орбитальные радиусы спутников
const double LUNA_ORBIT_RADIUS = 3.844e8;
const double IO_ORBIT_RADIUS = 4.217e8;
const double EUROPA_ORBIT_RADIUS = 6.709e8;
const double GANYMEDE_ORBIT_RADIUS = 1.0704e9;
const double CALLISTO_ORBIT_RADIUS = 1.8827e9;
const double TITAN_ORBIT_RADIUS = 1.22187e9;
const double RHEA_ORBIT_RADIUS = 5.2704e8;
const double TETHYS_ORBIT_RADIUS = 2.9467e8;
const double DIONE_ORBIT_RADIUS = 3.7742e8;
const double TITANIA_ORBIT_RADIUS = 4.363e8;
const double OBERON_ORBIT_RADIUS = 5.835e8;
const double TRITON_ORBIT_RADIUS = 3.5476e8;
const double PHOBOS_ORBIT_RADIUS = 9.376e6;
const double DEIMOS_ORBIT_RADIUS = 2.345e7;

// Радиусы
const float RADIUS_SUN = 695780000.0f * SCALE * RADIUS_SCALE;
const float RADIUS_MERCURY = 2439700.0f * SCALE * RADIUS_SCALE;
const float RADIUS_VENUS = 6051800.0f * SCALE * RADIUS_SCALE;
const float RADIUS_EARTH = 6371000.0f * SCALE * RADIUS_SCALE;
const float RADIUS_MARS = 3389500.0f * SCALE * RADIUS_SCALE;
const float RADIUS_JUPITER = 69911000.0f * SCALE * RADIUS_SCALE;
const float RADIUS_SATURN = 58232000.0f * SCALE * RADIUS_SCALE;
const float RADIUS_URAN = 25362000.0f * SCALE * RADIUS_SCALE;
const float RADIUS_NEPTUN = 24622000.0f * SCALE * RADIUS_SCALE;
const float RADIUS_LUNA = 1737400.0f * SCALE * RADIUS_SCALE;
const float RADIUS_IO = 1821600.0f * SCALE * RADIUS_SCALE;
const float RADIUS_EUROPA = 1560800.0f * SCALE * RADIUS_SCALE;
const float RADIUS_GANYMEDE = 2634100.0f * SCALE * RADIUS_SCALE;
const float RADIUS_CALLISTO = 2410300.0f * SCALE * RADIUS_SCALE;
const float RADIUS_TITAN = 2574700.0f * SCALE * RADIUS_SCALE;
const float RADIUS_RHEA = 763800.0f * SCALE * RADIUS_SCALE;
const float RADIUS_TETHYS = 531100.0f * SCALE * RADIUS_SCALE;
const float RADIUS_DIONE = 561400.0f * SCALE * RADIUS_SCALE;
const float RADIUS_TITANIA = 788400.0f * SCALE * RADIUS_SCALE;
const float RADIUS_OBERON = 761400.0f * SCALE * RADIUS_SCALE;
const float RADIUS_TRITON = 1353400.0f * SCALE * RADIUS_SCALE;
const float RADIUS_PHOBOS = 11080.0f * SCALE * RADIUS_SCALE;
const float RADIUS_DEIMOS = 6200.0f * SCALE * RADIUS_SCALE;

// Цвета
const float COLOR_SUN[3] = {1.0f, 1.0f, 0.0f};
const float COLOR_MERCURY[3] = {0.7f, 0.7f, 0.7f};
const float COLOR_VENUS[3] = {1.0f, 0.8f, 0.4f};
const float COLOR_EARTH[3] = {0.0f, 0.5f, 1.0f};
const float COLOR_MARS[3] = {1.0f, 0.3f, 0.1f};
const float COLOR_JUPITER[3] = {0.9f, 0.7f, 0.5f};
const float COLOR_SATURN[3] = {0.95f, 0.85f, 0.6f};
const float COLOR_URAN[3] = {0.15f, 0.35f, 0.6f};
const float COLOR_NEPTUN[3] = {0.1f, 0.2f, 0.9f};
const float COLOR_LUNA[3] = {0.8f, 0.8f, 0.8f};

// ==================== СТРУКТУРЫ ====================

struct TrailPoint {
    double x, y, z;
    float r, g, b;
};

struct Trail {
    std::vector<TrailPoint> points;
    int maxPoints;
    float width;
    float startAlpha;
    float endAlpha;
    
    Trail(int max = 500, float w = 2.0f, float startA = 0.1f, float endA = 1.0f) 
        : maxPoints(max), width(w), startAlpha(startA), endAlpha(endA) {}
    
    void addPoint(double x, double y, double z, float r, float g, float b) {
        TrailPoint point;
        point.x = x;
        point.y = y;
        point.z = z;
        point.r = r;
        point.g = g;
        point.b = b;
        points.push_back(point);
        
        if(points.size() > maxPoints) {
            points.erase(points.begin());
        }
    }
    
    void draw() {
        if(points.size() < 2) return;
        
        glLineWidth(width);
        glBegin(GL_LINE_STRIP);
        
        for(size_t i = 0; i < points.size(); ++i) {
            float t = (float)i / (points.size() - 1);
            float alpha = startAlpha + (endAlpha - startAlpha) * t;
            
            glColor3f(points[i].r * alpha, 
                      points[i].g * alpha, 
                      points[i].b * alpha);
            glVertex3f(float(points[i].x * SCALE), 
                      float(points[i].y * SCALE), 
                      float(points[i].z * SCALE));
        }
        
        glEnd();
        glLineWidth(1.0f);
    }
};

struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct Planet {
    std::vector<Vertex> ver;
    std::vector<unsigned int> ind;
    float radius;
    double mass;
    double x, y, z;           // Текущие координаты (абсолютные или относительные)
    double vx, vy, vz;        // Текущие скорости
    double ax, ay, az;        // Текущие ускорения
    std::string name;
    float a;
    float orbit_eccentricity;
    float orbit_tilt;
    bool is_satellite = false;
    std::string parent;
    Trail trail;
    
    // Для восстановления при переключении систем
    double saved_x, saved_y, saved_z;
    double saved_vx, saved_vy, saved_vz;
};

struct NonInertialFrame {
    bool is_active = false;
    Planet* reference_planet = nullptr;
    double omega;  // Угловая скорость вращения НСО вокруг Солнца
    double omega_x, omega_y, omega_z;  // Направление угловой скорости
    
    void activate(Planet* planet, const std::vector<Planet*>& all_planets) {
        if (!planet) {
            deactivate(all_planets);
            return;
        }
        
        // Сохраняем текущие абсолютные координаты всех планет
        for (Planet* p : all_planets) {
            p->saved_x = p->x;
            p->saved_y = p->y;
            p->saved_z = p->z;
            p->saved_vx = p->vx;
            p->saved_vy = p->vy;
            p->saved_vz = p->vz;
        }
        
        reference_planet = planet;
        is_active = true;
        
        // Вычисляем угловую скорость вращения планеты вокруг Солнца
        compute_angular_velocity();
        
        // Преобразуем в относительные координаты
        convert_to_relative(all_planets);
        
        std::cout << "Активирована НСО: " << planet->name << std::endl;
        std::cout << "Угловая скорость: " << omega << " рад/с" << std::endl;
    }
    
    void deactivate(const std::vector<Planet*>& all_planets) {
        if (!is_active) return;
        
        // Восстанавливаем абсолютные координаты
        for (Planet* p : all_planets) {
            p->x = p->saved_x;
            p->y = p->saved_y;
            p->z = p->saved_z;
            p->vx = p->saved_vx;
            p->vy = p->saved_vy;
            p->vz = p->saved_vz;
        }
        
        is_active = false;
        reference_planet = nullptr;
        omega = 0;
        omega_x = omega_y = omega_z = 0;
        
        std::cout << "Деактивирована НСО (возврат в инерциальную систему)" << std::endl;
    }
    
    void compute_angular_velocity() {
        if (!reference_planet) return;
        
        // Угловая скорость для круговой орбиты: ω = sqrt(G*M_sun / r^3)
        double r = sqrt(reference_planet->saved_x * reference_planet->saved_x + 
                       reference_planet->saved_y * reference_planet->saved_y + 
                       reference_planet->saved_z * reference_planet->saved_z);
        
        if (r > 1e-12) {
            omega = sqrt(G * MASS_SUN / (r * r * r));
            
            // Направление угловой скорости (перпендикулярно плоскости орбиты)
            // Используем векторное произведение r × v для определения направления
            double rx = reference_planet->saved_x;
            double ry = reference_planet->saved_y;
            double rz = reference_planet->saved_z;
            double vx = reference_planet->saved_vx;
            double vy = reference_planet->saved_vy;
            double vz = reference_planet->saved_vz;
            
            // r × v
            double cross_x = ry * vz - rz * vy;
            double cross_y = rz * vx - rx * vz;
            double cross_z = rx * vy - ry * vx;
            
            double cross_mag = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
            
            if (cross_mag > 1e-12) {
                omega_x = cross_x / cross_mag;
                omega_y = cross_y / cross_mag;
                omega_z = cross_z / cross_mag;
            } else {
                // По умолчанию вращение вокруг оси Y
                omega_x = 0;
                omega_y = 1;
                omega_z = 0;
            }
        }
    }
    
    void convert_to_relative(const std::vector<Planet*>& all_planets) {
        if (!reference_planet || !is_active) return;
        
        double ref_x = reference_planet->saved_x;
        double ref_y = reference_planet->saved_y;
        double ref_z = reference_planet->saved_z;
        double ref_vx = reference_planet->saved_vx;
        double ref_vy = reference_planet->saved_vy;
        double ref_vz = reference_planet->saved_vz;
        
        for (Planet* planet : all_planets) {
            if (planet == reference_planet) {
                // Референсная планета в центре
                planet->x = 0;
                planet->y = 0;
                planet->z = 0;
                planet->vx = 0;
                planet->vy = 0;
                planet->vz = 0;
            } else {
                // Относительные координаты
                planet->x = planet->saved_x - ref_x;
                planet->y = planet->saved_y - ref_y;
                planet->z = planet->saved_z - ref_z;
                
                // Относительные скорости (с учетом вращения НСО)
                planet->vx = planet->saved_vx - ref_vx;
                planet->vy = planet->saved_vy - ref_vy;
                planet->vz = planet->saved_vz - ref_vz;
            }
        }
    }
    
    // Вычисление сил инерции
    void compute_inertial_forces(double x, double y, double z,
                                 double vx, double vy, double vz,
                                 double& fx, double& fy, double& fz) {
        // Центробежная сила: F_cent = m * ω × (ω × r)
        // Где ω - угловая скорость НСО, r - радиус-вектор в НСО
        
        // ω × r
        double omega_cross_r_x = omega_y * z - omega_z * y;
        double omega_cross_r_y = omega_z * x - omega_x * z;
        double omega_cross_r_z = omega_x * y - omega_y * x;
        
        // ω × (ω × r)
        fx = omega_y * omega_cross_r_z - omega_z * omega_cross_r_y;
        fy = omega_z * omega_cross_r_x - omega_x * omega_cross_r_z;
        fz = omega_x * omega_cross_r_y - omega_y * omega_cross_r_x;
        
        // Умножаем на ω²
        double omega_sq = omega * omega;
        fx *= omega_sq;
        fy *= omega_sq;
        fz *= omega_sq;
        
        // Сила Кориолиса: F_cor = -2m * (ω × v)
        double cor_x = -2.0 * (omega_y * vz - omega_z * vy);
        double cor_y = -2.0 * (omega_z * vx - omega_x * vz);
        double cor_z = -2.0 * (omega_x * vy - omega_y * vx);
        
        // Умножаем на ω
        cor_x *= omega;
        cor_y *= omega;
        cor_z *= omega;
        
        // Суммарная сила инерции
        fx += cor_x;
        fy += cor_y;
        fz += cor_z;
    }
};

struct Camera {
    float posX, posY, posZ;
    float targetX, targetY, targetZ;
    float upX, upY, upZ;
    float distance;
    float angleX, angleY;
    
    Camera() {
        posX = 100.0f; posY = 150.0f; posZ = 100.0f;
        targetX = 0.0f; targetY = 0.0f; targetZ = 0.0f;
        upX = 0.0f; upY = 1.0f; upZ = 0.0f;
        distance = 200.0f;
        angleX = 0.0f;
        angleY = 0.3f;
    }
    
    void update() {
        posX = targetX + distance * cos(angleY) * sin(angleX);
        posY = targetY + distance * sin(angleY);
        posZ = targetZ + distance * cos(angleY) * cos(angleX);
    }
    
    void focus_on_planet(Planet* planet, const NonInertialFrame& frame) {
        if (!planet) return;
        
        if (frame.is_active && frame.reference_planet) {
            // В НСО используем относительные координаты
            targetX = float(planet->x * SCALE);
            targetY = float(planet->y * SCALE);
            targetZ = float(planet->z * SCALE);
        } else {
            // В инерциальной системе используем абсолютные координаты
            targetX = float(planet->x * SCALE);
            targetY = float(planet->y * SCALE);
            targetZ = float(planet->z * SCALE);
        }
        
        // Настраиваем расстояние
        if (planet->name == "Sun") {
            distance = 10.0f;
        } else if (planet->is_satellite) {
            distance = 0.5f;
        } else {
            distance = 3.0f;
        }
        
        angleX = 0.0f;
        angleY = 0.0f;
        update();
    }
};

struct Tracker {
    double prev_angle = 0;
    double prev_time = 0;
    double period = 0;
    int revolutions = 0;
    bool has_period = false;
    
    void update(double x, double z, double time) {
        time *= TIME_SCALE;
        double current_angle = atan2(z, x);
        
        if (prev_time == 0) {
            prev_angle = current_angle;
            prev_time = time;
            return;
        }
        
        double angle_diff = current_angle - prev_angle;
        if (angle_diff > M_PI) angle_diff -= 2 * M_PI;
        if (angle_diff < -M_PI) angle_diff += 2 * M_PI;
        
        if (fabs(angle_diff) > 0.1) {  // Фильтр шума
            revolutions++;
            
            if (revolutions >= 1) {
                period = time / revolutions;
                has_period = true;
            }
        }
        
        prev_angle = current_angle;
        prev_time = time;
    }
    
    double get_period() { return period; }
    bool ready() { return has_period; }
};

// ==================== ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ ====================
std::vector<Planet*> planets;
NonInertialFrame non_inertial_frame;
Camera camera;

// ==================== ФУНКЦИИ ====================

void create_sphere(Planet* planet, float radius, float r, float g, float b) {
    planet->radius = radius;
    int segments = 30;
    
    for (int i = 0; i <= segments; ++i) {
        float phi = M_PI * i / segments;
        for (int j = 0; j <= segments; ++j) {
            float theta = 2.0f * M_PI * j / segments;
            Vertex v;
            v.x = radius * sin(phi) * cos(theta);
            v.y = radius * cos(phi);
            v.z = radius * sin(phi) * sin(theta);
            v.r = r;
            v.g = g;
            v.b = b;
            planet->ver.push_back(v);
        }
    }
    
    for (int i = 0; i < segments; ++i) {
        for (int j = 0; j < segments; ++j) {
            int first = i * (segments + 1) + j;
            int second = first + 1;
            int third = first + (segments + 1);
            int fourth = third + 1;
            
            planet->ind.push_back(first);
            planet->ind.push_back(second);
            planet->ind.push_back(third);
            
            planet->ind.push_back(second);
            planet->ind.push_back(fourth);
            planet->ind.push_back(third);
        }
    }
}

void draw_planet(Planet* planet) {
    glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < planet->ind.size(); ++i) {
        Vertex* v = &planet->ver[planet->ind[i]];
        glColor3f(v->r, v->g, v->b);
        glVertex3f(v->x, v->y, v->z);
    }
    glEnd();
}

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

void render_planet(Planet* planet) {
    glPushMatrix();
    glTranslatef(float(planet->x * SCALE), 
                 float(planet->y * SCALE), 
                 float(planet->z * SCALE));
    draw_planet(planet);
    glPopMatrix();
}

void draw_coordinate_system() {
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(5.0f, 0.0f, 0.0f);
    
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 5.0f, 0.0f);
    
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 5.0f);
    
    glEnd();
    glLineWidth(1.0f);
}

// Обновление физики в инерциальной системе
void update_physics_inertial(Planet* planet, const std::vector<Planet*>& all_planets, double dt) {
    // Полушаговый метод (leapfrog)
    planet->vx += 0.5 * planet->ax * dt;
    planet->vy += 0.5 * planet->ay * dt;
    planet->vz += 0.5 * planet->az * dt;
    
    planet->x += planet->vx * dt;
    planet->y += planet->vy * dt;
    planet->z += planet->vz * dt;
    
    double total_ax = 0, total_ay = 0, total_az = 0;
    for (Planet* other : all_planets) {
        if (planet == other) continue;
        
        double dx = other->x - planet->x;
        double dy = other->y - planet->y;
        double dz = other->z - planet->z;
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-12) continue;
        
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

// Обновление физики в неинерциальной системе
void update_physics_non_inertial(Planet* planet, const std::vector<Planet*>& all_planets, double dt) {
    // Если это планета-референс, она неподвижна
    if (planet == non_inertial_frame.reference_planet) {
        planet->vx = planet->vy = planet->vz = 0;
        planet->ax = planet->ay = planet->az = 0;
        return;
    }
    
    // Полушаговый метод с учетом сил инерции
    planet->vx += 0.5 * planet->ax * dt;
    planet->vy += 0.5 * planet->ay * dt;
    planet->vz += 0.5 * planet->az * dt;
    
    planet->x += planet->vx * dt;
    planet->y += planet->vy * dt;
    planet->z += planet->vz * dt;
    
    // Вычисляем гравитационные ускорения
    double total_ax = 0, total_ay = 0, total_az = 0;
    for (Planet* other : all_planets) {
        if (planet == other) continue;
        
        double dx, dy, dz;
        if (other == non_inertial_frame.reference_planet) {
            // Референсная планета в центре
            dx = -planet->x;
            dy = -planet->y;
            dz = -planet->z;
        } else {
            dx = other->x - planet->x;
            dy = other->y - planet->y;
            dz = other->z - planet->z;
        }
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-12) continue;
        
        double r = sqrt(r2);
        double a = G * other->mass / r2;
        
        total_ax += a * dx / r;
        total_ay += a * dy / r;
        total_az += a * dz / r;
    }
    
    // Добавляем силы инерции
    double inertial_x, inertial_y, inertial_z;
    non_inertial_frame.compute_inertial_forces(planet->x, planet->y, planet->z,
                                              planet->vx, planet->vy, planet->vz,
                                              inertial_x, inertial_y, inertial_z);
    
    total_ax += inertial_x;
    total_ay += inertial_y;
    total_az += inertial_z;
    
    planet->vx += 0.5 * total_ax * dt;
    planet->vy += 0.5 * total_ay * dt;
    planet->vz += 0.5 * total_az * dt;
    
    planet->ax = total_ax;
    planet->ay = total_ay;
    planet->az = total_az;
}

void update_physics(Planet* planet, const std::vector<Planet*>& all_planets, double dt) {
    if (non_inertial_frame.is_active) {
        update_physics_non_inertial(planet, all_planets, dt);
    } else {
        update_physics_inertial(planet, all_planets, dt);
    }
}

// Инициализация планет
void initialize_planets() {
    // Вычисление начальных скоростей
    double VEL_MERCURY = sqrt(2 * ((-G * MASS_SUN / (2 * A_MERCURY)) + (G * MASS_SUN / PERIHELION_MERCURY)));
    double VEL_VENUS = sqrt(2 * ((-G * MASS_SUN / (2 * A_VENUS)) + (G * MASS_SUN / PERIHELION_VENUS)));
    double VEL_EARTH = sqrt(2 * ((-G * MASS_SUN / (2 * A_EARTH)) + (G * MASS_SUN / PERIHELION_EARTH)));
    double VEL_MARS = sqrt(2 * ((-G * MASS_SUN / (2 * A_MARS)) + (G * MASS_SUN / PERIHELION_MARS)));
    double VEL_JUPITER = sqrt(2 * ((-G * MASS_SUN / (2 * A_JUPITER)) + (G * MASS_SUN / PERIHELION_JUPITER)));
    double VEL_SATURN = sqrt(2 * ((-G * MASS_SUN / (2 * A_SATURN)) + (G * MASS_SUN / PERIHELION_SATURN)));
    double VEL_URAN = sqrt(2 * ((-G * MASS_SUN / (2 * A_URAN)) + (G * MASS_SUN / PERIHELION_URAN)));
    double VEL_NEPTUN = sqrt(2 * ((-G * MASS_SUN / (2 * A_NEPTUN)) + (G * MASS_SUN / PERIHELION_NEPTUN)));
    double VEL_LUNA = sqrt(G * MASS_EARTH / LUNA_ORBIT_RADIUS);
    
    // Солнце
    Planet* sun = new Planet();
    sun->name = "Sun";
    sun->mass = MASS_SUN;
    sun->x = sun->y = sun->z = 0;
    sun->vx = sun->vy = sun->vz = 0;
    sun->ax = sun->ay = sun->az = 0;
    create_sphere(sun, RADIUS_SUN, COLOR_SUN[0], COLOR_SUN[1], COLOR_SUN[2]);
    sun->trail = Trail(TRIAL_STEPS, 2.0f);
    planets.push_back(sun);
    
    // Меркурий
    Planet* mercury = new Planet();
    mercury->name = "Mercury";
    mercury->mass = MASS_MERCURY;
    mercury->x = PERIHELION_MERCURY;
    mercury->y = mercury->z = 0;
    mercury->vx = 0;
    mercury->vy = VEL_MERCURY * sin(TILT_MERCURY);
    mercury->vz = VEL_MERCURY * cos(TILT_MERCURY);
    mercury->ax = mercury->ay = mercury->az = 0;
    mercury->orbit_eccentricity = EX_MERCURY;
    mercury->a = A_MERCURY;
    mercury->orbit_tilt = TILT_MERCURY;
    mercury->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(mercury, RADIUS_MERCURY, COLOR_MERCURY[0], COLOR_MERCURY[1], COLOR_MERCURY[2]);
    planets.push_back(mercury);
    
    // Венера
    Planet* venus = new Planet();
    venus->name = "Venus";
    venus->mass = MASS_VENUS;
    venus->x = PERIHELION_VENUS;
    venus->y = venus->z = 0;
    venus->vx = 0;
    venus->vy = VEL_VENUS * sin(TILT_VENUS);
    venus->vz = VEL_VENUS * cos(TILT_VENUS);
    venus->ax = venus->ay = venus->az = 0;
    venus->orbit_eccentricity = EX_VENUS;
    venus->a = A_VENUS;
    venus->orbit_tilt = TILT_VENUS;
    venus->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(venus, RADIUS_VENUS, COLOR_VENUS[0], COLOR_VENUS[1], COLOR_VENUS[2]);
    planets.push_back(venus);
    
    // Земля
    Planet* earth = new Planet();
    earth->name = "Earth";
    earth->mass = MASS_EARTH;
    earth->x = PERIHELION_EARTH;
    earth->y = earth->z = 0;
    earth->vx = 0;
    earth->vy = VEL_EARTH * sin(TILT_EARTH);
    earth->vz = VEL_EARTH * cos(TILT_EARTH);
    earth->ax = earth->ay = earth->az = 0;
    earth->orbit_eccentricity = EX_EARTH;
    earth->a = A_EARTH;
    earth->orbit_tilt = TILT_EARTH;
    earth->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(earth, RADIUS_EARTH, COLOR_EARTH[0], COLOR_EARTH[1], COLOR_EARTH[2]);
    planets.push_back(earth);
    
    // Марс
    Planet* mars = new Planet();
    mars->name = "Mars";
    mars->mass = MASS_MARS;
    mars->x = PERIHELION_MARS;
    mars->y = mars->z = 0;
    mars->vx = 0;
    mars->vy = VEL_MARS * sin(TILT_MARS);
    mars->vz = VEL_MARS * cos(TILT_MARS);
    mars->ax = mars->ay = mars->az = 0;
    mars->orbit_eccentricity = EX_MARS;
    mars->a = A_MARS;
    mars->orbit_tilt = TILT_MARS;
    mars->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(mars, RADIUS_MARS, COLOR_MARS[0], COLOR_MARS[1], COLOR_MARS[2]);
    planets.push_back(mars);
    
    // Юпитер
    Planet* jupiter = new Planet();
    jupiter->name = "Jupiter";
    jupiter->mass = MASS_JUPITER;
    jupiter->x = PERIHELION_JUPITER;
    jupiter->y = jupiter->z = 0;
    jupiter->vx = 0;
    jupiter->vy = VEL_JUPITER * sin(TILT_JUPITER);
    jupiter->vz = VEL_JUPITER * cos(TILT_JUPITER);
    jupiter->ax = jupiter->ay = jupiter->az = 0;
    jupiter->orbit_eccentricity = EX_JUPITER;
    jupiter->a = A_JUPITER;
    jupiter->orbit_tilt = TILT_JUPITER;
    jupiter->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(jupiter, RADIUS_JUPITER, COLOR_JUPITER[0], COLOR_JUPITER[1], COLOR_JUPITER[2]);
    planets.push_back(jupiter);
    
    // Сатурн
    Planet* saturn = new Planet();
    saturn->name = "Saturn";
    saturn->mass = MASS_SATURN;
    saturn->x = PERIHELION_SATURN;
    saturn->y = saturn->z = 0;
    saturn->vx = 0;
    saturn->vy = VEL_SATURN * sin(TILT_SATURN);
    saturn->vz = VEL_SATURN * cos(TILT_SATURN);
    saturn->ax = saturn->ay = saturn->az = 0;
    saturn->orbit_eccentricity = EX_SATURN;
    saturn->a = A_SATURN;
    saturn->orbit_tilt = TILT_SATURN;
    saturn->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(saturn, RADIUS_SATURN, COLOR_SATURN[0], COLOR_SATURN[1], COLOR_SATURN[2]);
    planets.push_back(saturn);
    
    // Уран
    Planet* uran = new Planet();
    uran->name = "Uran";
    uran->mass = MASS_URAN;
    uran->x = PERIHELION_URAN;
    uran->y = uran->z = 0;
    uran->vx = 0;
    uran->vy = VEL_URAN * sin(TILT_URAN);
    uran->vz = VEL_URAN * cos(TILT_URAN);
    uran->ax = uran->ay = uran->az = 0;
    uran->orbit_eccentricity = EX_URAN;
    uran->a = A_URAN;
    uran->orbit_tilt = TILT_URAN;
    uran->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(uran, RADIUS_URAN, COLOR_URAN[0], COLOR_URAN[1], COLOR_URAN[2]);
    planets.push_back(uran);
    
    // Нептун
    Planet* neptun = new Planet();
    neptun->name = "Neptun";
    neptun->mass = MASS_NEPTUN;
    neptun->x = PERIHELION_NEPTUN;
    neptun->y = neptun->z = 0;
    neptun->vx = 0;
    neptun->vy = VEL_NEPTUN * sin(TILT_NEPTUN);
    neptun->vz = VEL_NEPTUN * cos(TILT_NEPTUN);
    neptun->ax = neptun->ay = neptun->az = 0;
    neptun->orbit_eccentricity = EX_NEPTUN;
    neptun->a = A_NEPTUN;
    neptun->orbit_tilt = TILT_NEPTUN;
    neptun->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(neptun, RADIUS_NEPTUN, COLOR_NEPTUN[0], COLOR_NEPTUN[1], COLOR_NEPTUN[2]);
    planets.push_back(neptun);
    
    // Луна
    Planet* luna = new Planet();
    luna->name = "Luna";
    luna->mass = MASS_LUNA;
    luna->x = PERIHELION_EARTH + LUNA_ORBIT_RADIUS;
    luna->y = luna->z = 0;
    luna->vx = 0;
    luna->vy = VEL_EARTH * sin(TILT_EARTH);
    luna->vz = VEL_EARTH * cos(TILT_EARTH) + VEL_LUNA;
    luna->ax = luna->ay = luna->az = 0;
    luna->is_satellite = true;
    luna->parent = "Earth";
    luna->trail = Trail(TRIAL_STEPS, 2.0f);
    create_sphere(luna, RADIUS_LUNA, COLOR_LUNA[0], COLOR_LUNA[1], COLOR_LUNA[2]);
    planets.push_back(luna);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        float cameraSpeed = 1.0f;
        
        // Управление камерой
        if (key == GLFW_KEY_UP) camera.targetY += cameraSpeed;
        if (key == GLFW_KEY_DOWN) camera.targetY -= cameraSpeed;
        if (key == GLFW_KEY_LEFT) camera.targetX -= cameraSpeed;
        if (key == GLFW_KEY_RIGHT) camera.targetX += cameraSpeed;
        if (key == GLFW_KEY_W) camera.targetZ -= cameraSpeed;
        if (key == GLFW_KEY_S) camera.targetZ += cameraSpeed;
        
        // Приближение/отдаление
        if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_KP_ADD) camera.distance *= 0.9f;
        if (key == GLFW_KEY_MINUS || key == GLFW_KEY_KP_SUBTRACT) camera.distance *= 1.1f;
        
        // Фокусировка на планетах (всегда работает, даже в НСО)
        if (key == GLFW_KEY_1) {
            for (Planet* p : planets) if (p->name == "Sun") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_2) {
            for (Planet* p : planets) if (p->name == "Mercury") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_3) {
            for (Planet* p : planets) if (p->name == "Venus") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_4) {
            for (Planet* p : planets) if (p->name == "Earth") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_5) {
            for (Planet* p : planets) if (p->name == "Mars") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_6) {
            for (Planet* p : planets) if (p->name == "Jupiter") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_7) {
            for (Planet* p : planets) if (p->name == "Saturn") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_8) {
            for (Planet* p : planets) if (p->name == "Uran") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_9) {
            for (Planet* p : planets) if (p->name == "Neptun") camera.focus_on_planet(p, non_inertial_frame);
        }
        if (key == GLFW_KEY_0) {
            for (Planet* p : planets) if (p->name == "Luna") camera.focus_on_planet(p, non_inertial_frame);
        }
        
        // Сброс камеры
        if (key == GLFW_KEY_R) {
            camera = Camera();
            camera.update();
        }
        
        // Фокус на центре
        if (key == GLFW_KEY_F) {
            camera.targetX = camera.targetY = camera.targetZ = 0;
            camera.distance = 200.0f;
            camera.angleX = camera.angleY = 0;
            camera.update();
        }
        
        // Активация НСО
        if (key == GLFW_KEY_F1) {
            for (Planet* p : planets) if (p->name == "Mercury") non_inertial_frame.activate(p, planets);
        }
        if (key == GLFW_KEY_F2) {
            for (Planet* p : planets) if (p->name == "Venus") non_inertial_frame.activate(p, planets);
        }
        if (key == GLFW_KEY_F3) {
            for (Planet* p : planets) if (p->name == "Earth") non_inertial_frame.activate(p, planets);
        }
        if (key == GLFW_KEY_F4) {
            for (Planet* p : planets) if (p->name == "Mars") non_inertial_frame.activate(p, planets);
        }
        if (key == GLFW_KEY_F5) {
            for (Planet* p : planets) if (p->name == "Jupiter") non_inertial_frame.activate(p, planets);
        }
        if (key == GLFW_KEY_F6) {
            for (Planet* p : planets) if (p->name == "Saturn") non_inertial_frame.activate(p, planets);
        }
        if (key == GLFW_KEY_F7) {
            for (Planet* p : planets) if (p->name == "Uran") non_inertial_frame.activate(p, planets);
        }
        if (key == GLFW_KEY_F8) {
            for (Planet* p : planets) if (p->name == "Neptun") non_inertial_frame.activate(p, planets);
        }
        
        // Возврат в инерциальную систему
        if (key == GLFW_KEY_F12) {
            non_inertial_frame.deactivate(planets);
        }
        
        camera.update();
    }
}

void display_info() {
    static int counter = 0;
    counter++;
    
    if (counter % 60 == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Система отсчета: " << (non_inertial_frame.is_active ? "НЕИНЕРЦИАЛЬНАЯ" : "ИНЕРЦИАЛЬНАЯ") << std::endl;
        
        if (non_inertial_frame.is_active && non_inertial_frame.reference_planet) {
            std::cout << "Центр НСО: " << non_inertial_frame.reference_planet->name << std::endl;
            std::cout << "Угловая скорость: " << non_inertial_frame.omega << " рад/с" << std::endl;
            
            // Для примера покажем относительные координаты Солнца
            for (Planet* p : planets) {
                if (p->name == "Sun") {
                    std::cout << "Солнце в НСО: (" << p->x << ", " << p->y << ", " << p->z << ")" << std::endl;
                    break;
                }
            }
        }
        std::cout << "========================================\n" << std::endl;
    }
}

int main() {
    // Инициализация GLFW
    if (!glfwInit()) {
        std::cout << "Ошибка инициализации GLFW" << std::endl;
        return -1;
    }
    
    // Создание окна
    GLFWwindow* window = glfwCreateWindow(1200, 800, "Солнечная система", NULL, NULL);
    if (!window) {
        std::cout << "Ошибка создания окна" << std::endl;
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent(window);
    
    // Инициализация GLEW
    if (glewInit() != GLEW_OK) {
        std::cout << "Ошибка инициализации GLEW" << std::endl;
        return -1;
    }
    
    // Установка callback
    glfwSetKeyCallback(window, key_callback);
    
    // Инициализация камеры
    camera.update();
    
    // Инициализация планет
    initialize_planets();
    
    // Настройка OpenGL
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.05f, 1.0f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Инициализация трекеров
    std::map<std::string, Tracker> trackers;
    for (Planet* p : planets) {
        trackers[p->name] = Tracker();
    }
    
    double last_time = glfwGetTime();
    
    // Главный цикл
    while (!glfwWindowShouldClose(window)) {
        double current_time = glfwGetTime();
        double delta_time = current_time - last_time;
        last_time = current_time;
        
        double dt = delta_time * TIME_SCALE;
        
        // Обновление физики
        for (Planet* planet : planets) {
            update_physics(planet, planets, dt);
        }
        
        // Добавление точек к трейлам
        static int frame_counter = 0;
        frame_counter++;
        if (frame_counter % 3 == 0) {
            for (Planet* planet : planets) {
                planet->trail.addPoint(planet->x, planet->y, planet->z,
                                      planet->ver[0].r, planet->ver[0].g, planet->ver[0].b);
            }
        }
        
        // Обновление трекеров
        for (Planet* planet : planets) {
            trackers[planet->name].update(planet->x, planet->z, current_time);
        }
        
        // Отрисовка
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Настройка viewport
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        
        // Настройка проекции
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)width / (float)height;
        gluPerspective(60.0f, aspect, 0.001f, 10000.0f);
        
        // Настройка вида
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(camera.posX, camera.posY, camera.posZ,
                  camera.targetX, camera.targetY, camera.targetZ,
                  camera.upX, camera.upY, camera.upZ);
        
        // Отрисовка системы координат
        draw_coordinate_system();
        
        // Отрисовка орбит (только в инерциальной системе)
        if (!non_inertial_frame.is_active) {
            glDisable(GL_DEPTH_TEST);
            glLineWidth(0.5f);
            
            for (Planet* planet : planets) {
                if (!planet->is_satellite && planet->name != "Sun") {
                    float color[3];
                    if (planet->name == "Mercury") { color[0]=0.5f; color[1]=0.5f; color[2]=0.5f; }
                    else if (planet->name == "Venus") { color[0]=1.0f; color[1]=0.8f; color[2]=0.4f; }
                    else if (planet->name == "Earth") { color[0]=0.3f; color[1]=0.3f; color[2]=0.8f; }
                    else if (planet->name == "Mars") { color[0]=0.8f; color[1]=0.3f; color[2]=0.2f; }
                    else if (planet->name == "Jupiter") { color[0]=1.0f; color[1]=0.4f; color[2]=0.5f; }
                    else if (planet->name == "Saturn") { color[0]=0.2f; color[1]=0.7f; color[2]=1.0f; }
                    else if (planet->name == "Uran") { color[0]=0.45f; color[1]=0.85f; color[2]=0.15f; }
                    else if (planet->name == "Neptun") { color[0]=0.12f; color[1]=0.95f; color[2]=1.0f; }
                    
                    draw_orbit(planet->a, planet->orbit_eccentricity, planet->orbit_tilt,
                              color[0], color[1], color[2]);
                }
            }
            glEnable(GL_DEPTH_TEST);
            glLineWidth(1.0f);
        }
        
        // Отрисовка трейлов
        glDisable(GL_DEPTH_TEST);
        for (Planet* planet : planets) {
            planet->trail.draw();
        }
        glEnable(GL_DEPTH_TEST);
        
        // Отрисовка планет
        for (Planet* planet : planets) {
            render_planet(planet);
        }
        
        // Вывод информации
        display_info();
        
        // Обмен буферов
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    // Очистка
    for (Planet* p : planets) delete p;
    glfwDestroyWindow(window);
    glfwTerminate();
    
    return 0;
}