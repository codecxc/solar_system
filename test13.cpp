#include<iostream>
#include<GL/glew.h>
#include<GLFW/glfw3.h>
#include<vector>
#include<cmath>
#include<unistd.h>
#include<algorithm>
#include<string>
#include<sstream>
#include<map>

double SIMULATION_SPEED_MULTIPLIER = 1.0;
const double AU = 1.496e11;
const double G = 6.67430e-11;
const double SCALE = 250.0 / AU;  // Уменьшил масштаб для лучшей видимости
const double TIME_STEP = 3600.0;  // Шаг времени 1 час (в секундах)

const float MAX_DISTANCE = 100.0f;
const float CLIP_NEAR = 0.1f;
const float CLIP_FAR = 10000.0f;

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

const float EX_SUN = 0.0f;
const float EX_MERCURY = 0.206f;
const float EX_VENUS = 0.007f;
const float EX_EARTH = 0.017f;
const float EX_MARS = 0.093f;
const float EX_JUPITER = 0.049f;
const float EX_SATURN = 0.057f;
const float EX_URAN = 0.046f;
const float EX_NEPTUN = 0.011f;

const float A_SUN = 0.0f;
const float A_MERCURY = 0.387 * AU;
const float A_VENUS = 0.723 * AU;
const float A_EARTH = 1.000 * AU;
const float A_MARS = 1.524 * AU;
const float A_JUPITER = 5.2044 * AU;
const float A_SATURN = 9.5826 * AU;
const float A_URAN = 19.21840 * AU;
const float A_NEPTUN = 30.11000 * AU;

const float TILT_MERCURY = 7.01f * M_PI / 180.0f;
const float TILT_VENUS = 3.39f * M_PI / 180.0f;
const float TILT_EARTH = 0.0f;
const float TILT_MARS = 1.85f * M_PI / 180.0f;
const float TILT_JUPITER = 1.31f * M_PI / 180.0f;
const float TILT_SATURN = 2.49f * M_PI / 180.0f;
const float TILT_URAN = 0.77f * M_PI / 180.0f;
const float TILT_NEPTUN = 1.77f * M_PI / 180.0f;

// Вычисляем перигелии правильно
const double PERIHELION_MERCURY = A_MERCURY * (1 - EX_MERCURY);
const double PERIHELION_VENUS = A_VENUS * (1 - EX_VENUS);
const double PERIHELION_EARTH = A_EARTH * (1 - EX_EARTH);
const double PERIHELION_MARS = A_MARS * (1 - EX_MARS);
const double PERIHELION_JUPITER = A_JUPITER * (1 - EX_JUPITER);
const double PERIHELION_SATURN = A_SATURN * (1 - EX_SATURN);
const double PERIHELION_URAN = A_URAN * (1 - EX_URAN);
const double PERIHELION_NEPTUN = A_NEPTUN * (1 - EX_NEPTUN);

// Радиусы для отображения (условные)
const float RADIUS_SUN = 10.0f;
const float RADIUS_MERCURY = 1.0f;
const float RADIUS_VENUS = 1.5f;
const float RADIUS_EARTH = 1.6f;
const float RADIUS_MARS = 1.4f;
const float RADIUS_JUPITER = 5.0f;
const float RADIUS_SATURN = 4.5f;
const float RADIUS_URAN = 3.5f;
const float RADIUS_NEPTUN = 3.5f;
const float RADIUS_LUNA = 0.5f;
const float RADIUS_SATELLITE = 1.0f;

const float COLOR_SUN[3] = {1.0f, 1.0f, 0.0f};
const float COLOR_MERCURY[3] = {0.7f, 0.7f, 0.7f};
const float COLOR_VENUS[3] = {1.0f, 0.8f, 0.4f};
const float COLOR_EARTH[3] = {0.0f, 0.5f, 1.0f};
const float COLOR_MARS[3] = {1.0f, 0.3f, 0.1f};
const float COLOR_JUPITER[3] = {0.9f, 0.7f, 0.5f};
const float COLOR_SATURN[3] = {0.95f, 0.85f, 0.6f};
const float COLOR_NEPTUN[3] = {0.1f, 0.2f, 0.0f};
const float COLOR_URAN[3] = {0.15f, 0.35f, 0.6f};
const float COLOR_LUNA[3] = {1.0f, 1.0f, 1.0f};
const float COLOR_SATELLITE[3] = {1.0f, 0.0f, 1.0f};

const double LUNA_ORBIT_RADIUS = 3.844e8;

bool showLagrangePoints = true;
bool lagrangeSatelliteActive = false;
int selectedLagrangeSystem = 0;
int selectedLagrangePoint = 1;

struct Vertex {
    float x, y, z;
    float r, g, b;
};

struct LagrangePoint {
    double x, y, z;
    int type;
    float color[3];
    std::string system;
    double distanceFromPrimary;
    double R;
    double mu;
};

struct Planet {
    std::vector<Vertex> ver;
    std::vector<unsigned int> ind;
    float radius;
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    std::string name;
    float a;
    float orbit_eccentricity;
    float orbit_tilt;
    bool is_satellite = false;
    bool is_lagrange_satellite = false;
    std::string parent;
    int lagrange_point_type = 0;
    double true_anomaly = 0;  // Для инициализации орбиты
};

// Вычисление скорости в заданной точке орбиты
double calculateOrbitalVelocity(double a, double e, double trueAnomaly, double centralMass) {
    double r = a * (1 - e * e) / (1 + e * cos(trueAnomaly));
    double energy = -G * centralMass / (2 * a);
    double v = sqrt(2 * (G * centralMass / r + energy));
    return v;
}

void calculateCenterOfMass(Planet* primary, Planet* secondary, double& com_x, double& com_y, double& com_z) {
    if (!primary || !secondary) {
        com_x = com_y = com_z = 0.0;
        return;
    }
    double totalMass = primary->mass + secondary->mass;
    com_x = (primary->mass * primary->x + secondary->mass * secondary->x) / totalMass;
    com_y = (primary->mass * primary->y + secondary->mass * secondary->y) / totalMass;
    com_z = (primary->mass * primary->z + secondary->mass * secondary->z) / totalMass;
}

std::vector<LagrangePoint> calculateLagrangePoints(Planet* primary, Planet* secondary) {
    std::vector<LagrangePoint> points;
    if (!primary || !secondary) return points;
    
    double alpha = secondary->mass / (primary->mass + secondary->mass);
    double delta_x = secondary->x - primary->x;
    double delta_y = secondary->y - primary->y;
    double delta_z = secondary->z - primary->z;
    double distance_R = std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
    
    double center_mass_x, center_mass_y, center_mass_z;
    calculateCenterOfMass(primary, secondary, center_mass_x, center_mass_y, center_mass_z);
    
    double vector_x = secondary->x - center_mass_x;
    double vector_y = secondary->y - center_mass_y;
    double vector_z = secondary->z - center_mass_z;
    double vector_length = std::sqrt(vector_x * vector_x + vector_y * vector_y + vector_z * vector_z);
    
    if (vector_length > 0) {
        vector_x /= vector_length;
        vector_y /= vector_length;
        vector_z /= vector_length;
    }
    
    double perpendicular_x = -vector_z;
    double perpendicular_y = 0;
    double perpendicular_z = vector_x;
    double perpendicular_length = std::sqrt(perpendicular_x * perpendicular_x + perpendicular_y * perpendicular_y + perpendicular_z * perpendicular_z);
    
    if (perpendicular_length > 0) {
        perpendicular_x /= perpendicular_length;
        perpendicular_y /= perpendicular_length;
        perpendicular_z /= perpendicular_length;
    }
    
    double alpha_root = std::cbrt(alpha / 3.0);
    double sin_60 = std::sqrt(3.0) / 2.0;
    double cos_60 = 0.5;
    
    double L1_offset = distance_R * (1.0 - alpha_root);
    double L2_offset = distance_R * (1.0 + alpha_root);
    double L3_offset = -distance_R * (1.0 + (5.0 / 12.0) * alpha);
    
    for (int point_type = 1; point_type <= 5; point_type++) {
        LagrangePoint point;
        point.type = point_type;
        point.system = primary->name + "-" + secondary->name;
        point.mu = alpha;
        point.R = distance_R;
        
        if (point_type == 1) {
            point.x = center_mass_x + L1_offset * vector_x;
            point.y = center_mass_y + L1_offset * vector_y;
            point.z = center_mass_z + L1_offset * vector_z;
            point.color[0] = 1.0f; point.color[1] = 0.0f; point.color[2] = 0.0f;
        } else if (point_type == 2) {
            point.x = center_mass_x + L2_offset * vector_x;
            point.y = center_mass_y + L2_offset * vector_y;
            point.z = center_mass_z + L2_offset * vector_z;
            point.color[0] = 0.0f; point.color[1] = 1.0f; point.color[2] = 0.0f;
        } else if (point_type == 3) {
            point.x = center_mass_x + L3_offset * vector_x;
            point.y = center_mass_y + L3_offset * vector_y;
            point.z = center_mass_z + L3_offset * vector_z;
            point.color[0] = 0.0f; point.color[1] = 0.0f; point.color[2] = 1.0f;
        } else if (point_type == 4) {
            point.x = primary->x + distance_R * (cos_60 * vector_x + sin_60 * perpendicular_x);
            point.y = primary->y + distance_R * (cos_60 * vector_y + sin_60 * perpendicular_y);
            point.z = primary->z + distance_R * (cos_60 * vector_z + sin_60 * perpendicular_z);
            point.color[0] = 1.0f; point.color[1] = 1.0f; point.color[2] = 0.0f;
        } else if (point_type == 5) {
            point.x = primary->x + distance_R * (cos_60 * vector_x - sin_60 * perpendicular_x);
            point.y = primary->y + distance_R * (cos_60 * vector_y - sin_60 * perpendicular_y);
            point.z = primary->z + distance_R * (cos_60 * vector_z - sin_60 * perpendicular_z);
            point.color[0] = 1.0f; point.color[1] = 0.0f; point.color[2] = 1.0f;
        }
        
        double dx = point.x - primary->x;
        double dy = point.y - primary->y;
        double dz = point.z - primary->z;
        point.distanceFromPrimary = std::sqrt(dx * dx + dy * dy + dz * dz);
        points.push_back(point);
    }
    return points;
}

void updateLagrangePoints(std::vector<LagrangePoint>& points, Planet* primary, Planet* secondary) {
    if (!primary || !secondary) return;
    
    double delta_x = secondary->x - primary->x;
    double delta_y = secondary->y - primary->y;
    double delta_z = secondary->z - primary->z;
    double current_R = std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
    
    double center_mass_x, center_mass_y, center_mass_z;
    calculateCenterOfMass(primary, secondary, center_mass_x, center_mass_y, center_mass_z);
    
    double vector_x = secondary->x - center_mass_x;
    double vector_y = secondary->y - center_mass_y;
    double vector_z = secondary->z - center_mass_z;
    double vector_length = std::sqrt(vector_x * vector_x + vector_y * vector_y + vector_z * vector_z);
    
    if (vector_length > 0) {
        vector_x /= vector_length;
        vector_y /= vector_length;
        vector_z /= vector_length;
    }
    
    double perpendicular_x = -vector_z;
    double perpendicular_y = 0;
    double perpendicular_z = vector_x;
    double perpendicular_length = std::sqrt(perpendicular_x * perpendicular_x + perpendicular_y * perpendicular_y + perpendicular_z * perpendicular_z);
    
    if (perpendicular_length > 0) {
        perpendicular_x /= perpendicular_length;
        perpendicular_y /= perpendicular_length;
        perpendicular_z /= perpendicular_length;
    }
    
    double sin_60 = std::sqrt(3.0) / 2.0;
    double cos_60 = 0.5;
    
    for (auto& point : points) {
        if (point.system != primary->name + "-" + secondary->name) continue;
        
        double alpha_root = std::cbrt(point.mu / 3.0);
        double L1_offset = current_R * (1.0 - alpha_root);
        double L2_offset = current_R * (1.0 + alpha_root);
        double L3_offset = -current_R * (1.0 + (5.0 / 12.0) * point.mu);
        
        if (point.type == 1) {
            point.x = center_mass_x + L1_offset * vector_x;
            point.y = center_mass_y + L1_offset * vector_y;
            point.z = center_mass_z + L1_offset * vector_z;
        } else if (point.type == 2) {
            point.x = center_mass_x + L2_offset * vector_x;
            point.y = center_mass_y + L2_offset * vector_y;
            point.z = center_mass_z + L2_offset * vector_z;
        } else if (point.type == 3) {
            point.x = center_mass_x + L3_offset * vector_x;
            point.y = center_mass_y + L3_offset * vector_y;
            point.z = center_mass_z + L3_offset * vector_z;
        } else if (point.type == 4) {
            point.x = primary->x + current_R * (cos_60 * vector_x + sin_60 * perpendicular_x);
            point.y = primary->y + current_R * (cos_60 * vector_y + sin_60 * perpendicular_y);
            point.z = primary->z + current_R * (cos_60 * vector_z + sin_60 * perpendicular_z);
        } else if (point.type == 5) {
            point.x = primary->x + current_R * (cos_60 * vector_x - sin_60 * perpendicular_x);
            point.y = primary->y + current_R * (cos_60 * vector_y - sin_60 * perpendicular_y);
            point.z = primary->z + current_R * (cos_60 * vector_z - sin_60 * perpendicular_z);
        }
    }
}

Planet* createLagrangeSatellite(const LagrangePoint& point, const std::string& systemName) {
    Planet* satellite = new Planet();
    satellite->name = "Satellite_L" + std::to_string(point.type);
    satellite->mass = 1000.0;
    satellite->x = point.x;
    satellite->y = point.y;
    satellite->z = point.z;
    satellite->vx = 0.0;
    satellite->vy = 0.0;
    satellite->vz = 0.0;
    satellite->ax = 0.0;
    satellite->ay = 0.0;
    satellite->az = 0.0;
    satellite->is_satellite = true;
    satellite->is_lagrange_satellite = true;
    satellite->lagrange_point_type = point.type;
    satellite->parent = systemName;
    satellite->radius = RADIUS_SATELLITE;
    return satellite;
}

void create_sphere(Planet* planet, float radius, float r, float g, float b) {
    planet->radius = radius;
    int s = 20;  // Уменьшил для производительности
    
    for (int i = 0; i <= s; ++i) {
        float phi = M_PI * i / s;
        for (int j = 0; j <= s; ++j) {
            float theta = 2.0f * M_PI * j / s;
            Vertex v;
            v.x = radius * sin(phi) * cos(theta);
            v.y = radius * cos(phi);
            v.z = radius * sin(phi) * sin(theta);
            v.r = r;
            v.g = g;
            v.b = b;
            (planet->ver).push_back(v);
        }
    }
    
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            int first = i * (s + 1) + j;
            int second = first + 1;
            int third = first + (s + 1);
            int fourth = third + 1;
            
            (planet->ind).push_back(first);
            (planet->ind).push_back(second);
            (planet->ind).push_back(third);
            
            (planet->ind).push_back(second);
            (planet->ind).push_back(fourth);
            (planet->ind).push_back(third);
        }
    }
}

void draw_planet(Planet* planet) {
    glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < (planet->ind).size(); ++i) {
        Vertex* v = &(planet->ver[planet->ind[i]]);
        glColor3f(v->r, v->g, v->b);
        glVertex3f(v->x, v->y, v->z);
    }
    glEnd();
}

void draw_orbit(float a, float eccentricity, float tilt, float r, float g, float b) {
    glColor3f(r, g, b);
    glBegin(GL_LINE_LOOP);
    
    for (int i = 0; i < 100; i++) {
        float angle = 2.0f * M_PI * i / 100.0f;
        float r_orbit = a * (1 - eccentricity * eccentricity) / (1 + eccentricity * cos(angle));
        r_orbit *= SCALE;
        
        float x = r_orbit * cos(angle);
        float z = r_orbit * sin(angle);
        float y = 0.0f;
        
        // Учет наклона орбиты
        float y_tilted = z * sin(tilt);
        float z_tilted = z * cos(tilt);
        
        glVertex3f(x, y_tilted, z_tilted);
    }
    
    glEnd();
}

// Метод Leapfrog для стабильного интегрирования
void update_physics(Planet* planet, std::vector<Planet*>* all_planets, double dt) {
    if (planet->is_lagrange_satellite) {
        return;
    }
    
    // Leapfrog метод (сохраняет энергию)
    // 1. Обновляем скорость на половину шага
    planet->vx += planet->ax * dt * 0.5;
    planet->vy += planet->ay * dt * 0.5;
    planet->vz += planet->az * dt * 0.5;
    
    // 2. Обновляем позицию
    planet->x += planet->vx * dt;
    planet->y += planet->vy * dt;
    planet->z += planet->vz * dt;
    
    // 3. Вычисляем новые ускорения
    double total_ax = 0.0, total_ay = 0.0, total_az = 0.0;
    for (Planet* other : *all_planets) {
        if (planet == other) continue;
        
        double dx = other->x - planet->x;
        double dy = other->y - planet->y;
        double dz = other->z - planet->z;
        double r2 = dx * dx + dy * dy + dz * dz;
        
        // Избегаем деления на ноль
        if (r2 < 1e10) continue;
        
        double r = sqrt(r2);
        double force = G * other->mass / (r2 * r);  // F/m = G*M/r^2 * (dx/r)
        total_ax += force * dx;
        total_ay += force * dy;
        total_az += force * dz;
    }
    
    // 4. Завершаем обновление скорости
    planet->ax = total_ax;
    planet->ay = total_ay;
    planet->az = total_az;
    
    planet->vx += planet->ax * dt * 0.5;
    planet->vy += planet->ay * dt * 0.5;
    planet->vz += planet->az * dt * 0.5;
}

void calculateLagrangeForces(Planet* satellite, std::vector<Planet*>* all_planets) {
    double total_ax = 0.0, total_ay = 0.0, total_az = 0.0;
    for (Planet* other : *all_planets) {
        if (satellite == other) continue;
        
        double dx = other->x - satellite->x;
        double dy = other->y - satellite->y;
        double dz = other->z - satellite->z;
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < 1e-12) continue;
        
        double r = sqrt(r2);
        double force = G * other->mass / (r2 * r);
        total_ax += force * dx;
        total_ay += force * dy;
        total_az += force * dz;
    }
    
    satellite->ax = total_ax;
    satellite->ay = total_ay;
    satellite->az = total_az;
}

void render_planet(Planet* planet) {
    glPushMatrix();
    glTranslatef(float(planet->x * SCALE), float(planet->y * SCALE), float(planet->z * SCALE));
    draw_planet(planet);
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

// Вычисление начальных ускорений
void computeInitialAccelerations(std::vector<Planet*>& planets) {
    for (Planet* planet : planets) {
        if (planet->name == "Sun") continue;
        
        double total_ax = 0.0, total_ay = 0.0, total_az = 0.0;
        for (Planet* other : planets) {
            if (planet == other) continue;
            
            double dx = other->x - planet->x;
            double dy = other->y - planet->y;
            double dz = other->z - planet->z;
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < 1e10) continue;
            
            double r = sqrt(r2);
            double force = G * other->mass / (r2 * r);
            total_ax += force * dx;
            total_ay += force * dy;
            total_az += force * dz;
        }
        
        planet->ax = total_ax;
        planet->ay = total_ay;
        planet->az = total_az;
    }
}

void initialize_planets(std::vector<Planet*>* planets) {
    // Солнце
    Planet* sun = new Planet();
    sun->name = "Sun";
    sun->mass = MASS_SUN;
    sun->x = 0.0; sun->y = 0.0; sun->z = 0.0;
    sun->vx = 0.0; sun->vy = 0.0; sun->vz = 0.0;
    sun->ax = 0.0; sun->ay = 0.0; sun->az = 0.0;
    sun->orbit_eccentricity = EX_SUN;
    sun->a = A_SUN;
    sun->orbit_tilt = 0.0f;
    create_sphere(sun, RADIUS_SUN, COLOR_SUN[0], COLOR_SUN[1], COLOR_SUN[2]);
    planets->push_back(sun);
    
    // Функция для инициализации планеты
    auto initPlanet = [&](const std::string& name, double mass, double a, double e, 
                         double tilt, double perihelion, float radius, const float* color, double trueAnomaly = 0.0) {
        Planet* planet = new Planet();
        planet->name = name;
        planet->mass = mass;
        planet->orbit_eccentricity = e;
        planet->a = a;
        planet->orbit_tilt = tilt;
        planet->true_anomaly = trueAnomaly;
        
        // Положение в орбитальной плоскости
        double r = a * (1 - e * e) / (1 + e * cos(trueAnomaly));
        planet->x = r * cos(trueAnomaly);
        planet->z = r * sin(trueAnomaly);
        planet->y = 0.0;
        
        // Учет наклона орбиты (вращение вокруг оси X)
        double y_temp = planet->y;
        double z_temp = planet->z;
        planet->y = y_temp * cos(tilt) - z_temp * sin(tilt);
        planet->z = y_temp * sin(tilt) + z_temp * cos(tilt);
        
        // Скорость в орбитальной плоскости
        double v = calculateOrbitalVelocity(a, e, trueAnomaly, MASS_SUN);
        double vx_orb = -v * sin(trueAnomaly);
        double vz_orb = v * cos(trueAnomaly);
        
        planet->vx = vx_orb;
        planet->vy = 0.0;
        planet->vz = vz_orb;
        
        // Учет наклона для скорости
        double vy_temp = planet->vy;
        double vz_temp = planet->vz;
        planet->vy = vy_temp * cos(tilt) - vz_temp * sin(tilt);
        planet->vz = vy_temp * sin(tilt) + vz_temp * cos(tilt);
        
        planet->ax = 0.0;
        planet->ay = 0.0;
        planet->az = 0.0;
        
        create_sphere(planet, radius, color[0], color[1], color[2]);
        planets->push_back(planet);
    };
    
    // Планеты с разными начальными фазами для лучшей визуализации
    initPlanet("Mercury", MASS_MERCURY, A_MERCURY, EX_MERCURY, TILT_MERCURY, 
               PERIHELION_MERCURY, RADIUS_MERCURY, COLOR_MERCURY, 0.0);
    initPlanet("Venus", MASS_VENUS, A_VENUS, EX_VENUS, TILT_VENUS, 
               PERIHELION_VENUS, RADIUS_VENUS, COLOR_VENUS, M_PI/3);
    initPlanet("Earth", MASS_EARTH, A_EARTH, EX_EARTH, TILT_EARTH, 
               PERIHELION_EARTH, RADIUS_EARTH, COLOR_EARTH, M_PI/2);
    initPlanet("Mars", MASS_MARS, A_MARS, EX_MARS, TILT_MARS, 
               PERIHELION_MARS, RADIUS_MARS, COLOR_MARS, M_PI);
    initPlanet("Jupiter", MASS_JUPITER, A_JUPITER, EX_JUPITER, TILT_JUPITER, 
               PERIHELION_JUPITER, RADIUS_JUPITER, COLOR_JUPITER, 3*M_PI/2);
    initPlanet("Saturn", MASS_SATURN, A_SATURN, EX_SATURN, TILT_SATURN, 
               PERIHELION_SATURN, RADIUS_SATURN, COLOR_SATURN, M_PI/6);
    initPlanet("Uran", MASS_URAN, A_URAN, EX_URAN, TILT_URAN, 
               PERIHELION_URAN, RADIUS_URAN, COLOR_URAN, M_PI/4);
    initPlanet("Neptun", MASS_NEPTUN, A_NEPTUN, EX_NEPTUN, TILT_NEPTUN, 
               PERIHELION_NEPTUN, RADIUS_NEPTUN, COLOR_NEPTUN, M_PI/2);
    
    // Луна (относительно Земли)
    Planet* earth = planets->at(3);  // Земля
    Planet* luna = new Planet();
    luna->name = "luna";
    luna->mass = MASS_LUNA;
    luna->is_satellite = true;
    luna->parent = "Earth";
    
    // Орбита Луны вокруг Земли
    double moon_distance = LUNA_ORBIT_RADIUS;
    double moon_v = sqrt(G * earth->mass / moon_distance);
    
    // Положение Луны относительно Земли
    luna->x = earth->x + moon_distance;
    luna->y = earth->y;
    luna->z = earth->z;
    
    // Скорость Луны (перпендикулярно радиус-вектору)
    luna->vx = earth->vx;
    luna->vy = earth->vy;
    luna->vz = earth->vz + moon_v;
    
    luna->ax = 0.0;
    luna->ay = 0.0;
    luna->az = 0.0;
    
    create_sphere(luna, RADIUS_LUNA, COLOR_LUNA[0], COLOR_LUNA[1], COLOR_LUNA[2]);
    planets->push_back(luna);
}

struct Camera {
    float posX, posY, posZ;
    float targetX, targetY, targetZ;
    float upX, upY, upZ;
    float distance;
    float angleX, angleY;
    int followingPlanet = -1;
    float followDistance = 50.0f;
    float followAngleX = 0.0f;
    float followAngleY = 0.3f;
    
    Camera() {
        posX = 30.0f; posY = 80.0f; posZ = 30.0f;
        targetX = 0.0f; targetY = 0.0f; targetZ = 0.0f;
        upX = 0.0f; upY = 1.0f; upZ = 0.0f;
        distance = 100.0f;
        angleX = 0.0f;
        angleY = 0.3f;
        update();
    }
    
    void update() {
        posX = targetX + distance * cos(angleY) * sin(angleX);
        posY = targetY + distance * sin(angleY);
        posZ = targetZ + distance * cos(angleY) * cos(angleX);
    }
    
    void updateFollowing(Planet* planet) {
        if (!planet) return;
        targetX = float(planet->x * SCALE);
        targetY = float(planet->y * SCALE);
        targetZ = float(planet->z * SCALE);
        angleX = followAngleX;
        angleY = followAngleY;
        distance = followDistance;
        update();
    }
};

Camera camera;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        float cameraSpeed = 0.5f;
        
        if (camera.followingPlanet >= 0) {
            if (key == GLFW_KEY_UP) camera.followAngleY += 0.05f;
            if (key == GLFW_KEY_DOWN) camera.followAngleY -= 0.05f;
            if (key == GLFW_KEY_LEFT) camera.followAngleX -= 0.05f;
            if (key == GLFW_KEY_RIGHT) camera.followAngleX += 0.05f;
            if (key == GLFW_KEY_W) camera.followDistance -= 1.0f;
            if (key == GLFW_KEY_S) camera.followDistance += 1.0f;
        } else {
            if (key == GLFW_KEY_UP) camera.targetY += cameraSpeed;
            if (key == GLFW_KEY_DOWN) camera.targetY -= cameraSpeed;
            if (key == GLFW_KEY_LEFT) camera.targetX -= cameraSpeed;
            if (key == GLFW_KEY_RIGHT) camera.targetX += cameraSpeed;
            if (key == GLFW_KEY_W) camera.targetZ -= cameraSpeed;
            if (key == GLFW_KEY_S) camera.targetZ += cameraSpeed;
        }
        
        if (key == GLFW_KEY_R) {
            camera.followingPlanet = -1;
            camera.posX = 30.0f; camera.posY = 80.0f; camera.posZ = 30.0f;
            camera.targetX = 0.0f; camera.targetY = 0.0f; camera.targetZ = 0.0f;
            camera.distance = 100.0f;
            camera.angleX = 0.0f;
            camera.angleY = 0.3f;
            camera.update();
        }
        
        if (key == GLFW_KEY_F) {
            camera.followingPlanet = -1;
            camera.posX = 0.0f; camera.posY = 0.0f; camera.posZ = 50.0f;
            camera.targetX = 0.0f; camera.targetY = 0.0f; camera.targetZ = 0.0f;
            camera.angleX = 0.0f; camera.angleY = 0.0f;
            camera.distance = 50.0f;
            camera.update();
        }
        
        if (key == GLFW_KEY_Z) {
            SIMULATION_SPEED_MULTIPLIER *= 2.0;
            std::cout << "Скорость симуляции: " << SIMULATION_SPEED_MULTIPLIER << "x" << std::endl;
        }
        
        if (key == GLFW_KEY_X) {
            SIMULATION_SPEED_MULTIPLIER /= 2.0;
            std::cout << "Скорость симуляции: " << SIMULATION_SPEED_MULTIPLIER << "x" << std::endl;
        }
        
        if (key == GLFW_KEY_C) {
            SIMULATION_SPEED_MULTIPLIER = (SIMULATION_SPEED_MULTIPLIER == 0.0) ? 1.0 : 0.0;
            std::cout << (SIMULATION_SPEED_MULTIPLIER == 0.0 ? "Пауза" : "Возобновлено") << std::endl;
        }
        
        if (key == GLFW_KEY_L) {
            showLagrangePoints = !showLagrangePoints;
        }
        
        if (key == GLFW_KEY_P) {
            lagrangeSatelliteActive = !lagrangeSatelliteActive;
        }
        
        if (key == GLFW_KEY_F1) {
            selectedLagrangeSystem = 0;
            std::cout << "Выбрана система: Солнце-Земля" << std::endl;
        } else if (key == GLFW_KEY_F2) {
            selectedLagrangeSystem = 1;
            std::cout << "Выбрана система: Земля-Луна" << std::endl;
        } else if (key == GLFW_KEY_F3) {
            selectedLagrangeSystem = 2;
            std::cout << "Выбрана система: Солнце-Юпитер" << std::endl;
        }
        
        if (key >= GLFW_KEY_1 && key <= GLFW_KEY_5) {
            selectedLagrangePoint = key - GLFW_KEY_1 + 1;
            std::cout << "Выбрана точка L" << selectedLagrangePoint << std::endl;
        }
        
        // Клавиши для слежения за планетами
        if (key == GLFW_KEY_F4) camera.followingPlanet = 0;
        else if (key == GLFW_KEY_F5) camera.followingPlanet = 1;
        else if (key == GLFW_KEY_F6) camera.followingPlanet = 2;
        else if (key == GLFW_KEY_F7) camera.followingPlanet = 3;
        else if (key == GLFW_KEY_F8) camera.followingPlanet = 4;
        else if (key == GLFW_KEY_F9) camera.followingPlanet = 5;
        else if (key == GLFW_KEY_F10) camera.followingPlanet = 6;
        else if (key == GLFW_KEY_F11) camera.followingPlanet = 7;
        else if (key == GLFW_KEY_F12) camera.followingPlanet = 8;
        
        camera.update();
    }
}

struct Tracker {
    double last_angle = 0;
    double last_time = 0;
    int revolutions = 0;
    double period_sum = 0;
    double last_period = 0;
    
    void update(double x, double z, double time) {
        double current_angle = atan2(z, x);
        
        if (last_time == 0) {
            last_angle = current_angle;
            last_time = time;
            return;
        }
        
        double angle_diff = current_angle - last_angle;
        if (angle_diff < -M_PI) angle_diff += 2 * M_PI;
        if (angle_diff > M_PI) angle_diff -= 2 * M_PI;
        
        static double accumulated_angle = 0;
        accumulated_angle += angle_diff;
        
        if (fabs(accumulated_angle) >= 2 * M_PI) {
            accumulated_angle -= copysign(2 * M_PI, accumulated_angle);
            revolutions++;
            double current_period = time - last_time;
            period_sum += current_period;
            last_period = current_period;
            last_time = time;
        }
        
        last_angle = current_angle;
    }
    
    double get_period() {
        if (revolutions == 0) return 0;
        return period_sum / revolutions;
    }
    
    bool ready() {
        return revolutions >= 1;
    }
};

void print_data(Planet* planet, double t) {
    std::cout << planet->name << ":" << std::endl;
    std::cout << "  x: " << planet->x/AU << " AU" << std::endl;
    std::cout << "  y: " << planet->y/AU << " AU" << std::endl;
    std::cout << "  z: " << planet->z/AU << " AU" << std::endl;
    std::cout << "  v: " << sqrt(planet->vx*planet->vx + planet->vy*planet->vy + planet->vz*planet->vz) << " m/s" << std::endl;
    std::cout << "  a: " << planet->a/AU << " AU" << std::endl;
    
    if (t > 0) {
        double period_days = t / (24 * 3600);
        double ratio = (t * t) / (planet->a * planet->a * planet->a);
        std::cout << "  Период: " << period_days << " дней" << std::endl;
        std::cout << "  T²/a³: " << ratio << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    
    GLFWwindow* window = glfwCreateWindow(1200, 800, "Солнечная система", NULL, NULL);
    if (!window) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent(window);
    
    if (glewInit() != GLEW_OK) {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return -1;
    }
    
    glfwSetKeyCallback(window, key_callback);
    
    std::vector<Planet*> planets;
    initialize_planets(&planets);
    
    // Вычисляем начальные ускорения
    computeInitialAccelerations(planets);
    
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.05f, 1.0f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    double last_time = glfwGetTime();
    double accumulator = 0.0;
    double sim_time = 0.0;
    
    Tracker mercury_tracker, venus_tracker, earth_tracker, mars_tracker;
    Tracker jupiter_tracker, saturn_tracker, uran_tracker, neptun_tracker;
    
    double period_mercury = 0, period_venus = 0, period_earth = 0, period_mars = 0;
    double period_jupiter = 0, period_saturn = 0, period_uran = 0, period_neptun = 0;
    
    std::vector<LagrangePoint> allLagrangePoints;
    Planet* lagrangeSatellite = nullptr;
    
    Planet* sun = nullptr;
    Planet* earth = nullptr;
    Planet* moon = nullptr;
    Planet* jupiter = nullptr;
    
    for (Planet* p : planets) {
        if (p->name == "Sun") sun = p;
        if (p->name == "Earth") earth = p;
        if (p->name == "luna") moon = p;
        if (p->name == "Jupiter") jupiter = p;
    }
    
    if (sun && earth) {
        std::vector<LagrangePoint> sunEarthPoints = calculateLagrangePoints(sun, earth);
        allLagrangePoints.insert(allLagrangePoints.end(), sunEarthPoints.begin(), sunEarthPoints.end());
    }
    if (earth && moon) {
        std::vector<LagrangePoint> earthMoonPoints = calculateLagrangePoints(earth, moon);
        allLagrangePoints.insert(allLagrangePoints.end(), earthMoonPoints.begin(), earthMoonPoints.end());
    }
    if (sun && jupiter) {
        std::vector<LagrangePoint> sunJupiterPoints = calculateLagrangePoints(sun, jupiter);
        allLagrangePoints.insert(allLagrangePoints.end(), sunJupiterPoints.begin(), sunJupiterPoints.end());
    }
    
    while (!glfwWindowShouldClose(window)) {
        double current_time = glfwGetTime();
        double delta_time = current_time - last_time;
        last_time = current_time;
        
        // Фиксированный шаг времени для физики
        accumulator += delta_time * SIMULATION_SPEED_MULTIPLIER;
        
        while (accumulator >= 0.01) {  // Шаг 0.01 секунды реального времени
            double dt = TIME_STEP * 0.01;  // 36 секунд симуляции за шаг
            
            for (Planet* planet : planets) {
                update_physics(planet, &planets, dt);
            }
            
            if (showLagrangePoints) {
                if (sun && earth) updateLagrangePoints(allLagrangePoints, sun, earth);
                if (earth && moon) updateLagrangePoints(allLagrangePoints, earth, moon);
                if (sun && jupiter) updateLagrangePoints(allLagrangePoints, sun, jupiter);
            }
            
            if (lagrangeSatelliteActive && !allLagrangePoints.empty()) {
                std::string targetSystem;
                if (selectedLagrangeSystem == 0) targetSystem = "Sun-Earth";
                else if (selectedLagrangeSystem == 1) targetSystem = "Earth-luna";
                else if (selectedLagrangeSystem == 2) targetSystem = "Sun-Jupiter";
                
                for (const auto& point : allLagrangePoints) {
                    if (point.system == targetSystem && point.type == selectedLagrangePoint) {
                        if (!lagrangeSatellite) {
                            lagrangeSatellite = createLagrangeSatellite(point, targetSystem);
                            create_sphere(lagrangeSatellite, RADIUS_SATELLITE, 
                                        COLOR_SATELLITE[0], COLOR_SATELLITE[1], COLOR_SATELLITE[2]);
                            planets.push_back(lagrangeSatellite);
                        } else {
                            lagrangeSatellite->x = point.x;
                            lagrangeSatellite->y = point.y;
                            lagrangeSatellite->z = point.z;
                        }
                        break;
                    }
                }
                
                if (lagrangeSatellite) {
                    calculateLagrangeForces(lagrangeSatellite, &planets);
                }
            } else if (lagrangeSatellite) {
                auto it = std::find(planets.begin(), planets.end(), lagrangeSatellite);
                if (it != planets.end()) {
                    planets.erase(it);
                    delete lagrangeSatellite;
                    lagrangeSatellite = nullptr;
                }
            }
            
            sim_time += dt;
            accumulator -= 0.01;
        }
        
        // Обновление трекеров периодов
        for (Planet* planet : planets) {
            if (planet->name == "Mercury") {
                mercury_tracker.update(planet->x, planet->z, sim_time);
                if (mercury_tracker.ready()) period_mercury = mercury_tracker.get_period();
            } else if (planet->name == "Venus") {
                venus_tracker.update(planet->x, planet->z, sim_time);
                if (venus_tracker.ready()) period_venus = venus_tracker.get_period();
            } else if (planet->name == "Earth") {
                earth_tracker.update(planet->x, planet->z, sim_time);
                if (earth_tracker.ready()) period_earth = earth_tracker.get_period();
            } else if (planet->name == "Mars") {
                mars_tracker.update(planet->x, planet->z, sim_time);
                if (mars_tracker.ready()) period_mars = mars_tracker.get_period();
            } else if (planet->name == "Jupiter") {
                jupiter_tracker.update(planet->x, planet->z, sim_time);
                if (jupiter_tracker.ready()) period_jupiter = jupiter_tracker.get_period();
            } else if (planet->name == "Saturn") {
                saturn_tracker.update(planet->x, planet->z, sim_time);
                if (saturn_tracker.ready()) period_saturn = saturn_tracker.get_period();
            } else if (planet->name == "Uran") {
                uran_tracker.update(planet->x, planet->z, sim_time);
                if (uran_tracker.ready()) period_uran = uran_tracker.get_period();
            } else if (planet->name == "Neptun") {
                neptun_tracker.update(planet->x, planet->z, sim_time);
                if (neptun_tracker.ready()) period_neptun = neptun_tracker.get_period();
            }
        }
        
        if (camera.followingPlanet >= 0 && camera.followingPlanet < planets.size()) {
            camera.updateFollowing(planets[camera.followingPlanet]);
        } else {
            camera.update();
        }
        
        // Отрисовка
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)width / (float)height;
        gluPerspective(60.0f, aspect, CLIP_NEAR, CLIP_FAR);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(camera.posX, camera.posY, camera.posZ,
                  camera.targetX, camera.targetY, camera.targetZ,
                  camera.upX, camera.upY, camera.upZ);
        
        draw_coordinate_system();
        
        glDisable(GL_DEPTH_TEST);
        glLineWidth(0.5f);
        
        for (Planet* planet : planets) {
            if (planet->name == "Mercury") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_MERCURY, 0.5f, 0.5f, 0.5f);
            else if (planet->name == "Venus") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_VENUS, 1.0f, 0.8f, 0.4f);
            else if (planet->name == "Earth") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_EARTH, 0.3f, 0.3f, 0.8f);
            else if (planet->name == "Mars") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_MARS, 0.8f, 0.3f, 0.2f);
            else if (planet->name == "Jupiter") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_JUPITER, 1.0f, 0.4f, 0.5f);
            else if (planet->name == "Saturn") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_SATURN, 0.2f, 0.7f, 1.0f);
            else if (planet->name == "Uran") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_URAN, 0.45f, 0.85f, 0.15f);
            else if (planet->name == "Neptun") draw_orbit(planet->a, planet->orbit_eccentricity, TILT_NEPTUN, 0.12f, 0.95f, 1.0f);
        }
        
        if (showLagrangePoints) {
            for (const auto& point : allLagrangePoints) {
                glPushMatrix();
                glTranslatef(float(point.x * SCALE), float(point.y * SCALE), float(point.z * SCALE));
                glColor3f(point.color[0], point.color[1], point.color[2]);
                glPointSize(5.0f);
                glBegin(GL_POINTS);
                glVertex3f(0, 0, 0);
                glEnd();
                glPopMatrix();
            }
        }
        
        glEnable(GL_DEPTH_TEST);
        glLineWidth(1.0f);
        
        for (Planet* planet : planets) {
            render_planet(planet);
        }
        
        static int frame_count = 0;
        if (frame_count++ % 60 == 0) {
            system("clear");
            std::cout << "=== СОЛНЕЧНАЯ СИСТЕМА ===" << std::endl;
            std::cout << "Время симуляции: " << sim_time / (365.25 * 24 * 3600) << " лет" << std::endl;
            std::cout << "Скорость: " << SIMULATION_SPEED_MULTIPLIER << "x" << std::endl;
            std::cout << std::endl;
            
            if (showLagrangePoints) {
                std::cout << "Точки Лагранжа: ВКЛ" << std::endl;
                if (lagrangeSatelliteActive && lagrangeSatellite) {
                    std::cout << "Спутник в точке L" << selectedLagrangePoint;
                    if (selectedLagrangeSystem == 0) std::cout << " (Солнце-Земля)";
                    else if (selectedLagrangeSystem == 1) std::cout << " (Земля-Луна)";
                    else if (selectedLagrangeSystem == 2) std::cout << " (Солнце-Юпитер)";
                    std::cout << std::endl;
                }
            }
            
            std::cout << std::endl << "=== ОРБИТАЛЬНЫЕ ПАРАМЕТРЫ ===" << std::endl;
            for (Planet* planet : planets) {
                if (planet->name == "Mercury") print_data(planet, period_mercury);
                else if (planet->name == "Venus") print_data(planet, period_venus);
                else if (planet->name == "Earth") print_data(planet, period_earth);
                else if (planet->name == "Mars") print_data(planet, period_mars);
                else if (planet->name == "Jupiter") print_data(planet, period_jupiter);
                else if (planet->name == "Saturn") print_data(planet, period_saturn);
                else if (planet->name == "Uran") print_data(planet, period_uran);
                else if (planet->name == "Neptun") print_data(planet, period_neptun);
            }
        }
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    for (Planet* p : planets) delete p;
    if (lagrangeSatellite) delete lagrangeSatellite;
    
    glfwDestroyWindow(window);
    glfwTerminate();
    
    std::cout << std::endl << "Симуляция завершена" << std::endl;
    return 0;
}