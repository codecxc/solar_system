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
const double AU=1.496e11;
const double G=6.67430e-11;
const double SCALE=5000.0/AU;
const double TIME_SCALE=86400;
const double RADIUS_SCALE=1.0;

const float MAX_DISTANCE=100.0f;
const float CLIP_NEAR=0.0001f;
const float CLIP_FAR=5000.0f;

const double MASS_SUN=1.9891e30;
const double MASS_MERCURY=3.302e23;
const double MASS_VENUS=4.868e24;
const double MASS_EARTH=5.9736e24;
const double MASS_MARS=6.4185e23;
const double MASS_JUPITER=1.8986e27;
const double MASS_SATURN=5.68460e26;
const double MASS_URAN=8.6832e25;
const double MASS_NEPTUN=1.02430e26;
const double MASS_LUNA=7.342e22;

double VEL_MERCURY,VEL_VENUS,VEL_EARTH,VEL_MARS,VEL_JUPITER,VEL_SATURN,VEL_URAN,VEL_NEPTUN;
double VEL_LUNA;

const float EX_SUN=0.0f;
const float EX_MERCURY=0.206f;
const float EX_VENUS=0.007f;
const float EX_EARTH=0.017f;
const float EX_MARS=0.093f;
const float EX_JUPITER=0.049f;
const float EX_SATURN=0.057f;
const float EX_URAN=0.046f;
const float EX_NEPTUN=0.011f;

const float A_SUN=0.0f;
const float A_MERCURY=0.387*AU;
const float A_VENUS=0.723*AU;
const float A_EARTH=1.000*AU;
const float A_MARS=1.524*AU;
const float A_JUPITER=5.2044*AU;
const float A_SATURN=9.5826*AU;
const float A_URAN=19.21840*AU;
const float A_NEPTUN=30.11000*AU;

const float TILT_MERCURY=7.01f*M_PI/180.0f;
const float TILT_VENUS=3.39f*M_PI/180.0f;
const float TILT_EARTH=0.0f;
const float TILT_MARS=1.85f*M_PI/180.0f;
const float TILT_JUPITER=1.31f*M_PI/180.0f;
const float TILT_SATURN=2.49f*M_PI/180.0f;
const float TILT_URAN=0.77f*M_PI/180.0f;
const float TILT_NEPTUN=1.77f*M_PI/180.0f;

const double PERIHELION_MERCURY=A_MERCURY*(1-EX_MERCURY);
const double PERIHELION_VENUS=A_VENUS*(1-EX_VENUS);
const double PERIHELION_EARTH=A_EARTH*(1-EX_EARTH);
const double PERIHELION_MARS=A_MARS*(1-EX_MARS);
const double PERIHELION_JUPITER=A_JUPITER*(1-EX_JUPITER);
const double PERIHELION_SATURN=A_SATURN*(1-EX_SATURN);
const double PERIHELION_URAN=A_URAN*(1-EX_URAN);
const double PERIHELION_NEPTUN=A_NEPTUN*(1-EX_NEPTUN);

const float RADIUS_SUN=695780000.0f*SCALE*RADIUS_SCALE;
const float RADIUS_MERCURY=2439700.0f*SCALE*RADIUS_SCALE;
const float RADIUS_VENUS=6051800.0f*SCALE*RADIUS_SCALE;
const float RADIUS_EARTH=6371000.0f*SCALE*RADIUS_SCALE;
const float RADIUS_MARS=3389500.0f*SCALE*RADIUS_SCALE;
const float RADIUS_JUPITER=69911000.0f*SCALE*RADIUS_SCALE;
const float RADIUS_SATURN=58232000.0f*SCALE*RADIUS_SCALE;
const float RADIUS_URAN=25362000.0f*SCALE*RADIUS_SCALE;
const float RADIUS_NEPTUN=24622000.0f*SCALE*RADIUS_SCALE;
const float RADIUS_LUNA=1737400.0f*SCALE*RADIUS_SCALE;
const float RADIUS_SATELLITE=50000000.0f*SCALE*RADIUS_SCALE;

const float COLOR_SUN[3]={1.0f,1.0f,0.0f};
const float COLOR_MERCURY[3]={0.7f,0.7f,0.7f};
const float COLOR_VENUS[3]={1.0f,0.8f,0.4f};
const float COLOR_EARTH[3]={0.0f,0.5f,1.0f};
const float COLOR_MARS[3]={1.0f,0.3f,0.1f};
const float COLOR_JUPITER[3]={0.9f,0.7f,0.5f};
const float COLOR_SATURN[3]={0.95f,0.85f,0.6f};
const float COLOR_NEPTUN[3]={0.1f,0.2f,0.0f};
const float COLOR_URAN[3]={0.15f,0.35f,0.6f};
const float COLOR_LUNA[3]={1.0f,1.0f,1.0f};
const float COLOR_SATELLITE[3]={1.0f,0.0f,1.0f};

const double LUNA_ORBIT_RADIUS=3.844e8+6371e3;

bool showLagrangePoints=true;
bool lagrangeSatelliteActive=false;
int selectedLagrangeSystem=0;
int selectedLagrangePoint=1;

struct Vertex {
	float x,y,z;
	float r,g,b;
};

struct LagrangePoint {
	double x,y,z;
	int type;
	float color[3];
	std::string system;
	double distanceFromPrimary;
	double R;
	double mu;
};

struct Planet {
	std::vector<Vertex>ver;
	std::vector<unsigned int>ind;
	float radius;
	double mass;
	double x,y,z;
	double vx,vy,vz;
	double ax,ay,az;
	std::string name;
	float a;
	float orbit_eccentricity;
	float orbit_tilt;
	bool is_satellite=false;
	bool is_lagrange_satellite=false;
	std::string parent;
	int lagrange_point_type=0;
};

void calculateCenterOfMass(Planet*primary,Planet*secondary,double&com_x,double&com_y,double&com_z) {
	if(!primary||!secondary) {
		com_x=com_y=com_z=0.0;
		return;
	}
	double totalMass=primary->mass+secondary->mass;
	com_x=(primary->mass*primary->x+secondary->mass*secondary->x)/totalMass;
	com_y=(primary->mass*primary->y+secondary->mass*secondary->y)/totalMass;
	com_z=(primary->mass*primary->z+secondary->mass*secondary->z)/totalMass;
}
std::vector<LagrangePoint> calculateLagrangePoints(Planet* primary, Planet* secondary) {
    std::vector<LagrangePoint> points;
    if (!primary || !secondary) return points;
    
    // параметр μ = m2/(m1+m2)
    double mu = secondary->mass / (primary->mass + secondary->mass);
    
    // Вектор от primary к secondary в плоскости XZ
    double dx = secondary->x - primary->x;
    double dz = secondary->z - primary->z;
    double dy = 0; 
    
    // Расстояние a между телами в плоскости XZ
    double a = std::sqrt(dx*dx + dz*dz);
    
    // центр масс
    double center_x = (primary->mass * primary->x + secondary->mass * secondary->x) / 
                      (primary->mass + secondary->mass);
    double center_y = (primary->mass * primary->y + secondary->mass * secondary->y) / 
                      (primary->mass + secondary->mass);
    double center_z = (primary->mass * primary->z + secondary->mass * secondary->z) / 
                      (primary->mass + secondary->mass);
    
    // единичный вектор от центра масс к secondary в плоскости XZ
    double vector_x = secondary->x - center_x;
    double vector_z = secondary->z - center_z;
    double vector_y = 0;  // В плоскости XZ
    
    double length = std::sqrt(vector_x*vector_x + vector_z*vector_z);
    if (length > 0) {
        vector_x /= length;
        vector_z /= length;
    }
    
    // Перпендикулярный вектор в плоскости XZ (поворот на 90°)
    double perp_x = -vector_z;
    double perp_y = 0;
    double perp_z = vector_x;
    
    // Нормализуем перпендикулярный вектор
    double perp_length = std::sqrt(perp_x*perp_x + perp_z*perp_z);
    if (perp_length > 0) {
        perp_x /= perp_length;
        perp_z /= perp_length;
    }
    
    // Параметр δ = (μ/3)^(1/3)
    double delta = std::cbrt(mu / 3.0);
    
    //вычисляем 5 точек Лагранжа
    for (int i = 1; i <= 5; i++) {
        LagrangePoint point;
        point.type = i;
        point.system = primary->name + "-" + secondary->name;
        point.mu = mu;
        point.R = a;
        
        switch(i) {
            case 1: // L1 между телами
                point.x = center_x + a * (1 - delta) * vector_x;
                point.y = center_y;  // В плоскости XZ, Y остается таким же
                point.z = center_z + a * (1 - delta) * vector_z;
                point.color[0] = 1.0f; point.color[1] = 0.0f; point.color[2] = 0.0f;
                break;
                
            case 2: // L2 за secondary
                point.x = center_x + a * (1 + delta) * vector_x;
                point.y = center_y;
                point.z = center_z + a * (1 + delta) * vector_z;
                point.color[0] = 0.0f; point.color[1] = 1.0f; point.color[2] = 0.0f;
                break;
                
            case 3: // L3 за primary
                point.x = center_x - a * (1.0 + (5.0 * mu / 12.0)) * vector_x;
                point.y = center_y;
                point.z = center_z - a * (1.0 + (5.0 * mu / 12.0)) * vector_z;
                point.color[0] = 0.0f; point.color[1] = 0.0f; point.color[2] = 1.0f;
                break;
                
            case 4: // L4 (равносторонний треугольник)
                // В плоскости XZ: x = a/2 - μa, z = a√3/2
                point.x = center_x + (a/2.0 - mu*a) * vector_x + (a * std::sqrt(3.0)/2.0) * perp_x;
                point.y = center_y;  // Y остается постоянным
                point.z = center_z + (a/2.0 - mu*a) * vector_z + (a * std::sqrt(3.0)/2.0) * perp_z;
                point.color[0] = 1.0f; point.color[1] = 1.0f; point.color[2] = 0.0f;
                break;
                
            case 5: // L5 (симметрично L4)
                point.x = center_x + (a/2.0 - mu*a) * vector_x - (a * std::sqrt(3.0)/2.0) * perp_x;
                point.y = center_y;
                point.z = center_z + (a/2.0 - mu*a) * vector_z - (a * std::sqrt(3.0)/2.0) * perp_z;
                point.color[0] = 1.0f; point.color[1] = 0.0f; point.color[2] = 1.0f;
                break;
        }
        
        // Расстояние до primary
        double dxp = point.x - primary->x;
        double dyp = point.y - primary->y;
        double dzp = point.z - primary->z;
        point.distanceFromPrimary = std::sqrt(dxp*dxp + dyp*dyp + dzp*dzp);
        
        points.push_back(point);
    }
    
    return points;
}

void updateLagrangePoints(std::vector<LagrangePoint>& points, Planet* primary, Planet* secondary) {
    if(!primary || !secondary) return;
    
    double delta_x = secondary->x - primary->x;
    double delta_z = secondary->z - primary->z;
    double current_R = std::sqrt(delta_x*delta_x + delta_z*delta_z);
    
    double center_mass_x, center_mass_y, center_mass_z;
    calculateCenterOfMass(primary, secondary, center_mass_x, center_mass_y, center_mass_z);

    double vector_x = secondary->x - center_mass_x;
    double vector_z = secondary->z - center_mass_z;
    double vector_y = 0;
    
    double vector_length = std::sqrt(vector_x*vector_x + vector_z*vector_z);
    if(vector_length > 0) {
        vector_x /= vector_length;
        vector_z /= vector_length;
    }
    
    double perp_x = -vector_z;
    double perp_y = 0;
    double perp_z = vector_x;
    
    double perp_length = std::sqrt(perp_x*perp_x + perp_z*perp_z);
    if(perp_length > 0) {
        perp_x /= perp_length;
        perp_z /= perp_length;
    }
    
    double sin_60 = std::sqrt(3.0) / 2.0;
    
    for(auto& point : points) {
        if(point.system != primary->name + "-" + secondary->name) continue;
        
        double delta = std::cbrt(point.mu / 3.0);
        
        switch(point.type) {
            case 1: // L1
                point.x = center_mass_x + current_R * (1 - delta) * vector_x;
                point.y = center_mass_y;
                point.z = center_mass_z + current_R * (1 - delta) * vector_z;
                break;
                
            case 2: // L2
                point.x = center_mass_x + current_R * (1 + delta) * vector_x;
                point.y = center_mass_y;
                point.z = center_mass_z + current_R * (1 + delta) * vector_z;
                break;
                
            case 3: // L3
                point.x = center_mass_x - current_R * (1.0 + (5.0 * point.mu / 12.0)) * vector_x;
                point.y = center_mass_y;
                point.z = center_mass_z - current_R * (1.0 + (5.0 * point.mu / 12.0)) * vector_z;
                break;
                
            case 4: // L4
                {
                    double x_offset = current_R * (0.5 - point.mu);
                    double z_offset = current_R * sin_60;
                    point.x = center_mass_x + x_offset * vector_x + z_offset * perp_x;
                    point.y = center_mass_y;
                    point.z = center_mass_z + x_offset * vector_z + z_offset * perp_z;
                }
                break;
                
            case 5: // L5
                {
                    double x_offset = current_R * (0.5 - point.mu);
                    double z_offset = current_R * sin_60;
                    point.x = center_mass_x + x_offset * vector_x - z_offset * perp_x;
                    point.y = center_mass_y;
                    point.z = center_mass_z + x_offset * vector_z - z_offset * perp_z;
                }
                break;
        }
        

        double dxp = point.x - primary->x;
        double dyp = point.y - primary->y;
        double dzp = point.z - primary->z;
        point.distanceFromPrimary = std::sqrt(dxp*dxp + dyp*dyp + dzp*dzp);
    }
}
Planet*createLagrangeSatellite(const LagrangePoint&point,const std::string&systemName) {
	Planet*satellite=new Planet();
	satellite->name="Satellite_L"+std::to_string(point.type);
	satellite->mass=1000.0;
	satellite->x=point.x;
	satellite->y=point.y;
	satellite->z=point.z;
	satellite->vx=0.0;
	satellite->vy=0.0;
	satellite->vz=0.0;
	satellite->ax=0.0;
	satellite->ay=0.0;
	satellite->az=0.0;
	satellite->is_satellite=true;
	satellite->is_lagrange_satellite=true;
	satellite->lagrange_point_type=point.type;
	satellite->parent=systemName;
	satellite->radius=2.0f;
	return satellite;
}

void vel_calc(std::vector<Planet*>*planets) {
	for(Planet*planet:*planets) {
		if(planet->name=="Mercury")VEL_MERCURY=sqrt(2*((-G*MASS_SUN/(2*A_MERCURY))+(G*MASS_SUN/PERIHELION_MERCURY)));
		else if(planet->name=="Venus")VEL_VENUS=sqrt(2*((-G*MASS_SUN/(2*A_VENUS))+(G*MASS_SUN/PERIHELION_VENUS)));
		else if(planet->name=="Earth")VEL_EARTH=sqrt(2*((-G*MASS_SUN/(2*A_EARTH))+(G*MASS_SUN/PERIHELION_EARTH)));
		else if(planet->name=="Mars")VEL_MARS=sqrt(2*((-G*MASS_SUN/(2*A_MARS))+(G*MASS_SUN/PERIHELION_MARS)));
		else if(planet->name=="Jupiter")VEL_JUPITER=sqrt(2*((-G*MASS_SUN/(2*A_JUPITER))+(G*MASS_SUN/PERIHELION_JUPITER)));
		else if(planet->name=="Saturn")VEL_SATURN=sqrt(2*((-G*MASS_SUN/(2*A_SATURN))+(G*MASS_SUN/PERIHELION_SATURN)));
		else if(planet->name=="Uran")VEL_URAN=sqrt(2*((-G*MASS_SUN/(2*A_URAN))+(G*MASS_SUN/PERIHELION_URAN)));
		else if(planet->name=="Neptun")VEL_NEPTUN=sqrt(2*((-G*MASS_SUN/(2*A_NEPTUN))+(G*MASS_SUN/PERIHELION_NEPTUN)));
		else if(planet->name=="luna")VEL_LUNA=sqrt(G*MASS_EARTH/LUNA_ORBIT_RADIUS);
	}
}

void create_sphere(Planet*planet,float radius,float r,float g,float b) {
	planet->radius=radius;
	int s=50;
	for(int i=0;i<=s;++i) {
		float phi=M_PI*i/s;
		for(int j=0;j<=s;++j) {
			float theta=2.0f*M_PI*j/s;
			Vertex v;
			v.x=radius*sin(phi)*cos(theta);
			v.y=radius*cos(phi);
			v.z=radius*sin(phi)*sin(theta);
			v.r=r;
			v.g=g;
			v.b=b;
			(planet->ver).push_back(v);
		}
	}
	for(int i=0;i<s;++i) {
		for(int j=0;j<s;++j) {
			int first=i*(s+1)+j;
			int second=first+1;
			int third=first+(s+1);
			int fourth=third+1;
			(planet->ind).push_back(first);
			(planet->ind).push_back(second);
			(planet->ind).push_back(third);
			(planet->ind).push_back(second);
			(planet->ind).push_back(fourth);
			(planet->ind).push_back(third);
		}
	}
}

void draw_planet(Planet*planet) {
	glBegin(GL_TRIANGLES);
	for(unsigned int i=0;i<(planet->ind).size();++i) {
		Vertex*v=&(planet->ver[planet->ind[i]]);
		glColor3f(v->r,v->g,v->b);
		glVertex3f(v->x,v->y,v->z);
	}
	glEnd();
}

void draw_orbit(float a,float eccentricity,float tilt,float r,float g,float b) {
	glColor3f(r,g,b);
	glBegin(GL_LINE_LOOP);
	for(int i=0;i<500;i++) {
		float angle=2.0f*M_PI*i/500.0f;
		float r_orbit=a*(1-eccentricity*eccentricity)/(1+eccentricity*cos(angle));
		r_orbit*=SCALE;
		float x=r_orbit*cos(angle);
		float z=r_orbit*sin(angle);
		float y=0.0f;
		float y_tilted=z*sin(tilt);
		float z_tilted=z*cos(tilt);
		glVertex3f(x,y_tilted,z_tilted);
	}
	glEnd();
}

void update_physics(Planet* planet, std::vector<Planet*>* all_planets, double dt) {
    if(planet->is_lagrange_satellite) return;
    
    if(planet->name == "Sun") return;
    
    planet->vx += planet->ax * dt * 0.5;
    planet->vy += planet->ay * dt * 0.5;
    planet->vz += planet->az * dt * 0.5;
    
    planet->x += planet->vx * dt;
    planet->y += planet->vy * dt;
    planet->z += planet->vz * dt;
    
    double total_ax = 0.0, total_ay = 0.0, total_az = 0.0;
    for(Planet* other : *all_planets) {
        if(planet == other) continue;
        
        double dx = other->x - planet->x;
        double dy = other->y - planet->y;
        double dz = other->z - planet->z;
        double r2 = dx*dx + dy*dy + dz*dz;
        double r = sqrt(r2);

        if(r < 1e6) continue;
        
        double force = G * other->mass / (r2 * r); // F/m = G*M/r^2 * (dx/r)
        total_ax += force * dx;
        total_ay += force * dy;
        total_az += force * dz;
    }
    
    planet->ax = total_ax;
    planet->ay = total_ay;
    planet->az = total_az;
    
    planet->vx += planet->ax * dt * 0.5;
    planet->vy += planet->ay * dt * 0.5;
    planet->vz += planet->az * dt * 0.5;
}

void calculateLagrangeForces(Planet*satellite,std::vector<Planet*>*all_planets) {
	double total_ax=0.0,total_ay=0.0,total_az=0.0;
	for(Planet*other:*all_planets) {
		if(satellite==other)continue;
		double dx=other->x-satellite->x;
		double dy=other->y-satellite->y;
		double dz=other->z-satellite->z;
		double r2=dx*dx+dy*dy+dz*dz;
		if(r2<1e-12)continue;
		double r=sqrt(r2);
		double a=G*other->mass/r2;
		total_ax+=a*dx/r;
		total_ay+=a*dy/r;
		total_az+=a*dz/r;
	}
	satellite->ax=total_ax;
	satellite->ay=total_ay;
	satellite->az=total_az;
}

void render_planet(Planet*planet) {
	glPushMatrix();
	glTranslatef(float(planet->x*SCALE),float(planet->y*SCALE),float(planet->z*SCALE));
	draw_planet(planet);
	glPopMatrix();
}

void draw_coordinate_system() {
	glLineWidth(2.0f);
	glBegin(GL_LINES);
	glColor3f(1.0f,0.0f,0.0f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(2.0f,0.0f,0.0f);
	glColor3f(0.0f,1.0f,0.0f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(0.0f,2.0f,0.0f);
	glColor3f(0.0f,0.0f,1.0f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(0.0f,0.0f,2.0f);
	glEnd();
	glLineWidth(1.0f);
}

void initialize_planets(std::vector<Planet*>*planets) {
	Planet*sun=new Planet();
	sun->name="Sun";
	sun->mass=MASS_SUN;
	sun->x=0.0;sun->y=0.0;sun->z=0.0;
	sun->vx=0.0;sun->vy=0.0;sun->vz=0.0;
	sun->ax=0.0;sun->ay=0.0;sun->az=0.0;
	sun->orbit_eccentricity=0.0f;
	sun->a=0.0f;
	sun->orbit_tilt=0.0f;
	create_sphere(sun,RADIUS_SUN,COLOR_SUN[0],COLOR_SUN[1],COLOR_SUN[2]);
	planets->push_back(sun);
	Planet*mercury=new Planet();
	mercury->name="Mercury";
	mercury->mass=MASS_MERCURY;
	mercury->x=PERIHELION_MERCURY;mercury->y=0.0;mercury->z=0.0;
	mercury->vx=0.0;mercury->vy=VEL_MERCURY*sin(TILT_MERCURY);mercury->vz=VEL_MERCURY;
	mercury->ax=0.0;mercury->ay=0.0;mercury->az=0.0;
	mercury->orbit_eccentricity=EX_MERCURY;
	mercury->a=A_MERCURY;
	mercury->orbit_tilt=TILT_MERCURY;
	create_sphere(mercury,RADIUS_MERCURY,COLOR_MERCURY[0],COLOR_MERCURY[1],COLOR_MERCURY[2]);
	planets->push_back(mercury);
	Planet*venus=new Planet();
	venus->name="Venus";
	venus->mass=MASS_VENUS;
	venus->x=PERIHELION_VENUS;venus->y=0.0;venus->z=0.0;
	venus->vx=0.0;venus->vy=VEL_VENUS*sin(TILT_VENUS);venus->vz=VEL_VENUS;
	venus->ax=0.0;venus->ay=0.0;venus->az=0.0;
	venus->orbit_eccentricity=EX_VENUS;
	venus->a=A_VENUS;
	venus->orbit_tilt=TILT_VENUS;
	create_sphere(venus,RADIUS_VENUS,COLOR_VENUS[0],COLOR_VENUS[1],COLOR_VENUS[2]);
	planets->push_back(venus);
	Planet*earth=new Planet();
	earth->name="Earth";
	earth->mass=MASS_EARTH;
	earth->x=PERIHELION_EARTH;earth->y=0.0;earth->z=0.0;
	earth->vx=0.0;earth->vy=VEL_EARTH*sin(TILT_EARTH);earth->vz=VEL_EARTH;
	earth->ax=0.0;earth->ay=0.0;earth->az=0.0;
	earth->orbit_eccentricity=EX_EARTH;
	earth->a=A_EARTH;
	earth->orbit_tilt=TILT_EARTH;
	create_sphere(earth,RADIUS_EARTH,COLOR_EARTH[0],COLOR_EARTH[1],COLOR_EARTH[2]);
	planets->push_back(earth);
	Planet*mars=new Planet();
	mars->name="Mars";
	mars->mass=MASS_MARS;
	mars->x=PERIHELION_MARS;mars->y=0.0;mars->z=0.0;
	mars->vx=0.0;mars->vy=VEL_MARS*sin(TILT_MARS);mars->vz=VEL_MARS;
	mars->ax=0.0;mars->ay=0.0;mars->az=0.0;
	mars->orbit_eccentricity=EX_MARS;
	mars->a=A_MARS;
	mars->orbit_tilt=TILT_MARS;
	create_sphere(mars,RADIUS_MARS,COLOR_MARS[0],COLOR_MARS[1],COLOR_MARS[2]);
	planets->push_back(mars);
	Planet*jupiter=new Planet();
	jupiter->name="Jupiter";
	jupiter->mass=MASS_JUPITER;
	jupiter->x=PERIHELION_JUPITER;jupiter->y=0.0;jupiter->z=0.0;
	jupiter->vx=0.0;jupiter->vy=VEL_JUPITER*sin(TILT_JUPITER);jupiter->vz=VEL_JUPITER;
	jupiter->ax=0.0;jupiter->ay=0.0;jupiter->az=0.0;
	jupiter->orbit_eccentricity=EX_JUPITER;
	jupiter->a=A_JUPITER;
	jupiter->orbit_tilt=TILT_JUPITER;
	create_sphere(jupiter,RADIUS_JUPITER,COLOR_JUPITER[0],COLOR_JUPITER[1],COLOR_JUPITER[2]);
	planets->push_back(jupiter);
	Planet*saturn=new Planet();
	saturn->name="Saturn";
	saturn->mass=MASS_SATURN;
	saturn->x=PERIHELION_SATURN;saturn->y=0.0;saturn->z=0.0;
	saturn->vx=0.0;saturn->vy=VEL_SATURN*sin(TILT_SATURN);saturn->vz=VEL_SATURN;
	saturn->ax=0.0;saturn->ay=0.0;saturn->az=0.0;
	saturn->orbit_eccentricity=EX_SATURN;
	saturn->a=A_SATURN;
	saturn->orbit_tilt=TILT_SATURN;
	create_sphere(saturn,RADIUS_SATURN,COLOR_SATURN[0],COLOR_SATURN[1],COLOR_SATURN[2]);
	planets->push_back(saturn);
	Planet*uran=new Planet();
	uran->name="Uran";
	uran->mass=MASS_URAN;
	uran->x=PERIHELION_URAN;uran->y=0.0;uran->z=0.0;
	uran->vx=0.0;uran->vy=VEL_URAN*sin(TILT_URAN);uran->vz=VEL_URAN;
	uran->ax=0.0;uran->ay=0.0;uran->az=0.0;
	uran->orbit_eccentricity=EX_URAN;
	uran->a=A_URAN;
	uran->orbit_tilt=TILT_URAN;
	create_sphere(uran,RADIUS_URAN,COLOR_URAN[0],COLOR_URAN[1],COLOR_URAN[2]);
	planets->push_back(uran);
	Planet*neptun=new Planet();
	neptun->name="Neptun";
	neptun->mass=MASS_NEPTUN;
	neptun->x=PERIHELION_NEPTUN;neptun->y=0.0;neptun->z=0.0;
	neptun->vx=0.0;neptun->vy=VEL_NEPTUN*sin(TILT_NEPTUN);neptun->vz=VEL_NEPTUN;
	neptun->ax=0.0;neptun->ay=0.0;neptun->az=0.0;
	neptun->orbit_eccentricity=EX_NEPTUN;
	neptun->a=A_NEPTUN;
	neptun->orbit_tilt=TILT_NEPTUN;
	create_sphere(neptun,RADIUS_NEPTUN,COLOR_NEPTUN[0],COLOR_NEPTUN[1],COLOR_NEPTUN[2]);
	planets->push_back(neptun);
	Planet*luna=new Planet();
	luna->name="luna";
	luna->mass=MASS_LUNA;
	luna->x=PERIHELION_EARTH+LUNA_ORBIT_RADIUS;luna->y=0.0;luna->z=0.0;
	luna->vx=0;luna->vy=0;luna->vz=VEL_EARTH+VEL_LUNA;
	luna->ax=0.0;luna->ay=0.0;luna->az=0.0;
	luna->is_satellite=true;
	luna->parent="Earth";
	create_sphere(luna,RADIUS_LUNA,COLOR_LUNA[0],COLOR_LUNA[1],COLOR_LUNA[2]);
	planets->push_back(luna);
}

struct Camera {
	float posX,posY,posZ;
	float targetX,targetY,targetZ;
	float upX,upY,upZ;
	float distance;
	float angleX,angleY;
	int followingPlanet=-1;
	float followDistance=0.5f;
	float followAngleX=0.0f;
	float followAngleY=0.3f;
	Camera() {
		posX=30.0f;posY=80.0f;posZ=30.0f;
		targetX=0.0f;targetY=0.0f;targetZ=0.0f;
		upX=0.0f;upY=1.0f;upZ=0.0f;
		distance=50.0f;
		angleX=0.0f;
		angleY=0.3f;
	}
	void update() {
		posX=targetX+distance*cos(angleY)*sin(angleX);
		posY=targetY+distance*sin(angleY);
		posZ=targetZ+distance*cos(angleY)*cos(angleX);
	}
	void updateFollowing(Planet*planet) {
		if(!planet)return;
		targetX=float(planet->x*SCALE);
		targetY=float(planet->y*SCALE);
		targetZ=float(planet->z*SCALE);
		angleX=followAngleX;
		angleY=followAngleY;
		distance=followDistance;
		update();
	}
};

Camera camera;

void key_callback(GLFWwindow*window,int key,int scancode,int action,int mods) {
	if(action==GLFW_PRESS||action==GLFW_REPEAT) {
		float cameraSpeed=0.5f;
		if(camera.followingPlanet>=0) {
			if(key==GLFW_KEY_UP)camera.followAngleY+=0.05f;
			if(key==GLFW_KEY_DOWN)camera.followAngleY-=0.05f;
			if(key==GLFW_KEY_LEFT)camera.followAngleX-=0.05f;
			if(key==GLFW_KEY_RIGHT)camera.followAngleX+=0.05f;
			if(key==GLFW_KEY_W)camera.followDistance-=1.0f;
			if(key==GLFW_KEY_S)camera.followDistance+=1.0f;
		}
		else {
			if(key==GLFW_KEY_UP)camera.targetY+=cameraSpeed;
			if(key==GLFW_KEY_DOWN)camera.targetY-=cameraSpeed;
			if(key==GLFW_KEY_LEFT)camera.targetX-=cameraSpeed;
			if(key==GLFW_KEY_RIGHT)camera.targetX+=cameraSpeed;
			if(key==GLFW_KEY_W)camera.targetZ-=cameraSpeed;
			if(key==GLFW_KEY_S)camera.targetZ+=cameraSpeed;
		}
		if(key==GLFW_KEY_R) {
			camera.followingPlanet=-1;
			camera.posX=30.0f;camera.posY=80.0f;camera.posZ=30.0f;
			camera.targetX=0.0f;camera.targetY=0.0f;camera.targetZ=0.0f;
			camera.distance=30.0f;
			camera.angleX=0.0f;
			camera.angleY=0.3f;
			camera.update();
		}
		if(key==GLFW_KEY_F) {
			camera.followingPlanet=-1;
			camera.posX=0.0f;camera.posY=0.0f;camera.posZ=0.0f;
			camera.targetX=0.0f;camera.targetY=0.0f;camera.targetZ=0.0f;
			camera.angleX=0.0f;camera.angleY=0.0f;
			camera.distance=1.0f;
			camera.update();
		}
		if(key==GLFW_KEY_V) {
			camera.followingPlanet=-1;
			camera.posX=0.0f;camera.posY=100.0f;camera.posZ=0.0f;
			camera.targetX=0.0f;camera.targetY=0.0f;camera.targetZ=0.0f;
			camera.upX=0.0f;camera.upY=0.0f;camera.upZ=-1.0f;
			camera.distance=70.0f;
			camera.angleX=0.0f;
			camera.angleY=-M_PI/2.0f;
			camera.update();
			std::cout<<"Вид сверху над Солнцем"<<std::endl;
		}
		if(key == GLFW_KEY_Z) {
    		SIMULATION_SPEED_MULTIPLIER *= 2.0;
    		std::cout << "Скорость симуляции: " << SIMULATION_SPEED_MULTIPLIER*TIME_SCALE/3600/24 << std::endl;
		}
		if(key == GLFW_KEY_X) {
    		SIMULATION_SPEED_MULTIPLIER /= 2.0;
    		std::cout << "Скорость симуляции: " << SIMULATION_SPEED_MULTIPLIER*TIME_SCALE/3600/24 << std::endl;
		}
		if(key == GLFW_KEY_C) {
    		SIMULATION_SPEED_MULTIPLIER = (SIMULATION_SPEED_MULTIPLIER == 0.0) ? 1.0 : 0.0;
    		std::cout << (SIMULATION_SPEED_MULTIPLIER == 0.0 ? "Пауза" : "Возобновлено") << std::endl;
		}
		if(key==GLFW_KEY_L) {
			showLagrangePoints=!showLagrangePoints;
		}
		if(key==GLFW_KEY_P) {
			lagrangeSatelliteActive=!lagrangeSatelliteActive;
		}
		if(key==GLFW_KEY_F1) {
			selectedLagrangeSystem=0;
			std::cout<<"Выбрана система: Солнце-Земля"<<std::endl;
		}
		else if(key==GLFW_KEY_F2) {
			selectedLagrangeSystem=1;
			std::cout<<"Выбрана система: Земля-Луна"<<std::endl;
		}
		else if(key==GLFW_KEY_F3) {
			selectedLagrangeSystem=2;
			std::cout<<"Выбрана система: Солнце-Юпитер"<<std::endl;
		}
		if(key>=GLFW_KEY_1&&key<=GLFW_KEY_5) {
			selectedLagrangePoint=key-GLFW_KEY_1+1;
			std::cout<<"Выбрана точка L"<<selectedLagrangePoint<<std::endl;
		}
		if(key==GLFW_KEY_F4) {
			camera.followingPlanet=0;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Солнцем"<<std::endl;
		}
		else if(key==GLFW_KEY_F5) {
			camera.followingPlanet=1;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Меркурием"<<std::endl;
		}
		else if(key==GLFW_KEY_F6) {
			camera.followingPlanet=2;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Венерой"<<std::endl;
		}
		else if(key==GLFW_KEY_F7) {
			camera.followingPlanet=3;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Землей"<<std::endl;
		}
		else if(key==GLFW_KEY_F8) {
			camera.followingPlanet=4;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Марсом"<<std::endl;
		}
		else if(key==GLFW_KEY_F9) {
			camera.followingPlanet=5;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Юпитером"<<std::endl;
		}
		else if(key==GLFW_KEY_F10) {
			camera.followingPlanet=6;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Сатурном"<<std::endl;
		}
		else if(key==GLFW_KEY_F11) {
			camera.followingPlanet=7;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Ураном"<<std::endl;
		}
		else if(key==GLFW_KEY_F12) {
			camera.followingPlanet=8;
			camera.followAngleX=0.0f;
			camera.followAngleY=0.3f;
			camera.followDistance=0.5f;
			std::cout<<"Слежение за Нептуном"<<std::endl;
		}
		camera.update();
	}
}

struct Tracker {
    double last_angle = 0;
    double last_time = 0;
    int revolutions = 0;
    double period_sum = 0;
    double last_period = 0;
    double accumulated_angle = 0; 
    
    void update(double x, double z, double time) {
        double current_angle = atan2(z, x);
        time*=TIME_SCALE;
        if(last_time == 0) {
            last_angle = current_angle;
            last_time = time;
            return;
        }
        
        double angle_diff = current_angle - last_angle;

        if(angle_diff < -M_PI) angle_diff += 2 * M_PI;
        if(angle_diff > M_PI) angle_diff -= 2 * M_PI;
        

        
        accumulated_angle += angle_diff;

        if(accumulated_angle >= 2 * M_PI) {
            accumulated_angle -= 2 * M_PI;
            revolutions++;
            double current_period = time - last_time;
            period_sum += current_period;
            last_period = current_period;
            last_time = time;
        }
        else if(accumulated_angle <= -2 * M_PI) {
            accumulated_angle += 2 * M_PI;
            revolutions--;
            double current_period = time - last_time;
            period_sum += current_period;
            last_period = current_period;
            last_time = time;
        }
        
        last_angle = current_angle;
    }
    
    double get_period() {
        if(revolutions == 0) return 0;
        return period_sum / abs(revolutions);
    }
    
    double get_last_period() {
        return last_period;
    }
    
    bool ready() {
        return abs(revolutions) >= 1;
    }
    
    void reset() {
        last_angle = 0;
        last_time = 0;
        revolutions = 0;
        period_sum = 0;
        last_period = 0;
        accumulated_angle = 0;
    }
};
void print_data(Planet* planet, double t) {
    std::cout<<planet->name<<":";
    std::cout<<"vx: "<<planet->vx<<", vy: "<<planet->vy<<", vz: "<<planet->vz<<"\n";
    std::cout<<"x: "<<planet->x<<", y: "<<planet->y<<", z: "<<planet->z<<"\n";
    std::cout<<"ax: "<<planet->ax<<", ay: "<<planet->ay<<", az: "<<planet->az<<"\n";
    std::cout<<"mass: "<<planet->mass<<", a: "<<planet->a<<"\n";
    std::cout<<"period: "<<t<<", t^2/a^3: "<<pow(t,2)/pow(planet->a,3)<<"\n";
    std::cout<<"\n";
}
int main() {
	if(!glfwInit()) {
		std::cout<<"Failed to initialize GLFW"<<std::endl;
		return -1;
	}
	GLFWwindow*window=glfwCreateWindow(1200,800,"Солнечная система",NULL,NULL);
	if(!window) {
		std::cout<<"Failed to create GLFW window"<<std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	if(glewInit()!=GLEW_OK) {
		std::cout<<"Failed to initialize GLEW"<<std::endl;
		return -1;
	}
	glfwSetKeyCallback(window,key_callback);
	camera.update();
	std::vector<Planet*>planets;
	initialize_planets(&planets);
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.0f,0.0f,0.05f,1.0f);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	vel_calc(&planets);
	for(Planet*planet:planets) {
		if(planet->name=="Mercury"){planet->vz=VEL_MERCURY*cos(TILT_MERCURY);planet->vy=VEL_MERCURY*sin(TILT_MERCURY);}
		else if(planet->name=="Venus"){planet->vz=VEL_VENUS*cos(TILT_VENUS);planet->vy=VEL_VENUS*sin(TILT_VENUS);}
		else if(planet->name=="Earth"){planet->vz=VEL_EARTH*cos(TILT_EARTH);planet->vy=VEL_EARTH*sin(TILT_EARTH);}
		else if(planet->name=="Mars"){planet->vz=VEL_MARS*cos(TILT_MARS);planet->vy=VEL_MARS*sin(TILT_MARS);}
		else if(planet->name=="Jupiter"){planet->vz=VEL_JUPITER*cos(TILT_JUPITER);planet->vy=VEL_JUPITER*sin(TILT_JUPITER);}
		else if(planet->name=="Saturn"){planet->vz=VEL_SATURN*cos(TILT_SATURN);planet->vy=VEL_SATURN*sin(TILT_SATURN);}
		else if(planet->name=="Uran"){planet->vz=VEL_URAN*cos(TILT_URAN);planet->vy=VEL_URAN*sin(TILT_URAN);}
		else if(planet->name=="Neptun"){planet->vz=VEL_NEPTUN*cos(TILT_NEPTUN);planet->vy=VEL_NEPTUN*sin(TILT_NEPTUN);}
		else if(planet->name=="luna"){planet->vz=(VEL_EARTH)*cos(TILT_EARTH)+VEL_LUNA;planet->vy=(VEL_EARTH)*sin(TILT_EARTH);}
	}
	double last_time=glfwGetTime();
	Tracker mercury_tracker;
	Tracker venus_tracker;
	Tracker earth_tracker;
	Tracker mars_tracker;
	Tracker jupiter_tracker;
	Tracker saturn_tracker;
	Tracker uran_tracker;
	Tracker neptun_tracker;
	double period_mercury=0,period_venus=0,period_earth=0,period_mars=0,period_jupiter=0,period_saturn=0,period_uran=0,period_neptun=0;
	std::vector<LagrangePoint>allLagrangePoints;
	Planet*lagrangeSatellite=nullptr;
	Planet*sun=nullptr;
	Planet*earth=nullptr;
	Planet*moon=nullptr;
	Planet*jupiter=nullptr;
	for(Planet*p:planets) {
		if(p->name=="Sun")sun=p;
		if(p->name=="Earth")earth=p;
		if(p->name=="luna")moon=p;
		if(p->name=="Jupiter")jupiter=p;
	}
	if(sun&&earth) {
		std::vector<LagrangePoint>sunEarthPoints=calculateLagrangePoints(sun,earth);
		for(auto&point:sunEarthPoints) {
			allLagrangePoints.push_back(point);
		}
	}
	if(earth&&moon) {
		std::vector<LagrangePoint>earthMoonPoints=calculateLagrangePoints(earth,moon);
		for(auto&point:earthMoonPoints) {
			allLagrangePoints.push_back(point);
		}
	}
	if(sun&&jupiter) {
		std::vector<LagrangePoint>sunJupiterPoints=calculateLagrangePoints(sun,jupiter);
		for(auto&point:sunJupiterPoints) {
			allLagrangePoints.push_back(point);
		}
	}
	while(!glfwWindowShouldClose(window)) {
		double current_time=glfwGetTime();
		double delta_time=current_time-last_time;
		last_time=current_time;
		double dt=delta_time*TIME_SCALE*SIMULATION_SPEED_MULTIPLIER;
		for(Planet*planet:planets) {
			update_physics(planet,&planets,dt);
		}
		if(showLagrangePoints) {
			if(sun&&earth) {
				updateLagrangePoints(allLagrangePoints,sun,earth);
			}
			if(earth&&moon) {
				updateLagrangePoints(allLagrangePoints,earth,moon);
			}
			if(sun&&jupiter) {
				updateLagrangePoints(allLagrangePoints,sun,jupiter);
			}
		}
		if(lagrangeSatelliteActive&&!allLagrangePoints.empty()) {
			std::string targetSystem;
			if(selectedLagrangeSystem==0)targetSystem="Sun-Earth";
			else if(selectedLagrangeSystem==1)targetSystem="Earth-luna";
			else if(selectedLagrangeSystem==2)targetSystem="Sun-Jupiter";
			for(const auto&point:allLagrangePoints) {
				if(point.system==targetSystem&&point.type==selectedLagrangePoint) {
					if(!lagrangeSatellite) {
						lagrangeSatellite=createLagrangeSatellite(point,targetSystem);
						create_sphere(lagrangeSatellite,RADIUS_SATELLITE,COLOR_SATELLITE[0],COLOR_SATELLITE[1],COLOR_SATELLITE[2]);
						planets.push_back(lagrangeSatellite);
					}
					else {
						lagrangeSatellite->x=point.x;
						lagrangeSatellite->y=point.y;
						lagrangeSatellite->z=point.z;
					}
					break;
				}
			}
			if(lagrangeSatellite) {
				calculateLagrangeForces(lagrangeSatellite,&planets);
			}
		}
		else if(lagrangeSatellite) {
			auto it=std::find(planets.begin(),planets.end(),lagrangeSatellite);
			if(it!=planets.end()) {
				planets.erase(it);
				delete lagrangeSatellite;
				lagrangeSatellite=nullptr;
			}
		}
		for(Planet*planet:planets) {
			if(planet->name=="Mercury") {
				mercury_tracker.update(planet->x,planet->z,current_time);
				if(mercury_tracker.ready()) {
					period_mercury=mercury_tracker.get_period();
					mercury_tracker.reset();
				}
			}
			if(planet->name=="Venus") {
				venus_tracker.update(planet->x,planet->z,current_time);
				if(venus_tracker.ready()) {
					period_venus=venus_tracker.get_period();
					venus_tracker.reset();
				}
			}
			if(planet->name=="Earth") {
				earth_tracker.update(planet->x,planet->z,current_time);
				if(earth_tracker.ready()) {
					period_earth=earth_tracker.get_period();
					earth_tracker.reset();
				}
			}
			if(planet->name=="Mars") {
				mars_tracker.update(planet->x,planet->z,current_time);
				if(mars_tracker.ready()) {
					period_mars=mars_tracker.get_period();
					mars_tracker.reset();
				}
			}
			if(planet->name=="Jupiter") {
				jupiter_tracker.update(planet->x,planet->z,current_time);
				if(jupiter_tracker.ready()) {
					period_jupiter=jupiter_tracker.get_period();
				}
			}
			if(planet->name=="Saturn") {
				saturn_tracker.update(planet->x,planet->z,current_time);
				if(saturn_tracker.ready()) {
					period_saturn=saturn_tracker.get_period();
				}
			}
			if(planet->name=="Uran") {
				uran_tracker.update(planet->x,planet->z,current_time);
				if(uran_tracker.ready()) {
					period_uran=uran_tracker.get_period();
				}
			}
			if(planet->name=="Neptun") {
				neptun_tracker.update(planet->x,planet->z,current_time);
				if(neptun_tracker.ready()) {
					period_neptun=neptun_tracker.get_period();
				}
			}
		}
		if(camera.followingPlanet>=0&&camera.followingPlanet<planets.size()) {
			camera.updateFollowing(planets[camera.followingPlanet]);
		}
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		int width,height;
		glfwGetFramebufferSize(window,&width,&height);
		glViewport(0,0,width,height);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		float aspect=(float)width/(float)height;
		gluPerspective(60.0f,aspect,CLIP_NEAR,CLIP_FAR);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(camera.posX,camera.posY,camera.posZ,camera.targetX,camera.targetY,camera.targetZ,camera.upX,camera.upY,camera.upZ);
		draw_coordinate_system();
		glDisable(GL_DEPTH_TEST);
		glLineWidth(0.5f);
		for(Planet*planet:planets) {
			if(planet->name=="Mercury")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_MERCURY,0.5f,0.5f,0.5f);
			else if(planet->name=="Venus")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_VENUS,1.0f,0.8f,0.4f);
			else if(planet->name=="Earth")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_EARTH,0.3f,0.3f,0.8f);
			else if(planet->name=="Mars")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_MARS,0.8f,0.3f,0.2f);
			else if(planet->name=="Jupiter")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_JUPITER,1.0f,0.4f,0.5f);
			else if(planet->name=="Saturn")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_SATURN,0.2f,0.7f,1.0f);
			else if(planet->name=="Uran")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_URAN,0.45f,0.85f,0.15f);
			else if(planet->name=="Neptun")draw_orbit(planet->a,planet->orbit_eccentricity,TILT_NEPTUN,0.12f,0.95f,1.0f);
		}
		if(showLagrangePoints) {
			for(const auto&point:allLagrangePoints) {
				glPushMatrix();
				glTranslatef(float(point.x*SCALE),float(point.y*SCALE),float(point.z*SCALE));
				glColor3f(point.color[0],point.color[1],point.color[2]);
				glPointSize(5.0f);
				glBegin(GL_POINTS);
				glVertex3f(0,0,0);
				glEnd();
				glPopMatrix();
			}
		}
		glEnable(GL_DEPTH_TEST);
		glLineWidth(1.0f);
		for(Planet*planet:planets) {
			render_planet(planet);
		}
		static int frame_count=0;
		if(frame_count++%60==0) {
			system("clear");
			std::cout<<"=== СОЛНЕЧНАЯ СИСТЕМА С ТОЧКАМИ ЛАГРАНЖА ==="<<std::endl;
			std::cout<<"Управление:"<<std::endl;
			std::cout<<"  L - показ/скрытие точек Лагранжа"<<std::endl;
			std::cout<<"  P - активация/деактивация спутника в точке Лагранжа"<<std::endl;
			std::cout<<"  F1,F2,F3 - выбор системы (Солнце-Земля, Земля-Луна, Солнце-Юпитер)"<<std::endl;
			std::cout<<"  1-5 - выбор точки Лагранжа (1-5)"<<std::endl;
			std::cout<<"  F4-F12 - слежение за планетой"<<std::endl;
			std::cout<<"  R - сброс камеры"<<std::endl;
			std::cout<<"  F - вид из центра"<<std::endl;
			std::cout<<std::endl;
			if(showLagrangePoints) {
				std::cout<<"Точки Лагранжа: ВКЛ"<<std::endl;
				if(lagrangeSatelliteActive&&lagrangeSatellite) {
					std::cout<<"Спутник активен в точке L"<<selectedLagrangePoint;
					std::cout<<" системы ";
					if(selectedLagrangeSystem==0)std::cout<<"Солнце-Земля";
					else if(selectedLagrangeSystem==1)std::cout<<"Земля-Луна";
					else if(selectedLagrangeSystem==2)std::cout<<"Солнце-Юпитер";
					std::cout<<std::endl;
					std::cout<<"Силы на спутнике: ";
					std::cout<<"ax="<<lagrangeSatellite->ax<<", ";
					std::cout<<"ay="<<lagrangeSatellite->ay<<", ";
					std::cout<<"az="<<lagrangeSatellite->az<<std::endl;
					double totalForce=sqrt(lagrangeSatellite->ax*lagrangeSatellite->ax+lagrangeSatellite->ay*lagrangeSatellite->ay+lagrangeSatellite->az*lagrangeSatellite->az);
					std::cout<<"Суммарная сила: "<<totalForce<<" m/s^2"<<std::endl;
					if(totalForce<1e-10) {
						std::cout<<"Спутник остается в покое (силы скомпенсированы)"<<std::endl;
					}
				}
				std::cout<<std::endl;
			}
			std::cout<<"=== ПЛАНЕТЫ ==="<<std::endl;
			for(Planet* planet:planets) {
                if(planet->name=="Mercury") {print_data(planet, period_mercury);}
                else if(planet->name=="Venus") {print_data(planet,period_venus);}
                else if(planet->name=="Earth") {print_data(planet,period_earth);}
                else if(planet->name=="Mars") {print_data(planet,period_mars);}
                else if(planet->name=="Jupiter") {print_data(planet,period_jupiter);}
                else if(planet->name=="Saturn") {print_data(planet,period_saturn);}
                else if(planet->name=="Uran") {print_data(planet,period_uran);}
                else if(planet->name=="Neptun") {print_data(planet,period_neptun);}
                else if(planet->name=="luna") {print_data(planet,0);}
            }
		}
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	for(Planet*p:planets)delete p;
	if(lagrangeSatellite)delete lagrangeSatellite;
	glfwDestroyWindow(window);
	glfwTerminate();
	std::cout<<std::endl<<"Симуляция завершена"<<std::endl;
	return 0;
}
