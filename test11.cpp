#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <algorithm>

const double AU=1.496e11;
const double G=6.67430e-11;
const double SCALE = 1000/AU;
const double TIME_SCALE=86400*10;
const double RADIUS_SCALE=1;
unsigned int trial_steps=50000;

const double MASS_SUN=1.9891e30;
const double MASS_MERCURY=3.302e23;
const double MASS_VENUS=4.868e24;
const double MASS_EARTH=5.9736e24;
const double MASS_MARS=6.4185e23;
const double MASS_JUPITER=1.8986e27;
const double MASS_SATURN=5.68460e26;
const double MASS_URAN=8.6832e25;
const double MASS_NEPTUN=1.02430e26;
const double MASS_MIPT_1=440075;
const double MASS_LUNA=7.342e22;
const double MASS_IO=8.9319e22;
const double MASS_EUROPA=4.7998e22;
const double MASS_GANYMEDE=1.4819e23;
const double MASS_CALLISTO=1.0759e23;
const double MASS_TITAN=1.3452e23;
const double MASS_RHEA=2.306e21;
const double MASS_TETHYS=6.174e20;
const double MASS_DIONE=1.095e21;
const double MASS_TITANIA=3.527e21;
const double MASS_OBERON=3.014e21;
const double MASS_TRITON=2.14e22;
const double MASS_PHOBOS=1.0659e16;
const double MASS_DEIMOS=1.4762e15;

const double MIPT_1_ORBIT_RADIUS=6.371e6+418200;
const double VEL_MIPT_1=sqrt(G*MASS_EARTH/MIPT_1_ORBIT_RADIUS);
const double LUNA_ORBIT_RADIUS=3.844e8;
const double IO_ORBIT_RADIUS=4.217e8;
const double EUROPA_ORBIT_RADIUS=6.709e8;
const double GANYMEDE_ORBIT_RADIUS=1.0704e9;
const double CALLISTO_ORBIT_RADIUS=1.8827e9;
const double TITAN_ORBIT_RADIUS=1.22187e9;
const double RHEA_ORBIT_RADIUS=5.2704e8;
const double TETHYS_ORBIT_RADIUS=2.9467e8;
const double DIONE_ORBIT_RADIUS=3.7742e8;
const double TITANIA_ORBIT_RADIUS=4.363e8;
const double OBERON_ORBIT_RADIUS=5.835e8;
const double TRITON_ORBIT_RADIUS=3.5476e8;
const double PHOBOS_ORBIT_RADIUS=9.376e6;
const double DEIMOS_ORBIT_RADIUS=2.345e7;

double VEL_MERCURY, VEL_VENUS, VEL_EARTH, VEL_MARS, VEL_JUPITER, VEL_SATURN, VEL_URAN, VEL_NEPTUN;
double VEL_LUNA, VEL_IO, VEL_EUROPA, VEL_GANYMEDE, VEL_CALLISTO, VEL_TITAN, VEL_RHEA, VEL_TETHYS, VEL_DIONE, VEL_TITANIA, VEL_OBERON, VEL_TRITON, VEL_PHOBOS, VEL_DEIMOS;

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

const float TILT_MERCURY=7.01f*M_PI / 180.0f;
const float TILT_VENUS=3.39f*M_PI / 180.0f;
const float TILT_EARTH=0.0f;
const float TILT_MARS=1.85f*M_PI / 180.0f;
const float TILT_JUPITER=1.31f*M_PI / 180.0f;
const float TILT_SATURN=2.49f*M_PI / 180.0f;
const float TILT_URAN=0.77f*M_PI / 180.0f;
const float TILT_NEPTUN=1.77f*M_PI / 180.0f;

const double PERIHELION_MERCURY=A_MERCURY*(1-EX_MERCURY);
const double PERIHELION_VENUS=A_VENUS*(1-EX_VENUS);
const double PERIHELION_EARTH=A_EARTH*(1-EX_EARTH);
const double PERIHELION_MARS=A_MARS*(1-EX_MARS);
const double PERIHELION_JUPITER=A_JUPITER*(1-EX_JUPITER);
const double PERIHELION_SATURN=A_SATURN*(1-EX_SATURN);
const double PERIHELION_URAN=A_URAN*(1-EX_URAN);
const double PERIHELION_NEPTUN=A_NEPTUN*(1-EX_NEPTUN);

/*const float RADIUS_SUN=0.15f;
const float RADIUS_MERCURY=0.0005f;
const float RADIUS_VENUS=0.0005f;
const float RADIUS_EARTH=0.0005f;
const float RADIUS_MARS=0.0005f;
const float RADIUS_JUPITER=0.0005f;
const float RADIUS_SATURN=0.0005f;
const float RADIUS_URAN=0.0005f;
const float RADIUS_NEPTUN=0.0005f;
const float RADIUS_MIPT_1=0.0005;
const float RADIUS_LUNA=0.0002f;
const float RADIUS_IO=0.0002f;
const float RADIUS_EUROPA=0.0002f;
const float RADIUS_GANYMEDE=0.0002f;
const float RADIUS_CALLISTO=0.0002f;
const float RADIUS_TITAN=0.0002f;
const float RADIUS_RHEA=0.0002f;
const float RADIUS_TETHYS=0.0002f;
const float RADIUS_DIONE=0.0002f;
const float RADIUS_TITANIA=0.0002f;
const float RADIUS_OBERON=0.0002f;
const float RADIUS_TRITON=0.0002f;
const float RADIUS_PHOBOS=0.0002f;
const float RADIUS_DEIMOS=0.0002f;*/


const float RADIUS_SUN = SCALE*695780*RADIUS_SCALE;
const float RADIUS_MERCURY = 2439.7f * SCALE*RADIUS_SCALE;
const float RADIUS_VENUS = 6051.8f * SCALE*RADIUS_SCALE;
const float RADIUS_EARTH = 6371.0f * SCALE*RADIUS_SCALE;
const float RADIUS_MARS = 3389.5f * SCALE*RADIUS_SCALE;
const float RADIUS_JUPITER = 69911.0f * SCALE*RADIUS_SCALE;
const float RADIUS_SATURN = 58232.0f * SCALE*RADIUS_SCALE;
const float RADIUS_URAN = 25362.0f * SCALE*RADIUS_SCALE;
const float RADIUS_NEPTUN = 24622.0f * SCALE*RADIUS_SCALE;
const float RADIUS_MIPT_1 = 0.0005f*RADIUS_SCALE;
const float RADIUS_LUNA = 1737.4f *SCALE*RADIUS_SCALE;
const float RADIUS_IO = 1821.6f * SCALE*RADIUS_SCALE;
const float RADIUS_EUROPA = 1560.8f * SCALE*RADIUS_SCALE;
const float RADIUS_GANYMEDE = 2634.1f * SCALE*RADIUS_SCALE;
const float RADIUS_CALLISTO = 2410.3f * SCALE*RADIUS_SCALE;
const float RADIUS_TITAN = 2574.7f * SCALE*RADIUS_SCALE;
const float RADIUS_RHEA = 763.8f * SCALE*RADIUS_SCALE;
const float RADIUS_TETHYS = 531.1f * SCALE*RADIUS_SCALE;
const float RADIUS_DIONE = 561.4f * SCALE*RADIUS_SCALE;
const float RADIUS_TITANIA = 788.4f *SCALE*RADIUS_SCALE;
const float RADIUS_OBERON = 761.4f *SCALE*RADIUS_SCALE;
const float RADIUS_TRITON = 1353.4f *SCALE*RADIUS_SCALE;
const float RADIUS_PHOBOS = 11.08f * SCALE*RADIUS_SCALE;
const float RADIUS_DEIMOS = 6.2f*SCALE*RADIUS_SCALE;



const float COLOR_SUN[3]={1.0f,1.0f,0.0f};
const float COLOR_MERCURY[3]={0.7f,0.7f,0.7f};
const float COLOR_VENUS[3]={1.0f,0.8f,0.4f};
const float COLOR_EARTH[3]={0.0f,0.5f,1.0f};
const float COLOR_MARS[3]={1.0f,0.3f,0.1f};
const float COLOR_JUPITER[3]={0.9f,0.7f,0.5f};
const float COLOR_SATURN[3]={0.95f,0.85f,0.6f};
const float COLOR_URAN[3]={0.15f,0.35f,0.6f};
const float COLOR_NEPTUN[3]={0.1f,0.2f,0.0f};
const float COLOR_MIPT_1[3]={0.5f,0.55f,0.7f};
const float COLOR_LUNA[3]={1.0f,1.0f,1.0f};
const float COLOR_IO[3]={1.0f,1.0f,1.0f};
const float COLOR_EUROPA[3]={1.0f,1.0f,1.0f};
const float COLOR_GANYMEDE[3]={1.0f,1.0f,1.0f};
const float COLOR_CALLISTO[3]={1.0f,1.0f,1.0f};
const float COLOR_TITAN[3]={1.0f,1.0f,1.0f};
const float COLOR_RHEA[3]={1.0f,1.0f,1.0f};
const float COLOR_TETHYS[3]={1.0f,1.0f,1.0f};
const float COLOR_DIONE[3]={1.0f,1.0f,1.0f};
const float COLOR_TITANIA[3]={1.0f,1.0f,1.0f};
const float COLOR_OBERON[3]={1.0f,1.0f,1.0f};
const float COLOR_TRITON[3]={1.0f,1.0f,1.0f};
const float COLOR_PHOBOS[3]={1.0f,1.0f,1.0f};
const float COLOR_DEIMOS[3]={1.0f,1.0f,1.0f};

struct TrailPoint {
    double x, y, z;
    float r, g, b;
    double time;
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
        point.time = glfwGetTime();
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
	float x,y,z;
	float r,g,b;
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
	bool is_sputnik=false;
	std::string parent;
    Trail trail;
};

void vel_calc(std::vector<Planet*>* planets) {
	for(Planet* planet:*planets) {
		if(planet->name=="Mercury") VEL_MERCURY=sqrt(2*((-G*MASS_SUN/(2*A_MERCURY))+(G*MASS_SUN/PERIHELION_MERCURY)));
        else if(planet->name=="Venus") VEL_VENUS=sqrt(2*((-G*MASS_SUN/(2*A_VENUS))+(G*MASS_SUN/PERIHELION_VENUS)));
        else if(planet->name=="Earth") VEL_EARTH=sqrt(2*((-G*MASS_SUN/(2*A_EARTH))+(G*MASS_SUN/PERIHELION_EARTH)));
        else if(planet->name=="Mars") VEL_MARS=sqrt(2*((-G*MASS_SUN/(2*A_MARS))+(G*MASS_SUN/PERIHELION_MARS)));
		else if(planet->name=="Jupiter") VEL_JUPITER=sqrt(2*((-G*MASS_SUN/(2*A_JUPITER))+(G*MASS_SUN/PERIHELION_JUPITER)));
		else if(planet->name=="Saturn") VEL_SATURN=sqrt(2*((-G*MASS_SUN/(2*A_SATURN))+(G*MASS_SUN/PERIHELION_SATURN)));
		else if(planet->name=="Uran") VEL_URAN=sqrt(2*((-G*MASS_SUN/(2*A_URAN))+(G*MASS_SUN/PERIHELION_URAN)));
		else if(planet->name=="Neptun") VEL_NEPTUN=sqrt(2*((-G*MASS_SUN/(2*A_NEPTUN))+(G*MASS_SUN/PERIHELION_NEPTUN)));
		else if(planet->name=="luna") VEL_LUNA=sqrt(G*MASS_EARTH/LUNA_ORBIT_RADIUS);
		else if(planet->name=="Io") VEL_IO=sqrt(G*MASS_JUPITER/IO_ORBIT_RADIUS);
		else if(planet->name=="Europa") VEL_EUROPA=sqrt(G*MASS_JUPITER/EUROPA_ORBIT_RADIUS);
		else if(planet->name=="Ganymede") VEL_GANYMEDE=sqrt(G*MASS_JUPITER/GANYMEDE_ORBIT_RADIUS);
		else if(planet->name=="Callisto") VEL_CALLISTO=sqrt(G*MASS_JUPITER/CALLISTO_ORBIT_RADIUS);
		else if(planet->name=="Titan") VEL_TITAN=sqrt(G*MASS_SATURN/TITAN_ORBIT_RADIUS);
		else if(planet->name=="Rhea") VEL_RHEA=sqrt(G*MASS_SATURN/RHEA_ORBIT_RADIUS);
		else if(planet->name=="Tethys") VEL_TETHYS=sqrt(G*MASS_SATURN/TETHYS_ORBIT_RADIUS);
		else if(planet->name=="Dione") VEL_DIONE=sqrt(G*MASS_SATURN/DIONE_ORBIT_RADIUS);
		else if(planet->name=="Titania") VEL_TITANIA=sqrt(G*MASS_URAN/TITANIA_ORBIT_RADIUS);
		else if(planet->name=="Oberon") VEL_OBERON=sqrt(G*MASS_URAN/OBERON_ORBIT_RADIUS);
		else if(planet->name=="Triton") VEL_TRITON=sqrt(G*MASS_NEPTUN/TRITON_ORBIT_RADIUS);
		else if(planet->name=="Phobos") VEL_PHOBOS=sqrt(G*MASS_MARS/PHOBOS_ORBIT_RADIUS);
		else if(planet->name=="Deimos") VEL_DEIMOS=sqrt(G*MASS_MARS/DEIMOS_ORBIT_RADIUS);
    }
}

void create_sphere(Planet* planet,float radius,float r,float g,float b) {
	planet->radius=radius;
	int s=30;
	
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

void draw_planet(Planet* planet) {
	glBegin(GL_TRIANGLES);
	for(unsigned int i=0;i<(planet->ind).size();++i) {
		Vertex* v=&(planet->ver[planet->ind[i]]);
		glColor3f(v->r,v->g,v->b);
		glVertex3f(v->x,v->y,v->z);
	}
	glEnd();
}

void draw_orbit(float a, float eccentricity, float tilt, float r, float g, float b) {
	glColor3f(r, g, b);
	glBegin(GL_LINE_LOOP);
	for(int i=0;i<200;i++) {
		float angle=2.0f*M_PI*i/200.0f;
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
    planet->vx+=0.5*planet->ax*dt;
    planet->vy+=0.5*planet->ay*dt;
    planet->vz+=0.5*planet->az*dt;
    planet->x+=planet->vx*dt;
    planet->y += planet->vy*dt;
    planet->z+=planet->vz*dt;
    double total_ax=0.0,total_ay = 0.0, total_az = 0.0;
    for(Planet* other:*all_planets) {
        if(planet==other) continue;
        
        double dx=other->x-planet->x;
        double dy=other->y-planet->y;
        double dz=other->z-planet->z;
        
        double r2=dx*dx+dy*dy+dz*dz;
        if(r2<1e-12) continue;
        
        double r=sqrt(r2);
        double a=G*other->mass/r2;
        
        total_ax+=a*dx/r;
        total_ay+=a*dy/r;
        total_az+=a*dz/r;
    }
    planet->vx+=0.5*total_ax*dt;
    planet->vy+=0.5*total_ay*dt;
    planet->vz+=0.5*total_az*dt;
    planet->ax=total_ax;
    planet->ay=total_ay;
    planet->az=total_az;
}

void render_planet(Planet* planet) {
	glPushMatrix();
	glTranslatef(float(planet->x*SCALE), float(planet->y*SCALE), float(planet->z*SCALE));
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

void initialize_planets(std::vector<Planet*>* planets) {
	Planet* sun=new Planet();
	sun->name="Sun";
	sun->mass=MASS_SUN;
	sun->x=0.0; sun->y=0.0; sun->z=0.0;
	sun->vx=0.0; sun->vy=0.0; sun->vz=0.0;
	sun->ax=0.0; sun->ay=0.0; sun->az=0.0;
	sun->orbit_eccentricity=0.0f;
	sun->a=0.0f;
	sun->orbit_tilt=0.0f;
	create_sphere(sun,RADIUS_SUN,COLOR_SUN[0],COLOR_SUN[1],COLOR_SUN[2]);
	planets->push_back(sun);
	
	Planet* mercury=new Planet();
	mercury->name="Mercury";
	mercury->mass=MASS_MERCURY;
	mercury->x=PERIHELION_MERCURY; mercury->y=0.0; mercury->z=0.0;
	mercury->vx=0.0; mercury->vy=VEL_MERCURY*sin(TILT_MERCURY); mercury->vz=VEL_MERCURY;
	mercury->ax=0.0; mercury->ay=0.0; mercury->az=0.0;
	mercury->orbit_eccentricity=EX_MERCURY;
    mercury->a=A_MERCURY;
    mercury->orbit_tilt=TILT_MERCURY;
    mercury->trail = Trail(trial_steps, 2.0f); 
	create_sphere(mercury,RADIUS_MERCURY,COLOR_MERCURY[0],COLOR_MERCURY[1],COLOR_MERCURY[2]);
	planets->push_back(mercury);
	
	Planet* venus=new Planet();
	venus->name="Venus";
	venus->mass=MASS_VENUS;
	venus->x=PERIHELION_VENUS; venus->y=0.0; venus->z=0.0;
	venus->vx=0.0; venus->vy=VEL_VENUS*sin(TILT_VENUS); venus->vz=VEL_VENUS;
	venus->ax=0.0; venus->ay=0.0; venus->az=0.0;
	venus->orbit_eccentricity=EX_VENUS;
    venus->a=A_VENUS;
    venus->orbit_tilt=TILT_VENUS;
    venus->trail = Trail(trial_steps, 2.0f); 
	create_sphere(venus,RADIUS_VENUS,COLOR_VENUS[0],COLOR_VENUS[1],COLOR_VENUS[2]);
	planets->push_back(venus);
	
	Planet* earth=new Planet();
	earth->name="Earth";
	earth->mass=MASS_EARTH;
	earth->x=PERIHELION_EARTH; earth->y=0.0; earth->z=0.0;
	earth->vx=0.0; earth->vy=VEL_EARTH*sin(TILT_EARTH); earth->vz=VEL_EARTH;
	earth->ax=0.0; earth->ay=0.0; earth->az=0.0;
	earth->orbit_eccentricity=EX_EARTH;
    earth->a=A_EARTH;
    earth->orbit_tilt=TILT_EARTH;
    earth->trail = Trail(trial_steps, 2.0f); 
	create_sphere(earth,RADIUS_EARTH,COLOR_EARTH[0],COLOR_EARTH[1],COLOR_EARTH[2]);
	planets->push_back(earth);
	
	Planet* mars=new Planet();
	mars->name="Mars";
	mars->mass=MASS_MARS;
	mars->x=PERIHELION_MARS; mars->y=0.0; mars->z=0.0;
	mars->vx=0.0; mars->vy=VEL_MARS*sin(TILT_MARS); mars->vz=VEL_MARS;
	mars->ax=0.0; mars->ay=0.0; mars->az=0.0;
	mars->orbit_eccentricity=EX_MARS;
    mars->a=A_MARS;
    mars->orbit_tilt=TILT_MARS;
    mars->trail = Trail(trial_steps, 2.0f); 
	create_sphere(mars,RADIUS_MARS,COLOR_MARS[0],COLOR_MARS[1],COLOR_MARS[2]);
	planets->push_back(mars);

	Planet* jupiter=new Planet();
    jupiter->name="Jupiter";
    jupiter->mass=MASS_JUPITER;
    jupiter->x=PERIHELION_JUPITER; jupiter->y=0.0; jupiter->z=0.0;
    jupiter->vx=0.0; jupiter->vy=VEL_JUPITER*sin(TILT_JUPITER); jupiter->vz=VEL_JUPITER;
    jupiter->ax=0.0; jupiter->ay=0.0; jupiter->az=0.0;
    jupiter->orbit_eccentricity=EX_JUPITER;
    jupiter->a=A_JUPITER;
    jupiter->orbit_tilt=TILT_JUPITER;
    jupiter->trail = Trail(trial_steps, 2.0f); 
    create_sphere(jupiter,RADIUS_JUPITER,COLOR_JUPITER[0],COLOR_JUPITER[1],COLOR_JUPITER[2]);
    planets->push_back(jupiter);

	Planet* saturn=new Planet();
	saturn->name="Saturn";
	saturn->mass=MASS_SATURN;
	saturn->x=PERIHELION_SATURN; saturn->y=0.0; saturn->z=0.0;
	saturn->vx=0.0; saturn->vy=VEL_SATURN*sin(TILT_SATURN); saturn->vz=VEL_SATURN;
	saturn->ax=0.0; saturn->ay=0.0; saturn->az=0.0;
	saturn->orbit_eccentricity=EX_SATURN;
	saturn->a=A_SATURN;
	saturn->orbit_tilt=TILT_SATURN;
    saturn->trail = Trail(trial_steps, 2.0f); 
	create_sphere(saturn,RADIUS_SATURN,COLOR_SATURN[0],COLOR_SATURN[1],COLOR_SATURN[2]);
	planets->push_back(saturn);

	Planet* uran=new Planet();
	uran->name="Uran";
	uran->mass=MASS_URAN;
	uran->x=PERIHELION_URAN; uran->y=0.0; uran->z=0.0;
	uran->vx=0.0; uran->vy=VEL_URAN*sin(TILT_URAN); uran->vz=VEL_URAN;
	uran->ax=0.0; uran->ay=0.0; uran->az=0.0;
	uran->orbit_eccentricity=EX_URAN;
	uran->a=A_URAN;
	uran->orbit_tilt=TILT_URAN;
    uran->trail = Trail(trial_steps, 2.0f); 
	create_sphere(uran,RADIUS_URAN,COLOR_URAN[0],COLOR_URAN[1],COLOR_URAN[2]);
	planets->push_back(uran);

	Planet* neptun=new Planet();
	neptun->name="Neptun";
	neptun->mass=MASS_NEPTUN;
	neptun->x=PERIHELION_NEPTUN; neptun->y=0.0; neptun->z=0.0;
	neptun->vx=0.0; neptun->vy=VEL_NEPTUN*sin(TILT_NEPTUN); neptun->vz=VEL_NEPTUN;
	neptun->ax=0.0; neptun->ay=0.0; neptun->az=0.0;
	neptun->orbit_eccentricity=EX_NEPTUN;
	neptun->a=A_NEPTUN;
	neptun->orbit_tilt=TILT_NEPTUN;
    neptun->trail = Trail(trial_steps, 2.0f); 
	create_sphere(neptun,RADIUS_NEPTUN,COLOR_NEPTUN[0],COLOR_NEPTUN[1],COLOR_NEPTUN[2]);
	planets->push_back(neptun);

	Planet* luna=new Planet();
    luna->name="luna";
    luna->mass=MASS_LUNA;
    luna->x=PERIHELION_EARTH+LUNA_ORBIT_RADIUS; luna->y=0.0; luna->z=0.0;
    luna->vx=0; luna->vy=0; luna->vz=VEL_EARTH+VEL_LUNA;
    luna->ax=0.0; luna->ay=0.0; luna->az=0.0;
	luna->is_sputnik=true;
	luna->parent="Earth";
    luna->trail = Trail(trial_steps, 2.0f); 
    create_sphere(luna,RADIUS_LUNA,COLOR_LUNA[0],COLOR_LUNA[1],COLOR_LUNA[2]);
    planets->push_back(luna);

	Planet* io=new Planet();
    io->name="Io";
    io->mass=MASS_IO;
    io->x=PERIHELION_JUPITER+IO_ORBIT_RADIUS; io->y=0.0; io->z=0.0;
    io->vx=0; io->vy=0; io->vz=VEL_JUPITER+VEL_IO;
    io->ax=0.0; io->ay=0.0; io->az=0.0;
	io->is_sputnik=true;
	io->parent="Jupiter";
    io->trail = Trail(trial_steps, 2.0f); 
    create_sphere(io,RADIUS_IO,COLOR_IO[0],COLOR_IO[1],COLOR_IO[2]);
    planets->push_back(io);

	Planet* europa=new Planet();
    europa->name="Europa";
    europa->mass=MASS_EUROPA;
    europa->x=PERIHELION_JUPITER+EUROPA_ORBIT_RADIUS; europa->y=0.0; europa->z=0.0;
    europa->vx=0; europa->vy=0; europa->vz=VEL_JUPITER+VEL_EUROPA;
    europa->ax=0.0; europa->ay=0.0; europa->az=0.0;
	europa->is_sputnik=true;
	europa->parent="Jupiter";
    europa->trail = Trail(trial_steps, 2.0f); 
    create_sphere(europa,RADIUS_EUROPA,COLOR_EUROPA[0],COLOR_EUROPA[1],COLOR_EUROPA[2]);
    planets->push_back(europa);

	Planet* ganymede=new Planet();
    ganymede->name="Ganymede";
    ganymede->mass=MASS_GANYMEDE;
    ganymede->x=PERIHELION_JUPITER+GANYMEDE_ORBIT_RADIUS; ganymede->y=0.0; ganymede->z=0.0;
    ganymede->vx=0; ganymede->vy=0; ganymede->vz=VEL_JUPITER+VEL_GANYMEDE;
    ganymede->ax=0.0; ganymede->ay=0.0; ganymede->az=0.0;
	ganymede->is_sputnik=true;
	ganymede->parent="Jupiter";
    ganymede->trail = Trail(trial_steps, 2.0f); 
    create_sphere(ganymede,RADIUS_GANYMEDE,COLOR_GANYMEDE[0],COLOR_GANYMEDE[1],COLOR_GANYMEDE[2]);
    planets->push_back(ganymede);

	Planet* callisto=new Planet();
    callisto->name="Callisto";
    callisto->mass=MASS_CALLISTO;
    callisto->x=PERIHELION_JUPITER+CALLISTO_ORBIT_RADIUS; callisto->y=0.0; callisto->z=0.0;
    callisto->vx=0; callisto->vy=0; callisto->vz=VEL_JUPITER+VEL_CALLISTO;
    callisto->ax=0.0; callisto->ay=0.0; callisto->az=0.0;
	callisto->is_sputnik=true;
	callisto->parent="Jupiter";
    callisto->trail = Trail(trial_steps, 2.0f); 
    create_sphere(callisto,RADIUS_CALLISTO,COLOR_CALLISTO[0],COLOR_CALLISTO[1],COLOR_CALLISTO[2]);
    planets->push_back(callisto);

	Planet* titan=new Planet();
    titan->name="Titan";
    titan->mass=MASS_TITAN;
    titan->x=PERIHELION_SATURN+TITAN_ORBIT_RADIUS; titan->y=0.0; titan->z=0.0;
    titan->vx=0; titan->vy=0; titan->vz=VEL_SATURN+VEL_TITAN;
    titan->ax=0.0; titan->ay=0.0; titan->az=0.0;
	titan->is_sputnik=true;
	titan->parent="Saturn";
    titan->trail = Trail(trial_steps, 2.0f); 
    create_sphere(titan,RADIUS_TITAN,COLOR_TITAN[0],COLOR_TITAN[1],COLOR_TITAN[2]);
    planets->push_back(titan);

	Planet* rhea=new Planet();
    rhea->name="Rhea";
    rhea->mass=MASS_RHEA;
    rhea->x=PERIHELION_SATURN+RHEA_ORBIT_RADIUS; rhea->y=0.0; rhea->z=0.0;
    rhea->vx=0; rhea->vy=0; rhea->vz=VEL_SATURN+VEL_RHEA;
    rhea->ax=0.0; rhea->ay=0.0; rhea->az=0.0;
	rhea->is_sputnik=true;
	rhea->parent="Saturn";
    rhea->trail = Trail(trial_steps, 2.0f); 
    create_sphere(rhea,RADIUS_RHEA,COLOR_RHEA[0],COLOR_RHEA[1],COLOR_RHEA[2]);
    planets->push_back(rhea);

	Planet* tethys=new Planet();
    tethys->name="Tethys";
    tethys->mass=MASS_TETHYS;
    tethys->x=PERIHELION_SATURN+TETHYS_ORBIT_RADIUS; tethys->y=0.0; tethys->z=0.0;
    tethys->vx=0; tethys->vy=0; tethys->vz=VEL_SATURN+VEL_TETHYS;
    tethys->ax=0.0; tethys->ay=0.0; tethys->az=0.0;
	tethys->is_sputnik=true;
	tethys->parent="Saturn";
    tethys->trail = Trail(trial_steps, 2.0f); 
    create_sphere(tethys,RADIUS_TETHYS,COLOR_TETHYS[0],COLOR_TETHYS[1],COLOR_TETHYS[2]);
    planets->push_back(tethys);

	Planet* dione=new Planet();
    dione->name="Dione";
    dione->mass=MASS_DIONE;
    dione->x=PERIHELION_SATURN+DIONE_ORBIT_RADIUS; dione->y=0.0; dione->z=0.0;
    dione->vx=0; dione->vy=0; dione->vz=VEL_SATURN+VEL_DIONE;
    dione->ax=0.0; dione->ay=0.0; dione->az=0.0;
	dione->is_sputnik=true;
	dione->parent="Saturn";
    dione->trail = Trail(trial_steps, 2.0f); 
    create_sphere(dione,RADIUS_DIONE,COLOR_DIONE[0],COLOR_DIONE[1],COLOR_DIONE[2]);
    planets->push_back(dione);

	Planet* titania=new Planet();
    titania->name="Titania";
    titania->mass=MASS_TITANIA;
    titania->x=PERIHELION_URAN+TITANIA_ORBIT_RADIUS; titania->y=0.0; titania->z=0.0;
    titania->vx=0; titania->vy=0; titania->vz=VEL_URAN+VEL_TITANIA;
    titania->ax=0.0; titania->ay=0.0; titania->az=0.0;
	titania->is_sputnik=true;
	titania->parent="Uran";
    titania->trail = Trail(trial_steps, 2.0f); 
    create_sphere(titania,RADIUS_TITANIA,COLOR_TITANIA[0],COLOR_TITANIA[1],COLOR_TITANIA[2]);
    planets->push_back(titania);

	Planet* oberon=new Planet();
    oberon->name="Oberon";
    oberon->mass=MASS_OBERON;
    oberon->x=PERIHELION_URAN+OBERON_ORBIT_RADIUS; oberon->y=0.0; oberon->z=0.0;
    oberon->vx=0; oberon->vy=0; oberon->vz=VEL_URAN+VEL_OBERON;
    oberon->ax=0.0; oberon->ay=0.0; oberon->az=0.0;
	oberon->is_sputnik=true;
	oberon->parent="Uran";
    oberon->trail = Trail(trial_steps, 2.0f); 
    create_sphere(oberon,RADIUS_OBERON,COLOR_OBERON[0],COLOR_OBERON[1],COLOR_OBERON[2]);
    planets->push_back(oberon);

	Planet* triton=new Planet();
    triton->name="Triton";
    triton->mass=MASS_TRITON;
    triton->x=PERIHELION_NEPTUN+TRITON_ORBIT_RADIUS; triton->y=0.0; triton->z=0.0;
    triton->vx=0; triton->vy=0; triton->vz=VEL_NEPTUN+VEL_TRITON;
    triton->ax=0.0; triton->ay=0.0; triton->az=0.0;
	triton->is_sputnik=true;
	triton->parent="Neptun";
    triton->trail = Trail(trial_steps, 2.0f); 
    create_sphere(triton,RADIUS_TRITON,COLOR_TRITON[0],COLOR_TRITON[1],COLOR_TRITON[2]);
    planets->push_back(triton);

	Planet* phobos=new Planet();
    phobos->name="Phobos";
    phobos->mass=MASS_PHOBOS;
    phobos->x=PERIHELION_MARS+PHOBOS_ORBIT_RADIUS; phobos->y=0.0; phobos->z=0.0;
    phobos->vx=0; phobos->vy=0; phobos->vz=VEL_MARS+VEL_PHOBOS;
    phobos->ax=0.0; phobos->ay=0.0; phobos->az=0.0;
	phobos->is_sputnik=true;
	phobos->parent="Mars";
    phobos->trail = Trail(trial_steps, 2.0f); 
    create_sphere(phobos,RADIUS_PHOBOS,COLOR_PHOBOS[0],COLOR_PHOBOS[1],COLOR_PHOBOS[2]);
    planets->push_back(phobos);

	Planet* deimos=new Planet();
    deimos->name="Deimos";
    deimos->mass=MASS_DEIMOS;
    deimos->x=PERIHELION_MARS+DEIMOS_ORBIT_RADIUS; deimos->y=0.0; deimos->z=0.0;
    deimos->vx=0; deimos->vy=0; deimos->vz=VEL_MARS+VEL_DEIMOS;
    deimos->ax=0.0; deimos->ay=0.0; deimos->az=0.0;
	deimos->is_sputnik=true;
	deimos->parent="Mars";
    deimos->trail = Trail(trial_steps, 2.0f); 
    create_sphere(deimos,RADIUS_DEIMOS,COLOR_DEIMOS[0],COLOR_DEIMOS[1],COLOR_DEIMOS[2]);
    planets->push_back(deimos);
}

std::vector<Planet*> global_planets;

struct Camera {
    float posX, posY, posZ;
    float targetX, targetY, targetZ;
    float upX, upY, upZ;
    float distance;
    float angleX, angleY;
    float nearPlane;
    float farPlane;
    
    Camera() {
        posX = 30.0f; posY = 80.0f; posZ = 30.0f;
        targetX = 0.0f; targetY = 0.0f; targetZ = 0.0f;
        upX = 0.0f; upY = 1.0f; upZ = 0.0f;
        distance = 50.0f;
        angleX = 0.0f;
        angleY = 0.3f;
        nearPlane = 0.1f;
        farPlane = 10000.0f;
    }
    
    void update() {
        posX = targetX + distance * cos(angleY) * sin(angleX);
        posY = targetY + distance * sin(angleY);
        posZ = targetZ + distance * cos(angleY) * cos(angleX);
        
        nearPlane = std::max(0.001f, distance * 0.01f);
        farPlane = std::max(1000.0f, distance * 100.0f);
    }
    
    void focus_on_planet(const std::string& planetName) {
        for(Planet* planet : global_planets) {
            if(planet->name == planetName) {
                targetX = float(planet->x * SCALE);
                targetY = float(planet->y * SCALE);
                targetZ = float(planet->z * SCALE);
                
                if(planet->name == "Sun") {
                    distance = 1.0f;
                } else if(planet->is_sputnik) {
                    distance = 0.01f;
                } else {
                    distance = 0.01f;
                }
                
                angleX = 0.0f;
                angleY = 0.0f;
                update();
                return;
            }
        }
    }
};

Camera camera;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        float cameraSpeed = 0.01f;
        
        if (key == GLFW_KEY_UP) camera.targetY += cameraSpeed;
        if (key == GLFW_KEY_DOWN) camera.targetY -= cameraSpeed;
        if (key == GLFW_KEY_LEFT) camera.targetX -= cameraSpeed;
        if (key == GLFW_KEY_RIGHT) camera.targetX += cameraSpeed;
        if (key == GLFW_KEY_W) camera.targetZ -= cameraSpeed;
        if (key == GLFW_KEY_S) camera.targetZ += cameraSpeed;
        
        if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_KP_ADD) {
            camera.distance *= 0.9f;
            camera.update();
        }
        if (key == GLFW_KEY_MINUS || key == GLFW_KEY_KP_SUBTRACT) {
            camera.distance *= 1.1f;
            camera.update();
        }
        
        if (key == GLFW_KEY_1) camera.focus_on_planet("Sun");
        if (key == GLFW_KEY_2) camera.focus_on_planet("Mercury");
        if (key == GLFW_KEY_3) camera.focus_on_planet("Venus");
        if (key == GLFW_KEY_4) camera.focus_on_planet("Earth");
        if (key == GLFW_KEY_5) camera.focus_on_planet("Mars");
        if (key == GLFW_KEY_6) camera.focus_on_planet("Jupiter");
        if (key == GLFW_KEY_7) camera.focus_on_planet("Saturn");
        if (key == GLFW_KEY_8) camera.focus_on_planet("Uran");
        if (key == GLFW_KEY_9) camera.focus_on_planet("Neptun");
        if (key == GLFW_KEY_0) camera.focus_on_planet("luna");
        
        if (key == GLFW_KEY_R) {
            camera = Camera();
            camera.update();
        }
        
        if (key == GLFW_KEY_F) {
            camera.targetX = 0.0f;
            camera.targetY = 0.0f;
            camera.targetZ = 0.0f;
            camera.distance = 30.0f;
            camera.angleX = 0.0f;
            camera.angleY = 0.3f;
            camera.update();
        }
        
        camera.update();
    }
}

struct Tracker {
    double p_angle=0;
    double p_time=0;
    double period_start=0;
    int num_oborots=0;
    double current_period=0;
    bool has_period=false;
    double sum_angle=0;

    void update(double x, double z, double time) {
        time*=TIME_SCALE;
        double current_angle=atan2(z, x);
        
        if (p_time==0) {
            p_angle=current_angle;
            p_time=time;
            return;
        }
        double angle_diff=current_angle-p_angle;
        sum_angle+=fabs(angle_diff);
        
        if (sum_angle>=(4*M_PI)) {
            sum_angle-=4*M_PI;
            num_oborots++;

            if (num_oborots==1) {
                period_start=time;
                has_period=true;
            }
        }
        p_angle=current_angle;
        p_time=time;
    }
    double get_period() {
        return period_start;
    }
    bool ready() {
        return has_period;
    }
    double get_angle() {
        return sum_angle;
    }
};
struct SatelliteTracker {
    double p_angle = 0;
    double p_time = 0;
    double period_start = 0;
    int num_oborots = 0;
    double current_period = 0;
    bool has_period = false;
    double sum_angle = 0;

    void update(double rel_x, double rel_z, double time) {
        time *= TIME_SCALE;
        double current_angle = atan2(rel_z, rel_x);
        
        if (p_time == 0) {
            p_angle = current_angle;
            p_time = time;
            return;
        }
        
        double angle_diff = current_angle - p_angle;
        if (angle_diff > M_PI) angle_diff -= 2 * M_PI;
        if (angle_diff < -M_PI) angle_diff += 2 * M_PI;
        
        sum_angle += fabs(angle_diff);
        
        if (sum_angle >= (2 * M_PI)) {
            sum_angle -= 2 * M_PI;
            num_oborots++;

            if (num_oborots == 1) {
                period_start = time;
                has_period = true;
            }
        }
        
        p_angle = current_angle;
        p_time = time;
    }
    
    double get_period() {
        return period_start;
    }
    
    bool ready() {
        return has_period;
    }
    
    double get_angle() {
        return sum_angle;
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
    
    GLFWwindow* window=glfwCreateWindow(1200,800,"Солнечная система",NULL,NULL);
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

    glfwSetKeyCallback(window, key_callback);
    
    camera.update();
    
    std::vector<Planet*>planets;
    initialize_planets(&planets);
    global_planets = planets;
    glfwSetWindowUserPointer(window, &planets);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f,0.0f,0.05f,1.0f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    vel_calc(&planets);

    for(Planet* planet:planets) {
        if(planet->name=="Mercury") {planet->vz=VEL_MERCURY*cos(TILT_MERCURY);planet->vy=VEL_MERCURY*sin(TILT_MERCURY);}
        else if(planet->name=="Venus") {planet->vz=VEL_VENUS*cos(TILT_VENUS);planet->vy=VEL_VENUS*sin(TILT_VENUS);}
        else if(planet->name=="Earth") {planet->vz=VEL_EARTH*cos(TILT_EARTH);planet->vy=VEL_EARTH*sin(TILT_EARTH);}
        else if(planet->name=="Mars") {planet->vz=VEL_MARS*cos(TILT_MARS);planet->vy=VEL_MARS*sin(TILT_MARS);}
        else if(planet->name=="Jupiter") {planet->vz=VEL_JUPITER*cos(TILT_JUPITER);planet->vy=VEL_JUPITER*sin(TILT_JUPITER);}
        else if(planet->name=="Saturn") {planet->vz=VEL_SATURN*cos(TILT_SATURN);planet->vy=VEL_SATURN*sin(TILT_SATURN);}
        else if(planet->name=="Uran") {planet->vz=VEL_URAN*cos(TILT_URAN);planet->vy=VEL_URAN*sin(TILT_URAN);}
        else if(planet->name=="Neptun") {planet->vz=VEL_NEPTUN*cos(TILT_NEPTUN);planet->vy=VEL_NEPTUN*sin(TILT_NEPTUN);}
        else if(planet->name=="luna") {planet->vz=(VEL_EARTH)*cos(TILT_EARTH)+VEL_LUNA;planet->vy=(VEL_EARTH)*sin(TILT_EARTH);}
        else if(planet->name=="Io") {planet->vz=(VEL_JUPITER)*cos(TILT_JUPITER)+VEL_IO;planet->vy=(VEL_JUPITER)*sin(TILT_JUPITER);}
        else if(planet->name=="Europa") {planet->vz=(VEL_JUPITER)*cos(TILT_JUPITER)+VEL_EUROPA;planet->vy=(VEL_JUPITER)*sin(TILT_JUPITER);}
        else if(planet->name=="Ganymede") {planet->vz=(VEL_JUPITER)*cos(TILT_JUPITER)+VEL_GANYMEDE;planet->vy=(VEL_JUPITER)*sin(TILT_JUPITER);}
        else if(planet->name=="Callisto") {planet->vz=(VEL_JUPITER)*cos(TILT_JUPITER)+VEL_CALLISTO;planet->vy=(VEL_JUPITER)*sin(TILT_JUPITER);}
        else if(planet->name=="Titan") {planet->vz=(VEL_SATURN)*cos(TILT_SATURN)+VEL_TITAN;planet->vy=(VEL_SATURN)*sin(TILT_SATURN);}
        else if(planet->name=="Rhea") {planet->vz=(VEL_SATURN)*cos(TILT_SATURN)+VEL_RHEA;planet->vy=(VEL_SATURN)*sin(TILT_SATURN);}
        else if(planet->name=="Tethys") {planet->vz=(VEL_SATURN)*cos(TILT_SATURN)+VEL_TETHYS;planet->vy=(VEL_SATURN)*sin(TILT_SATURN);}
        else if(planet->name=="Dione") {planet->vz=(VEL_SATURN)*cos(TILT_SATURN)+VEL_DIONE;planet->vy=(VEL_SATURN)*sin(TILT_SATURN);}
        else if(planet->name=="Titania") {planet->vz=(VEL_URAN)*cos(TILT_URAN)+VEL_TITANIA;planet->vy=(VEL_URAN)*sin(TILT_URAN);}
        else if(planet->name=="Oberon") {planet->vz=(VEL_URAN)*cos(TILT_URAN)+VEL_OBERON;planet->vy=(VEL_URAN)*sin(TILT_URAN);}
        else if(planet->name=="Triton") {planet->vz=(VEL_NEPTUN)*cos(TILT_NEPTUN)+VEL_TRITON;planet->vy=(VEL_NEPTUN)*sin(TILT_NEPTUN);}
        else if(planet->name=="Phobos") {planet->vz=(VEL_MARS)*cos(TILT_MARS)+VEL_PHOBOS;planet->vy=(VEL_MARS)*sin(TILT_MARS);}
        else if(planet->name=="Deimos") {planet->vz=(VEL_MARS)*cos(TILT_MARS)+VEL_DEIMOS;planet->vy=(VEL_MARS)*sin(TILT_MARS);}
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

    SatelliteTracker luna_tracker;
SatelliteTracker io_tracker;
SatelliteTracker europa_tracker;
SatelliteTracker ganymede_tracker;
SatelliteTracker callisto_tracker;
SatelliteTracker titan_tracker;
SatelliteTracker rhea_tracker;
SatelliteTracker tethys_tracker;
SatelliteTracker dione_tracker;
SatelliteTracker titania_tracker;
SatelliteTracker oberon_tracker;
SatelliteTracker triton_tracker;
SatelliteTracker phobos_tracker;
SatelliteTracker deimos_tracker;
    double period_mercury=0,period_venus=0,period_earth=0,period_mars=0,period_jupiter=0,period_saturn=0,period_uran=0,period_neptun=0;
    double period_luna=0, period_io=0, period_europa=0, period_ganymede=0, period_callisto=0;
double period_titan=0, period_rhea=0, period_tethys=0, period_dione=0;
double period_titania=0, period_oberon=0, period_triton=0;
double period_phobos=0, period_deimos=0;
    for(Planet* planet:planets) {
        if(planet->name=="Mercury") {print_data(planet, period_mercury);}
        else if(planet->name=="Venus") {print_data(planet,period_venus);}
        else if(planet->name=="Earth") {print_data(planet,period_earth);}
        else if(planet->name=="Mars") {print_data(planet,period_mars);}
        else if(planet->name=="Jupiter") {print_data(planet,period_jupiter);}
        else if(planet->name=="Saturn") {print_data(planet,period_saturn);}
        else if(planet->name=="Uran") {print_data(planet,period_uran);}
        else if(planet->name=="Neptun") {print_data(planet,period_neptun);}
        else if(planet->name=="Sun") {print_data(planet,0);}
        else if(planet->name=="luna") {print_data(planet,0);}
        else if(planet->name=="Io") {print_data(planet,0);}
        else if(planet->name=="Europa") {print_data(planet,0);}
        else if(planet->name=="Ganymede") {print_data(planet,0);}
        else if(planet->name=="Callisto") {print_data(planet,0);}
        else if(planet->name=="Titan") {print_data(planet,0);}
        else if(planet->name=="Rhea") {print_data(planet,0);}
        else if(planet->name=="Tethys") {print_data(planet,0);}
        else if(planet->name=="Dione") {print_data(planet,0);}
        else if(planet->name=="Titania") {print_data(planet,0);}
        else if(planet->name=="Oberon") {print_data(planet,0);}
        else if(planet->name=="Triton") {print_data(planet,0);}
        else if(planet->name=="Phobos") {print_data(planet,0);}
        else if(planet->name=="Deimos") {print_data(planet,0);}
    }
    while(!glfwWindowShouldClose(window)) {
        double current_time=glfwGetTime();
        double delta_time=current_time-last_time;
        last_time=current_time;
        
        double dt=delta_time*TIME_SCALE;
        for(Planet* planet:planets) {
            update_physics(planet,&planets,dt);
        }
        static int frame_counter = 0;
        frame_counter++;
        for(Planet* planet:planets) {

                if(frame_counter % 3 == 0) {
                    if(planet->name == "luna") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_LUNA[0], COLOR_LUNA[1], COLOR_LUNA[2]);
                    else if(planet->name == "Io") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_IO[0], COLOR_IO[1], COLOR_IO[2]);
                    else if(planet->name == "Europa") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_EUROPA[0], COLOR_EUROPA[1], COLOR_EUROPA[2]);
                    else if(planet->name == "Ganymede") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_GANYMEDE[0], COLOR_GANYMEDE[1], COLOR_GANYMEDE[2]);
                    else if(planet->name == "Callisto") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_CALLISTO[0], COLOR_CALLISTO[1], COLOR_CALLISTO[2]);
                    else if(planet->name == "Titan") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_TITAN[0], COLOR_TITAN[1], COLOR_TITAN[2]);
                    else if(planet->name == "Rhea") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_RHEA[0], COLOR_RHEA[1], COLOR_RHEA[2]);
                    else if(planet->name == "Tethys") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_TETHYS[0], COLOR_TETHYS[1], COLOR_TETHYS[2]);
                    else if(planet->name == "Dione") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_DIONE[0], COLOR_DIONE[1], COLOR_DIONE[2]);
                    else if(planet->name == "Titania") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_TITANIA[0], COLOR_TITANIA[1], COLOR_TITANIA[2]);
                    else if(planet->name == "Oberon") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_OBERON[0], COLOR_OBERON[1], COLOR_OBERON[2]);
                    else if(planet->name == "Triton") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_TRITON[0], COLOR_TRITON[1], COLOR_TRITON[2]);
                    else if(planet->name == "Phobos") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_PHOBOS[0], COLOR_PHOBOS[1], COLOR_PHOBOS[2]);
                    else if(planet->name == "Deimos") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_DEIMOS[0], COLOR_DEIMOS[1], COLOR_DEIMOS[2]);
                    else if(planet->name == "Mercury") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_MERCURY[0], COLOR_MERCURY[1], COLOR_MERCURY[2]);
                    else if(planet->name == "Venus") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_VENUS[0], COLOR_VENUS[1], COLOR_VENUS[2]);
                    else if(planet->name == "Earth") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_EARTH[0], COLOR_EARTH[1], COLOR_EARTH[2]);
                    else if(planet->name == "Mars") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_MARS[0], COLOR_MARS[1], COLOR_MARS[2]);
                    else if(planet->name == "Jupiter") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_JUPITER[0], COLOR_JUPITER[1], COLOR_JUPITER[2]);
                    else if(planet->name == "Saturn") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_SATURN[0], COLOR_SATURN[1], COLOR_SATURN[2]);
                    else if(planet->name == "Uran") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_URAN[0], COLOR_URAN[1], COLOR_URAN[2]);
                    else if(planet->name == "Neptun") planet->trail.addPoint(planet->x, planet->y, planet->z, COLOR_NEPTUN[0], COLOR_NEPTUN[1], COLOR_NEPTUN[2]);
                }
            
        }

for(Planet* planet : planets) {
    if(planet->name == "Mercury") {
        mercury_tracker.update(planet->x, planet->z, current_time);
        if(mercury_tracker.ready()) {
            period_mercury = mercury_tracker.get_period();
        }
    }
    else if(planet->name == "Venus") {
        venus_tracker.update(planet->x, planet->z, current_time);
        if(venus_tracker.ready()) {
            period_venus = venus_tracker.get_period();
        }
    }
    else if(planet->name == "Earth") {
        earth_tracker.update(planet->x, planet->z, current_time);
        if(earth_tracker.ready()) {
            period_earth = earth_tracker.get_period();
        }
    }
    else if(planet->name == "Mars") {
        mars_tracker.update(planet->x, planet->z, current_time);
        if(mars_tracker.ready()) {
            period_mars = mars_tracker.get_period();
        }
    }
    else if(planet->name == "Jupiter") {
        jupiter_tracker.update(planet->x, planet->z, current_time);
        if(jupiter_tracker.ready()) {
            period_jupiter = jupiter_tracker.get_period();
        }
    }
    else if(planet->name == "Saturn") {
        saturn_tracker.update(planet->x, planet->z, current_time);
        if(saturn_tracker.ready()) {
            period_saturn = saturn_tracker.get_period();
        }
    }
    else if(planet->name == "Uran") {
        uran_tracker.update(planet->x, planet->z, current_time);
        if(uran_tracker.ready()) {
            period_uran = uran_tracker.get_period();
        }
    }
    else if(planet->name == "Neptun") {
        neptun_tracker.update(planet->x, planet->z, current_time);
        if(neptun_tracker.ready()) {
            period_neptun = neptun_tracker.get_period();
        }
    }
    
    // Спутники
    if(planet->is_sputnik) {
        Planet* parent_planet = nullptr;
        for(Planet* p : planets) {
            if(p->name == planet->parent) {
                parent_planet = p;
                break;
            }
        }
        
        if(parent_planet) {

            double rel_x = planet->x - parent_planet->x;
            double rel_z = planet->z - parent_planet->z;
            
            if(planet->name == "luna") {
                luna_tracker.update(rel_x, rel_z, current_time);
                if(luna_tracker.ready()) {
                    period_luna = luna_tracker.get_period();
                }
            }
            else if(planet->name == "Io") {
                io_tracker.update(rel_x, rel_z, current_time);
                if(io_tracker.ready()) {
                    period_io = io_tracker.get_period();
                }
            }
            else if(planet->name == "Europa") {
                europa_tracker.update(rel_x, rel_z, current_time);
                if(europa_tracker.ready()) {
                    period_europa = europa_tracker.get_period();
                }
            }
            else if(planet->name == "Ganymede") {
                ganymede_tracker.update(rel_x, rel_z, current_time);
                if(ganymede_tracker.ready()) {
                    period_ganymede = ganymede_tracker.get_period();
                }
            }
            else if(planet->name == "Callisto") {
                callisto_tracker.update(rel_x, rel_z, current_time);
                if(callisto_tracker.ready()) {
                    period_callisto = callisto_tracker.get_period();
                }
            }
            else if(planet->name == "Titan") {
                titan_tracker.update(rel_x, rel_z, current_time);
                if(titan_tracker.ready()) {
                    period_titan = titan_tracker.get_period();
                }
            }
            else if(planet->name == "Rhea") {
                rhea_tracker.update(rel_x, rel_z, current_time);
                if(rhea_tracker.ready()) {
                    period_rhea = rhea_tracker.get_period();
                }
            }
            else if(planet->name == "Tethys") {
                tethys_tracker.update(rel_x, rel_z, current_time);
                if(tethys_tracker.ready()) {
                    period_tethys = tethys_tracker.get_period();
                }
            }
            else if(planet->name == "Dione") {
                dione_tracker.update(rel_x, rel_z, current_time);
                if(dione_tracker.ready()) {
                    period_dione = dione_tracker.get_period();
                }
            }
            else if(planet->name == "Titania") {
                titania_tracker.update(rel_x, rel_z, current_time);
                if(titania_tracker.ready()) {
                    period_titania = titania_tracker.get_period();
                }
            }
            else if(planet->name == "Oberon") {
                oberon_tracker.update(rel_x, rel_z, current_time);
                if(oberon_tracker.ready()) {
                    period_oberon = oberon_tracker.get_period();
                }
            }
            else if(planet->name == "Triton") {
                triton_tracker.update(rel_x, rel_z, current_time);
                if(triton_tracker.ready()) {
                    period_triton = triton_tracker.get_period();
                }
            }
            else if(planet->name == "Phobos") {
                phobos_tracker.update(rel_x, rel_z, current_time);
                if(phobos_tracker.ready()) {
                    period_phobos = phobos_tracker.get_period();
                }
            }
            else if(planet->name == "Deimos") {
                deimos_tracker.update(rel_x, rel_z, current_time);
                if(deimos_tracker.ready()) {
                    period_deimos = deimos_tracker.get_period();
                }
            }
        }
    }
}

        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        
        int width,height;
        glfwGetFramebufferSize(window,&width,&height);
        glViewport(0,0,width,height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect=(float)width/(float)height;
        float nearPlane = camera.nearPlane;
        float farPlane = camera.farPlane;
        float fov = 90.0f;
        gluPerspective(fov, aspect, 0.000001f, 1000000000.0f);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(camera.posX, camera.posY, camera.posZ, camera.targetX, camera.targetY, camera.targetZ, camera.upX, camera.upY, camera.upZ);
        draw_coordinate_system();
        
        glDisable(GL_DEPTH_TEST);
        glLineWidth(0.5f);
        
        for(Planet* planet:planets) {
            if(planet->name=="Mercury") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_MERCURY,0.5f,0.5f,0.5f);
            else if(planet->name=="Venus") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_VENUS,1.0f,0.8f,0.4f);
            else if(planet->name=="Earth") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_EARTH,0.3f,0.3f,0.8f);
            else if(planet->name=="Mars") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_MARS,0.8f,0.3f,0.2f);
            else if(planet->name=="Jupiter") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_JUPITER,1.0f,0.4f,0.5f);
            else if(planet->name=="Saturn") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_SATURN,0.2f,0.7f,1.0f);
            else if(planet->name=="Uran") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_URAN,0.45f,0.85f,0.15f);
            else if(planet->name=="Neptun") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_NEPTUN,0.12f,0.95f,1.0f);
        }
        
        for(Planet* planet:planets) {
            if(planet->name=="Mercury") {print_data(planet, period_mercury);}
            else if(planet->name=="Venus") {print_data(planet,period_venus);}
            else if(planet->name=="Earth") {print_data(planet,period_earth);}
            else if(planet->name=="Mars") {print_data(planet,period_mars);}
            else if(planet->name=="Jupiter") {print_data(planet,period_jupiter);}
            else if(planet->name=="Saturn") {print_data(planet,period_saturn);}
            else if(planet->name=="Uran") {print_data(planet,period_uran);}
            else if(planet->name=="Neptun") {print_data(planet,period_neptun);}
            else if(planet->name=="Sun") {print_data(planet,0);}
            else if(planet->name=="luna") {print_data(planet,period_luna);}
            else if(planet->name=="Io") {print_data(planet,period_io);}
            else if(planet->name=="Europa") {print_data(planet,period_europa);}
            else if(planet->name=="Ganymede") {print_data(planet,period_ganymede);}
            else if(planet->name=="Callisto") {print_data(planet,period_callisto);}
            else if(planet->name=="Titan") {print_data(planet,period_titan);}
            else if(planet->name=="Rhea") {print_data(planet,period_rhea);}
            else if(planet->name=="Tethys") {print_data(planet,period_tethys);}
            else if(planet->name=="Dione") {print_data(planet,period_dione);}
            else if(planet->name=="Titania") {print_data(planet,period_titania);}
            else if(planet->name=="Oberon") {print_data(planet,period_oberon);}
            else if(planet->name=="Triton") {print_data(planet,period_triton);}
            else if(planet->name=="Phobos") {print_data(planet,period_phobos);}
            else if(planet->name=="Deimos") {print_data(planet,period_deimos);}
        }
        glEnable(GL_DEPTH_TEST);
        glLineWidth(1.0f);
        
        glDisable(GL_DEPTH_TEST); 
        for(Planet* planet:planets) {

                planet->trail.draw();
            
        }
        
        glEnable(GL_DEPTH_TEST);
        for(Planet* planet:planets) {
            render_planet(planet);
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    for(Planet* p:planets) delete p;
    glfwDestroyWindow(window);
    glfwTerminate();
    std::cout<<std::endl<<"система энд"<<std::endl;
    return 0;
}
