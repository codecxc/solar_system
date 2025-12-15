#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <algorithm>
// const
const double AU=1.496e11;
const double G=6.67430e-11;
const double SCALE=1.0/AU;
const double TIME_SCALE=86400*10;

const double MASS_SUN=1.9891e30;
const double MASS_MERCURY=3.302e23;
const double MASS_VENUS=4.868e24;
const double MASS_EARTH=5.9736e24;
//const double MASS_EARTH=5.972e40;
const double MASS_MARS=6.4185e23;
const double MASS_JUPITER=1.8986e27;
const double MASS_SATURN=5.68460e26;
const double MASS_URAN=8.6832e25;
const double MASS_NEPTUN=1.02430e26;
const double MASS_MIPT_1=440075;


const double MIPT_1_ORBIT_RADIUS=6.371e6+418200;
const double VEL_MIPT_1=sqrt(G*MASS_EARTH/MIPT_1_ORBIT_RADIUS);

/*const double ORBIT_MERCURY=0.387*AU;
const double ORBIT_VENUS=0.723*AU;
const double ORBIT_EARTH=1.0*AU;
const double ORBIT_MARS=1.524*AU;
const double ORBIT_JUPITER=5.203*AU;
const double ORBIT_SATURN=9.537*AU;

const double VEL_MERCURY=0.0;
const double VEL_VENUS=0.0;
const double VEL_EARTH=0.0;
const double VEL_MARS=0.0;
const double VEL_JUPITER=0.0;
const double VEL_SATURN=0.0;
*/

double VEL_MERCURY, VEL_VENUS, VEL_EARTH, VEL_MARS, VEL_JUPITER, VEL_SATURN, VEL_URAN, VEL_NEPTUN;


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


const float RADIUS_SUN=0.15f;
const float RADIUS_MERCURY=0.03f;
const float RADIUS_VENUS=0.04f;
const float RADIUS_EARTH=0.05f;
const float RADIUS_MARS=0.04f;
const float RADIUS_JUPITER=0.12f;
const float RADIUS_SATURN=0.1f;
const float RADIUS_URAN=0.08f;
const float RADIUS_NEPTUN=0.07f;
const float RADIUS_MIPT_1=0.02;

const float COLOR_SUN[3]={1.0f,1.0f,0.0f};
const float COLOR_MERCURY[3]={0.7f,0.7f,0.7f};
const float COLOR_VENUS[3]={1.0f,0.8f,0.4f};
const float COLOR_EARTH[3]={0.0f,0.5f,1.0f};
const float COLOR_MARS[3]={1.0f,0.3f,0.1f};
const float COLOR_JUPITER[3]={0.9f,0.7f,0.5f};
const float COLOR_SATURN[3]={0.95f,0.85f,0.6f};
const float COLOR_NEPTUN[3]={0.1f,0.2f,0.0f};
const float COLOR_URAN[3]={0.15f,0.35f,0.6f};
const float COLOR_MIPT_1[3]={0.5f,0.55f,0.7f};

const double MASS_ASTEROID=1.0e19;
const float RADIUS_ASTEROID=0.05f;
const float COLOR_ASTEROID[3]={0.8f,0.4f,0.1f};

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
	float a; // большая полуось
	float orbit_eccentricity; // эксцентриситет
	float orbit_tilt; // угол наклона орбиты
};

/*
struct Sput {
        std::vector<Vertex>vertices;
        std::vector<unsigned int>indices;
        float visual_radius;
        double mass;
        double x,y,z;
        double vx,vy,vz;
        double ax,ay,az;
        std::string name;
};
*/

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
        }
}


// сфера
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


// отрисовка сферы
void draw_planet(Planet* planet) {
	glBegin(GL_TRIANGLES);
	for(unsigned int i=0;i<(planet->ind).size();++i) {
		Vertex* v=&(planet->ver[planet->ind[i]]);
		glColor3f(v->r,v->g,v->b);
		glVertex3f(v->x,v->y,v->z);
	}
	glEnd();
}


// отрисовка орбиты
/*void draw_orbit(float radius,float r,float g,float b) {
	glColor3f(r,g,b);
	glBegin(GL_LINE_LOOP);
	for(int i=0;i<100;i++) {
		float angle=2.0f*M_PI*i/100.0f;
		glVertex3f(radius*cos(angle),0.0f,radius*sin(angle));
	}
	glEnd();
}*/


void draw_orbit(float a, float eccentricity, float tilt, float r, float g, float b) {
	glColor3f(r, g, b);
	glBegin(GL_LINE_LOOP);
	for(int i=0;i<200;i++) {
		float angle=2.0f*M_PI*i/200.0f;
		float r_orbit=a*(1-eccentricity*eccentricity)/(1+eccentricity*cos(angle)); // r=p/(1+ecos(q))
		
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

// ФИЗИКА
/*void update_physics(Planet* planet,const std::vector<Planet*>& all_planets,double dt) {
	double total_ax=0.0,total_ay=0.0,total_az=0.0;
	
	for(const auto& other:all_planets) {
		if(planet==other) continue;
		
		double dx=other->x-planet->x;
		double dy=other->y-planet->y;
		double dz=other->z-planet->z;
		
		double r2=dx*dx+dy*dy+dz*dz;
		if(r2<1e-12) r2=1e-12;
		
		double r=sqrt(r2);
		double a=G*other->mass/r2;
		
		total_ax+=a*dx/r;
		total_ay+=a*dy/r;
		total_az+=a*dz/r;
		if(planet->name=="MIPT_1") {
			std::cout<<"\n"<<other->name<<": "<<r<<"\n";
		}
	}
	
	double old_ax=planet->ax;
	double old_ay=planet->ay;
	double old_az=planet->az;
	
/*	planet->x+=planet->vx*dt+0.5*old_ax*dt*dt;
	planet->y+=planet->vy*dt+0.5*old_ay*dt*dt;
	planet->z+=planet->vz*dt+0.5*old_az*dt*dt;
	
	planet->ax=total_ax;
	planet->ay=total_ay;
	planet->az=total_az;
	
	planet->vx+=0.5*(old_ax+planet->ax)*dt;
	planet->vy+=0.5*(old_ay+planet->ay)*dt;
	planet->vz+=0.5*(old_az+planet->az)*dt;

	planet->vx+=total_ax*dt;
    	planet->vy+=total_ay*dt;
    	planet->vz+=total_az*dt;
    	planet->x+=planet->vx*dt;
    	planet->y+=planet->vy*dt;
    	planet->z+=planet->vz*dt;

    	planet->ax=total_ax;
    	planet->ay=total_ay;
    	planet->az=total_az;
} */
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





// рендер
void render_planet(Planet* planet) {
	glPushMatrix();
	glTranslatef(float(planet->x*SCALE), float(planet->y*SCALE), float(planet->z*SCALE));
	draw_planet(planet);
	glPopMatrix();
}


// отрисовка координат
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
















// инит планет
void initialize_planets(std::vector<Planet*>* planets) {
//	for(Planet* p:*planets) delete p;
//	planets->clear();
//	for(auto p:sputs) delete p;
//      sputs.clear();
	Planet* sun=new Planet();
	sun->name="Sun";
	sun->mass=MASS_SUN;
	sun->x=0.0; sun->y=0.0; sun->z=0.0;
	sun->vx=0.0; sun->vy=0.0; sun->vz=0.0;
	sun->ax=0.0; sun->ay=0.0; sun->az=0.0;
//	sun->orbit_radius=0.0f;
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
//	mercury->orbit_radius=static_cast<float>(ORBIT_MERCURY*SCALE);
	mercury->orbit_eccentricity=EX_MERCURY;
        mercury->a=A_MERCURY;
        mercury->orbit_tilt=TILT_MERCURY;
	create_sphere(mercury,RADIUS_MERCURY,COLOR_MERCURY[0],COLOR_MERCURY[1],COLOR_MERCURY[2]);
	planets->push_back(mercury);
	
	Planet* venus=new Planet();
	venus->name="Venus";
	venus->mass=MASS_VENUS;
	venus->x=PERIHELION_VENUS; venus->y=0.0; venus->z=0.0;
	venus->vx=0.0; venus->vy=VEL_VENUS*sin(TILT_VENUS); venus->vz=VEL_VENUS;
	venus->ax=0.0; venus->ay=0.0; venus->az=0.0;
//	venus->orbit_radius=static_cast<float>(ORBIT_VENUS*SCALE);
	venus->orbit_eccentricity=EX_VENUS;
        venus->a=A_VENUS;
        venus->orbit_tilt=TILT_VENUS;
	create_sphere(venus,RADIUS_VENUS,COLOR_VENUS[0],COLOR_VENUS[1],COLOR_VENUS[2]);
	planets->push_back(venus);
	
	Planet* earth=new Planet();
	earth->name="Earth";
	earth->mass=MASS_EARTH;
	earth->x=PERIHELION_EARTH; earth->y=0.0; earth->z=0.0;
	earth->vx=0.0; earth->vy=VEL_EARTH*sin(TILT_EARTH); earth->vz=VEL_EARTH;
	earth->ax=0.0; earth->ay=0.0; earth->az=0.0;
//	earth->orbit_radius=static_cast<float>(ORBIT_EARTH*SCALE);
	earth->orbit_eccentricity=EX_EARTH;
        earth->a=A_EARTH;
        earth->orbit_tilt=TILT_EARTH;
	create_sphere(earth,RADIUS_EARTH,COLOR_EARTH[0],COLOR_EARTH[1],COLOR_EARTH[2]);
	planets->push_back(earth);
	
	Planet* mars=new Planet();
	mars->name="Mars";
	mars->mass=MASS_MARS;
	mars->x=PERIHELION_MARS; mars->y=0.0; mars->z=0.0;
	mars->vx=0.0; mars->vy=VEL_MARS*sin(TILT_MARS); mars->vz=VEL_MARS;
	mars->ax=0.0; mars->ay=0.0; mars->az=0.0;
//	mars->orbit_radius=static_cast<float>(ORBIT_MARS*SCALE);
	mars->orbit_eccentricity=EX_MARS;
        mars->a=A_MARS;
        mars->orbit_tilt=TILT_MARS;
	create_sphere(mars,RADIUS_MARS,COLOR_MARS[0],COLOR_MARS[1],COLOR_MARS[2]);
	planets->push_back(mars);


	Planet* jupiter=new Planet();
        jupiter->name="Jupiter";
        jupiter->mass=MASS_JUPITER;
        jupiter->x=PERIHELION_JUPITER; jupiter->y=0.0; jupiter->z=0.0;
        jupiter->vx=0.0; jupiter->vy=VEL_JUPITER*sin(TILT_JUPITER); jupiter->vz=VEL_JUPITER;
        jupiter->ax=0.0; jupiter->ay=0.0; jupiter->az=0.0;
//      mars->orbit_radius=static_cast<float>(ORBIT_MARS*SCALE);
        jupiter->orbit_eccentricity=EX_JUPITER;
        jupiter->a=A_JUPITER;
        jupiter->orbit_tilt=TILT_JUPITER;
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
	create_sphere(neptun,RADIUS_NEPTUN,COLOR_NEPTUN[0],COLOR_NEPTUN[1],COLOR_NEPTUN[2]);
	planets->push_back(neptun);
	/*Planet* jupiter=new Planet();
	jupiter->name="Jupiter";
	jupiter->mass=MASS_JUPITER;
	jupiter->x=ORBIT_JUPITER; jupiter->y=0.0; jupiter->z=0.0;
	jupiter->vx=0.0; jupiter->vy=0.0; jupiter->vz=VEL_JUPITER;
	jupiter->ax=0.0; jupiter->ay=0.0; jupiter->az=0.0;
	jupiter->orbit_radius=static_cast<float>(ORBIT_JUPITER*SCALE);
	create_sphere(*jupiter,RADIUS_JUPITER,COLOR_JUPITER[0],COLOR_JUPITER[1],COLOR_JUPITER[2]);
	planets.push_back(jupiter);
	
	Planet* saturn=new Planet();
	saturn->name="Saturn";
	saturn->mass=MASS_SATURN;
	saturn->x=ORBIT_SATURN; saturn->y=0.0; saturn->z=0.0;
	saturn->vx=0.0; saturn->vy=0.0; saturn->vz=VEL_SATURN;
	saturn->ax=0.0; saturn->ay=0.0; saturn->az=0.0;
	saturn->orbit_radius=static_cast<float>(ORBIT_SATURN*SCALE);
	create_sphere(*saturn,RADIUS_SATURN,COLOR_SATURN[0],COLOR_SATURN[1],COLOR_SATURN[2]);
	planets.push_back(saturn);*/

	Planet* mipt_1=new Planet();
        mipt_1->name="MIPT_1";
        mipt_1->mass=MASS_MIPT_1;
        mipt_1->x=PERIHELION_EARTH+MIPT_1_ORBIT_RADIUS; mipt_1->y=0.0; mipt_1->z=0.0;
        mipt_1->vx=0; mipt_1->vy=0; mipt_1->vz=0;
        mipt_1->ax=0.0; mipt_1->ay=0.0; mipt_1->az=0.0;
//      mars->orbit_radius=static_cast<float>(ORBIT_MARS*SCALE);
        create_sphere(mipt_1,RADIUS_MIPT_1,COLOR_MIPT_1[0],COLOR_MIPT_1[1],COLOR_MIPT_1[2]);
        planets->push_back(mipt_1);
}







// камера


struct Camera {
    	float posX, posY, posZ;
    	float targetX, targetY, targetZ;
    	float upX, upY, upZ;
    	float distance;
    	float angleX, angleY;
        Camera() {
        	posX = 30.0f; posY = 80.0f; posZ = 30.0f;
        	targetX = 0.0f; targetY = 0.0f; targetZ = 0.0f;
        	upX = 0.0f; upY = 1.0f; upZ = 0.0f;
        	distance = 50.0f;
        	angleX = 0.0f;
        	angleY = 0.3f;
    	}
    	void update() {
        	posX = targetX + distance * cos(angleY) * sin(angleX);
        	posY = targetY + distance * sin(angleY);
        	posZ = targetZ + distance * cos(angleY) * cos(angleX);
    	}
};

Camera camera;



void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    	if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        	float cameraSpeed = 0.2f;
        
        	if (key == GLFW_KEY_UP) camera.targetY += cameraSpeed;
        	if (key == GLFW_KEY_DOWN) camera.targetY -= cameraSpeed;
        	if (key == GLFW_KEY_LEFT) camera.targetX -= cameraSpeed;
        	if (key == GLFW_KEY_RIGHT) camera.targetX += cameraSpeed;
        	if (key == GLFW_KEY_W) camera.targetZ -= cameraSpeed;
        	if (key == GLFW_KEY_S) camera.targetZ += cameraSpeed;
        
        	if (key == GLFW_KEY_R) {
            		camera.posX = 30.0f; camera.posY = 80.0f; camera.posZ = 30.0f;
            		camera.targetX = 0.0f; camera.targetY = 0.0f; camera.targetZ = 0.0f;
            		camera.distance = 30.0f;
            		camera.angleX = 0.0f;
            		camera.angleY = 0.3f;
            		camera.update();
        	}
        
        	if (key == GLFW_KEY_F) {
			camera.posX=0.0f; 
            		camera.posY=0.0f; 
            		camera.posZ= 0.0f;
            		camera.targetX = 0.0f;
            		camera.targetY = 0.0f;
            		camera.targetZ = 0.0f;
			camera.angleX = 0.0f;
            		camera.angleY = 0.0f;
            		camera.distance = 1.0f;
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
        // 0 радиан против x
        // π/2 радиан против z
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
		else if(planet->name=="MIPT_1") {planet->vz=(VEL_MIPT_1+VEL_EARTH)*cos(TILT_EARTH);planet->vy=(VEL_MIPT_1+VEL_EARTH)*sin(TILT_EARTH);}
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
		else if(planet->name=="MIPT_1") {print_data(planet,0);}
        }
	while(!glfwWindowShouldClose(window)) {
		double current_time=glfwGetTime();
		double delta_time=current_time-last_time;
		last_time=current_time;
		
		double dt=delta_time*TIME_SCALE;
		for(Planet* planet:planets) {
			update_physics(planet,&planets,dt);
		}


		for(Planet* planet:planets) {
    			if(planet->name=="Mercury") {
        			mercury_tracker.update(planet->x, planet->z, current_time);
        			if(mercury_tracker.ready()) {

            				period_mercury=mercury_tracker.get_period();
        			}
    			}
			if(planet->name=="Venus") {
                                venus_tracker.update(planet->x, planet->z, current_time);

                                if(venus_tracker.ready()) {

                                        period_venus=venus_tracker.get_period();
                                }
                        }

			if(planet->name=="Earth") {
                                earth_tracker.update(planet->x, planet->z, current_time);

                                if(earth_tracker.ready()) {

                                        period_earth=earth_tracker.get_period();
                                }
                        }
			if(planet->name=="Mars") {
                                mars_tracker.update(planet->x, planet->z, current_time);

                                if(mars_tracker.ready()) {

                                        period_mars=mars_tracker.get_period();
                                }
                        }
			if(planet->name=="Jupiter") {
                                jupiter_tracker.update(planet->x, planet->z, current_time);

                                if(jupiter_tracker.ready()) {

                                        period_jupiter=jupiter_tracker.get_period();
                                }
                        }
			if(planet->name=="Saturn") {
                                saturn_tracker.update(planet->x, planet->z, current_time);

                                if(saturn_tracker.ready()) {

                                        period_saturn=saturn_tracker.get_period();
                                }
                        }
			if(planet->name=="Uran") {
                                uran_tracker.update(planet->x, planet->z, current_time);

                                if(uran_tracker.ready()) {

                                        period_uran=uran_tracker.get_period();
                                }
                        }
			if(planet->name=="Neptun") {
                                neptun_tracker.update(planet->x, planet->z, current_time);

                                if(neptun_tracker.ready()) {

                                        period_neptun=neptun_tracker.get_period();
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
		gluPerspective(60.0f,aspect,0.1f,100.0f);
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
//		gluLookAt(30.0f,20.0f,30.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f);
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
		}
		glEnable(GL_DEPTH_TEST);
		glLineWidth(1.0f);
		
		for(Planet* planet:planets) {
			render_planet(planet);
		}
		glfwSwapBuffers(window);
		glfwPollEvents();
//		usleep(20000);
	}
	
	for(Planet* p:planets) delete p;
	glfwDestroyWindow(window);
	glfwTerminate();
	std::cout<<std::endl<<"система энд"<<std::endl;
	return 0;
}
