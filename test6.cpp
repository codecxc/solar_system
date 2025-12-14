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
const double SCALE=1.0/AU*15.0;
const double TIME_SCALE=86400.0*0.001;

const double MASS_SUN=1.989e30;
const double MASS_MERCURY=3.285e23;
const double MASS_VENUS=4.867e24;
const double MASS_EARTH=5.972e24;
const double MASS_MARS=6.39e23;
const double MASS_JUPITER=1.898e27;
const double MASS_SATURN=5.683e26;

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

double VEL_MERCURY, VEL_VENUS, VEL_EARTH, VEL_MARS;
const float EX_SUN=0.0f;
const float EX_MERCURY=0.206f;
const float EX_VENUS=0.007f;
const float EX_EARTH=0.017f;
const float EX_MARS=0.093f;
const float EX_JUPITER=0.049f;
const float EX_SATURN=0.057f;


const float A_SUN=0.0f;
const float A_MERCURY=0.387*AU;
const float A_VENUS=0.723*AU;
const float A_EARTH=1.000*AU;
const float A_MARS=1.524*AU;
const float A_JUPITER=5.2044*AU;
const float A_SATURN=9.5826*AU;

const float TILT_MERCURY=7.01f*M_PI / 180.0f;
const float TILT_VENUS=3.39f*M_PI / 180.0f;
const float TILT_EARTH=0.0f;
const float TILT_MARS=1.85f*M_PI / 180.0f;
const float TILT_JUPITER=1.31f*M_PI / 180.0f;
const float TILT_SATURN=2.49f*M_PI / 180.0f;

const double PERIHELION_MERCURY=A_MERCURY*(1-EX_MERCURY);
const double PERIHELION_VENUS=A_VENUS*(1-EX_VENUS);
const double PERIHELION_EARTH=A_EARTH*(1-EX_EARTH);
const double PERIHELION_MARS=A_MARS*(1-EX_MARS);


const float RADIUS_SUN=0.4f;
const float RADIUS_MERCURY=0.08f;
const float RADIUS_VENUS=0.12f;
const float RADIUS_EARTH=0.15f;
const float RADIUS_MARS=0.10f;
const float RADIUS_JUPITER=0.35f;
const float RADIUS_SATURN=0.3f;

const float COLOR_SUN[3]={1.0f,1.0f,0.0f};
const float COLOR_MERCURY[3]={0.7f,0.7f,0.7f};
const float COLOR_VENUS[3]={1.0f,0.8f,0.4f};
const float COLOR_EARTH[3]={0.0f,0.5f,1.0f};
const float COLOR_MARS[3]={1.0f,0.3f,0.1f};
const float COLOR_JUPITER[3]={0.9f,0.7f,0.5f};
const float COLOR_SATURN[3]={0.95f,0.85f,0.6f};

const double MASS_ASTEROID=1.0e19;
const float RADIUS_ASTEROID=0.05f;
const float COLOR_ASTEROID[3]={0.8f,0.4f,0.1f};

struct Vertex {
	float x,y,z;
	float r,g,b;
};

struct Planet {
	std::vector<Vertex>vertices;
	std::vector<unsigned int>indices;
	float visual_radius;
	double mass;
	double x,y,z;
	double vx,vy,vz;
	double ax,ay,az;
	std::string name;
	float a; // большая полуось
	float orbit_eccentricity; // эксцентриситет
	float orbit_tilt; // угол наклона орбиты




    // Добавь только эти три строки в конец структуры:
    std::vector<float> traj_x;  // x координаты траектории
    std::vector<float> traj_y;  // y координаты траектории  
    std::vector<float> traj_z;  // z координаты траектории
};
// Добавь эту функцию после всех других функций
void update_trajectory(Planet* planet) {
    // Добавляем текущие координаты в траекторию
    planet->traj_x.push_back(static_cast<float>(planet->x * SCALE));
    planet->traj_y.push_back(static_cast<float>(planet->y * SCALE));
    planet->traj_z.push_back(static_cast<float>(planet->z * SCALE));
    
    // Ограничиваем длину траектории (500 точек)
    if(planet->traj_x.size() > 500) {
        planet->traj_x.erase(planet->traj_x.begin());
        planet->traj_y.erase(planet->traj_y.begin());
        planet->traj_z.erase(planet->traj_z.begin());
    }
}

// Добавь эту функцию после update_trajectory
void draw_trajectory(Planet* planet) {
    // Нужно минимум 2 точки для отрисовки линии
    if(planet->traj_x.size() < 2) return;
    
    // Выбираем цвет в зависимости от планеты
    if(planet->name == "Sun") {
        glColor3f(1.0f, 0.0f, 0.0f);  // красный
    } else if(planet->name == "Earth") {
        glColor3f(0.0f, 0.0f, 1.0f);  // синий
    } else if(planet->name == "Asteroid") {
        glColor3f(0.0f, 1.0f, 0.0f);  // зеленый
    } else {
        return;  // для остальных планет не рисуем
    }
    
    // Рисуем линию по точкам траектории
    glBegin(GL_LINE_STRIP);
    for(size_t i = 0; i < planet->traj_x.size(); i++) {
        glVertex3f(planet->traj_x[i], planet->traj_y[i], planet->traj_z[i]);
    }
    glEnd();
}


// Добавь эту функцию после draw_trajectory
Planet* create_asteroid() {
    Planet* asteroid = new Planet();
    asteroid->name = "Asteroid";
    asteroid->mass = 1.0e40;  // маленькая масса
    asteroid->x = 1.2 * AU;   // немного дальше Земли
    asteroid->y = 0.0;
    asteroid->z = 0.0;
    
    // Скорость чуть больше земной
    double earth_v = 29784.8;
    asteroid->vx = 0.0;
    asteroid->vy = 0.0;
    asteroid->vz = earth_v * 1.05;
    
    asteroid->ax = 0.0; 
    asteroid->ay = 0.0; 
    asteroid->az = 0.0;
    
    // Создаем сферу (цвет коричневый)
    asteroid->visual_radius = 0.08f;
    int sectors = 30;
    int stacks = 30;
    
    for(int i = 0; i <= stacks; ++i) {
        float phi = M_PI * i / stacks;
        for(int j = 0; j <= sectors; ++j) {
            float theta = 2.0f * M_PI * j / sectors;
            Vertex v;
            v.x = asteroid->visual_radius * sin(phi) * cos(theta);
            v.y = asteroid->visual_radius * cos(phi);
            v.z = asteroid->visual_radius * sin(phi) * sin(theta);
            v.r = 0.8f;  // коричневый
            v.g = 0.4f;
            v.b = 0.1f;
            asteroid->vertices.push_back(v);
        }
    }
    
    for(int i = 0; i < stacks; ++i) {
        for(int j = 0; j < sectors; ++j) {
            int first = i * (sectors + 1) + j;
            int second = first + 1;
            int third = first + (sectors + 1);
            int fourth = third + 1;
            
            asteroid->indices.push_back(first);
            asteroid->indices.push_back(second);
            asteroid->indices.push_back(third);
            asteroid->indices.push_back(second);
            asteroid->indices.push_back(fourth);
            asteroid->indices.push_back(third);
        }
    }
    
    return asteroid;
}
void vel_calculate(const std::vector<Planet*>& planets) {
	for(const auto& planet:planets) {
		if(planet->name=="Mercury") VEL_MERCURY=sqrt(2*((-G*MASS_SUN/(2*A_MERCURY))+(G*MASS_SUN/PERIHELION_MERCURY)));
                else if(planet->name=="Venus") VEL_VENUS=sqrt(2*((-G*MASS_SUN/(2*A_VENUS))+(G*MASS_SUN/PERIHELION_VENUS)));
                else if(planet->name=="Earth") VEL_EARTH=sqrt(2*((-G*MASS_SUN/(2*A_EARTH))+(G*MASS_SUN/PERIHELION_EARTH)));
                else if(planet->name=="Mars") VEL_MARS=sqrt(2*((-G*MASS_SUN/(2*A_MARS))+(G*MASS_SUN/PERIHELION_MARS)));
        }
}


// сфера
void create_sphere(Planet& planet,float visual_radius,float r,float g,float b) {
	planet.visual_radius=visual_radius;
	int sectors=30;
	int stacks=30;
	
	for(int i=0;i<=stacks;++i) {
		float phi=M_PI*i/stacks;
		for(int j=0;j<=sectors;++j) {
			float theta=2.0f*M_PI*j/sectors;
			Vertex v;
			v.x=visual_radius*sin(phi)*cos(theta);
			v.y=visual_radius*cos(phi);
			v.z=visual_radius*sin(phi)*sin(theta);
			v.r=r;
			v.g=g;
			v.b=b;
			planet.vertices.push_back(v);
		}
	}
	
	for(int i=0;i<stacks;++i) {
		for(int j=0;j<sectors;++j) {
			int first=i*(sectors+1)+j;
			int second=first+1;
			int third=first+(sectors+1);
			int fourth=third+1;
			
			planet.indices.push_back(first);
			planet.indices.push_back(second);
			planet.indices.push_back(third);
			
			planet.indices.push_back(second);
			planet.indices.push_back(fourth);
			planet.indices.push_back(third);
		}
	}
}


// отрисовка сферы
void draw_planet(const Planet& planet) {
	glBegin(GL_TRIANGLES);
	for(size_t i=0;i<planet.indices.size();++i) {
		const Vertex& v=planet.vertices[planet.indices[i]];
		glColor3f(v.r,v.g,v.b);
		glVertex3f(v.x,v.y,v.z);
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

// ФИЗИКА
void update_physics(Planet* planet,const std::vector<Planet*>& all_planets,double dt) {
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
*/
	planet->vx+=total_ax*dt;
    	planet->vy+=total_ay*dt;
    	planet->vz+=total_az*dt;
    	planet->x+=planet->vx*dt;
    	planet->y+=planet->vy*dt;
    	planet->z+=planet->vz*dt;

    	planet->ax=total_ax;
    	planet->ay=total_ay;
    	planet->az=total_az;


    update_trajectory(planet);
} 







// рендер
void render_planet(const Planet* planet) {
	glPushMatrix();
	glTranslatef(static_cast<float>(planet->x*SCALE),
				 static_cast<float>(planet->y*SCALE),
				 static_cast<float>(planet->z*SCALE));
	draw_planet(*planet);
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
void initialize_planets(std::vector<Planet*>& planets) {
	for(auto p:planets) delete p;
	planets.clear();
	
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
	create_sphere(*sun,RADIUS_SUN,COLOR_SUN[0],COLOR_SUN[1],COLOR_SUN[2]);
	planets.push_back(sun);
	
	Planet* mercury=new Planet();
	mercury->name="Mercury";
	mercury->mass=MASS_MERCURY;
	mercury->x=PERIHELION_MERCURY; mercury->y=0.0; mercury->z=0.0;
	mercury->vx=0.0; mercury->vy=0.0; mercury->vz=VEL_MERCURY;
	mercury->ax=0.0; mercury->ay=0.0; mercury->az=0.0;
//	mercury->orbit_radius=static_cast<float>(ORBIT_MERCURY*SCALE);
	mercury->orbit_eccentricity=EX_MERCURY;
        mercury->a=A_MERCURY;
        mercury->orbit_tilt=TILT_MERCURY;
	create_sphere(*mercury,RADIUS_MERCURY,COLOR_MERCURY[0],COLOR_MERCURY[1],COLOR_MERCURY[2]);
	planets.push_back(mercury);
	
	Planet* venus=new Planet();
	venus->name="Venus";
	venus->mass=MASS_VENUS;
	venus->x=PERIHELION_VENUS; venus->y=0.0; venus->z=0.0;
	venus->vx=0.0; venus->vy=0.0; venus->vz=VEL_VENUS;
	venus->ax=0.0; venus->ay=0.0; venus->az=0.0;
//	venus->orbit_radius=static_cast<float>(ORBIT_VENUS*SCALE);
	venus->orbit_eccentricity=EX_VENUS;
        venus->a=A_VENUS;
        venus->orbit_tilt=TILT_VENUS;
	create_sphere(*venus,RADIUS_VENUS,COLOR_VENUS[0],COLOR_VENUS[1],COLOR_VENUS[2]);
	planets.push_back(venus);
	
	Planet* earth=new Planet();
	earth->name="Earth";
	earth->mass=MASS_EARTH;
	earth->x=PERIHELION_EARTH; earth->y=0.0; earth->z=0.0;
	earth->vx=0.0; earth->vy=0.0; earth->vz=VEL_EARTH;
	earth->ax=0.0; earth->ay=0.0; earth->az=0.0;
//	earth->orbit_radius=static_cast<float>(ORBIT_EARTH*SCALE);
	earth->orbit_eccentricity=EX_EARTH;
        earth->a=A_EARTH;
        earth->orbit_tilt=TILT_EARTH;
	create_sphere(*earth,RADIUS_EARTH,COLOR_EARTH[0],COLOR_EARTH[1],COLOR_EARTH[2]);
	planets.push_back(earth);
	
	Planet* mars=new Planet();
	mars->name="Mars";
	mars->mass=MASS_MARS;
	mars->x=PERIHELION_MARS; mars->y=0.0; mars->z=0.0;
	mars->vx=0.0; mars->vy=0.0; mars->vz=VEL_MARS;
	mars->ax=0.0; mars->ay=0.0; mars->az=0.0;
//	mars->orbit_radius=static_cast<float>(ORBIT_MARS*SCALE);
	mars->orbit_eccentricity=EX_MARS;
        mars->a=A_MARS;
        mars->orbit_tilt=TILT_MARS;
	create_sphere(*mars,RADIUS_MARS,COLOR_MARS[0],COLOR_MARS[1],COLOR_MARS[2]);
	planets.push_back(mars);
	
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
}

int main() {
	if(!glfwInit()) {
		std::cerr<<"Failed to initialize GLFW"<<std::endl;
		return -1;
	}
	
	GLFWwindow* window=glfwCreateWindow(1200,800,"Solar System Simulation",NULL,NULL);
	if(!window) {
		std::cerr<<"Failed to create GLFW window"<<std::endl;
		glfwTerminate();
		return -1;
	}
	
	glfwMakeContextCurrent(window);
	
	if(glewInit()!=GLEW_OK) {
		std::cerr<<"Failed to initialize GLEW"<<std::endl;
		return -1;
	}
	
	std::vector<Planet*>planets;
	initialize_planets(planets);

    Planet* asteroid = create_asteroid();
    planets.push_back(asteroid);


	glEnable(GL_DEPTH_TEST);
	glClearColor(0.0f,0.0f,0.05f,1.0f);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	vel_calculate(planets);
	for(auto planet:planets) {
		if(planet->name=="Mercury") planet->vz=VEL_MERCURY;
		else if(planet->name=="Venus") planet->vz=VEL_VENUS;
		else if(planet->name=="Earth") planet->vz=VEL_EARTH;
		else if(planet->name=="Mars") planet->vz=VEL_MARS;
	}
	double last_time=glfwGetTime();
	
	while(!glfwWindowShouldClose(window)) {
		double current_time=glfwGetTime();
		double delta_time=current_time-last_time;
		last_time=current_time;
		
		double dt=delta_time*TIME_SCALE;
		for(auto& planet:planets) {
			update_physics(planet,planets,dt);
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
		gluLookAt(30.0f,20.0f,30.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f);	
		draw_coordinate_system();







    glDisable(GL_DEPTH_TEST);
    glLineWidth(2.0f);
    for(auto planet : planets) {
        if(planet->name == "Sun" || planet->name == "Earth" || planet->name == "Asteroid") {
            draw_trajectory(planet);
        }
    }
    glLineWidth(1.0f);
    glEnable(GL_DEPTH_TEST);



glDisable(GL_DEPTH_TEST);


		glLineWidth(0.5f);
		
		


		

		for(const auto& planet:planets) {
			if(planet->name=="Mercury") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_MERCURY,0.5f,0.5f,0.5f);
			else if(planet->name=="Venus") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_VENUS,1.0f,0.8f,0.4f);
			else if(planet->name=="Earth") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_EARTH,0.3f,0.3f,0.8f);
			else if(planet->name=="Mars") draw_orbit(planet->a,planet->orbit_eccentricity,TILT_MARS,0.8f,0.3f,0.2f);
			
		}
		glEnable(GL_DEPTH_TEST);
		glLineWidth(1.0f);
		
		for(const auto& planet:planets) {
			render_planet(planet);
		}
		glfwSwapBuffers(window);
		glfwPollEvents();
		
		usleep(16666);
	}
	
	for(auto p:planets) delete p;
	glfwDestroyWindow(window);
	glfwTerminate();
	std::cout<<std::endl<<"система энд"<<std::endl;
	return 0;
}
