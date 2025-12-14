#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <unistd.h>

struct Vertex {
	float x,y,z;
	float r,g,b;
};

struct Planet {
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
	float visual_radius;
	double mass;
	double x,y,z;
	double vx,vy,vz;
	double ax,ay,az;
};

void create_sphere(Planet& planet,float visual_radius,float r,float g,float b) {
	planet.visual_radius=visual_radius;
	int s=30;
	for(int i=0;i<=s;++i) {
		float phi=M_PI*i/s;
		for(int j=0;j<=s;++j) {
			float theta=2.0f*M_PI*j/s;
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
	for(int i=0;i<s;++i) {
		for(int j=0;j<s;++j) {
			int first=i*(s+1)+j;
			int second=first+1;
			int third=first+(s+1);
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

void draw_planet(const Planet& planet) {
	glBegin(GL_TRIANGLES);
	for(size_t i=0;i<planet.indices.size();++i) {
		const Vertex& v=planet.vertices[planet.indices[i]];
		glColor3f(v.r,v.g,v.b);
		glVertex3f(v.x,v.y,v.z);
	}
	glEnd();
}

void draw_coordinate() {
	glBegin(GL_LINES);
	glColor3f(1.0f,0.0f,0.0f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(1.0f,0.0f,0.0f);
	glColor3f(0.0f,1.0f,0.0f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(0.0f,1.0f,0.0f);
	glColor3f(0.0f,0.0f,1.0f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(0.0f,0.0f,1.0f);
	glEnd();
}
const double G=6.67*pow(10,-11);
void update_physics(Planet& planet,const Planet& sun,double dt) {
	double dx=sun.x-planet.x;
	double dy=sun.y-planet.y;
	double dz=sun.z-planet.z;
	double r2=dx*dx+dy*dy+dz*dz;
	if(r2<1e-6) r2=1e-6;
	double r=sqrt(r2);
	double a=G*sun.mass/r2;
	planet.ax=a*dx/r;
	planet.ay=a*dy/r;
	planet.az=a*dz/r;
	planet.vx+=planet.ax*dt;
	planet.vy+=planet.ay*dt;
	planet.vz+=planet.az*dt;
	planet.x+=planet.vx*dt;
	planet.y+=planet.vy*dt;
	planet.z+=planet.vz*dt;
}

void render_planet(const Planet& planet) {
	glPushMatrix();
	glTranslatef(static_cast<float>(planet.x),static_cast<float>(planet.y),static_cast<float>(planet.z));
	draw_planet(planet);
	glPopMatrix();
}

int main() {
	if(!glfwInit()) {
		return -1;
	}
	GLFWwindow* window=glfwCreateWindow(800,600,"Solar System",NULL,NULL);
	if(!window) {
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	if(glewInit()!=GLEW_OK) {
		return -1;
	}
	Planet sun,earth;
	sun.mass=1.98*pow(10,30);
	sun.x=0.0;
	sun.y=0.0;
	sun.z=0.0;
	sun.vx=0.0;
	sun.vy=0.0;
	sun.vz=0.0;
	sun.ax=0.0;
	sun.ay=0.0;
	sun.az=0.0;
	create_sphere(sun,0.6f,1.0f,1.0f,0.0f);
	earth.mass=6*pow(10,24);
	earth.x=1.5;
	earth.y=0.0;
	earth.z=0.0;
	earth.vx=0.0;
	earth.vy=0.0;
	earth.vz=sqrt(G*1.98*pow(10,30)/(149*pow(10,18)));
	earth.ax=0.0;
	earth.ay=0.0;
	earth.az=0.0;
	create_sphere(earth,0.2,0.0f,0.5f,1.0f);
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.0f,0.0f,0.1f,1.0f);
	double last_time=glfwGetTime();
	double fixed_dt=0.016;
	while(!glfwWindowShouldClose(window)) {
		double current_time=glfwGetTime();
		double delta_time=current_time-last_time;
		last_time=current_time;
		update_physics(earth,sun,fixed_dt);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		int width,height;
		glfwGetFramebufferSize(window,&width,&height);
		glViewport(0,0,width,height);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		float aspect=(float)width/(float)height;
		gluPerspective(45.0f,aspect,0.1f,100.0f);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(0.0f,5.0f,8.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f);
		draw_coordinate();
		render_planet(sun);
		render_planet(earth);
		glColor3f(0.3f,0.3f,0.3f);
		glBegin(GL_LINE_LOOP);
		for(int i=0;i<100;i++) {
			float angle=2.0f*M_PI*i/100.0f;
			glVertex3f(2.0f*cos(angle),0.0f,2.0f*sin(angle));
		}
		glEnd();
		std::cout<<"Pos:("<<earth.x<<","<<earth.z<<")|Vel:("<<earth.vx<<","<<earth.vz<<")|Acc:("<<earth.ax<<","<<earth.az<<")\n";
		glfwSwapBuffers(window);
		glfwPollEvents();
		if(glfwGetKey(window,GLFW_KEY_ESCAPE)==GLFW_PRESS) {
			break;
		}
		usleep(160000);
	}
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
