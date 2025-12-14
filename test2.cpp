#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <unistd.h>
struct Vertex {
    	double x, y, z;
    	float r, g, b;
};

struct Planet {
    	std::vector<Vertex> vertices;
    	std::vector<unsigned int> indices;
    	double radius;
    	double mass;
    	double x, y, z;
    	double vx, vy, vz;
    	double ax, ay, az;
};

void create_sphere(Planet& planet, double radius, float r, float g, float b) {
    	planet.radius=radius;
    	int s=30;
    	for (int i=0;i<=s;++i) {
        	double phi=M_PI*i/s;
        	for (int j=0;j<=s;++j) {
            		double theta=2.0f*M_PI*j/s;
            		Vertex v;
            		v.x=radius*sin(phi)*cos(theta);
            		v.y=radius*cos(phi);
            		v.z=radius*sin(phi)*sin(theta);
            		v.r=r;
            		v.g=g;
            		v.b=b;
            		planet.vertices.push_back(v);
        	}
    	}
    	for (int i=0;i<s;++i) {
        	for (int j=0;j<s;++j) {
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

void update_physics(Planet& planet, const Planet& sun, double dt) {
    const double G = 1.0f;
    double dx = sun.x - planet.x;
    double dy = sun.y - planet.y;
    double dz = sun.z - planet.z;
    
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r2);
    
    if (r < 0.01f) r = 0.01f;
    
    float a = G * sun.mass / r2;
    
    planet.ax = a * dx / r;
    planet.ay = a * dy / r;
    planet.az = a * dz / r;

    planet.vx += planet.ax * dt;
    planet.vy += planet.ay * dt;
    planet.vz += planet.az * dt;    
    planet.x += planet.vx * dt;
    planet.y += planet.vy * dt;
    planet.z += planet.vz * dt;
}

void render_planet(const Planet& planet) {
    glPushMatrix();
    glTranslatef(planet.x, planet.y, planet.z);
    draw_planet(planet);
    glPopMatrix();
}

int main() {
    	if (!glfwInit()) {
        	return -1;
    	}
    	GLFWwindow* window = glfwCreateWindow(800, 600, "Solar System", NULL, NULL);
    	if (!window) {
        	glfwTerminate();
        	return -1;
    	}
    	glfwMakeContextCurrent(window);
    	if (glewInit() != GLEW_OK) {
        	return -1;
    	}
    	Planet sun, earth;
    	sun.mass = 100.0f;
    	sun.x = 0.0f;
    	sun.y = 0.0f;
    	sun.z = 0.0f;
    	sun.vx = 0.0f;
    	sun.vy = 0.0f;
    	sun.vz = 0.0f;
    	create_sphere(sun, 0.5f, 1.0f, 1.0f, 0.0f);

    	earth.mass = 1.0f;
    	earth.x = 1.0f;
    	earth.y = 0.0f;
    	earth.z = 0.0f;
    	earth.vx = 0.0f;
    	earth.vy = 0.01f;
    	earth.vz = 0.0f;
    	create_sphere(earth, 0.2f, 0.0f, 0.5f, 1.0f);
    	glEnable(GL_DEPTH_TEST);
    	glClearColor(0.0f, 0.0f, 0.1f, 1.0f);
	float last_time = glfwGetTime();
    	while (!glfwWindowShouldClose(window)) {
        	float current_time=glfwGetTime();
        	float delta_time=current_time-last_time;
		last_time=current_time;
        	float fixed_dt=0.016f;
	        update_physics(earth,sun,fixed_dt * 2.0f);
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
        gluLookAt(0.0f, 5.0f, 8.0f,   // Позиция камеры
                  0.0f, 0.0f, 0.0f,   // Смотрим на центр
                  0.0f, 1.0f, 0.0f);  // Вектор "вверх"
        
        	draw_coordinate();
        
        	render_planet(sun);
        
        	render_planet(earth);
        
        	glColor3f(0.3f, 0.3f, 0.3f);
        	glBegin(GL_LINE_LOOP);
        	for (int i = 0; i < 100; i++) {
            		float angle = 2.0f * M_PI * i / 100.0f;
            		glVertex3f(3.0f * cos(angle), 0.0f, 3.0f * sin(angle));
        	}
        	glEnd();
        	glfwSwapBuffers(window);
        	glfwPollEvents();
        
        	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            		break;
        	}
		std::cout<<earth.ax<<" "<<earth.ay<<" "<<earth.vx<<" "<<earth.vy<<" "<<earth.x<<" "<<earth.y<<"\n";
		usleep(5000000);
    	}
    	glfwDestroyWindow(window);
    	glfwTerminate();
    	return 0;
}
