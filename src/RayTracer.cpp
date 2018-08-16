/*========================================================================
* COSC 363  Computer Graphics (2018)
* Ray tracer 
* See Lab07.pdf for details.
*=========================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "TextureBMP.h"
#include "Ray.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Triangle.h"
#include <GL/glut.h>
#include <math.h>

using namespace std;

const float WIDTH = 20.0;  
const float HEIGHT = 20.0;
const float EDIST = 40.0;
const int NUMDIV = 600;
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;
TextureBMP skytexture;
TextureBMP marbletexture;

vector<SceneObject*> sceneObjects;  //A global list containing pointers to objects in the scene


//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0.18, 0.56, 0.87);
	glm::vec3 light1(200, 200, -3);
	glm::vec3 ambientCol(0.2);   //Ambient color of light
	glm::vec3 eye(0., 0., 0.);

    ray.closestPt(sceneObjects);		//Compute the closest point of intersetion of objects with the ray

    if(ray.xindex == -1) return backgroundCol;      //If there is no intersection return background colour
 
    glm::vec3 materialCol = sceneObjects[ray.xindex]->getColor(); //else return object's colour
    glm::vec3 normalVector = sceneObjects[ray.xindex]->normal(ray.xpt);
    glm::vec3 lightVector1 = light1 - ray.xpt;
    lightVector1 = normalize(lightVector1);
    glm::vec3 reflVector1 = reflect(-lightVector1, normalVector);
    
    Ray shadow1(ray.xpt, lightVector1);
    shadow1.closestPt(sceneObjects);
    
    float ldotn1 = glm::dot(normalVector, lightVector1);
    float lightDist1 = glm::distance(light1, ray.xpt);
    
    glm::vec3 viewVec = -ray.dir;
    
    double rdotv1 = glm::dot(reflVector1, viewVec);
    
    glm::vec3 specCol1;
    
    if (rdotv1 < 0.0) {
		specCol1 = glm::vec3(0., 0., 0.);
	} else {
		specCol1 = glm::vec3(1., 1., 1.);
	}
	
	float specularTerm1 = pow(rdotv1, 10);
	
	glm::vec3 spec1 = specCol1 * specularTerm1;
	
	//floor plane texture
	if (ray.xindex == 1) {
		glm::vec3 color1 = glm::vec3(0.6,0.6,0);
		glm::vec3 color2 = glm::vec3(0.3,0.3,0.6);
		if ((-fmod(ray.xpt.z, 8) < 4.0)) {
			if (ray.xpt.x >= 0.0) {
				if ((fmod(ray.xpt.x, 8) < 4.0)) {
					materialCol = color1;
				} else {
					materialCol = color2;
				}
			} else {
				if ((-fmod(ray.xpt.x, 8) > 4.0)) {
					materialCol = color1;
				} else {
					materialCol = color2;
				}
			}
		} else {
			if (ray.xpt.x >= 0.0) {
				if ((fmod(ray.xpt.x, 8) > 4.0)) {
					materialCol = color1;
				} else {
					materialCol = color2;
				}
			} else {
				if ((-fmod(ray.xpt.x, 8) < 4.0)) {
					materialCol = color1;
				} else {
					materialCol = color2;
				}
			}
		}
	}
	
	//procedural texturing of sphere
	if (ray.xindex == 3) {
		glm::vec3 spherecolor1 = glm::vec3(0.8,0,0);
		glm::vec3 spherecolor2 = glm::vec3(0.8,0.4,0.0);
		glm::vec3 spherecolor3 = glm::vec3(0.8,0.8,0);
		float function = fmod((ray.xpt.y + ray.xpt.x), 3); 
		if (function <= 1.0) {
			materialCol = spherecolor1;
		} else if (function <= 2.0 && function > 1.0){
			materialCol = spherecolor2;
		} else {
			materialCol = spherecolor3;
		}
	}
	
	//the variable holding colour value to be returned
	glm::vec3 colorSum = ambientCol * materialCol; 
	
	//the material's specular component
	glm::vec3 matSpec1 = ldotn1 * materialCol + spec1;


	if (ldotn1 < 0.0) {
		//point is facing away from light source
		colorSum = ambientCol * materialCol;
	} else if ((shadow1.xindex != -1 && shadow1.xdist < lightDist1)) {
		//point is occluded and is in shadow
		
		//soften the shadows of the transparent cube and sphere
		if (shadow1.xindex == 4 || (shadow1.xindex >= 7 && shadow1.xindex <= 12)) {
			colorSum = ambientCol * materialCol + (matSpec1 * 0.5f);
			//make cube shadow colour more similar to the cube colour
			if (shadow1.xindex >= 7 && shadow1.xindex <= 12) {
				colorSum.r += 0.2;
				colorSum.b += 0.1;
			}
		} else {
			colorSum = ambientCol * materialCol;
		}
	} else {
		//point is in the path of the light source
		colorSum = ambientCol * materialCol + matSpec1;
	}
	
	
	//transparent box
	if (ray.xindex >= 7 && ray.xindex <= 12) {
		Ray transparentRay(ray.xpt, ray.dir); 
		glm::vec3 transparentColor = trace(transparentRay, step+1);
		colorSum += 0.3f*colorSum + 0.7f*transparentColor;
	}
	
	//background texture
	if (ray.xindex == 0) {
		float s = (ray.xpt.x - -100)/(100 - -100);
		float t = (ray.xpt.y - -60)/(100 - -60);
		colorSum = skytexture.getColorAt(s, t);
	}
	
	//cylinder texture
	if (ray.xindex == 5) {
		float s = (ray.xpt.x - 10)/(30 - 10);
		float t = (ray.xpt.y - -20)/(15 - -20);
		colorSum = marbletexture.getColorAt(s, t);
	}
	
	//sphere refraction
	if (ray.xindex == 4 and step < MAX_STEPS) {
		float eta = 1.0/1.03;
		glm::vec3 g = glm::refract(ray.dir, normalVector, eta);
		Ray refractionRay(ray.xpt, g);
		refractionRay.closestPt(sceneObjects);
		if (refractionRay.xindex != -1) {
			glm::vec3 m = sceneObjects[refractionRay.xindex]->normal(refractionRay.xpt);
			glm::vec3 h = glm::refract(g, -m, 1.0f/eta);
			Ray exitray(refractionRay.xpt, h);
			colorSum = colorSum*0.3f + (trace(exitray, step+1) + spec1)*0.7f;
		}	
	}
    
	
	//reflections
	if (((ray.xindex >= 1 && ray.xindex <= 13) || (ray.xindex >= 13 && ray.xindex <= 16)) && step < MAX_STEPS) {
		float ratio = 0.5f;
		if (ray.xindex == 4) {
			ratio = 0.2f;
		} else if (ray.xindex >= 13) {
			ratio = 0.3f;
		} else if (ray.xindex == 5) {
			ratio = 0.2f;
		}
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
		Ray reflectedRay(ray.xpt, reflectedDir);
		glm::vec3 reflectedCol = trace(reflectedRay, step+1);
		
		if (ray.xindex >= 7 && ray.xindex <= 12) {
			ratio = 0.2f;
			colorSum = ((1.0f - ratio) * colorSum) + (ratio * reflectedCol);
		} else {
			colorSum = colorSum + (ratio * reflectedCol);
		}
	}
	return colorSum;
}

/**
 * the anti-aliasing function, if the adaptive flag is set, it computes the adaptive supersample of a pixel
 * otherwise it only computes the supersample
*/
glm::vec3 antiAlias(glm::vec3 eye, glm::vec3 dir, float xcellsize, float ycellsize, int step, bool adaptive) {
	float red = 0;
	float green = 0;
	float blue = 0;
	float quarterx = xcellsize / 4;
	float quartery = ycellsize / 4;
	float halfx = xcellsize / 2;
	float halfy = ycellsize / 2;
	glm::vec3 samples[4];
	
	//the colour value of the first subpixel (top left)
	glm::vec3 sample1dir = glm::vec3(dir.x + quarterx, dir.y + quartery, dir.z);
	Ray ray1 = Ray(eye, sample1dir);
	ray1.normalize();
	glm::vec3 trace1color = trace(ray1, 1);
	samples[0] = trace1color;
	
	//the colour value of the second subpixel (top right)
	glm::vec3 sample2dir = glm::vec3(dir.x + 3*quarterx, dir.y + quartery, dir.z);
	Ray ray2 = Ray(eye, sample2dir);
	ray2.normalize();
	glm::vec3 trace2color = trace(ray2, 1);
	samples[1] = trace2color;

	//the colour value of the third subpixel (bottom left)	
	glm::vec3 sample3dir = glm::vec3(dir.x + quarterx, dir.y + 3*quartery, dir.z);
	Ray ray3 = Ray(eye, sample3dir);
	ray3.normalize();
	glm::vec3 trace3color = trace(ray3, 1);
	samples[2] = trace3color;

	//the colour value of the fourth subpixel (bottom right)	
	glm::vec3 sample4dir = glm::vec3(dir.x + 3*quarterx, dir.y + 3*quartery, dir.z);
	Ray ray4 = Ray(eye, sample4dir);
	ray4.normalize();
	glm::vec3 trace4color = trace(ray4, 1);
	samples[3] = trace4color;
	
	
	for (int i = 0; i < 4; i++) {
		if (step < MAX_STEPS && adaptive) {
			//if the adaptive flag is set, then check for differences between the four pixels
			bool differs = false;
			for (int j = 0; j < 4; j++) {
				if (i != j) {
					if (fabs(samples[i].r - samples[j].r) >= 0.3 || fabs(samples[i].g - samples[j].g) >= 0.3 || fabs(samples[i].b - samples[j].b) >= 0.3) {
						differs = true;
					}
				}
			}
			if (differs) {
				//if there is a pixel that differs from the others then recursively call the function on that subpixel
				glm::vec3 newdir = dir;
				if (i == 1) {
					newdir = glm::vec3(dir.x + halfx, dir.y, dir.z);
				} else if (i == 2) {
					newdir = glm::vec3(dir.x, dir.y + halfy, dir.z);
				} else if (i == 3) {
					newdir = glm::vec3(dir.x + halfx, dir.y + halfy, dir.z);
				}
				samples[i] = antiAlias(eye, newdir, xcellsize/4, ycellsize/4, step + 1, true);
			}
		}
		//after sampling is finished (adaptive or otherwise) then compute the average of the four subpixels and return that value
		red += samples[i].r;
		green += samples[i].g;
		blue += samples[i].b;
	}
	
	return glm::vec3(red/4, green/4, blue/4);
}


//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX-XMIN)/NUMDIV;  //cell width
	float cellY = (YMAX-YMIN)/NUMDIV;  //cell height

	glm::vec3 eye(0., 0., 0.);  //The eye position (source of primary rays) is the origin

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  

	for(int i = 0; i < NUMDIV; i++)  	
	{
		xp = XMIN + i*cellX;
		for(int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j*cellY;

		    glm::vec3 dir(xp + 0.5*cellX, yp+ 0.5*cellY, -EDIST);	//direction of the primary ray glm::vec3 dir(xp+0.5*cellX, yp+0.5*cellY, -EDIST);

		    Ray ray = Ray(eye, dir);		//primary ray definition
			ray.normalize();				
		    //glm::vec3 col = trace (ray, 1); 	//trace the ray and get its colour value (uncomment to disable anti-aliasing)
		    glm::vec3 col = antiAlias(eye, dir, cellX, cellY, 0, true);    //supersample and get the colour value of that function (comment out to disable anti-aliasing)

			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				
			glVertex2f(xp+cellX, yp);
			glVertex2f(xp+cellX, yp+cellY);
			glVertex2f(xp, yp+cellY);
        }
    }

    glEnd();
    glFlush();
}


//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glClearColor(0, 0, 0, 1);
    
    skytexture = TextureBMP("sky.bmp");
    marbletexture = TextureBMP("marble.bmp");


	//-- Create a pointer to a sphere object
	Sphere *sphere1 = new Sphere(glm::vec3(0.0, -5.0, -140.0), 15.0, glm::vec3(0.5, 0.5, 0.8));
	Sphere *sphere2 = new Sphere(glm::vec3(16.0, 22.0, -148.0), 7.0, glm::vec3(0, 1, 0));
	Sphere *glassSphere = new Sphere(glm::vec3(-11.0, -12, -120.0), 8.0, glm::vec3(0.3, 0.5, 0.9));
	Plane *floorplane = new Plane(glm::vec3(-60., -20, -40), //A
							glm::vec3(60., -20, -40), //B
							glm::vec3(60., -20, -200), //C
							glm::vec3(-60., -20, -200), //D
							glm::vec3(0.5, 0.5, 0)); //color
	Plane *backgroundplane = new Plane(glm::vec3(-100., -60, -300), //A
							glm::vec3(100., -60, -300), //B
							glm::vec3(100., 100, -300), //C
							glm::vec3(-100., 100, -300), //D
							glm::vec3(0.5, 0.7, 0.95)); //color
							
	//box point definitions
	glm::vec3 bpt1 = glm::vec3(5., -10., -110.);
	glm::vec3 bpt2 = glm::vec3(15., -10., -110.);
	glm::vec3 bpt3 = glm::vec3(5., -10., -100.);
	glm::vec3 bpt4 = glm::vec3(15., -10., -100.);
	glm::vec3 bpt5 = glm::vec3(5., -20., -110.);
	glm::vec3 bpt6 = glm::vec3(15., -20., -110.);
	glm::vec3 bpt7 = glm::vec3(5., -20., -100.);
	glm::vec3 bpt8 = glm::vec3(15., -20., -100.);
	glm::vec3 boxcolour = glm::vec3(0.52,0.21,0.71);
	
	//tetrahedron point definitions
	glm::vec3 tpt1 = glm::vec3(12, -20, -130.);
	glm::vec3 tpt2 = glm::vec3(28, -20, -130.);
	glm::vec3 tpt3 = glm::vec3(20, -20, -114.);
	glm::vec3 tpt4 = glm::vec3(20, -4, -122.);
	glm::vec3 tetrahedroncolour = glm::vec3(1,0,0);
	
	Plane *boxface1 = new Plane(bpt5, bpt6, bpt8, bpt7, boxcolour);
	Plane *boxface2 = new Plane(bpt5, bpt7, bpt3, bpt1, boxcolour);
	Plane *boxface3 = new Plane(bpt7, bpt8, bpt4, bpt3, boxcolour);
	Plane *boxface4 = new Plane(bpt8, bpt6, bpt2, bpt4, boxcolour);
	Plane *boxface5 = new Plane(bpt6, bpt5, bpt1, bpt2, boxcolour);
	Plane *boxface6 = new Plane(bpt3, bpt4, bpt2, bpt1, boxcolour);
	
	Triangle *tetrahedronface1 = new Triangle(tpt1, tpt3, tpt4, tetrahedroncolour);
	Triangle *tetrahedronface2 = new Triangle(tpt3, tpt2, tpt4, tetrahedroncolour);
	Triangle *tetrahedronface3 = new Triangle(tpt2, tpt1, tpt4, tetrahedroncolour);
	Triangle *tetrahedronface4 = new Triangle(tpt1, tpt2, tpt3, tetrahedroncolour);
	
	Cylinder *cyl = new Cylinder(glm::vec3(20.0, -20.0, -150.0), 7.0, 35.0, glm::vec3(1, 0, 0));
	
	Cone *cone = new Cone(glm::vec3(-22.0, -20.0, -140.0), 10.0, 30.0, glm::vec3(0, 1, 0));
							
							
	

	//--Add the above to the list of scene objects.
	sceneObjects.push_back(backgroundplane); 
	sceneObjects.push_back(floorplane); 
	sceneObjects.push_back(sphere1); 
	sceneObjects.push_back(sphere2); 
 	sceneObjects.push_back(glassSphere);
 	sceneObjects.push_back(cyl); 
 	sceneObjects.push_back(cone); 
	sceneObjects.push_back(boxface1);
	sceneObjects.push_back(boxface2); 
	sceneObjects.push_back(boxface3); 
	sceneObjects.push_back(boxface4); 
	sceneObjects.push_back(boxface5); 
	sceneObjects.push_back(boxface6); 
	sceneObjects.push_back(tetrahedronface1);
	sceneObjects.push_back(tetrahedronface2); 
	sceneObjects.push_back(tetrahedronface3); 
	sceneObjects.push_back(tetrahedronface4);

}



int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracer");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
