/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include <math.h>

/**
* Cylinder's intersection method.  The input is a ray (pos, dir). 
*/
float Cylinder::intersect(glm::vec3 posn, glm::vec3 dir)
{
    float a = (dir.x * dir.x) + (dir.z * dir.z);
    float b = 2 * (dir.x * (posn.x - center.x) + dir.z * (posn.z - center.z));
    float c = ((posn.x - center.x) * (posn.x - center.x)) + ((posn.z - center.z) * (posn.z - center.z)) - (radius * radius);
    
    float delta = b*b - (4*a*c);
    if (fabs(delta) < 0.001) {
		return -1.0;
	}
	
	float t1 = (-b - sqrt(delta)) / (2*a);
	float t2 = (-b + sqrt(delta)) / (2*a);
	
	if (t1 > t2) {
		return t2;
	} else {
		float pt1height = posn.y + t1 * dir.y;
		float pt2height = posn.y + t2 * dir.y;
		if (pt1height > center.y) {
			if (pt1height > center.y + height) {
				if (pt2height <= center.y + height) {
					return t2;
				} else {
					return -1.0;
				}
			} else {
				return t1;
			}
		}
				
	}
	
	return -1.0;
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the Cylinder.
*/
glm::vec3 Cylinder::normal(glm::vec3 p)
{
    glm::vec3 n = glm::vec3((p.x - center.x) / radius, 0, (p.z - center.z) / radius);
    n = glm::normalize(n);
    return n;
}
