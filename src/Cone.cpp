/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cone class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cone.h"
#include <math.h>

/**
* Cone's intersection method.  The input is a ray (pos, dir). 
*/
float Cone::intersect(glm::vec3 posn, glm::vec3 dir)
{
    //components of ray equation
    float xpart = posn.x - center.x;
    float zpart = posn.z - center.z;
    float hpart = height - posn.y + center.y;
    float rpart = radius / height;
	//squared components
    float xsqr = xpart * xpart;
    float zsqr = zpart * zpart;
    float hsqr = hpart * hpart;
    float rsqr = rpart * rpart;
        
    float a = (dir.x * dir.x) + (dir.z * dir.z) - (rsqr * (dir.y * dir.y));
    float b = (2 * xpart * dir.x) + (2 * zpart * dir.z) + (2 * rsqr * hpart * dir.y);
    float c = xsqr + zsqr - (rsqr * hsqr);
    
    float delta = b * b - 4 * (a * c);
    
	if (fabs(delta) < 0.001 || delta < 0.0) {
		return -1.0;
	} 
    
    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);
    
    float intersection;
    if (t1 > t2) {
		intersection = t2;
	} else {
		intersection = t1;
	}
    
    //check the point of intersection is within the height of the cone
    float intersection_height = posn.y + intersection * dir.y;
    if ((intersection_height >= center.y) && (intersection_height <= center.y + height)) {
		return intersection;
	} else {
		return -1.0;
	}
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the Cone.
*/
glm::vec3 Cone::normal(glm::vec3 p)
{
	float theta = atan(radius/height);
	float alpha = atan((p.x - center.x) / (p.z - center.z));
    glm::vec3 n = glm::vec3((sin(alpha) * cos(theta)), sin(theta), (cos(alpha) * cos(theta)));
    return n;
}
