#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal()
{
	myHalfedge *h1 = adjacent_halfedge;
	myPoint3D *v1 = h1->source->point;
	myHalfedge *h2 = h1->next;
	myPoint3D *v2 = h2->source->point;
	myHalfedge *h3 = h2->next;
	myPoint3D *v3 = h3->source->point;
	myVector3D v1v2(v2->X - v1->X, v2->Y - v1->Y, v2->Z - v1->Z);
	myVector3D v2v3(v3->X - v2->X, v3->Y - v2->Y, v3->Z - v2->Z);
	*normal = v1v2.crossproduct(v2v3);
	normal->normalize();
}
