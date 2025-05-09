#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"
#include <set>

myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,1.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal() {
    if (!originof) return;

    normal->clear();
    std::set<myFace*> visitedFaces;
    myHalfedge* start = originof;
    myHalfedge* current = start;

    // Traverse all incident faces
    do {
        // Check if adjacent_face is not null and valid
        if (!current || !current->adjacent_face || !current->adjacent_face->normal)
            break;

        if (current->adjacent_face && current->adjacent_face->normal &&
            visitedFaces.find(current->adjacent_face) == visitedFaces.end()) {

            // Safely access the face normal
            myVector3D* faceNormal = current->adjacent_face->normal;
            normal->dX += faceNormal->dX;
            normal->dY += faceNormal->dY;
            normal->dZ += faceNormal->dZ;

            visitedFaces.insert(current->adjacent_face);
        }

        // Check if the twin is valid before accessing it
        if (!current->twin || !current->twin->next) break;

        current = current->twin->next;

    } while (current != start && current != nullptr);

    // In case no faces were found (degenerate), leave default normal
    if (!visitedFaces.empty()) {
        normal->normalize();
    }
}


