#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <vector>
#include <string>

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void checkMesh();
	void checkSameFace();
	void checkIsTriangle();
	void check8EdgesMax();
	bool readFile(std::string filename);
	void computeNormals();
	void normalize();
	void simplify();

	void removeVertex(myVertex* v);
	void removeFace(myFace* f);
	void removeHalfedge(myHalfedge* h);
	myHalfedge* findShortestEdge();

	void subdivisionCatmullClark();

	void splitFaceTRIS(myFace *, myPoint3D *);

	void splitEdge(myHalfedge *, myPoint3D *);
	void splitFaceQUADS(myFace *, myPoint3D *);

	void triangulate();
	bool triangulate(myFace *);

	void clear();

	myMesh(void);
	~myMesh(void);
};

