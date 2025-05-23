#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	vector<myHalfedge *>::iterator it;
	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		if ((*it)->twin == NULL)
			break;

		myFace* f = new myFace();
	}
	if (it != halfedges.end())
		cout << "Error! Not all edges have their twins!\n";
	else cout << "Each edge has a twin!\n";

	checkSameFace();
	checkIsTriangle();
	check8EdgesMax();
}

void myMesh::checkSameFace() 
{
	bool allHalfEdgesHaveSameFace  = true;
	for (myFace* face : faces)
	{
		myHalfedge* start = face->adjacent_halfedge;
		myHalfedge* current = start;
		myFace* currentFace = start->adjacent_face;

		do
		{
			current = current->next;
			if (current->adjacent_face != currentFace) {
				allHalfEdgesHaveSameFace = false;
				cout << "Error : A halfedge has the wrong face\n";
			}
		} while (current != start);
	}
	if (allHalfEdgesHaveSameFace) {
		cout << "All halfedge have the right face\n";
	}
}

void myMesh::checkIsTriangle() 
{
	bool isTriangle = true;
	for (myFace* face : faces)
	{
		int edgeCount = 0;
		myHalfedge* start = face->adjacent_halfedge;
		if (start->prev != start->next->next) {
			cout << "Error : This face is not a triangle\n";
			isTriangle = false;
		}
	}
	if (isTriangle == true) {
		cout << "All faces are triangles\n";
	}
	else {
		cout << "Not all faces are triangles\n";
	}
}

void myMesh::check8EdgesMax()
{
	bool isTriangle = true;
	for (myFace* face : faces)
	{
		int edgeCount = 0;
		myHalfedge* start = face->adjacent_halfedge;
		myHalfedge* current = start;

		do
		{
			edgeCount++;
			current = current->next;
		} while (current != start);

		if (edgeCount > 8)
		{
			isTriangle = false;
			cout << "Error : A face has more than 8 edges\n";
			break;
		}
	}
}

bool myMesh::readFile(std::string filename)
{
	std::ifstream fin(filename);
	if (!fin.is_open()) {
		std::cerr << "Unable to open file: " << filename << std::endl;
		return false;
	}

	std::string line, prefix, token;
	std::map<std::pair<int, int>, myHalfedge*> twin_map;

	while (std::getline(fin, line)) {
		std::stringstream ss(line);
		ss >> prefix;

		if (prefix == "v") {
			float x, y, z;
			ss >> x >> y >> z;
			myVertex* v = new myVertex();
			v->point = new myPoint3D(x, y, z);
			vertices.push_back(v);
		}
		else if (prefix == "f") {
			std::vector<int> faceids;
			while (ss >> token) {
				size_t slash = token.find("/");
				int id = std::stoi((slash == std::string::npos) ? token : token.substr(0, slash));
				faceids.push_back(id - 1); // OBJ is 1-based
			}

			if (faceids.size() < 3)
				continue;

			myFace* face = new myFace();
			faces.push_back(face);

			std::vector<myHalfedge*> hedges;
			for (size_t i = 0; i < faceids.size(); ++i) {
				myHalfedge* he = new myHalfedge();
				he->adjacent_face = face;
				he->source = vertices[faceids[i]];
				hedges.push_back(he);
				halfedges.push_back(he);
			}

			// Link 'next' pointers in circular fashion
			for (size_t i = 0; i < hedges.size(); ++i) {
				hedges[i]->next = hedges[(i + 1) % hedges.size()];
			}

			face->adjacent_halfedge = hedges[0];

			// Build twin links
			for (size_t i = 0; i < hedges.size(); ++i) {
				int from = faceids[i];
				int to = faceids[(i + 1) % faceids.size()];
				auto key = std::make_pair(to, from); // search reverse

				auto it = twin_map.find(key);
				if (it != twin_map.end()) {
					hedges[i]->twin = it->second;
					it->second->twin = hedges[i];
				}
				else {
					twin_map[std::make_pair(from, to)] = hedges[i];
				}
			}

			// Assign originof if not already set
			for (size_t i = 0; i < faceids.size(); ++i) {
				if (vertices[faceids[i]]->originof == nullptr) {
					vertices[faceids[i]]->originof = hedges[i];
				}
			}
		}
	}

	checkMesh();   // Optional mesh validation
	normalize();   // Optional rescaling/centering

	return true;
}




void myMesh::computeNormals()
{
	for (myFace* face : faces) {
		face->computeNormal();
	}

	for (myVertex* vertex : vertices) {
		vertex->computeNormal();
	}
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}


void myMesh::splitFaceTRIS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{

	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}


void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
}


void myMesh::triangulate() {
	std::vector<myFace*> original_faces = faces;
	for (myFace* f : original_faces) {
		triangulate(f);
	}
}

bool myMesh::triangulate(myFace* f) {
	myHalfedge* start = f->adjacent_halfedge;

	if (start->next->next->next == start)
		return false;

	std::vector<myVertex*> verts;
	myHalfedge* curr = start;
	do {
		verts.push_back(curr->source);
		curr = curr->next;
	} while (curr != start);

	std::vector<myHalfedge*> new_halfedges;
	for (int i = 1; i < (int)verts.size() - 1; ++i) {
		myHalfedge* e0 = new myHalfedge();
		myHalfedge* e1 = new myHalfedge();
		myHalfedge* e2 = new myHalfedge();

		myFace* tri = new myFace();
		tri->adjacent_halfedge = e0;

		e0->source = verts[0];
		e1->source = verts[i];
		e2->source = verts[i + 1];

		e0->next = e1; e1->next = e2; e2->next = e0;
		e0->prev = e2; e1->prev = e0; e2->prev = e1;

		e0->adjacent_face = tri;
		e1->adjacent_face = tri;
		e2->adjacent_face = tri;

		new_halfedges.push_back(e0);
		new_halfedges.push_back(e1);
		new_halfedges.push_back(e2);

		halfedges.push_back(e0);
		halfedges.push_back(e1);
		halfedges.push_back(e2);
		faces.push_back(tri);
	}

	myHalfedge* h = start;
	do {
		myHalfedge* next = h->next;
		halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), h), halfedges.end());
		delete h;
		h = next;
	} while (h != start);

	faces.erase(std::remove(faces.begin(), faces.end(), f), faces.end());
	delete f;

	for (myVertex* v : verts) {
		for (myHalfedge* h : new_halfedges) {
			if (h->source == v) {
				v->originof = h;
				break;
			}
		}
	}

	return true;
}

void myMesh::simplify() {
	myHalfedge* halfedge = findShortestEdge();
	myHalfedge* nexth = halfedge->next;
	myHalfedge* twinh = halfedge->twin;
	for (auto h : this->halfedges) {
		if (!h->twin) {
			std::cout << "This halfege has no twin";
		}
		else {
			std::cout << "HAS TWIN";
		}
	}
	int twins = 0;
	for (auto h : halfedges) {
		if (h->twin) twins++;
	}
	std::cout << "Total halfedges: " << halfedges.size() << "\n";
	std::cout << "Halfedges with twins: " << twins << "\n";
	if (!halfedge || !halfedge->twin) {
		std::cout << "Skipping simplification: no valid twin found.\n";
		return;
	}
	if (twinh == nullptr) {
		std::cout << "twinh is nullptr";
	}
	/*
	myHalfedge* start = halfedge;
	myHalfedge* current = halfedge;
	myHalfedge* currenttwin = twinh;
	myHalfedge* starttwin = twinh;
	*/
	myHalfedge* nexttwin = nexth->twin;
	myHalfedge* twinnexttwin = twinh->next->twin;
	

	nexth->twin = nexth->next->twin;
	nexth->next->twin = nexttwin;

	twinh->next->twin = twinh->next->next->twin;
	twinh->next->next->twin = twinnexttwin;

	halfedge->source = nexth->source;
	
	removeVertex(nexth->source);

	removeFace(halfedge->adjacent_face);
	removeFace(twinh->adjacent_face);

	removeHalfedge(nexttwin->next);
	removeHalfedge(nexttwin);
	removeHalfedge(twinh);
	removeHalfedge(nexth->next);
	removeHalfedge(nexth);
	removeHalfedge(halfedge);

	
	
	/*
	while (current != start) {
		current = current->next;
		removeHalfedge(current);
	}
	

	while (currenttwin != starttwin) {
		currenttwin = currenttwin->next;
		removeHalfedge(currenttwin);
	}
	*/

}

void myMesh::removeVertex(myVertex* v) {
	auto it = std::find(vertices.begin(), vertices.end(), v);
	vertices.erase(it);
	delete v;
}

void myMesh::removeFace(myFace* f) {
	auto it = std::find(faces.begin(), faces.end(), f);
	faces.erase(it);
	delete f;
}

void myMesh::removeHalfedge(myHalfedge* h) {
	auto it = std::find(halfedges.begin(), halfedges.end(), h);
	halfedges.erase(it);
	delete h;
}

myHalfedge* myMesh::findShortestEdge() {
	double shortestLength = 0;
	myHalfedge* shortestHalfedge = halfedges[0];
	for (myHalfedge* halfedge : halfedges) {
		myHalfedge* nexth = halfedge->next;
		myPoint3D p1 = *halfedge->source->point;
		myPoint3D p2 = *nexth->source->point;
		myVector3D vect(p2.X - p1.X, p2.Y - p1.Y, p2.Z - p1.Z);
		if (vect.length() < shortestLength || shortestLength == 0) {
			shortestLength = vect.length();
			shortestHalfedge = halfedge;
		}
	}
	cout << shortestLength;
	return shortestHalfedge;
}