#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"
#include <unordered_set>
#include <algorithm>

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

			// Set 'next' and 'prev' pointers
			for (size_t i = 0; i < hedges.size(); ++i) {
				hedges[i]->next = hedges[(i + 1) % hedges.size()];
				hedges[i]->prev = hedges[(i + hedges.size() - 1) % hedges.size()];
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
	std::vector<myVertex*> new_vertices;
	std::vector<myFace*> new_faces;
	std::vector<myHalfedge*> new_halfedges;

	std::map<myFace*, myVertex*> face_points;
	std::map<myHalfedge*, myVertex*> edge_points;
	std::map<myVertex*, myVertex*> updated_vertices;

	for (myFace* f : this->faces)
	{
		myPoint3D* center_p = new myPoint3D(0, 0, 0);
		int vertex_count = 0;
		myHalfedge* h_start = f->adjacent_halfedge;
		myHalfedge* h_current = h_start;
		do
		{
			*center_p = *center_p + *(h_current->source->point);
			vertex_count++;
			h_current = h_current->next;
		} while (h_current != h_start);

		*center_p = *center_p / static_cast<float>(vertex_count);

		myVertex* v_face = new myVertex();
		v_face->point = center_p;
		face_points[f] = v_face;
		new_vertices.push_back(v_face);
	}

	for (myHalfedge* h : this->halfedges)
	{
		if (edge_points.find(h) == edge_points.end())
		{
			myPoint3D* edge_p = new myPoint3D(0, 0, 0);
			*edge_p = *edge_p + *(h->source->point);
			*edge_p = *edge_p + *(h->twin->source->point);
			*edge_p = *edge_p + *(face_points[h->adjacent_face]->point);
			*edge_p = *edge_p + *(face_points[h->twin->adjacent_face]->point);
			*edge_p = *edge_p / 4.0f;

			myVertex* v_edge = new myVertex();
			v_edge->point = edge_p;
			edge_points[h] = v_edge;
			edge_points[h->twin] = v_edge;
			new_vertices.push_back(v_edge);
		}
	}

	for (myVertex* v_orig : this->vertices)
	{
		myPoint3D avg_face_points(0, 0, 0);
		myPoint3D avg_edge_midpoints(0, 0, 0);
		int n = 0;

		myHalfedge* h_start = v_orig->originof;
		myHalfedge* h_current = h_start;
		do
		{
			avg_face_points = avg_face_points + *(face_points[h_current->adjacent_face]->point);

			myPoint3D edge_midpoint = (*(h_current->source->point) + *(h_current->twin->source->point)) / 2.0f;
			avg_edge_midpoints = avg_edge_midpoints + edge_midpoint;

			n++;
			h_current = h_current->twin->next;
		} while (h_current != h_start);

		avg_face_points = avg_face_points / static_cast<float>(n);
		avg_edge_midpoints = avg_edge_midpoints / static_cast<float>(n);

		myPoint3D* new_pos = new myPoint3D(
			(avg_face_points.X + 2 * avg_edge_midpoints.X + v_orig->point->X * (n - 3)) / n,
			(avg_face_points.Y + 2 * avg_edge_midpoints.Y + v_orig->point->Y * (n - 3)) / n,
			(avg_face_points.Z + 2 * avg_edge_midpoints.Z + v_orig->point->Z * (n - 3)) / n
		);

		myVertex* v_updated = new myVertex();
		v_updated->point = new_pos;
		updated_vertices[v_orig] = v_updated;
		new_vertices.push_back(v_updated);
	}

	for (myFace* f_orig : this->faces)
	{
		myHalfedge* h_start = f_orig->adjacent_halfedge;
		myHalfedge* h_current = h_start;
		do
		{
			myVertex* v1 = face_points[f_orig];
			myVertex* v2 = edge_points[h_current];
			myVertex* v3 = updated_vertices[h_current->source];
			myVertex* v4 = edge_points[h_current->prev];

			myFace* new_q = new myFace();
			myHalfedge* h1 = new myHalfedge(); h1->source = v1; h1->adjacent_face = new_q;
			myHalfedge* h2 = new myHalfedge(); h2->source = v2; h2->adjacent_face = new_q;
			myHalfedge* h3 = new myHalfedge(); h3->source = v3; h3->adjacent_face = new_q;
			myHalfedge* h4 = new myHalfedge(); h4->source = v4; h4->adjacent_face = new_q;

			h1->next = h2; h2->next = h3; h3->next = h4; h4->next = h1;
			h1->prev = h4; h2->prev = h1; h3->prev = h2; h4->prev = h3;

			v1->originof = h1;
			v2->originof = h2;
			v3->originof = h3;
			v4->originof = h4;

			new_q->adjacent_halfedge = h1;

			new_faces.push_back(new_q);
			new_halfedges.push_back(h1); new_halfedges.push_back(h2);
			new_halfedges.push_back(h3); new_halfedges.push_back(h4);

			h_current = h_current->next;
		} while (h_current != h_start);
	}
	
	clear();

	this->vertices = new_vertices;
	this->faces = new_faces;
	this->halfedges = new_halfedges;

	std::map<std::pair<myVertex*, myVertex*>, myHalfedge*> twin_finder;
	for (myHalfedge* h : this->halfedges)
	{
		myVertex* v_start = h->source;
		myVertex* v_end = h->next->source;
		auto it = twin_finder.find({ v_end, v_start });
		if (it != twin_finder.end())
		{
			h->twin = it->second;
			it->second->twin = h;
			twin_finder.erase(it);
		}
		else
		{
			twin_finder[{v_start, v_end}] = h;
		}
	}
}


void myMesh::triangulate() {
	std::vector<myFace*> original_faces = faces;
	for (myFace* f : original_faces) {
		triangulate(f);
	}
}

bool myMesh::triangulate(myFace* f)
{
	
	myHalfedge* start = f->adjacent_halfedge;

	if (start->next->next->next == start)
		return false;

	myVertex* commonSource = start->source;
	myHalfedge* current = start;
	myFace* currentFace = f;

	while (current->next->next->next != current) {
		myHalfedge* he1 = new myHalfedge();
		myHalfedge* he2 = new myHalfedge();

		he1->twin = he2;
		he2->twin = he1;

		he1->source = current->next->next->source;
		he2->source = commonSource;

		myFace* newFace = new myFace();
		newFace->adjacent_halfedge = he2;

		myHalfedge* he2Next = current->next->next;
		myHalfedge* he2Prev = current->prev;

		he1->next = current;
		he1->prev = current->next;
		current->next->next = he1;
		current->prev = he1;

		he2->next = he2Next;
		he2->prev = he2Prev;
		he2Next->prev = he2;
		he2Prev->next = he2;

		myHalfedge* cur = he2;
		do {
			cur->adjacent_face = newFace;
			cur = cur->next;
		} while (cur != he2);

		he1->adjacent_face = currentFace;

		faces.push_back(newFace);
		halfedges.push_back(he1);
		halfedges.push_back(he2);

		current = he2;
	}

	return true;
}

void myMesh::simplify() {
	myHalfedge* he = findShortestEdge();
	if (!he || !he->twin) return;

	myHalfedge* heTwin = he->twin;
	myVertex* v1 = he->source;
	myVertex* v2 = he->next->source;

	if (!v1 || !v2 || v1 == v2) {
		std::cerr << "Triangulez pour aller plus loin" << std::endl;
		return;
	}

	myPoint3D* mid = new myPoint3D(
		(v1->point->X + v2->point->X) / 2,
		(v1->point->Y + v2->point->Y) / 2,
		(v1->point->Z + v2->point->Z) / 2
	);
	myVertex* newVertex = new myVertex();
	newVertex->point = mid;
	vertices.push_back(newVertex);

	for (myHalfedge* h : halfedges) {
		if (h && h->source && (h->source == v1 || h->source == v2)) {
			h->source = newVertex;
		}
	}

	std::vector<myHalfedge*> heToRemove;

	auto collectHalfedgesFromFace = [&](myFace* f) {
		vertices.erase(std::remove(vertices.begin(), vertices.end(), nullptr), vertices.end());
		faces.erase(std::remove(faces.begin(), faces.end(), nullptr), faces.end());
		halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), nullptr), halfedges.end());
		if (!f || !f->adjacent_halfedge) return;
		myHalfedge* start = f->adjacent_halfedge;
		myHalfedge* h = start;
		std::unordered_set<myHalfedge*> visited;
		int count = 0;

		const int max_edges = 1000;
		do {
			if (!h || visited.count(h) || count++ > max_edges) {
				std::cerr << "Triangulez pour aller plus loin" << std::endl;
				return;
			}
			visited.insert(h);
			heToRemove.push_back(h);
			h = h->next;
		} while (h != start);
		};

	collectHalfedgesFromFace(he->adjacent_face);
	collectHalfedgesFromFace(heTwin->adjacent_face);

	auto faceSize = [](myFace* f) {
		if (!f || !f->adjacent_halfedge) return 0;
		myHalfedge* start = f->adjacent_halfedge;
		myHalfedge* h = start;
		int count = 0;
		do {
			h = h->next;
			count++;
		} while (h && h != start && count < 100);
		return count;
		};

	int f1size = faceSize(he->adjacent_face);
	int f2size = faceSize(heTwin->adjacent_face);

	if (f1size != 3 || f2size != 3) {
		std::cerr << "Triangulez pour aller plus loin" << std::endl;
		return;
	}


	removeFace(he->adjacent_face);
	removeFace(heTwin->adjacent_face);

	for (myHalfedge* h : heToRemove) {
		removeHalfedge(h);
	}

	removeVertex(v1);
	removeVertex(v2);

	std::cout << "Simplification réussie." << std::endl;
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