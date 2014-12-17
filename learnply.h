/*

Data structures for learnply

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_H__
#define __LEARNPLY_H__

#include<stdio.h>
#include<math.h>
#include<vector>
#include "ply.h"
#include "icVector.H"

using namespace std;

const double EPS = 1.0e-6;
const double PI=3.1415926535898;

/* forward declarations */
class Triangle;
class Corner;

class Vertex {
public:
  double x,y,z;
  int index;

  int ntris;
  Triangle **tris;
  int max_tris;

  vector<Corner*> cList;
  vector<double> normCurvList;
  icVector3 localX,localY;
  double l,m,n;//curvature tensor [l m;m n]
  icVector2 majEigv, minEigv;//local eigenvectors for curvature tensor
  icVector3 gloMajEigv, gloMinEigv;//global eigenvectors for curvature tensor
  double heat;
  int minMaxSad;//-1 min, 0 max, >=1 saddle
  icMatrix3x3 globalTensor;

	icVector3 normal;
  void *other_props;

  double AMixed;//mixed area for curvature computation
  double GaussCurv;
  double meanCurv;

public:
  Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
};

class Edge {
public:
  int index;
  Vertex *verts[2];
  int ntris;
  Triangle **tris;
	double length;
	bool sliEdge;//if true, silhouette edge.
	//bool containSliP;//if true, this edge contains one silhouette point
	Vertex *newVert;//new division point on edge

};

class Triangle {
public:
  int index;
  int nverts;
  Vertex *verts[3];
  Edge *edges[3];

	double angle[3];
	float area;

	bool front;//if true, front face;else, back face.

	icVector3 normal;
  void *other_props;

  bool containSliP;//if true, this edge contains one silhouette point
  icVector3 silPoints[2];
  int silPointIndex;

  bool obtuse;

};

class VertexList {

public:

  int num_verts;
  int max_verts;
  Vertex **verts;

  VertexList (int max) {
    max_verts = max;
    verts = new Vertex *[max_verts];
    num_verts = 0;
  }
	void finalize() {
		if (this != NULL) {
			free(verts);
			free(this);
		}
	}
	void append(Vertex *v)
	{
		int i;

		/* first make sure there is enough room for new vertex */

		if (num_verts >= max_verts) {
			max_verts += 10;
			Vertex **tlist = new Vertex *[max_verts];
			for (i = 0; i < num_verts; i++)
				tlist[i] = verts[i];
			delete (verts);
			verts = tlist;
		}

		/* add new vertex to list */

		verts[num_verts] = v;
		num_verts++;
	}

};

class TriangleList {
public:

  int num_tris;
  int max_tris;
  Triangle **tris;

  TriangleList (int max) {
    max_tris = max;
    tris = new Triangle *[max_tris];
    num_tris = 0;
  }
	void append(Triangle *t)
	{
		int i;

		/* first make sure there is enough room for new triangle */

		if (num_tris >= max_tris) {
			max_tris += 10;
			Triangle **tlist = new Triangle *[max_tris];
			for (i = 0; i < num_tris; i++)
				tlist[i] = tris[i];
			delete (tris);
			tris = tlist;
		}

		/* add new triangle to list */

		tris[num_tris] = t;
		num_tris++;
	}
};

class EdgeList {

public:

  int num_edges;
  int max_edges;
  Edge **edges;

  EdgeList (int max) {
    max_edges = max;
    edges = new Edge *[max_edges];
    num_edges = 0;
  }
	void append(Edge *e)
	{
		int i;

		/* first make sure there is enough room for new edge */

		if (num_edges >= max_edges) {
			max_edges += 10;
			Edge **tlist = new Edge *[max_edges];
			for (i = 0; i < num_edges; i++)
				tlist[i] = edges[i];
			delete (edges);
			edges = tlist;
		}

		/* add new edge to list */

		edges[num_edges] = e;
		num_edges++;
	}
};

class Corner{
public:
	int index;
	Vertex *v;
	Edge *e;
	Triangle *t;
	Corner *p, *n, *o;
	double cos, sin, cot;
    double angle;
	Corner(Triangle *tri, int nVertex)
	{
		t = tri;
		e = tri->edges[(nVertex+1)%3];
		v = tri->verts[nVertex];
		v->cList.push_back(this);
		n = NULL;
		p = NULL;
		o = NULL;
		cos = (pow(tri->edges[nVertex]->length,2) + pow(tri->edges[(nVertex+2)%3]->length,2)
			- pow(tri->edges[(nVertex+1)%3]->length,2))/(2*(tri->edges[nVertex]->length)*(tri->edges[(nVertex+2)%3]->length));
		angle = acos(cos);
		sin = sqrt(1 - cos*cos);
		cot = tan(PI/2.0 - angle);
		//cout<<"cos:	"<<cos<<"	angle:	"<<angle<<"	sin:	"<<sin<<"cot:	"<<cot<<endl;
		//if(cos>1.0 || cos<-1.0)
		//	cout<<tri->index<<pow(tri->edges[nVertex]->length,2)<<"
		//	pow(tri->edges[(nVertex+2)%2]->length,2)<<pow(tri->edges[(nVertex+1)%2]->length,2)<<endl;
		//cout<<"acos:	"<<acos(-0.5)<<endl;
		//cout<<"angle: "<<angle<<endl;
		//minMaxSad = -2;
	}


};

//class CornerList{
//public:
//	int num_corners;
//	int max_corners;
//	Corner** corners;
//
//};

class Polyhedron {
public:

	int index;

  Triangle **tlist;  /* list of triangles */
  int ntris;
  int max_tris;

  Vertex **vlist;    /* list of vertices */
  int nverts;
  int max_verts;

  Edge **elist;      /* list of edges */
  int nedges;
  int max_edges;

  Corner **clist;   /* list of corners of hwk2 */
  int ncorners;
  int max_corners;   

	icVector3 center;
	double radius;
	double area;

	int seed;

  PlyOtherProp *vert_other,*face_other;

	void average_normals();

  void create_edge(Vertex *, Vertex *);
  void create_edges();
  int face_to_vertex_ref(Triangle *, Vertex *);
  void order_vertex_to_tri_ptrs(Vertex *);
  void vertex_to_tri_ptrs();
  Triangle *find_common_edge(Triangle *, Vertex *, Vertex *);
  Triangle *other_triangle(Edge *, Triangle *);

	void calc_bounding_sphere();
	void calc_face_normals_and_area();
	void calc_edge_length();

	Polyhedron();
  Polyhedron(FILE *);
  void write_file(FILE *);

  void create_pointers();

	// initialization and finalization
	void initialize();
	void finalize();
/************************************************************************************/
  void createCList();//create cornerlist
  void findOppositeCorners();//find opposite corners for all corners in the cornerlist

  void markFBFace(icVector3& viewPoint);//mark every face as front or back.
  void markSliEdge();//mark every edge as silhouette edge or not

  void findSilPointOnEdge(icVector3& viewPoint);//find the silhouette points on all the edges

  void loopDivide();//Regular Loop subdivision

  void subDivEdgeL(double threshL);//Irregular subdivision based on edge length

/************************************************************************************/

};

#endif /* __LEARNPLY_H__ */