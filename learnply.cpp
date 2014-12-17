/*

Functions for learnply

Eugene Zhang, 2005
*/
#include<windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "glut.h"
#include <gl.h>			// Header File For The OpenGL32 Library
#include <glu.h>			// Header File For The GLu32 Library
#include <glaux.h>		// Header File For The Glaux Library
#include <string>
#include <fstream>
#include <iostream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include "trackball.h"
#include "tmatrix.h"
#include <iomanip>
#include <cmath>
#include <time.h>
#include <algorithm>
#include "nr.h"

//using namespace std;

static PlyFile *in_ply;

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width=1024;
const int win_height=1024;

double radius_factor = 0.9;

int display_mode = 0; 
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;	
int ACSIZE = 1; // for antialiasing
int view_mode=0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

struct jitter_struct{
	double x;
	double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125}, 
						  {0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375}, 
						  {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625}, 
						  {0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}, };

Polyhedron *poly;

void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);

/********************************************************/

icVector3 viewPoint(0.0,0.0,0.0);
void drawSilFace(Polyhedron *poly);
	//draw silhouette based on the normal of face
void drawSilVert(Polyhedron *poly);
	//draw silhouette based on the normal of vertex
void checkerboard(Polyhedron *this_poly,float L);

void smoothWgtUpdated(int scheme, double dt, int times, Polyhedron* poly);
void smoothWgtNotUpdated(int scheme, double dt, int times, Polyhedron* poly);
void heatDiffusion(int maxVertId, int minVertId, Polyhedron* poly);
void visHeat(Polyhedron* poly);
void identifyMMS(Polyhedron* poly);
void drawTexture(Polyhedron* poly);
void textureSyn(int maxVertId, int minVertId, int maxVertId2, Polyhedron* poly);
void computeAMixed(Polyhedron* poly);
void computeMeanCurv(Polyhedron* poly);
void computeGaussCurv(Polyhedron* poly);
void visCurv(Polyhedron* poly, int curvType);
void smoothCurv(int scheme, double dt, int times, int curvType, Polyhedron* poly);
void computeNormCurv(Polyhedron* poly);
void computeLocalFrm(Polyhedron* poly);
void computeTensor(Polyhedron* poly);
void computeGlobalTensor(Polyhedron* poly);
void smoothGlobalTensor(int scheme, double dt, int times, Polyhedron* poly);
void globTensorToLocal(Polyhedron* poly);
void computeTensorEig(Polyhedron* poly);
void drawCurvTensor(Polyhedron* poly);

double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*10)/CLOCKS_PER_SEC;
	return diffms;
}

Vec_INT *ija_p;
Vec_DP *sa_p;
int hitNum = 0;
int triSel[3];

/********************************************************/

/******************************************************************************
Main program.
******************************************************************************/
int operationNum = -1; //tell program to do which operation
float paraL = 0.1;//parameter L of checkerboard
bool checkerBrd = false;

int main(int argc, char *argv[])
{
  char *progname;
  int num = 1;
  FILE *this_file;

  progname = argv[0];

	//this_file = fopen("../tempmodels/icosahedron.ply", "r");
	//this_file = fopen("../tempmodels/cube.ply", "r");
   // this_file = fopen("../tempmodels/bunny.ply", "r");
	//poly = new Polyhedron (this_file);
	//fclose(this_file);
  	cout<<"We have the following ply models:"<<endl;
	cout<<"bunny"<<endl;
	cout<<"dodecahedron"<<endl;
	cout<<"dragon"<<endl;
	cout<<"feline"<<endl;
	cout<<"dragon"<<endl;
	cout<<"happy"<<endl;
	cout<<"hexahedron"<<endl;
	cout<<"icosahedron"<<endl;
	cout<<"octahedron"<<endl;
	cout<<"sphere"<<endl;
	cout<<"tetrahedron"<<endl;
	cout<<"torus"<<endl;
	cout<<"cube"<<endl;
	cout<<"Please input the name of the ply model you want to use(only these 11 models listed above):"<<endl;
	string plyName;
	cin>>plyName;
	plyName +=".ply";
	string pathName = "../tempmodels/";

	string fileName = pathName + plyName;
	char fileNameChar[128];
	strcpy(fileNameChar,fileName.c_str());

	this_file = fopen(fileNameChar, "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);
	mat_ident( rotmat );	

	poly->initialize(); // initialize everything
	
	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();
	poly->createCList();
	poly->findOppositeCorners();

	//if(hitNum==2)
	//{
	//int maxId = poly->tlist[triSel[0]]->verts[0]->index;
	//int minId = poly->tlist[triSel[1]]->verts[1]->index;
	//heatDiffusion(maxId,minId,1.0,poly);
	//}
	//heatDiffusion(5, 0, 1.0, poly);

	//double threshL = 1.0;//threshold of irregular subdivision based on edge length

	cout<<"We have these operations and their ID below"<<endl;
	cout<<"ID  operation"<<endl;
	cout<<"0 Display the polyhedron"<<endl;
	cout<<"1  Smooth"<<endl;
	cout<<"2  Heat Diffusion(you need to select two triangles to"<<
		" provide the max and min points), identify local max, min and saddle, print out the value of M "<<endl;
	cout<<"3  Texture Synthesis(you need to select three triangles"<<
		" to provide F and G for the perpendicular levelsets"<<endl;
	cout<<"4	Checkerboard"<<endl;
	cout<<"5	Gaussian Curvature"<<endl;
	cout<<"6	Mean Curvature"<<endl;
	cout<<"7	Curvature Tensor"<<endl;
	cin>>operationNum;

	/************************************************************************/
	//Hwk2

	if(operationNum==1)
	{
		cout<<"Please select one weighting scheme:"<<endl;
		cout<<"0	uniform scheme"<<endl;
		cout<<"1	cord scheme with recomputed weight"<<endl;
		cout<<"2	mean curvature scheme with recomputed weight"<<endl;
		cout<<"3	mean values with recomputed weight"<<endl;
		cout<<"4	cord scheme without recomputed weight"<<endl;
		cout<<"5	mean curvature without recomputed weight"<<endl;
		cout<<"6	mean values without recomputed weight"<<endl;
		int wghScheme;
		cin>>wghScheme;
		cout<<"How many times do you want to run the smoothing?"<<endl;
		int times;
		cin>>times;
		//for(int i=0;i<7;++i)
		//{
		clock_t begin=clock();
		if(wghScheme<4)
		smoothWgtUpdated(wghScheme, 0.5,times,poly);
		else
		smoothWgtNotUpdated(wghScheme-3, 0.5,times,poly);
		clock_t end=clock();
		cout << "Time elapsed for doing smooth: " << double(diffclock(end,begin)) << " ms"<< endl;
		//++wghScheme;
		//poly->finalize();
		//poly->finalize();
		//FILE* that_file = fopen("../tempmodels/bunny.ply","r");
		//poly = new Polyhedron (that_file);
		//poly->initialize(); // initialize everything
	
		//poly->calc_bounding_sphere();
		//poly->calc_face_normals_and_area();
		//poly->average_normals();
		//poly->createCList();
		//poly->findOppositeCorners();
		//fclose(that_file);
		//}
		cout<<"You can either display the smoothed object or do the checkerboard."<<endl;
		cout<<"1	display the object"<<endl;
		cout<<"2	checkerboard"<<endl;
		int checker;
		cin>>checker;
		if(checker==2)
			checkerBrd = true;
	}
	else if(operationNum==5)
	{
		//smoothWgtUpdated(1, 0.1,5,poly);
		//poly->calc_bounding_sphere();
	 //   poly->calc_face_normals_and_area();
	 //   poly->average_normals();
	 //   poly->createCList();
	 //   poly->findOppositeCorners();
		computeGaussCurv(poly);
		smoothCurv(2, 0.1, 500, 1, poly);
	}
	else if(operationNum==6)
	{
		//smoothWgtUpdated(1, 0.1,5,poly);
		//poly->calc_bounding_sphere();
	 //   poly->calc_face_normals_and_area();
	 //   poly->average_normals();
	 //   poly->createCList();
	 //   poly->findOppositeCorners();
		computeMeanCurv(poly);
		smoothCurv(2, 0.05, 500, 2, poly);
	}
	else if(operationNum==7)
		computeTensorEig(poly);


	//viewPoint.set(poly->center.x,poly->center.y,poly->center.z+10.0);
	//if(operationNum==2)
	//{
	//	poly->markFBFace(viewPoint);
	//	poly->markSliEdge();
	//}

	//if(operationNum==3)
	//	poly->findSilPointOnEdge(viewPoint);

	//if(operationNum==4)
	//{
	//	poly->loopDivide();
	//	poly->initialize(); // initialize everything
	//	int EulerChar = poly->nverts - poly->nedges + poly->ntris;
	//    cout<<endl<<"Euler Characteristics after regular subdivision is:"<<EulerChar<<endl;
	//	cout<<"Please input the parameter(float type) L for checkerboard"<<endl;
	//	cin>>paraL;
	//}

	//if(operationNum==5)
	//{
	//	cout<<"Please input the parameter(double type) threshL for irregular subdivision"<<endl;
	//	cin>>threshL;
	//	poly->subDivEdgeL(threshL);
	//	poly->initialize(); // initialize everything
	//	int EulerChar = poly->nverts - poly->nedges + poly->ntris;
	//    cout<<endl<<"Euler Characteristics after irregular subdivision is:"<<EulerChar<<endl;
	//	cout<<"Please input the parameter(float type) L for checkerboard"<<endl;
	//	cin>>paraL;
	//}
	//poly->initialize(); // initialize everything
	
	//poly->calc_bounding_sphere();
	//poly->calc_face_normals_and_area();
	//poly->average_normals();
	//poly->subDivEdgeL(0.01);

	

	/************************************************************************/

	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition (20, 20);
	glutInitWindowSize (win_width, win_height); 
	glutCreateWindow ("Geometric Modeling");
	init ();
	glutKeyboardFunc (keyboard);
	glutDisplayFunc(display); 
	glutMotionFunc (motion);
	glutMouseFunc (mouse);
	glutMainLoop(); 
	poly->finalize();  // finalize everything

  return 0;    /* ANSI C requires main to return int. */
}

void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0){
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0/3){
		col[0] = 1.0;
		col[1] = 1.0-percentage*3.0;
		col[2] = 1.0-percentage*3.0;
	}
	else if (percentage <= 2.0/3){
		col[0] = 1.0;
		col[1] = percentage*3.0-1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0/3){
		col[0] = 3.0-percentage*3.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/

Polyhedron::Polyhedron(FILE *file)
{
  int i,j;
  int elem_count;
  char *elem_name;

  /*** Read in the original PLY object ***/
  in_ply = read_ply (file);

  for (i = 0; i < in_ply->num_elem_types; i++) {

    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply (in_ply, i, &elem_count);

    if (equal_strings ("vertex", elem_name)) {

      /* create a vertex list to hold all the vertices */
      //nverts = max_verts = elem_count;
	  nverts = elem_count;
	  max_verts = 1000000;
      vlist = new Vertex *[max_verts];

      /* set up for getting vertex elements */

      setup_property_ply (in_ply, &vert_props[0]);
      setup_property_ply (in_ply, &vert_props[1]);
      setup_property_ply (in_ply, &vert_props[2]);
      vert_other = get_other_properties_ply (in_ply, 
					     offsetof(Vertex_io,other_props));

      /* grab all the vertex elements */
      for (j = 0; j < nverts; j++) {
        Vertex_io vert;
        get_element_ply (in_ply, (void *) &vert);

        /* copy info from the "vert" structure */
        vlist[j] = new Vertex (vert.x, vert.y, vert.z);
        vlist[j]->other_props = vert.other_props;
      }
    }
    else if (equal_strings ("face", elem_name)) {

      /* create a list to hold all the face elements */
      //ntris = max_tris = elem_count;
	  ntris =  elem_count;
	  max_tris = 1000000;
      tlist = new Triangle *[max_tris];

      /* set up for getting face elements */
      setup_property_ply (in_ply, &face_props[0]);
      face_other = get_other_properties_ply (in_ply, offsetof(Face_io,other_props));

      /* grab all the face elements */
      for (j = 0; j < elem_count; j++) {
        Face_io face;
        get_element_ply (in_ply, (void *) &face);

        if (face.nverts != 3) {
          fprintf (stderr, "Face has %d vertices (should be three).\n",
                   face.nverts);
          exit (-1);
        }

        /* copy info from the "face" structure */
        tlist[j] = new Triangle;
        tlist[j]->nverts = 3;
        tlist[j]->verts[0] = (Vertex *) face.verts[0];
        tlist[j]->verts[1] = (Vertex *) face.verts[1];
        tlist[j]->verts[2] = (Vertex *) face.verts[2];
        tlist[j]->other_props = face.other_props;
      }
    }
    else
      get_other_element_ply (in_ply);
  }

  /* close the file */
  close_ply (in_ply);

  /* fix up vertex pointers in triangles */
  for (i = 0; i < ntris; i++) {
    tlist[i]->verts[0] = vlist[(int) tlist[i]->verts[0]];
    tlist[i]->verts[1] = vlist[(int) tlist[i]->verts[1]];
    tlist[i]->verts[2] = vlist[(int) tlist[i]->verts[2]];
  }

  /* get rid of triangles that use the same vertex more than once */

  for (i = ntris-1; i >= 0; i--) {

    Triangle *tri = tlist[i];
    Vertex *v0 = tri->verts[0];
    Vertex *v1 = tri->verts[1];
    Vertex *v2 = tri->verts[2];

    if (v0 == v1 || v1 == v2 || v2 == v0) {
      free (tlist[i]);
      ntris--;
      tlist[i] = tlist[ntris];
    }
  }
}


/******************************************************************************
Write out a polyhedron to a file.
******************************************************************************/

void Polyhedron::write_file(FILE *file)
{
  int i;
  PlyFile *ply;
  char **elist;
  int num_elem_types;

  /*** Write out the transformed PLY object ***/

  elist = get_element_list_ply (in_ply, &num_elem_types);
  ply = write_ply (file, num_elem_types, elist, in_ply->file_type);

  /* describe what properties go into the vertex elements */

  describe_element_ply (ply, "vertex", nverts);
  describe_property_ply (ply, &vert_props[0]);
  describe_property_ply (ply, &vert_props[1]);
  describe_property_ply (ply, &vert_props[2]);
//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

  describe_element_ply (ply, "face", ntris);
  describe_property_ply (ply, &face_props[0]);

//  describe_other_properties_ply (ply, face_other,
//                                offsetof(Face_io,other_props));

//  describe_other_elements_ply (ply, in_ply->other_elems);

  copy_comments_ply (ply, in_ply);
	char mm[1024];
	sprintf(mm, "modified by learnply");
//  append_comment_ply (ply, "modified by simvizply %f");
	  append_comment_ply (ply, mm);
  copy_obj_info_ply (ply, in_ply);

  header_complete_ply (ply);

  /* set up and write the vertex elements */
  put_element_setup_ply (ply, "vertex");
  for (i = 0; i < nverts; i++) {
    Vertex_io vert;

    /* copy info to the "vert" structure */
    vert.x = vlist[i]->x;
    vert.y = vlist[i]->y;
    vert.z = vlist[i]->z;
    vert.other_props = vlist[i]->other_props;

    put_element_ply (ply, (void *) &vert);
  }

  /* index all the vertices */
  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  /* set up and write the face elements */
  put_element_setup_ply (ply, "face");

  Face_io face;
  face.verts = new int[3];
  
  for (i = 0; i < ntris; i++) {

    /* copy info to the "face" structure */
    face.nverts = 3;
    face.verts[0] = tlist[i]->verts[0]->index;
    face.verts[1] = tlist[i]->verts[1]->index;
    face.verts[2] = tlist[i]->verts[2]->index;
    face.other_props = tlist[i]->other_props;

    put_element_ply (ply, (void *) &face);
  }
  put_other_elements_ply (ply);

  close_ply (ply);
  free_ply (ply);
}

void Polyhedron::initialize(){
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	seed = -1;
}

void Polyhedron::finalize(){
	int i;

	for (i=0; i<ntris; i++){
		//free(tlist[i]->other_props);
		free(tlist[i]);
	}
	for (i=0; i<nedges; i++) {
		free(elist[i]->tris);
		free(elist[i]);
	}
	for (i=0; i<nverts; i++) {
		free(vlist[i]->tris);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}
	for (i=0; i<ncorners; i++) {
		free(clist[i]);
	}
	free(tlist);
	free(elist);
	free(vlist);
	free(clist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
  f1    - face that we're looking to share with
  v1,v2 - two vertices of f1 that define edge

Exit:
  return the matching face, or NULL if there is no such face
******************************************************************************/

Triangle *Polyhedron::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f2;
  Triangle *adjacent = NULL;

  /* look through all faces of the first vertex */

  for (i = 0; i < v1->ntris; i++) {
    f2 = v1->tris[i];
    if (f2 == f1)
      continue;
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < f2->nverts; j++) {

      /* look for a match */
      if (f2->verts[j] == v2) {

#if 0
	/* watch out for triple edges */

        if (adjacent != NULL) {

	  fprintf (stderr, "model has triple edges\n");

	  fprintf (stderr, "face 1: ");
	  for (k = 0; k < f1->nverts; k++)
	    fprintf (stderr, "%d ", f1->iverts[k]);
	  fprintf (stderr, "\nface 2: ");
	  for (k = 0; k < f2->nverts; k++)
	    fprintf (stderr, "%d ", f2->iverts[k]);
	  fprintf (stderr, "\nface 3: ");
	  for (k = 0; k < adjacent->nverts; k++)
	    fprintf (stderr, "%d ", adjacent->iverts[k]);
	  fprintf (stderr, "\n");

	}

	/* if we've got a match, remember this face */
        adjacent = f2;
#endif

#if 1
	/* if we've got a match, return this face */
        return (f2);
#endif

      }
    }
  }

  return (adjacent);
}


/******************************************************************************
Create an edge.

Entry:
  v1,v2 - two vertices of f1 that define edge
******************************************************************************/

void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f;

  /* make sure there is enough room for a new edge */

  if (nedges >= max_edges) {

    max_edges += 100;
    Edge **list = new Edge *[max_edges];

    /* copy the old list to the new one */
    for (i = 0; i < nedges; i++)
      list[i] = elist[i];

    /* replace list */
    free (elist);
    elist = list;
  }

  /* create the edge */

  elist[nedges] = new Edge;
  Edge *e = elist[nedges];
  e->index = nedges;
  e->verts[0] = v1;
  e->verts[1] = v2;
  nedges++;

  /* count all triangles that will share the edge, and do this */
  /* by looking through all faces of the first vertex */

  e->ntris = 0;

  for (i = 0; i < v1->ntris; i++) {
    f = v1->tris[i];
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++) {
      /* look for a match */
      if (f->verts[j] == v2) {
        e->ntris++;
        break;
      }
    }
  }

  /* make room for the face pointers (at least two) */
  if (e->ntris < 2)
    e->tris = new Triangle *[2];
  else
    e->tris = new Triangle *[e->ntris];

  /* create pointers from edges to faces and vice-versa */

  e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

  for (i = 0; i < v1->ntris; i++) {

    f = v1->tris[i];

    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++)
      if (f->verts[j] == v2) {

        e->tris[e->ntris] = f;
        e->ntris++;

        if (f->verts[(j+1)%3] == v1)
          f->edges[j] = e;
        else if (f->verts[(j+2)%3] == v1)
          f->edges[(j+2)%3] = e;
        else {
          fprintf (stderr, "Non-recoverable inconsistancy in create_edge()\n");
          exit (-1);
        }

        break;  /* we'll only find one instance of v2 */
      }

  }
}


/******************************************************************************
Create edges.
******************************************************************************/

void Polyhedron::create_edges()
{
  int i,j;
  Triangle *f;
  Vertex *v1,*v2;
  double count = 0;

  /* count up how many edges we may require */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      Triangle *result = find_common_edge (f, v1, v2);
      if (result)
        count += 0.5;
      else
        count += 1;
    }
  }

  /*
  printf ("counted %f edges\n", count);
  */

  /* create space for edge list */

  max_edges = (int) (count + 10);  /* leave some room for expansion */
  elist = new Edge *[max_edges];
  nedges = 0;

  /* zero out all the pointers from faces to edges */

  for (i = 0; i < ntris; i++)
    for (j = 0; j < 3; j++)
      tlist[i]->edges[j] = NULL;

  /* create all the edges by examining all the triangles */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < 3; j++) {
      /* skip over edges that we've already created */
      if (f->edges[j])
        continue;
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      create_edge (v1, v2);
    }
  }
}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Polyhedron::vertex_to_tri_ptrs()
{
  int i,j;
  Triangle *f;
  Vertex *v;

  /* zero the count of number of pointers to faces */

  for (i = 0; i < nverts; i++)
    vlist[i]->max_tris = 0;

  /* first just count all the face pointers needed for each vertex */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++)
      f->verts[j]->max_tris++;
  }

  /* allocate memory for face pointers of vertices */

  for (i = 0; i < nverts; i++) {
    vlist[i]->tris = (Triangle **)
		      malloc (sizeof (Triangle *) * vlist[i]->max_tris);
    vlist[i]->ntris = 0;
  }

  /* now actually create the face pointers */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v = f->verts[j];
      v->tris[v->ntris] = f;
      v->ntris++;
    }
  }
}


/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/

Triangle *Polyhedron::other_triangle(Edge *edge, Triangle *tri)
{
  /* search for any other triangle */

  for (int i = 0; i < edge->ntris; i++)
    if (edge->tris[i] != tri)
      return (edge->tris[i]);

  /* there is no such other triangle if we get here */
  return (NULL);
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
  v - vertex whose face list is to be ordered
******************************************************************************/

void Polyhedron::order_vertex_to_tri_ptrs(Vertex *v)
{
  int i,j;
  Triangle *f;
  Triangle *fnext;
  int nf;
  int vindex;
  int boundary;
  int count;

  nf = v->ntris;
  f = v->tris[0];

  /* go backwards (clockwise) around faces that surround a vertex */
  /* to find out if we reach a boundary */

  boundary = 0;

  for (i = 1; i <= nf; i++) {

    /* find reference to v in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[j] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #1\n");
      exit (-1);
    }

    /* corresponding face is the previous one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* see if we've reached a boundary, and if so then place the */
    /* current face in the first position of the vertice's face list */

    if (fnext == NULL) {
      /* find reference to f in v */
      for (j = 0; j < v->ntris; j++)
        if (v->tris[j] == f) {
	  v->tris[j] = v->tris[0];
	  v->tris[0] = f;
	  break;
	}
      boundary = 1;
      break;
    }

    f = fnext;
  }

  /* now walk around the faces in the forward direction and place */
  /* them in order */

  f = v->tris[0];
  count = 0;

  for (i = 1; i < nf; i++) {

    /* find reference to vertex in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[(j+1) % f->nverts] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #2\n");
      exit (-1);
    }

    /* corresponding face is next one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* break out of loop if we've reached a boundary */
    count = i;
    if (fnext == NULL) {
      break;
    }

    /* swap the next face into its proper place in the face list */
    for (j = 0; j < v->ntris; j++)
      if (v->tris[j] == fnext) {
	v->tris[j] = v->tris[i];
	v->tris[i] = fnext;
	break;
      }

    f = fnext;
  }
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
  f - face whose vertex list is to be searched
  v - vertex to return reference to

Exit:
  returns index in face's list, or -1 if vertex not found
******************************************************************************/

int Polyhedron::face_to_vertex_ref(Triangle *f, Vertex *v)
{
  int j;
  int vindex = -1;

  for (j = 0; j < f->nverts; j++)
    if (f->verts[j] == v) {
      vindex = j;
      break;
    }

  return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/

void Polyhedron::create_pointers()
{
  int i;

  /* index the vertices and triangles */

  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  for (i = 0; i < ntris; i++) 
    tlist[i]->index = i;

  /* create pointers from vertices to triangles */
  vertex_to_tri_ptrs();

  /* make edges */
  create_edges();


  /* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
    order_vertex_to_tri_ptrs(vlist[i]);
		
	}
  /* index the edges */

  for (i = 0; i < nedges; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
    elist[i]->index = i;
	}

}

void Polyhedron::calc_bounding_sphere()
{
  unsigned int i;
  icVector3 min, max;

  for (i=0; i<nverts; i++) {
    if (i==0)  {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
    }
    else {
			if (vlist[i]->x < min.entry[0])
			  min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
			  max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
			  min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
			  max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
			  min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
			  max.entry[2] = vlist[i]->z;
		}
  }
  center = (min + max) * 0.5;
  radius = length(center - min);
}

void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i=0; i<nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1-v2);
	}
}

void Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2;
  Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (i=0; i<ntris; i++){
		for (j=0; j<3; j++)
			length[j] = tlist[i]->edges[j]->length;
		double temp_s = (length[0] + length[1] + length[2])/2.0;
		tlist[i]->area = sqrt(temp_s*(temp_s-length[0])*(temp_s-length[1])*(temp_s-length[2]));

		area += tlist[i]->area;
		temp_t = tlist[i];
		v1.set(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		v2.set(vlist[tlist[i]->verts[1]->index]->x, vlist[tlist[i]->verts[1]->index]->y, vlist[tlist[i]->verts[1]->index]->z);
		v0.set(vlist[tlist[i]->verts[2]->index]->x, vlist[tlist[i]->verts[2]->index]->y, vlist[tlist[i]->verts[2]->index]->z);
		tlist[i]->normal = cross(v0-v1, v2-v1);
		normalize(tlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i=0; i<ntris; i++){
		icVector3 cent(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		signedvolume += dot(test-cent, tlist[i]->normal)*tlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume<0) 
		orientation = 0;
	else {
		orientation = 1;
		for (i=0; i<ntris; i++)
			tlist[i]->normal *= -1.0;
	}
}

void sort(unsigned int *A, unsigned int *B, unsigned int *C, unsigned int sid, unsigned int eid){
  unsigned int i;
	unsigned int *tempA, *tempB, *tempC;
	unsigned int current1, current2, current0;

  if (sid>=eid)
		return;
	sort(A, B, C, sid, (sid+eid)/2);
	sort(A, B, C, (sid+eid)/2+1, eid);
	tempA = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempB = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempC = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	for (i=0; i<eid-sid+1; i++){
		tempA[i] = A[i+sid];
		tempB[i] = B[i+sid];
		tempC[i] = C[i+sid];
	}
	current1 = sid;
	current2 = (sid+eid)/2+1;
	current0 = sid;
	while ((current1<=(sid+eid)/2) && (current2<=eid)){
		if (tempA[current1-sid] < tempA[current2-sid]) {
			A[current0] = tempA[current1-sid];
			B[current0] = tempB[current1-sid];
			C[current0] = tempC[current1-sid];
			current1++;		
		}
		else if (tempA[current1-sid] > tempA[current2-sid]){
			A[current0] = tempA[current2-sid];
			B[current0] = tempB[current2-sid];
			C[current0] = tempC[current2-sid];
			current2++;		
		}
		else {
			if (tempB[current1-sid] < tempB[current2-sid]) {
				A[current0] = tempA[current1-sid];
				B[current0] = tempB[current1-sid];
				C[current0] = tempC[current1-sid];
				current1++;		
			} else {
				A[current0] = tempA[current2-sid];
				B[current0] = tempB[current2-sid];
				C[current0] = tempC[current2-sid];
				current2++;		
			}
		}
		current0++;
	}
	if (current1<=(sid+eid)/2){
		for (i=current1; i<=(sid+eid)/2; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}
	if (current2<=eid){
		for (i=current2; i<=eid; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}

	free(tempA);
	free(tempB);
	free(tempC);
}

void init(void) {
  /* select clearing color */ 

  glClearColor (0.0, 0.0, 0.0, 0.0);  // background
  glShadeModel (GL_FLAT);
  //glPolygonMode(GL_FRONT, GL_FILL);
  glPolygonMode(GL_FRONT, GL_LINE);

  glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
	// may need it
  glPixelStorei(GL_PACK_ALIGNMENT,1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0) 
		glFrontFace(GL_CW);
	else 
		glFrontFace(GL_CCW);
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

  /* set escape key to exit */
  switch (key) {
    case 27:
			poly->finalize();  // finalize_everything
      exit(0);
      break;

		case '0':
			display_mode = 0;
			display();
			break;

		case '1':
			display_mode = 0;
			display();
			break;

		case '2':
			display_mode = 0;
			display();
			break;

		case '3':
			display_mode = 3;
			display();
			break;

		case '4':
			display_mode = 4;
			display();
			break;

		case '5':
			display_mode = 5;
			display();
			break;

		case '6':
			display_mode = 6;
			display();
			break;

		case '7':
			display_mode = 7;
			display();
			break;

		case '8':
			display_mode = 8;
			display();
			break;

		case '9':
			display_mode = 9;
			display();
			break;

		case 'x':
			switch(ACSIZE){
			case 1:
				ACSIZE = 16;
				break;

			case 16:
				ACSIZE = 1;
				break;

			default:
				ACSIZE = 1;
				break;
			}
			fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
			display();
			break;

		case '|':
			this_file = fopen("rotmat.txt", "w");
			for (i=0; i<4; i++) 
				fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
			fclose(this_file);
			break;

		case '^':
			this_file = fopen("rotmat.txt", "r");
			for (i=0; i<4; i++) 
				fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
			fclose(this_file);
			display();
			break;

	}
}

Polyhedron::Polyhedron()
{
	nverts = nedges = ntris = 0;
	max_verts = max_tris = 50;

	vlist = new Vertex *[max_verts];
	tlist = new Triangle *[max_tris];		
}


void multmatrix(const Matrix m)
{ 
  int i,j, index = 0;

  GLfloat mat[16];

  for ( i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      mat[index++] = m[i][j];

  glMultMatrixf (mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


  glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix( rotmat );

	glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

void motion(int x, int y) {
	float r[4];
	float xsize, ysize, s, t;

	switch(mouse_mode){
	case -1:

		xsize = (float) win_width;
		ysize = (float) win_height;
	
		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat( rotmat, rvec );
		trackball( r, s_old, t_old, s, t );
		add_quats( r, rvec, rvec );
		quat_to_mat( rvec, rotmat );

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth=1.0e+20, current_depth;
	int seed_id=-1; 
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint *) buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;
		
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	printf("triangle id = %d\n", seed_id);

	if(hitNum<3)
		triSel[hitNum++] = seed_id;
	return seed_id;
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		switch(mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float) win_width;
				float ysize = (float) win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1;  // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	} else if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLuint selectBuf[win_width];
		  GLint hits;
		  GLint viewport[4];

		  glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
		  (void) glRenderMode(GL_SELECT);

		  glInitNames();
		  glPushName(0);

		  glMatrixMode(GL_PROJECTION);
	    glPushMatrix();
			glLoadIdentity();
/*  create 5x5 pixel picking region near cursor location */
	    gluPickMatrix((GLdouble) x, (GLdouble) (viewport[3] - y),
                 1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix ();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
	    glPopMatrix();
		  glFlush();

	    hits = glRenderMode(GL_RENDER);
		  poly->seed = processHits(hits, selectBuf);
		  //heatDiffusion(5, 0, 1.0, poly);
		  //visHeat(poly);
			display();
		}
	}
}

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i=0; i<poly->ntris; i++) {
		Triangle *temp_t=poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};
		
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
   
		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (j=0; j<3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

  glEnable (GL_POLYGON_OFFSET_FILL);
  glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	for (i=0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
      glLoadName(i+1);

		Triangle *temp_t=this_poly->tlist[i];

		switch (display_mode) {
		case 0:
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			} else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i==this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 6:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 10:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
		
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}
	}
	glDisable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
}

bool heatDone = false;

void display(void)
{
  GLint viewport[4];
  int jitter;

  glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
  glGetIntegerv (GL_VIEWPORT, viewport);
 
  glClear(GL_ACCUM_BUFFER_BIT);
  for (jitter = 0; jitter < ACSIZE; jitter++) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  set_view(GL_RENDER, poly);
    glPushMatrix ();
		switch(ACSIZE){
		case 1:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;

		case 16:
			glTranslatef (ji16[jitter].x*2.0/viewport[2], ji16[jitter].y*2.0/viewport[3], 0.0);
			break;

		default:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;
		}
		//gluLookAt(viewPoint.x, viewPoint.y, viewPoint.z, poly->center.x, poly->center.y, poly->center.z, 0, 1, 0);
		set_scene(GL_RENDER, poly);
		//if(hitNum==2)
		//{
		//	int maxId = poly->tlist[triSel[0]]->verts[0]->index;
		//	int minId = poly->tlist[triSel[1]]->verts[1]->index;
		//	if(!heatDone)
		//	{
		//		heatDiffusion(maxId,minId,poly);
		//		heatDone = !heatDone;
		//		identifyMMS(poly);
		//	}
		//	visHeat( poly);
		//}

		//if(hitNum==3)
		//{
		//	int maxId = poly->tlist[triSel[0]]->verts[0]->index;
		//	int minId = poly->tlist[triSel[1]]->verts[0]->index;
		//	int maxId2 = poly->tlist[triSel[2]]->verts[0]->index;
		//	if(!heatDone)
		//	{
		//		textureSyn(maxId,minId,maxId2,poly);
		//		heatDone = !heatDone;
		//		//identifyMMS(poly);
		//	}
		//	drawTexture(poly);
		//}
		//else
		//display_shape(GL_RENDER, poly);
		
		switch(operationNum)
		{
		case 0:
			display_shape(GL_RENDER, poly);
			break;
		case 1:
			if(checkerBrd)
				checkerboard(poly,paraL);
			else
				display_shape(GL_RENDER, poly);
			break;
		case 2:
			//gluLookAt(viewPoint.x, viewPoint.y, viewPoint.z-19.0, poly->center.x, poly->center.y, poly->center.z, 0, 1, 0);
			if(hitNum==2)
			{
				int maxId = poly->tlist[triSel[0]]->verts[0]->index;
				int minId = poly->tlist[triSel[1]]->verts[1]->index;
				if(!heatDone)
				{
					heatDiffusion(maxId,minId,poly);
					heatDone = !heatDone;
					identifyMMS(poly);
				}
				visHeat( poly);
			}
			else
				display_shape(GL_RENDER, poly);
			break;
		case 3:
		if(hitNum==3)
			{
				int maxId = poly->tlist[triSel[0]]->verts[0]->index;
				int minId = poly->tlist[triSel[1]]->verts[0]->index;
				int maxId2 = poly->tlist[triSel[2]]->verts[0]->index;
				if(!heatDone)
				{
					textureSyn(maxId,minId,maxId2,poly);
					heatDone = !heatDone;
					//identifyMMS(poly);
				}
				drawTexture(poly);
			}
		else
			display_shape(GL_RENDER, poly);
			break;
		case 4:
			checkerboard(poly,paraL);
			break;
		case 5:
			visCurv(poly,1);
			break;
		case 6:
			visCurv(poly,2);
			break;
		case 7:
			drawCurvTensor(poly);
			display_shape(GL_RENDER, poly);
			break;
		case 8:
			gluLookAt(viewPoint.x, viewPoint.y, viewPoint.z, poly->center.x, poly->center.y, poly->center.z, 0, 0,1);
			display_shape(GL_RENDER, poly);
			break;
		}
    glPopMatrix ();
    glAccum(GL_ACCUM, 1.0/ACSIZE);
  }
  glAccum (GL_RETURN, 1.0);
  glFlush();
  glutSwapBuffers();
 	glFinish();
}

void Polyhedron::average_normals()
{
	int i, j;

	for (i=0; i<nverts; i++) {
		vlist[i]->normal = icVector3(0.0);
		for (j=0; j<vlist[i]->ntris; j++) 
			vlist[i]->normal += vlist[i]->tris[j]->normal;
		normalize(vlist[i]->normal);
	}
}
/************************************************************************/

void Polyhedron::createCList()
{
	  ncorners = max_corners = 3*ntris;
      clist = new Corner *[ncorners];

	  for(int i=0;i<ntris;++i)
		  for(int j=0;j<3;++j)
		  {
		  clist[i*3+j] = new Corner(tlist[i],j);
		  clist[i*3+j]->index = i*3+j;
		  }

	  for(int i=0;i<ntris;++i)
		  for(int j=0;j<3;++j)
		  {

			clist[i*3+j]->n = clist[i*3 + (j+1)%3];
			clist[i*3+j]->p = clist[i*3 + (j+2)%3];
		  }

}

void Polyhedron::findOppositeCorners()
{
	unsigned int *corners, *vlist1, *vlist2, prev_v1, prev_v2;
	int *oppo;
	int i;

	oppo = (int *)malloc(sizeof(int)*ntris*3);
	corners = (unsigned int *)malloc(sizeof(unsigned int)*ntris*3);
	vlist1 = (unsigned int *)malloc(sizeof(unsigned int)*ntris*3);
	vlist2 = (unsigned int *)malloc(sizeof(unsigned int)*ntris*3);

	for (i=0; i<ntris; i++) {
		Triangle *temp_t = tlist[i];
		corners[3*i] = 3*i;
		corners[3*i+1] = 3*i+1;
		corners[3*i+2] = 3*i+2;
		if (temp_t->verts[1]->index > temp_t->verts[2]->index){
			vlist1[3*i] = temp_t->verts[1]->index;
			vlist2[3*i] = temp_t->verts[2]->index;
		} else {
			vlist1[3*i] = temp_t->verts[2]->index;
			vlist2[3*i] = temp_t->verts[1]->index;
		}
		if (temp_t->verts[2]->index > temp_t->verts[0]->index){
			vlist1[3*i+1] = temp_t->verts[2]->index;
			vlist2[3*i+1] = temp_t->verts[0]->index;
		} else {
			vlist1[3*i+1] = temp_t->verts[0]->index;
			vlist2[3*i+1] = temp_t->verts[2]->index;
		}
		if (temp_t->verts[0]->index > temp_t->verts[1]->index){
			vlist1[3*i+2] = temp_t->verts[0]->index;
			vlist2[3*i+2] = temp_t->verts[1]->index;
		} else {
			vlist1[3*i+2] = temp_t->verts[1]->index;
			vlist2[3*i+2] = temp_t->verts[0]->index;
		}
	}


	sort(vlist1, vlist2, corners, 0, ntris*3-1);

	prev_v1 = vlist1[0];  
	prev_v2 = vlist2[0];
	for (i=0; i<ntris*3; i++)
		oppo[i] = -1;
	for (i=1; i<ntris*3; i++){
		if ((prev_v1 == vlist1[i]) && (prev_v2 == vlist2[i])){
			oppo[corners[i]] = corners[i-1];
			oppo[corners[i-1]] = corners[i];
			++i;
		} 
		prev_v1 = vlist1[i];
		prev_v2 = vlist2[i];
	}
	free(corners);
	free(vlist1);
	free(vlist2);
	for (i=0; i<ncorners; i++) {
		if (oppo[i] >= 0)
			clist[i]->o = clist[oppo[i]];
		else
			clist[i]->o = NULL;
	}
	free(oppo);
}

void Polyhedron::markFBFace(icVector3& viewPoint)
{
	icVector3 ray = viewPoint;
	icVector3 triCenter(0,0,0);
	
	for(int i=0;i<ntris;++i)
	{
		triCenter += icVector3(tlist[i]->verts[0]->x,tlist[i]->verts[0]->y,tlist[i]->verts[0]->z);
		triCenter += icVector3(tlist[i]->verts[1]->x,tlist[i]->verts[1]->y,tlist[i]->verts[1]->z);
		triCenter += icVector3(tlist[i]->verts[2]->x,tlist[2]->verts[0]->y,tlist[i]->verts[2]->z);
		triCenter /=3.0;
		ray -= triCenter;

		if(dot(ray,tlist[i]->normal)>=0)
			tlist[i]->front = true;
		else 
			tlist[i]->front = false;
	}
}

void Polyhedron::markSliEdge()
{
	for(int i=0;i<nedges;++i)
	{
		if(elist[i]->ntris==2)
		{
			if(elist[i]->tris[0]->front!=elist[i]->tris[1]->front)
				elist[i]->sliEdge = true;
			else
				elist[i]->sliEdge = false;
		}
		else 
			elist[i]->sliEdge = false;
	}
}

void drawSilFace(Polyhedron *poly)
	//draw silhouette based on the normal of face
{ 
    glDisable(GL_TEXTURE_2D); 
    glEnable(GL_BLEND); 
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
    glColor3f( 0.0, 0.0, 0.0); 

    glLineWidth(1.0); 
	
	for(int i=0;i<poly->nedges;++i)
	{
		if(poly->elist[i]->sliEdge)
		{
			glBegin(GL_LINES); 
			glVertex3d(poly->elist[i]->verts[0]->x, poly->elist[i]->verts[0]->y, poly->elist[i]->verts[0]->z); 
			glVertex3d(poly->elist[i]->verts[1]->x, poly->elist[i]->verts[1]->y, poly->elist[i]->verts[1]->z); 
			glEnd(); 
		}
	}

    glDisable(GL_BLEND); 
    glEnable(GL_TEXTURE_2D); 
}  

void Polyhedron::findSilPointOnEdge(icVector3& viewPoint)
{
	double dotPro1, dotPro2;
	icVector3 ray(0.0,0.0,0.0);

	for(int i=0;i<ntris;++i)
	{
		tlist[i]->containSliP = false;
		tlist[i]->silPointIndex = 0;
	}

	for(int i=0;i<nedges;++i)
	{
		ray = viewPoint - icVector3(elist[i]->verts[0]->x,elist[i]->verts[0]->y,elist[i]->verts[0]->z);
		dotPro1 = dot(elist[i]->verts[0]->normal,ray);
		ray = viewPoint - icVector3(elist[i]->verts[1]->x,elist[i]->verts[1]->y,elist[i]->verts[1]->z);
		dotPro2 = dot(elist[i]->verts[1]->normal,ray);
		//icVector3 silPoint;
		if(dotPro1*dotPro2<=0)
		{
			double ratio = dotPro2/(dotPro2-dotPro1);
			double x = elist[i]->verts[0]->x*(1.0-ratio) + elist[i]->verts[1]->x*ratio;
			double y = elist[i]->verts[0]->y*(1.0-ratio) + elist[i]->verts[1]->y*ratio;
			double z = elist[i]->verts[0]->z*(1.0-ratio) + elist[i]->verts[1]->z*ratio;
			for(int j=0;j<elist[i]->ntris;++j)
			{
				Triangle *tempTri = elist[i]->tris[j];
				tempTri->containSliP = true;
				tempTri->silPoints[tempTri->silPointIndex++].set(x,y,z);
			}
		}
	}
}

void drawSilVert(Polyhedron *poly)
{ 
    glDisable(GL_TEXTURE_2D); 
    glEnable(GL_BLEND); 
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
    glColor3f( 0.0, 0.0, 0.0); 

    glLineWidth(1.0); 
	
	for(int i=0;i<poly->ntris;++i)
	{
		if(poly->tlist[i]->containSliP)
		{
			glBegin(GL_LINES); 
			glVertex3d(poly->tlist[i]->silPoints[0].x, poly->tlist[i]->silPoints[0].y, poly->tlist[i]->silPoints[0].z); 
			glVertex3d(poly->tlist[i]->silPoints[1].x, poly->tlist[i]->silPoints[1].y, poly->tlist[i]->silPoints[1].z);
			glEnd(); 
		}
	}

    glDisable(GL_BLEND); 
    glEnable(GL_TEXTURE_2D); 
} 

void Polyhedron::loopDivide()
{
	for(int i=0;i<nedges;++i)
	{
		double x = 0.5*(elist[i]->verts[0]->x + elist[i]->verts[1]->x);
		double y = 0.5*(elist[i]->verts[0]->y + elist[i]->verts[1]->y);
		double z = 0.5*(elist[i]->verts[0]->z + elist[i]->verts[1]->z);
		vlist[nverts] = new Vertex(x,y,z);
		elist[i]->newVert = vlist[nverts];
		//vlist[nverts]->index = nverts;
		++nverts;
	}

	int OrigNTris = ntris;
	for(int i=0;i<OrigNTris;++i)
	{

		tlist[ntris] = new Triangle;
        tlist[ntris]->nverts = 3;
		tlist[ntris]->verts[0] = tlist[i]->edges[2]->newVert;
		tlist[ntris]->verts[1] = tlist[i]->verts[0];
        tlist[ntris++]->verts[2] = tlist[i]->edges[0]->newVert;
		

		tlist[ntris] = new Triangle;
        tlist[ntris]->nverts = 3;
		tlist[ntris]->verts[0] = tlist[i]->edges[0]->newVert;
		tlist[ntris]->verts[1] = tlist[i]->verts[1];
        tlist[ntris++]->verts[2] = tlist[i]->edges[1]->newVert;

		tlist[ntris] = new Triangle;
        tlist[ntris]->nverts = 3;
		tlist[ntris]->verts[0] = tlist[i]->edges[1]->newVert;
		tlist[ntris]->verts[1] = tlist[i]->verts[2];
        tlist[ntris++]->verts[2] = tlist[i]->edges[2]->newVert;

		tlist[i]->verts[0] = tlist[i]->edges[0]->newVert;
		tlist[i]->verts[1] = tlist[i]->edges[1]->newVert;
		tlist[i]->verts[2] = tlist[i]->edges[2]->newVert;
	}
}

void Polyhedron::subDivEdgeL(double threshL)
{
	for(int i=0;i<nedges;++i)
	{
		elist[i]->newVert = NULL;
		if(elist[i]->length > threshL)
		{
		double x = 0.5*(elist[i]->verts[0]->x + elist[i]->verts[1]->x);
		double y = 0.5*(elist[i]->verts[0]->y + elist[i]->verts[1]->y);
		double z = 0.5*(elist[i]->verts[0]->z + elist[i]->verts[1]->z);
		vlist[nverts] = new Vertex(x,y,z);
		elist[i]->newVert = vlist[nverts];
		//vlist[nverts]->index = nverts;
		++nverts;
		}
	}

	int OrigNTris = ntris;
	for(int i=0;i<OrigNTris;++i)
	{
		int divEIndex[3];
		int divENum = 0;
		for(int j=0;j<3;++j)
		{
			if(tlist[i]->edges[j]->newVert != NULL)
				divEIndex[divENum++] = j;
		}

		switch(divENum)
		{
		case 1:
			//int eNum = divEIndex[0];
			tlist[ntris] = new Triangle;
			tlist[ntris]->nverts = 3;
			tlist[ntris]->verts[0] = tlist[i]->verts[(divEIndex[0]+2)%3];
			tlist[ntris]->verts[1] = tlist[i]->verts[((divEIndex[0]+2)%3+1)%3];
			tlist[ntris++]->verts[2] = tlist[i]->edges[divEIndex[0]]->newVert;

			tlist[i]->verts[0] = tlist[i]->edges[divEIndex[0]]->newVert;
			tlist[i]->verts[1] = tlist[i]->verts[(divEIndex[0]+1)%3];
			tlist[i]->verts[2] = tlist[i]->verts[(divEIndex[0]+2)%3];
			break;

		case 2:
			//int eNum0 = divEIndex[0];
			//int eNum1 = divEIndex[1];

			if(divEIndex[0]==0 && divEIndex[1]==1)
			{
				tlist[ntris] = new Triangle;
				tlist[ntris]->nverts = 3;
				tlist[ntris]->verts[0] = tlist[i]->edges[0]->newVert;
				tlist[ntris]->verts[1] = tlist[i]->verts[1];
				tlist[ntris++]->verts[2] = tlist[i]->edges[1]->newVert;

				tlist[ntris] = new Triangle;
				tlist[ntris]->nverts = 3;
				tlist[ntris]->verts[0] = tlist[i]->edges[0]->newVert;
				tlist[ntris]->verts[1] = tlist[i]->edges[1]->newVert;
				tlist[ntris++]->verts[2] = tlist[i]->verts[2];

				tlist[i]->verts[0] = tlist[i]->edges[0]->newVert;
				tlist[i]->verts[1] = tlist[i]->verts[2];
				tlist[i]->verts[2] = tlist[i]->verts[0];
			}

			if(divEIndex[0]==1 && divEIndex[1]==2)
			{
				tlist[ntris] = new Triangle;
				tlist[ntris]->nverts = 3;
				tlist[ntris]->verts[0] = tlist[i]->edges[1]->newVert;
				tlist[ntris]->verts[1] = tlist[i]->verts[2];
				tlist[ntris++]->verts[2] = tlist[i]->edges[2]->newVert;

				tlist[ntris] = new Triangle;
				tlist[ntris]->nverts = 3;
				tlist[ntris]->verts[0] = tlist[i]->edges[1]->newVert;
				tlist[ntris]->verts[1] = tlist[i]->edges[2]->newVert;
				tlist[ntris++]->verts[2] = tlist[i]->verts[0];

				tlist[i]->verts[0] = tlist[i]->edges[1]->newVert;
				tlist[i]->verts[1] = tlist[i]->verts[0];
				tlist[i]->verts[2] = tlist[i]->verts[1];
			}

			if(divEIndex[0]==0 && divEIndex[1]==2)
			{
				tlist[ntris] = new Triangle;
				tlist[ntris]->nverts = 3;
				tlist[ntris]->verts[0] = tlist[i]->edges[0]->newVert;
				tlist[ntris]->verts[1] = tlist[i]->verts[2];
				tlist[ntris++]->verts[2] = tlist[i]->edges[2]->newVert;

				tlist[ntris] = new Triangle;
				tlist[ntris]->nverts = 3;
				tlist[ntris]->verts[0] = tlist[i]->edges[0]->newVert;
				tlist[ntris]->verts[1] = tlist[i]->edges[2]->newVert;
				tlist[ntris++]->verts[2] = tlist[i]->verts[0];

				tlist[i]->verts[0] = tlist[i]->edges[0]->newVert;
				tlist[i]->verts[1] = tlist[i]->verts[1];
				tlist[i]->verts[2] = tlist[i]->verts[2];
			}
			break;

		case 3:
			tlist[ntris] = new Triangle;
			tlist[ntris]->nverts = 3;
			tlist[ntris]->verts[0] = tlist[i]->edges[2]->newVert;
			tlist[ntris]->verts[1] = tlist[i]->verts[0];
			tlist[ntris++]->verts[2] = tlist[i]->edges[0]->newVert;
		

			tlist[ntris] = new Triangle;
			tlist[ntris]->nverts = 3;
			tlist[ntris]->verts[0] = tlist[i]->edges[0]->newVert;
			tlist[ntris]->verts[1] = tlist[i]->verts[1];
			tlist[ntris++]->verts[2] = tlist[i]->edges[1]->newVert;

			tlist[ntris] = new Triangle;
			tlist[ntris]->nverts = 3;
			tlist[ntris]->verts[0] = tlist[i]->edges[1]->newVert;
			tlist[ntris]->verts[1] = tlist[i]->verts[2];
			tlist[ntris++]->verts[2] = tlist[i]->edges[2]->newVert;

			tlist[i]->verts[0] = tlist[i]->edges[0]->newVert;
			tlist[i]->verts[1] = tlist[i]->edges[1]->newVert;
			tlist[i]->verts[2] = tlist[i]->edges[2]->newVert;
			break;
		}
	}
}

void checkerboard(Polyhedron *this_poly,float L=1.0)
{
	GLfloat mat_diffuse[4];

    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	Triangle *tri;
	//unsigned char x, y, z;
	//Vertex *verts;
	Vertex *temp_v;
	//ofstream checkerboard("../checkerboard.txt");

	for(unsigned int i=0;i<this_poly->ntris;++i)
	{
		tri = this_poly->tlist[i];
		int r,g,b;

		glBegin(GL_POLYGON);
		for(int j=0;j<tri->nverts;++j)
		{
		temp_v = tri->verts[j];
		r = ((int)floor(temp_v->x/L)+1)%2;
		g = ((int)floor(temp_v->y/L)+1)%2;
	    b = ((int)floor(temp_v->z/L)+1)%2;
		mat_diffuse[0] = r+0.0;
		mat_diffuse[1] = g+0.0;
		mat_diffuse[2] = b+0.0;
		mat_diffuse[3] = 1.0;

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glColor3i(r,g,b);
		//glColor3ub(255-x,255-y,255-z);
		//checkerboard<<"i: "<<i<<" r: "<<r<<" g: "<<g<<" b: "<<b<<endl;
		glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
		glVertex3d(temp_v->x,temp_v->y,temp_v->z);
		}
		glEnd();

	}
    //checkerboard.close();
	return;
}

/************************************************************************/
//Proj 3

void smoothWgtUpdated(int scheme, double dt, int times, Polyhedron* poly)
{
	if(scheme == 0)//uniform scheme 
	{
		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				updateX += tempN->x - tempV->x;
				updateY += tempN->y - tempV->y;
				updateZ += tempN->z - tempV->z;
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double weight = dt*(1.0/tempV->cList.size());
			tempV->x += weight*tempUpdateVList[i*3+0];
			tempV->y += weight*tempUpdateVList[i*3+1];
			tempV->z += weight*tempUpdateVList[i*3+2];
		}

		}

	}
	else if(scheme == 1)//cord scheme
	{
		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;

			for(int j=0;j<tempV->cList.size();++j)
				totalE += 1.0/tempV->cList[j]->p->e->length;

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}
		}

	}
	else if(scheme == 2)//mean curvature flow
	{
		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;
			vector<double> eValue;
			for(int j=0;j<tempV->cList.size();++j)
			{
				double tempE = (tempV->cList[j]->p->cot
							   + tempV->cList[j]->p->o->cot)/2.0;
				totalE += tempE;
				eValue.push_back(tempE);
			}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = eValue[j]/totalE;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}

		}
	}
	else if(scheme == 3)//mean values coordinates
	{
		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;
			vector<double> eValue;
			for(int j=0;j<tempV->cList.size();++j)
			{
				double tempE = (tan(tempV->cList[j]->angle/2.0)
							   + tan(tempV->cList[j]->p->o->p->angle/2.0))/2.0;
				totalE += tempE;
				eValue.push_back(tempE);
			}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = eValue[j]/totalE;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}
		}
	}

	return;
}

void smoothWgtNotUpdated(int scheme, double dt, int times, Polyhedron* poly)
{
	if(scheme == 1)//cord scheme
	{
		vector<double> weightList;

		if(times>0)
		{
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;

			for(int j=0;j<tempV->cList.size();++j)
				totalE += 1.0/tempV->cList[j]->p->e->length;

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
				weightList.push_back(weight);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}
		}

		for(int t=1;t<times;++t)
		{

		vector<double>::iterator it = weightList.begin();
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;

			//for(int j=0;j<tempV->cList.size();++j)
			//	totalE += 1.0/tempV->cList[j]->p->e->length;

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
				double weight = *it;
				++it;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}
		}

	}
	else if(scheme == 2)//mean curvature flow
	{
		vector<double> weightList;
		
		if(times>0)
		{
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;
			vector<double> eValue;
			for(int j=0;j<tempV->cList.size();++j)
			{
				double tempE = (tempV->cList[j]->p->cot
							   + tempV->cList[j]->p->o->cot)/2.0;
				totalE += tempE;
				eValue.push_back(tempE);
			}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = eValue[j]/totalE;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
				weightList.push_back(weight);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}

		}
		for(int t=1;t<times;++t)
		{
		vector<double>::iterator it = weightList.begin();
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;
			//vector<double> eValue;
			//for(int j=0;j<tempV->cList.size();++j)
			//{
			//	double tempE = (tempV->cList[j]->p->cos/tempV->cList[j]->p->sin
			//				   + tempV->cList[j]->p->o->cos/tempV->cList[j]->p->o->sin)/2.0;
			//	totalE += tempE;
			//	eValue.push_back(tempE);
			//}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = eValue[j]/totalE;
				double weight = *it;
				++it;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}

		}
	}
	else if(scheme == 3)//mean values coordinates
	{
		vector<double> weightList;
		
		if(times>0)
		{
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;
			vector<double> eValue;
			for(int j=0;j<tempV->cList.size();++j)
			{
				double tempE = (tan(tempV->cList[j]->angle/2.0)
							   + tan(tempV->cList[j]->p->o->p->angle/2.0))/2.0;
				totalE += tempE;
				eValue.push_back(tempE);
			}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = eValue[j]/totalE;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
				weightList.push_back(weight);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}

		}
		for(int t=1;t<times;++t)
		{
		vector<double>::iterator it = weightList.begin();
		vector<double> tempUpdateVList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double updateX=0.0;
			double updateY=0.0;
			double updateZ=0.0;
			double totalE = 0.0;
			//vector<double> eValue;
			//for(int j=0;j<tempV->cList.size();++j)
			//{
			//	double tempE = (tempV->cList[j]->p->cos/tempV->cList[j]->p->sin
			//				   + tempV->cList[j]->p->o->cos/tempV->cList[j]->p->o->sin)/2.0;
			//	totalE += tempE;
			//	eValue.push_back(tempE);
			//}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = eValue[j]/totalE;
				double weight = *it;
				++it;
				updateX += weight*(tempN->x - tempV->x);
				updateY += weight*(tempN->y - tempV->y);
				updateZ += weight*(tempN->z - tempV->z);
			}
			tempUpdateVList.push_back(updateX);
			tempUpdateVList.push_back(updateY);
			tempUpdateVList.push_back(updateZ);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->x += dt*tempUpdateVList[i*3+0];
			poly->vlist[i]->y += dt*tempUpdateVList[i*3+1];
			poly->vlist[i]->z += dt*tempUpdateVList[i*3+2];
		}

		}
		}

	return;
}

bool myCompare(Corner* c1, Corner* c2)
{
	return(c1->n->v->index < c2->n->v->index);
}

void heatDiffusion(int maxVertId, int minVertId, Polyhedron* poly)
{
	vector<int> vIja;
	vector<double> vSa;
	//double *heatDiff = (double*)malloc((poly->nverts)*sizeof(double));
	//Vec_DP heatDiff(poly->nverts);
	 Vec_DP b(poly->nverts-2),bcmp(poly->nverts-2),x(poly->nverts-2);
	 const int ITOL=1,ITMAX=75;
     const DP TOL=1.0e-9;
     int iter;
     DP err;

   for (int i=0;i<poly->nverts-2;i++) {
		x[i]=0.0;
		b[i]=0.0;
	}

	for (int j=0;j<=(poly->nverts-2);j++)
		{
			vSa.push_back(1.0);
			vIja.push_back(-1.0);
	}
	vIja[0]=poly->nverts-1;
	int k = poly->nverts-2;
	int jaIndex = 0;

	int maxIndex,minIndex;
	if(maxVertId>minVertId)
	{
		maxIndex = maxVertId;
		minIndex = minVertId;
	}
	else
	{
		maxIndex = minVertId;
		minIndex = maxVertId;
	}

	for(int i=0;i<poly->nverts;++i)
	{
		if(i==maxVertId || i==minVertId)
		//{
			//vIja[i+1]=k+1;
			//--i;
			continue;
		//}
		Vertex* tempV = poly->vlist[i];
		double totalE = 0.0;

		sort(tempV->cList.begin(),tempV->cList.end(),myCompare);

		for(int j=0;j<tempV->cList.size();++j)
			totalE += 1.0/tempV->cList[j]->p->e->length;

		for(int j=0;j<tempV->cList.size();++j)
		{
			Vertex* tempN = tempV->cList[j]->n->v;
			double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
			if(tempN->index==maxVertId || tempN->index==minVertId)
			{
				if(tempN->index==maxVertId)
				{
					int index = i;
					if(index>maxIndex)
						index -= 2;
					else if(index>minIndex)
						--index;
					b[index] = weight;
				}
			}
			else
			{
			vSa.push_back(weight*(-1.0));
			int index = tempN->index;
			if(index>maxIndex)
				index -= 2;
			else if(index>minIndex)
				--index;
			vIja.push_back(index);
			++k;
			}
		}

		vIja[++jaIndex]=k+1;
	}

	 ija_p=new Vec_INT(vIja);
     sa_p=new Vec_DP(vSa);

	// Vec_DP b(poly->nverts),bcmp(poly->nverts),x(poly->nverts);
	// const int ITOL=1,ITMAX=75;
 //    const DP TOL=1.0e-9;
 //    int iter;
 //    DP err;

 //  for (int i=0;i<poly->nverts;i++) {
	//	x[i]=0.0;
	//	b[i]=1.0e-15;
	//}
	//b[maxVertId]=1.0;
	//b[minVertId] = 0.0;
	NR::linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	/* cout << scientific << setprecision(6);
        cout << "Estimated error: " << setw(16) << err << endl;
        cout << "Iterations needed: " << setw(6) << iter;
        cout << endl << "Solution vector:" << endl;
        for (int ii=0;ii<NP/5;ii++) {
          for (i=5*ii;i<5*(ii+1);i++) cout << setw(15) << x[i];
          cout << endl;
        }
        for (i=0;i<(NP % 5);i++)
          cout << setw(12) << x[5*(NP/5)+i];
        cout << endl;*/
        //NR::sprsax(*sa_p,*ija_p,x,bcmp);
        //cout << endl << "press RETURN to continue..." << endl;
        //cin.get();
        //cout << "test of solution vector:" << endl;
        //cout << setw(9) << "a*x" << setw(13) << "b" << endl;
        //for (int i=0;i<poly->nverts;i++)
        //  cout << setw(14) << bcmp[i] << setw(15) << b[i] << endl;

	for(int i=0;i<poly->nverts;++i)
		{
			if(i==maxVertId || i==minVertId)
			//{
				//vIja[i+1]=k+1;
				//--i;
				continue;
			//}
			int index = i;
			if(index>maxIndex)
				index -= 2;
			else if(index>minIndex)
				--index;
			poly->vlist[i]->heat = x[index];
	        //cout<<x[index]<<endl;
	}

	 poly->vlist[maxVertId]->heat = 1.0;
	 poly->vlist[minVertId]->heat = 0.0;

	 delete sa_p;
     delete ija_p;

     return;
}

void visHeat(Polyhedron* poly)
{
	GLfloat mat_diffuse[4];

    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	Triangle *tri;
	//unsigned char x, y, z;
	//Vertex *verts;
	Vertex *temp_v;
	//ofstream normMap("../normMap.txt");

	for(unsigned int i=0;i<poly->ntris;++i)
	{
		tri = poly->tlist[i];
		//float r = abs(tri->normal.entry[0]);
		//float g = abs(tri->normal.entry[1]);
	 //   float b = abs(tri->normal.entry[2]);
		//mat_diffuse[0] = r;
		//mat_diffuse[1] = g;
		//mat_diffuse[2] = b;
		//mat_diffuse[3] = 1.0;

		//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		//glColor3f(r,g,b);
		//glColor3ub(255-x,255-y,255-z);
		//normMap<<"i: "<<i<<" r: "<<r<<" g: "<<g<<" b: "<<b<<endl;

		glBegin(GL_POLYGON);
		for(int j=0;j<tri->nverts;++j)
		{
		temp_v = tri->verts[j];

		//float r = temp_v->heat;
		//float g = 1-r;
	 //   float b = 1-r;

		float t = temp_v->heat;

	float r = 1.;
	float g = 0.0;
	float b = 1.  -  6. * ( t - (5./6.) );

        if( t <= (5./6.) )
        {
                r = 6. * ( t - (4./6.) );
                g = 0.;
                b = 1.;
        }

        if( t <= (4./6.) )
        {
                r = 0.;
                g = 1.  -  6. * ( t - (3./6.) );
                b = 1.;
        }

        if( t <= (3./6.) )
        {
                r = 0.;
                g = 1.;
                b = 6. * ( t - (2./6.) );
        }

        if( t <= (2./6.) )
        {
                r = 1.  -  6. * ( t - (1./6.) );
                g = 1.;
                b = 0.;
        }

        if( t <= (1./6.) )
        {
                r = 1.;
				b = 1.;
                g = 1.-6. * t;
        }

		mat_diffuse[0] = r;
		mat_diffuse[1] = g;
		mat_diffuse[2] = b;
		mat_diffuse[3] = 1.0;

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glColor3f(r,g,b);

		glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
		glVertex3d(temp_v->x,temp_v->y,temp_v->z);
		}
		glEnd();

	}
    
	return;

	}

	void identifyMMS(Polyhedron* poly)
	{
		int minNum = 0;
		int maxNum = 0;
		int sadNum = 0;
		for(int i=0;i<poly->nverts;++i)
		{
			int plusToMinus = 0;
			for(int j=0;j<poly->vlist[i]->cList.size();++j)
				if((poly->vlist[i]->cList[j]->n->v->heat > poly->vlist[i]->cList[j]->v->heat)&&
					(poly->vlist[i]->cList[j]->v->heat >= poly->vlist[i]->cList[j]->o->n->v->heat))
					++plusToMinus;

			if(plusToMinus==0)
			{
				if(poly->vlist[i]->cList[0]->n->v->heat > poly->vlist[i]->cList[0]->v->heat)
					{
						poly->vlist[i]->minMaxSad = -1;
						++minNum;
				}
				else
					{
						poly->vlist[i]->minMaxSad = 0;
						++maxNum;
				}
			}
			else
			{
				poly->vlist[i]->minMaxSad = plusToMinus-1;
				sadNum = sadNum + plusToMinus -1;
			}
		}

		int M = maxNum + minNum - sadNum;
		cout<<"M is: "<<M<<endl;

	}


GLuint	texture[1];			// Storage For One Texture ( NEW )

AUX_RGBImageRec *LoadBMP(char *Filename)				// Loads A Bitmap Image
{
	FILE *File=NULL;									// File Handle

	if (!Filename)										// Make Sure A Filename Was Given
	{
		return NULL;									// If Not Return NULL
	}

	File=fopen(Filename,"r");							// Check To See If The File Exists

	if (File)											// Does The File Exist?
	{
		fclose(File);									// Close The Handle
		return auxDIBImageLoad(Filename);				// Load The Bitmap And Return A Pointer
	}

	return NULL;										// If Load Failed Return NULL
}

vector<double> s,t;
void textureSyn(int maxVertId, int minVertId, int maxVertId2, Polyhedron* poly)
{
	/*icVector3 vecMaxToMin;
	double x = poly->vlist[maxVertId]->x - poly->vlist[minVertId]->x;
	double y = poly->vlist[maxVertId]->y - poly->vlist[minVertId]->y;
	double z = poly->vlist[maxVertId]->z - poly->vlist[minVertId]->z;
	vecMaxToMin.set(x,y,z);
	double len1 = length(vecMaxToMin);
	normalize(vecMaxToMin);
	int maxVertId2 = -1;

	for(int i=0;i<poly->nverts;++i)
	{
		icVector3 vecMaxToMin2;

		double x2 = poly->vlist[i]->x - poly->vlist[minVertId]->x;
		double y2 = poly->vlist[i]->y - poly->vlist[minVertId]->y;
		double z2 = poly->vlist[i]->z - poly->vlist[minVertId]->z;
		vecMaxToMin2.set(x2,y2,z2);
		double len2 = length(vecMaxToMin2);
		normalize(vecMaxToMin2);
		double dotPro = dot(vecMaxToMin,vecMaxToMin2);
		if(dotPro<0.0001)
			if((len2>0.5*len1) && (len2<1.1*len1))
				{
					maxVertId2 = i;
					cout<<"maxVertId2: "<<maxVertId2<<endl;
					break;
			}
	}*/

	/*vector<double> s,t;*/
	heatDiffusion(maxVertId,minVertId,poly);
	for(int i=0;i<poly->nverts;++i)
		s.push_back(poly->vlist[i]->heat);

	heatDiffusion(maxVertId2,minVertId,poly);
	for(int i=0;i<poly->nverts;++i)
		t.push_back(poly->vlist[i]->heat);

	return;
}

void drawTexture(Polyhedron* poly)
{
	glEnable(GL_TEXTURE_2D);							// Enable Texture Mapping ( NEW )
	glShadeModel(GL_SMOOTH);							// Enable Smooth Shading
	//glClearColor(0.0f, 0.0f, 0.0f, 0.5f);				// Black Background
	//glClearDepth(1.0f);									// Depth Buffer Setup
	//glEnable(GL_DEPTH_TEST);							// Enables Depth Testing

	AUX_RGBImageRec *TextureImage[1];					// Create Storage Space For The Texture

	memset(TextureImage,0,sizeof(void *)*1);           	// Set The Pointer To NULL

	// Load The Bitmap, Check For Errors, If Bitmap's Not Found Quit
	if (TextureImage[0]=LoadBMP("../texture/try.bmp"))
	{									
		cout<<"Load successfully!"<<endl;

		glGenTextures(1, &texture[0]);					// Create The Texture

		// Typical Texture Generation Using Data From The Bitmap
		glBindTexture(GL_TEXTURE_2D, texture[0]);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, TextureImage[0]->sizeX, TextureImage[0]->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureImage[0]->data);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	}

	//glBindTexture(GL_TEXTURE_2D, texture[0]);

	Triangle *tri;
	//unsigned char x, y, z;
	//Vertex *verts;
	Vertex *temp_v;
	//ofstream normMap("../normMap.txt");

	for(int i=0;i<poly->ntris;++i)
	{
		tri = poly->tlist[i];

		glBegin(GL_POLYGON);
		for(int j=0;j<tri->nverts;++j)
		{
			temp_v = tri->verts[j];
			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
			glTexCoord2d(s[temp_v->index], t[temp_v->index]);
			glVertex3d(temp_v->x,temp_v->y,temp_v->z);
		}
		glEnd();
	}
	
	if (TextureImage[0])									// If Texture Exists
	{
		if (TextureImage[0]->data)							// If Texture Image Exists
		{
			free(TextureImage[0]->data);					// Free The Texture Image Memory
		}

		free(TextureImage[0]);								// Free The Image Structure
	}

	return;
}

void computeAMixed(Polyhedron* poly)
{
	for(int i=0;i<poly->ntris;++i)
		poly->tlist[i]->obtuse = false;
	for(int i=0;i<poly->ncorners;++i)
		if(poly->clist[i]->angle > PI/2.0)
			poly->clist[i]->t->obtuse = true;
	for(int i=0;i<poly->nverts;++i)
		{
			poly->vlist[i]->AMixed = 0.0;
			for(int j=0;j<poly->vlist[i]->cList.size();++j)
			{
				if(!poly->vlist[i]->cList[j]->t->obtuse)
					poly->vlist[i]->AMixed += 0.125*(poly->vlist[i]->cList[j]->n->cot*pow(poly->vlist[i]->cList[j]->n->e->length,2)
					+ poly->vlist[i]->cList[j]->p->cot*pow(poly->vlist[i]->cList[j]->p->e->length,2));
				else 
					if(poly->vlist[i]->cList[j]->angle>PI/2.0)
						poly->vlist[i]->AMixed += 0.5*poly->vlist[i]->cList[j]->t->area;
					else 
						poly->vlist[i]->AMixed += 0.25*poly->vlist[i]->cList[j]->t->area;
				//if(!poly->vlist[i]->cList[j]->t->obtuse)
				//{
				//	cout<<"1:  "<<0.125*(poly->vlist[i]->cList[j]->n->cot*pow(poly->vlist[i]->cList[j]->n->e->length,2)
				//	+ poly->vlist[i]->cList[j]->p->cot*pow(poly->vlist[i]->cList[j]->p->e->length,2))<<endl;
				//	cout<<"cot:	"<<poly->vlist[i]->cList[j]->n->cot<<endl;
				//	cout<<"pow:	"<<pow(poly->vlist[i]->cList[j]->n->e->length,2)<<endl;
				//}
				
				//else 
				//	if(poly->vlist[i]->cList[j]->angle>PI/2.0)
				//		cout<<"2:  "<<0.5*poly->vlist[i]->cList[j]->t->area<<endl;
				//	else 
				//		cout<<"3:  "<<0.25*poly->vlist[i]->cList[j]->t->area<<endl;
			}
	}

	return;
}

void computeGaussCurv(Polyhedron* poly)
{
	computeAMixed(poly);
	double deficitAngle = 0.0;
	for(int i=0;i<poly->nverts;++i)
	{
			poly->vlist[i]->GaussCurv = 2.0*PI;
			for(int j=0;j<poly->vlist[i]->cList.size();++j)
				poly->vlist[i]->GaussCurv -= poly->vlist[i]->cList[j]->angle;
			deficitAngle += poly->vlist[i]->GaussCurv;
			poly->vlist[i]->GaussCurv /= poly->vlist[i]->AMixed;
			//cout<<"C "<<poly->vlist[i]->GaussCurv<<endl;
			//cout<<"A"<<poly->vlist[i]->AMixed<<endl;
	}
	cout<<"Deficit angle:	"<<deficitAngle<<endl;
	cout<<"Deficit angle divided by 2*PI:	"<<deficitAngle/(2.0*PI)<<endl;
	return;
}

void computeMeanCurv(Polyhedron* poly)
{
	computeAMixed(poly);
	for(int i=0;i<poly->nverts;++i)
	{
		poly->vlist[i]->meanCurv = 0.0;
		for(int j=0;j<poly->vlist[i]->cList.size();++j)
		{
			//icVector3 Xi(poly->vlist[i]->x,poly->vlist[i]->y,poly->vlist[i]->z);
			//icVector3 Xj(poly->vlist[i]->cList[j]->n->v->x,poly->vlist[i]->cList[j]->n->v->y,poly->vlist[i]->cList[j]->n->v->z);
			double x = poly->vlist[i]->x - poly->vlist[i]->cList[j]->n->v->x;
			double y = poly->vlist[i]->y - poly->vlist[i]->cList[j]->n->v->y;
			double z = poly->vlist[i]->z - poly->vlist[i]->cList[j]->n->v->z;
			icVector3 Xij(x,y,z);
			double dotPro = dot(Xij,poly->vlist[i]->normal);
			poly->vlist[i]->meanCurv += dotPro*(poly->vlist[i]->cList[j]->p->cot 
				+ poly->vlist[i]->cList[j]->p->o->cot);
		}
		poly->vlist[i]->meanCurv /= 4.0*poly->vlist[i]->AMixed;
		//cout<<"C "<<poly->vlist[i]->meanCurv<<endl;
	}
	return;
}
bool compareGaussCurv (Vertex* a, Vertex* b) 
{ 
	return (a->GaussCurv<b->GaussCurv);
}

bool compareMeanCurv (Vertex* a, Vertex* b) 
{ 
	return (a->meanCurv<b->meanCurv);
}


void visCurv(Polyhedron* poly, int curvType)//if curvType is 0, Gaussian curv, if 2 mean curv
{
	GLfloat mat_diffuse[4];

    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	Triangle *tri;
	//unsigned char x, y, z;
	//Vertex *verts;
	Vertex *temp_v;
	//ofstream normMap("../normMap.txt");
	float max, min;
	if(curvType == 1)
	{
		max = poly->vlist[0]->GaussCurv;
		min = poly->vlist[0]->GaussCurv;
	}
	else if(curvType == 2)
	{
		max = poly->vlist[0]->meanCurv;
		min = poly->vlist[0]->meanCurv;
	}
	for(int i=1;i<poly->nverts;++i)
	{
			if(curvType == 1)
		{
			if(max<poly->vlist[i]->GaussCurv)
				max = poly->vlist[i]->GaussCurv;
			if(min>poly->vlist[i]->GaussCurv)
				min = poly->vlist[i]->GaussCurv;
		}
		else if(curvType == 2)
		{
			if(max<poly->vlist[i]->meanCurv)
				max = poly->vlist[i]->meanCurv;
			if(min>poly->vlist[i]->meanCurv)
				min = poly->vlist[i]->meanCurv;
	}
	}
	//cout<<"max:	"<<max<<"	min:	"<<min<<endl;
	//vector<Vertex*> vertList;
	//for(int i=0;i<poly->nverts;++i)
	//{
	//		if(curvType == 1)
	//			vertList.push_back(poly->vlist[i]);
	//	else if(curvType == 2)
	//			vertList.push_back(poly->vlist[i]);
	//}
	//if(curvType == 1)
	//	sort(vertList.begin(),vertList.end(),compareGaussCurv);
	//else if(curvType == 2)
	//	sort(vertList.begin(),vertList.end(),compareMeanCurv);
	//for(int i=0;i<poly->nverts;++i)
	//	cout<<vertList[i]->GaussCurv<<endl;
	//cout<<"max:	"<<vertList.back()->GaussCurv<<"	min:	"<<vertList.front()->GaussCurv<<endl;
	for(int i=0;i<poly->ntris;++i)
	{
		tri = poly->tlist[i];

		glBegin(GL_POLYGON);
		for(int j=0;j<tri->nverts;++j)
		{
		temp_v = tri->verts[j];
		
		float t;
		if(curvType == 1)
		{
			t = temp_v->GaussCurv;
		}
		else if(curvType == 2)
		{
			t = temp_v->meanCurv;
		}
		else return;
		//cout<<"t "<<t<<endl;
		t = (t-min)/(max-min);
		//unsigned int clr = (unsigned int)t*2^10;
		////cout<<"t "<<t<<endl;
		//unsigned char x, y, z;
		//z = clr>>16;
		//y = (clr&0xFF00)>>8;
		//x = (clr&0xFF);
		//float r = (255.0-(float)x)/255.0;
		//float g = (255.0-(float)y)/255.0;
	 //   float b = (255.0-(float)z)/255.0;
		float r = 6. * ( t - (5./6.) );
		float g = 0.0;
		float b = 0.0;

        if( t <= (5./6.) )
        {
                r = 0.0;
                g = 6. * ( t - (4./6.) );
                b = 0.;
        }

        if( t <= (4./6.) )
        {
                r = 0.;
                g = 0.0;
                b = 6. * ( t - (3./6.) );
        }

        if( t <= (3./6.) )
        {
                r = 0.0;
                g = 6. * ( t - (2./6.) );
                b = 6. * ( t - (2./6.) );
        }

        if( t <= (2./6.) )
        {
                r = 1.  -  6. * ( t - (1./6.) );
                g = 1.  -  6. * ( t - (1./6.) );
                b = 0.;
        }

        if( t <= (1./6.) )
        {
                r = 1.-6. * t;
				b = 1.-6. * t;
                g = 1.-6. * t;
        }

		mat_diffuse[0] = r;
		mat_diffuse[1] = g;
		mat_diffuse[2] = b;
		mat_diffuse[3] = 1.0;

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glColor3f(r,g,b);

		glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
		glVertex3d(temp_v->x,temp_v->y,temp_v->z);
		}
		glEnd();

	}
    
	return;

	}

void smoothCurv(int scheme, double dt, int times, int curvType, Polyhedron* poly)
{//curType 1 Gausscurv 2 meancurv
	if(scheme == 0)//uniform scheme 
	{
		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
		    double updateCurv = 0.0;
			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				if(curvType == 1)
				updateCurv += tempN->GaussCurv - tempV->GaussCurv;
				else if(curvType == 2)
				updateCurv += tempN->meanCurv - tempV->meanCurv;
				else return;
			}
			tempUpdateList.push_back(updateCurv);
		}
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double weight = dt*(1.0/tempV->cList.size());
			tempV->GaussCurv += weight*tempUpdateList[i];
			if(curvType == 1)
				tempV->GaussCurv += weight*tempUpdateList[i];
			else if(curvType == 2)
				tempV->meanCurv += weight*tempUpdateList[i];
			else return;

		}

		}

	}
	else if(scheme == 1)//cord scheme
	{
		vector<double> weightList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double totalE = 0.0;
			double updateCurv = 0.0;

			for(int j=0;j<tempV->cList.size();++j)
				totalE += 1.0/tempV->cList[j]->p->e->length;

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
				weightList.push_back(weight);
	
			}
		}

		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateList;
		int wgtInd = 0;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			//double totalE = 0.0;
			double updateCurv = 0.0;

			//for(int j=0;j<tempV->cList.size();++j)
			//	totalE += 1.0/tempV->cList[j]->p->e->length;

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
				if(curvType == 1)
				updateCurv += weightList[wgtInd++]*(tempN->GaussCurv - tempV->GaussCurv);
				else if(curvType == 2)
				updateCurv += weightList[wgtInd++]*(tempN->meanCurv - tempV->meanCurv);
				else return;
			}
			tempUpdateList.push_back(updateCurv);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			if(curvType == 1)
				poly->vlist[i]->GaussCurv += dt*tempUpdateList[i];
			else if(curvType == 2)
				poly->vlist[i]->meanCurv += dt*tempUpdateList[i];
			else return;
			
		}
		}

	}
	else if(scheme == 2)//mean curvature flow
	{
		vector<double> weightList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double totalE = 0.0;
			vector<double> eValue;
			double updateCurv = 0.0;

			for(int j=0;j<tempV->cList.size();++j)
			{
				double tempE = (tempV->cList[j]->p->cot
							   + tempV->cList[j]->p->o->cot)/2.0;
				totalE += tempE;
				eValue.push_back(tempE);
			}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = eValue[j]/totalE;
				weightList.push_back(weight);
			}
		}

		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateList;
		int wgtInd = 0;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			//double totalE = 0.0;
			double updateCurv = 0.0;
			//vector<double> eValue;
			//for(int j=0;j<tempV->cList.size();++j)
			//{
			//	double tempE = (tempV->cList[j]->p->cot
			//				   + tempV->cList[j]->p->o->cot)/2.0;
			//	totalE += tempE;
			//	eValue.push_back(tempE);
			//}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = eValue[j]/totalE;
				//if(curvType == 1)
				//	updateCurv += weight*(tempN->GaussCurv - tempV->GaussCurv);
				//else if(curvType == 2)
				//	updateCurv += weight*(tempN->meanCurv - tempV->meanCurv);
				//else return;
				//cout<<"weight		"<<weight<<endl;
				//cout<<"weightList[wgtInd++]	 "<<weightList[wgtInd++]<<endl;

				if(curvType == 1)
					updateCurv += weightList[wgtInd++]*(tempN->GaussCurv - tempV->GaussCurv);
				else if(curvType == 2)
					updateCurv += weightList[wgtInd++]*(tempN->meanCurv - tempV->meanCurv);
				else return;
			}
			  tempUpdateList.push_back(updateCurv);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			if(curvType == 1)
				poly->vlist[i]->GaussCurv += dt*tempUpdateList[i];
			else if(curvType == 2)
				poly->vlist[i]->meanCurv += dt*tempUpdateList[i];
			else return;
		}

		}
	}
	else if(scheme == 3)//mean values coordinates
	{
		vector<double> weightList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double totalE = 0.0;
			vector<double> eValue;
			double updateCurv = 0.0;

			for(int j=0;j<tempV->cList.size();++j)
			{
				double tempE = (tan(tempV->cList[j]->angle/2.0)
							   + tan(tempV->cList[j]->p->o->p->angle/2.0))/2.0;
				totalE += tempE;
				eValue.push_back(tempE);
			}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				double weight = eValue[j]/totalE;
				weightList.push_back(weight);
			}
		}

		for(int t=0;t<times;++t)
		{
		vector<double> tempUpdateList;
		int wgtInd = 0;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			//double totalE = 0.0;
			double updateCurv = 0.0;
			//vector<double> eValue;
			//for(int j=0;j<tempV->cList.size();++j)
			//{
			//	double tempE = (tan(tempV->cList[j]->angle/2.0)
			//				   + tan(tempV->cList[j]->p->o->p->angle/2.0))/2.0;
			//	totalE += tempE;
			//	eValue.push_back(tempE);
			//}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = eValue[j]/totalE;
				if(curvType == 1)
					updateCurv += weightList[wgtInd++]*(tempN->GaussCurv - tempV->GaussCurv);
				else if(curvType == 2)
					updateCurv += weightList[wgtInd++]*(tempN->meanCurv - tempV->meanCurv);
				else return;
			}
			  tempUpdateList.push_back(updateCurv);
		}

		for(int i=0;i<poly->nverts;++i)
		{
			if(curvType == 1)
				poly->vlist[i]->GaussCurv += dt*tempUpdateList[i];
			else if(curvType == 2)
				poly->vlist[i]->meanCurv += dt*tempUpdateList[i];
			else return;
		}
		}
	}

	return;
}

void computeNormCurv(Polyhedron* poly)
{
	for(int i=0;i<poly->nverts;++i)
	{
		for(int j=0;j<poly->vlist[i]->cList.size();++j)
		{
			double x = poly->vlist[i]->x - poly->vlist[i]->cList[j]->n->v->x;
			double y = poly->vlist[i]->y - poly->vlist[i]->cList[j]->n->v->y;
			double z = poly->vlist[i]->z - poly->vlist[i]->cList[j]->n->v->z;
			icVector3 edge(x,y,z);
			double eNorm = x*x + y*y + z*z;
			edge /= eNorm;
			double normCurv = 2*dot(edge,poly->vlist[i]->normal);
			poly->vlist[i]->normCurvList.push_back(normCurv);
		}
	}
	return;
}


void computeLocalFrm(Polyhedron* poly)
{
	for(int i=0;i<poly->nverts;++i)
	{
		double x = poly->vlist[i]->cList[0]->n->v->x - poly->vlist[i]->x;
		double y = poly->vlist[i]->cList[0]->n->v->y - poly->vlist[i]->y;
		double z = poly->vlist[i]->cList[0]->n->v->z - poly->vlist[i]->z;
		icVector3 edge(x,y,z);
		poly->vlist[i]->localY = cross(edge,poly->vlist[i]->normal);
		normalize(poly->vlist[i]->localY);
		poly->vlist[i]->localX = cross(poly->vlist[i]->localY,poly->vlist[i]->normal);
		normalize(poly->vlist[i]->localX);
	}
	return;
}
void computeTensor(Polyhedron* poly)
{
		for(int i=0;i<poly->nverts;++i)
	{
		vector<double> matrixA;
		for(int j=0;j<poly->vlist[i]->cList.size();++j)
		{
			double x = poly->vlist[i]->cList[0]->n->v->x - poly->vlist[i]->x;
			double y = poly->vlist[i]->cList[0]->n->v->y - poly->vlist[i]->y;
			double z = poly->vlist[i]->cList[0]->n->v->z - poly->vlist[i]->z;
			icVector3 edge(x,y,z);
			double a = dot(edge,poly->vlist[i]->localX);
			double b = dot(edge,poly->vlist[i]->localY);
			matrixA.push_back(a*a);
			matrixA.push_back(a*b*2.0);
			matrixA.push_back(b*b);
		}
		//overdertermined system Ax=k. A'Ax=A'k x=inv(A'A)A'k
		//const double zero = 0.0;
		//Mat_DP proATranA(zero,3,3),k(zero,3,1);
		//for(int m=0;m<3;++m)
		//	for(int j=0;j<3;++j)
		//	{
		//		int rowNum = matrixA.size()/3;
		//		for(int k=0;k<rowNum;++k)
		//			proATranA[m][j] += matrixA[m+k*3]*matrixA[j+k*3];
		//	}
		//NR::gaussj(proATranA,k);//proATranA is returned as the inverse 
		icMatrix3x3 proATranA;
		for(int m=0;m<3;++m)
			for(int j=0;j<3;++j)
			{
				int rowNum = matrixA.size()/3;
				for(int k=0;k<rowNum;++k)
					proATranA.entry[m][j] += matrixA[m+k*3]*matrixA[j+k*3];
			}
		icMatrix3x3 proATranAInv = inverse(proATranA);
		double proATranK[3];
		for(int m=0;m<3;++m)
		{
			int rowNum = matrixA.size()/3;
			proATranK[m] = 0.0;
			for(int k=0;k<rowNum;++k)
				proATranK[m] += matrixA[m+k*3]*poly->vlist[i]->normCurvList[k];
		}
		double t[3];
		for(int m=0;m<3;++m)
		{
			t[m] = 0.0;
			for(int k=0;k<3;++k)
				t[m] += proATranAInv.entry[m][k]*proATranK[k];
		}
		poly->vlist[i]->l = t[0];
		poly->vlist[i]->m = t[1];
		poly->vlist[i]->n = t[2];
	}
		return;
}

inline icMatrix3x3 transfMtx(icVector3 newX,icVector3 newY,icVector3 newZ,icVector3 oldX,icVector3 oldY,icVector3 oldZ)
{
	icMatrix3x3 tmp;
	
	tmp.entry[0][0] = dot(newX,oldX);
	tmp.entry[1][0] = dot(newX,oldY);
	tmp.entry[2][0] = dot(newX,oldZ);

	tmp.entry[0][1] = dot(newY,oldX);
	tmp.entry[1][1] = dot(newY,oldY);
	tmp.entry[2][1] = dot(newY,oldZ);

	tmp.entry[0][2] = dot(newZ,oldX);
	tmp.entry[1][2] = dot(newZ,oldY);
	tmp.entry[2][2] = dot(newZ,oldZ);

	return tmp;
}

icMatrix3x3 chgBasMtx(icMatrix3x3 transfMtx, icMatrix3x3 oldMtx)
{
	icMatrix3x3 transfMtxInv = inverse(transfMtx);
	icMatrix3x3 pro1 = multiply(transfMtxInv,oldMtx);
	icMatrix3x3 finalPro = multiply(pro1,transfMtx);
	return finalPro;
}

void computeGlobalTensor(Polyhedron* poly)
{
	icVector3 globalX(1.0,0.0,0.0),globalY(0.0,1.0,0.0),globalZ(0.0,0.0,1.0);
	for(int i=0;i<poly->nverts;++i)
	{
		Vertex* tempV = poly->vlist[i];

		icMatrix3x3 trLToG = transfMtx(globalX,globalY,globalZ,tempV->localX,tempV->localY,tempV->normal);
		icMatrix3x3 localTensor(0.0);
		localTensor.entry[0][0] = tempV->l;
		localTensor.entry[0][1] = tempV->m;
		localTensor.entry[1][0] = tempV->m;
		localTensor.entry[1][1] = tempV->n;
		tempV->globalTensor = chgBasMtx(trLToG,localTensor);
	}
	return;
}

void smoothGlobalTensor(int scheme, double dt, int times, Polyhedron* poly)//with updated weight
{
	//icVector3 globalX(1.0,0.0,0.0),globalY(0.0,1.0,0.0),globalZ(0.0,0.0,1.0);
	if(scheme == 0)//uniform scheme 
	{
		for(int t=0;t<times;++t)
		{
		vector<icMatrix3x3> tempUpdateList;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			icMatrix3x3 updateTensor(0.0);
			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				updateTensor += tempN->globalTensor - tempV->globalTensor;
				
			}
			tempUpdateList.push_back(updateTensor);
		}
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			double weight = dt*(1.0/tempV->cList.size());
			tempV->globalTensor += tempUpdateList[i]*weight;
		}

		}

	}
	else if(scheme == 1)//cord scheme
	{
		vector<double> weightList;
		if(times>0)
		{
			vector<icMatrix3x3> tempUpdateList;
			for(int i=0;i<poly->nverts;++i)
			{
				Vertex* tempV = poly->vlist[i];
				double totalE = 0.0;
				icMatrix3x3 updateTensor(0.0);

				for(int j=0;j<tempV->cList.size();++j)
					totalE += 1.0/tempV->cList[j]->p->e->length;

				for(int j=0;j<tempV->cList.size();++j)
				{
					Vertex* tempN = tempV->cList[j]->n->v;
					double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
					weightList.push_back(weight);
					updateTensor += (tempN->globalTensor - tempV->globalTensor)*weight;
				}
				tempUpdateList.push_back(updateTensor);
			}

			for(int i=0;i<poly->nverts;++i)
				poly->vlist[i]->globalTensor += tempUpdateList[i]*dt;
		}
		for(int t=1;t<times;++t)
		{
		vector<icMatrix3x3> tempUpdateList;
		int wgtInd = 0;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			//double totalE = 0.0;
			icMatrix3x3 updateTensor(0.0);

			//for(int j=0;j<tempV->cList.size();++j)
			//	totalE += 1.0/tempV->cList[j]->p->e->length;

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = (1.0/tempV->cList[j]->p->e->length)/totalE;
				updateTensor += (tempN->globalTensor - tempV->globalTensor)*weightList[wgtInd++];
			}
			tempUpdateList.push_back(updateTensor);
		}

		for(int i=0;i<poly->nverts;++i)
			poly->vlist[i]->globalTensor += tempUpdateList[i]*dt;
		}

	}
	else if(scheme == 2)//mean curvature flow
	{
		vector<double> weightList;
		if(times>0)
		{
			vector<icMatrix3x3> tempUpdateList;
			for(int i=0;i<poly->nverts;++i)
			{
				Vertex* tempV = poly->vlist[i];
				double totalE = 0.0;
				vector<double> eValue;
				icMatrix3x3 updateTensor(0.0);
				for(int j=0;j<tempV->cList.size();++j)
				{
					double tempE = (tempV->cList[j]->p->cot
								   + tempV->cList[j]->p->o->cot)/2.0;
					totalE += tempE;
					eValue.push_back(tempE);
				}

				for(int j=0;j<tempV->cList.size();++j)
				{
					Vertex* tempN = tempV->cList[j]->n->v;
					double weight = eValue[j]/totalE;
					weightList.push_back(weight);
					updateTensor += (tempN->globalTensor - tempV->globalTensor)*weight;
				}
				  tempUpdateList.push_back(updateTensor);
			}

			for(int i=0;i<poly->nverts;++i)
				poly->vlist[i]->globalTensor += tempUpdateList[i]*dt;
		}

		for(int t=1;t<times;++t)
		{
		vector<icMatrix3x3> tempUpdateList;
		int wgtInd = 0;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			//double totalE = 0.0;
			//vector<double> eValue;
			icMatrix3x3 updateTensor(0.0);
			//for(int j=0;j<tempV->cList.size();++j)
			//{
			//	double tempE = (tempV->cList[j]->p->cot
			//				   + tempV->cList[j]->p->o->cot)/2.0;
			//	totalE += tempE;
			//	eValue.push_back(tempE);
			//}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = eValue[j]/totalE;
				updateTensor += (tempN->globalTensor - tempV->globalTensor)*weightList[wgtInd++];
			}
			  tempUpdateList.push_back(updateTensor);
		}

		for(int i=0;i<poly->nverts;++i)
			poly->vlist[i]->globalTensor += tempUpdateList[i]*dt;

		}
	}
	else if(scheme == 3)//mean values coordinates
	{
		vector<double> weightList;
		if(times>0)
		{
			vector<icMatrix3x3> tempUpdateList;
			for(int i=0;i<poly->nverts;++i)
			{
				Vertex* tempV = poly->vlist[i];
				double totalE = 0.0;
				icMatrix3x3 updateTensor(0.0);
				vector<double> eValue;
				for(int j=0;j<tempV->cList.size();++j)
				{
					double tempE = (tan(tempV->cList[j]->angle/2.0)
								   + tan(tempV->cList[j]->p->o->p->angle/2.0))/2.0;
					totalE += tempE;
					eValue.push_back(tempE);
				}

				for(int j=0;j<tempV->cList.size();++j)
				{
					Vertex* tempN = tempV->cList[j]->n->v;
					double weight = eValue[j]/totalE;
					weightList.push_back(weight);
					updateTensor += (tempN->globalTensor - tempV->globalTensor)*weight;
				}
				  tempUpdateList.push_back(updateTensor);
			}

			for(int i=0;i<poly->nverts;++i)
				poly->vlist[i]->globalTensor += tempUpdateList[i]*dt;
		}

		for(int t=1;t<times;++t)
		{
		vector<icMatrix3x3> tempUpdateList;
		int wgtInd = 0;
		for(int i=0;i<poly->nverts;++i)
		{
			Vertex* tempV = poly->vlist[i];
			//double totalE = 0.0;
			icMatrix3x3 updateTensor(0.0);
			//vector<double> eValue;
			//for(int j=0;j<tempV->cList.size();++j)
			//{
			//	double tempE = (tan(tempV->cList[j]->angle/2.0)
			//				   + tan(tempV->cList[j]->p->o->p->angle/2.0))/2.0;
			//	totalE += tempE;
			//	eValue.push_back(tempE);
			//}

			for(int j=0;j<tempV->cList.size();++j)
			{
				Vertex* tempN = tempV->cList[j]->n->v;
				//double weight = eValue[j]/totalE;
				updateTensor += (tempN->globalTensor - tempV->globalTensor)*weightList[wgtInd++];
			}
			  tempUpdateList.push_back(updateTensor);
		}

		for(int i=0;i<poly->nverts;++i)
			poly->vlist[i]->globalTensor += tempUpdateList[i]*dt;
		}
	}

	return;
}

void globTensorToLocal(Polyhedron* poly)
{
	icVector3 globalX(1.0,0.0,0.0),globalY(0.0,1.0,0.0),globalZ(0.0,0.0,1.0);
	for(int i=0;i<poly->nverts;++i)
	{
		Vertex* tempV = poly->vlist[i];
		icMatrix3x3 trGToL = transfMtx(tempV->localX,tempV->localY,tempV->normal,globalX,globalY,globalZ);
		icMatrix3x3 localTensor = chgBasMtx(trGToL,tempV->globalTensor);
		tempV->l = localTensor.entry[0][0];
		tempV->m = (localTensor.entry[0][1]+localTensor.entry[1][0])/2.0;
		tempV->n = localTensor.entry[1][1];
	}
	return;
}

void LocToGloEigV(Polyhedron* poly)
{
	icVector3 globalX(1.0,0.0,0.0),globalY(0.0,1.0,0.0),globalZ(0.0,0.0,1.0);
	for(int i=0;i<poly->nverts;++i)
	{
		Vertex* tempV = poly->vlist[i];
		icMatrix3x3 trLToG = transfMtx(globalX,globalY,globalZ,tempV->localX,tempV->localY,tempV->normal);
		icMatrix3x3 trLToGInv = inverse(trLToG);
		icVector3 lcoMajT(tempV->majEigv.entry[0],tempV->majEigv.entry[1],0.0);
		tempV->gloMajEigv = trLToGInv*lcoMajT;
		icVector3 lcoMinT(tempV->minEigv.entry[0],tempV->minEigv.entry[1],0.0);
		tempV->gloMinEigv = trLToGInv*lcoMinT;
	}
	return;
}
void computeTensorEig(Polyhedron* poly)
{
	computeNormCurv(poly);
	computeLocalFrm(poly);
	computeTensor(poly);
	computeGlobalTensor(poly);
	cout<<"Please input the smooth scheme:"<<endl;
	cout<<"0	uniform"<<endl;
	cout<<"1	cord"<<endl;
	cout<<"2	mean curvature"<<endl;
	cout<<"3	mean value"<<endl;
	int scheme;
	cin>>scheme;
	cout<<"How many times do you want to do the smooth?"<<endl;
	int times;
	cin>>times;
	smoothGlobalTensor(scheme, 0.3, times,poly);
	globTensorToLocal(poly);

	for(int i=0;i<poly->nverts;++i)
	{
		Mat_DP tensor(2,2), eigMatrix(2,2);
		Vec_DP eigVal(2);
		tensor[0][0]=poly->vlist[i]->l;
		tensor[0][1]=tensor[1][0]=poly->vlist[i]->m;
		tensor[1][1]=poly->vlist[i]->n;
		int nrot;
		NR::jacobi(tensor,eigVal,eigMatrix,nrot);
		int minind,maxind;
		double mineig, maxeig;
		  if(eigVal[0]<eigVal[1])
		  {
			  mineig = eigVal[0];
			  minind = 0;
			  maxeig = eigVal[1];
			  maxind = 1;
		  }
		  else
		  {
			  mineig = eigVal[1];
			  minind = 1;
			  maxeig = eigVal[0];
			  maxind = 0;
		  }
		poly->vlist[i]->majEigv.set(eigMatrix[maxind][0],eigMatrix[maxind][1]);
		normalize(poly->vlist[i]->majEigv);
		poly->vlist[i]->minEigv.set(eigMatrix[minind][0],eigMatrix[minind][1]);
		normalize(poly->vlist[i]->minEigv);
	}
	LocToGloEigV(poly);
	return;
}

void drawLine(icVector3 &x, icVector3 &y, float width, float r, float g, float b )
{ 
    glDisable(GL_TEXTURE_2D); 
    glEnable(GL_BLEND); 
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
    glColor3f( r, g, b); 

    glLineWidth(width); 
    glBegin(GL_LINES); 
    glVertex3d( x.x, x.y, x.z); 
    glVertex3d( y.x, y.y, y.z); 
    glEnd(); 

    glDisable(GL_BLEND); 
    glEnable(GL_TEXTURE_2D); 
	return;
}  

void drawCurvTensor(Polyhedron* poly)
{
	for(int i=0;i<poly->nverts;++i)
	{
		icVector3 point1(poly->vlist[i]->x,poly->vlist[i]->y,poly->vlist[i]->z);
		icVector3 point2 = point1 + poly->vlist[i]->gloMajEigv*0.02;
		icVector3 point3 = point1 - poly->vlist[i]->gloMajEigv*0.02;
		icVector3 point4 = point1 + poly->vlist[i]->gloMinEigv*0.02;
		icVector3 point5 = point1 - poly->vlist[i]->gloMinEigv*0.02;
		point2 += poly->vlist[i]->normal*0.005;
		point3 += poly->vlist[i]->normal*0.005;
		point4 += poly->vlist[i]->normal*0.005;
		point5 += poly->vlist[i]->normal*0.005;
		drawLine(point2,point3,0.6,1.0,0.0,0.0);//red major pricipal curvature
		drawLine(point4,point5,0.6,0.0,0.0,1.0);//green minor
	}
	return;
}

