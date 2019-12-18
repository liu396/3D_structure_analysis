#include <iostream>
#include <vector>
#include <cmath>
#include "mkl.h"
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <algorithm> 

using namespace std;


vector<vector<double> > InternalAtomPosition;
vector<vector<int> > InternalAtomInfo;
vector<vector<double> > SurfaceAtomPosition;
vector<vector<int> > SurfaceAtomInfo;
vector<vector<double> > AtomPosition;//x,y,z position
vector<vector<int> > AtomInfo;//type and id
vector<vector<double> > AtomOther;//other information

int L[3];
double h[3][3];
int *head;
int *lscl;
int *shead;
int *slscl;
int *inhead;
int *inlscl;
double rc_body=3.0;
double rc_surface=6.0;
double rc_alpha[3];
double fit_threshold=1.5;
vector<int> surface_atoms;

void initialize_vector(vector<double> & v);
void print_vector(vector<double> & to_be_print);
void write_surface_cfg();
void read_CFG(ifstream &file);
void read_internal_cfg(ifstream &file);
void read_internal_cfg(ifstream &file);
void read_surface_cfg(ifstream &file);
void find_nearest_neighbor(vector<int> &neighbor_list,int num, double r,vector<vector<double> > &AtomPosition,int *head,int *lscl);
extern void print_matrix(char* desc, MKL_INT m, MKL_INT n, double *a, MKL_INT lda);
extern void print_int_vector(char* desc, MKL_INT n, MKL_INT *a);
extern void print_ge_vector(char *name, MKL_INT n,double *b);
vector<vector<double> > new_points(vector<vector<double> > &old_points, double vertical_axe[3], double new_normal[3]);
vector<double> fitting_surface(vector<vector<double> >&point_list);
vector<double> fitting_qudratic_surface(vector<vector<double> > &point_list);
void neighbors_for_surface(vector<int> &neighbor_list,int num,double r);
vector<int> neighbor_expansion(int num,int steps,double r,vector<vector<double> > &AtomPosition, int *head, int *slscl);
vector<vector<double> > neighbor_expansion_v2(vector<double> &position,int steps,double r, vector<vector<double> > &AtomPosition, int *head, int *lscl);
vector<double> find_normal(vector<double> &position,double rc);
double diameter(vector<double> &position,double movement);
vector<vector<double> > find_atoms_around(vector<double> &position, double rc,vector<vector<double> > &AtomPosition,int* head, int* lscl);
vector<vector<double> > find_atoms_around_v2(vector<double> &position, double rc,vector<vector<double> > &AtomPosition,int* head, int* lscl);
vector<double> bivariate_surface_gradient(vector<double>cur_pos,vector<double> &parameters);
int consistency_of_normal(vector<double> &cur_pos, vector<double> &normal_direction, double movement);