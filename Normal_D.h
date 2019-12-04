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
double rc_surface=5.0;
double rc_alpha[3];
vector<int> surface_atoms;


void read_CFG(ifstream &file);
void read_internal_cfg(ifstream &file);
void read_surface_cfg(ifstream &file);
double diameter(vector<double> &position,double movement);
vector<double> find_normal(vector<double> &position,double rc);
vector<vector<double> > find_atoms_around(vector<double> &position, double rc,vector<vector<double> > &AtomPosition,int* head, int* lscl);
