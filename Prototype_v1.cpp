#include "Prototype_v1.h"

using namespace std;

int main(int argc, const char *argv[]){
	int i,j,a;
	int layers;
	int mc[3];
	int c;
	double Sv=0.0359;
	//double Sv=0.00342983;
	
	ifstream cfgfile,surfacecfgfile,internalcfgfile;
    vector<vector<double> > point_list;
    vector<double> parameters;
    double vec[3]={0.0,0.0,0.0};
    double norm_factor;
    vector<int> neighbor_list;
    double k_h,k_k,k_1,k_2,phix,phiy,phiz,phixx,phiyy,phizz,phixy,phiyz,phixz;
    double x,y,z;
    vector<double> normal{3,0};
    vector<vector<double> > new_neighbors;
	double distance;
    ofstream normalfile;
	ofstream cur_output;
	double up_limit=0.2, low_limit=-0.2;
	double interval,k1_min_value,k1_max_value,k2_min_value,k2_max_value;
	double large,small;
	int total_count=0;
	int bin=40;
	int *cur_data;
	int m,n,failed=0;
	vector<double> cur_pos(3,0);
	vector<double> analytical_normal;
	
	large=3;
	small=1;
	k1_min_value=-large*Sv;
	k1_max_value=small*Sv;
	k2_min_value=-small*Sv;
	k2_max_value=large*Sv;
	
    
	//=============read cfg file=============================================
	cfgfile.open("original.cfg");
	surfacecfgfile.open("original_surface.cfg");
	internalcfgfile.open("original_internal.cfg");
	read_CFG(cfgfile);
	read_surface_cfg(surfacecfgfile);
	read_internal_cfg(internalcfgfile);
    cout<<"There are "<<AtomPosition.size()<<" saved\n";
    cout<<"Box vectors\n";
    cout<<h[0][0]<<' '<<h[1][0]<<' '<<h[2][0]<<'\n';
    cout<<h[0][1]<<' '<<h[1][1]<<' '<<h[2][1]<<'\n';
    cout<<h[0][2]<<' '<<h[1][2]<<' '<<h[2][2]<<'\n';
    
    cout<<"Total box divived in to x,y,z:" <<L[0]<<' '<<L[1]<<' '<<L[2]<<' '<<'\n';
    cout<<"Single box size: "<<rc_alpha[0]<<' '<<rc_alpha[1]<<' '<<rc_alpha[2]<<'\n';
	//=======================================================================
	
//    find_nearest_neighbor(neighbor_list,0,4.2,AtomPosition,L,rc_alpha,head,lscl,h);
//    for(i=0;i<neighbor_list.size();i++){
//        cout<<AtomInfo[neighbor_list[i]][0]<<'\n';
//    }
//    exit(0);
    //=======================================================================
	//================Create Surface atom list===============================
		//     cout<<"Finding Surface Atoms\n";
		//     for(i=0;i<AtomPosition.size();i++){
		//         find_nearest_neighbor(neighbor_list, i, 4.0, AtomPosition,head,lscl);
		//         //if (i%10000==0) cout<<"Done "<<i<<' '<<neighbor_list.size()<<"neighbors\n";
		// //cout<<"neignbor number: "<<neighbor_list.size()<<'\n';
		//         if (neighbor_list.size()<9){
		//             surface_atoms.push_back(i);
		//         }
		//     }
		//     cout<<"Surface List Created\n";
		//     cout<<"There are in total "<<surface_atoms.size()<<" Surface atoms\n";
	//=======================================================================
	//=========================Create Surface Atom Linked List===============//
    // shead=new int[L[0]*L[1]*L[2]];
    // slscl=new int[AtomPosition.size()];
    // for(i=0;i<L[0]*L[1]*L[2];i++) shead[i]=-1;
    // for(i=0;i<AtomPosition.size();i++) slscl[i]=-1;
    // for(i=0;i<surface_atoms.size();i++){
    //     for(a=0;a<3;a++) mc[a]=AtomPosition[surface_atoms[i]][a]/rc_alpha[a];
    //     c=mc[0]*L[1]*L[2]+mc[1]*L[2]+mc[2];
    //     slscl[surface_atoms[i]]=shead[c];
    //     shead[c]=surface_atoms[i];
    // }
	// for(i=0;i<L[0]*L[1]*L[2];i++){
	// 	if(shead[i]!=-1) {
	// 		printf("shead[%d]=%d\n",i,shead[i]);
	// 		j=shead[i];
	// 		while(j!=-1){
	// 			printf("it contains atom : %d\n",slscl[j]);
	// 			j=slscl[j];
	// 		}
	// 	}
	// }
	
	//cout<<"Surface atom linked list created\n";
	
	//=======================================================================//
	//write_surface_cfg();
	//cout<<"Surface atoms written in cfg file.\n";
	//=====find normal on the surface and write it and its normal direction in a cfg file
	//normal=find_normal(AtomPosition[surface_atoms[40]],6.5);
    //find_nearest_neighbor(neighbor_list,surface_atoms[40],6.5,AtomPosition,L,rc_alpha,head,lscl,h);
    //cout<<"Finding neighbors of atom "<<AtomInfo[surface_atoms[40]][0]<<'\n';
    //cout<<"Its neighbors: \n";
    // for (i=0;i<neighbor_list.size();i++){
    //     cout<<AtomInfo[neighbor_list[i]][0]<<'\n';
    //     normal[0]=normal[0]+AtomPosition[surface_atoms[40]][0]-AtomPosition[neighbor_list[i]][0];
    //     normal[1]=normal[1]+AtomPosition[surface_atoms[40]][1]-AtomPosition[neighbor_list[i]][1];
    //     normal[2]=normal[2]+AtomPosition[surface_atoms[40]][2]-AtomPosition[neighbor_list[i]][2];
    // }
    //normalfile.open("normal.cfg");
    //normalfile<<"Number of particles = "<<2<<'\n';
	//     normalfile<<"A = 1.0 Angstrom (basic length-scale)\n";
	//     for (i=0;i<3;i++){
	//         for(j=0;j<3;j++){
	//            normalfile<<"H0("<<j+1<<','<<i+1<<") = "<<h[j][i]<<" A\n";
	//        }
	//     }
	//     normalfile<<".NO_VELOCITY.\n";
	//     normalfile<<"entry_count = 5\n"; //Information Just Type and Grain ID
	//     //What is element type, type 1 for what? Type 2 for what?
	//     normalfile<<"auxiliary[0] = type\n";
	//     normalfile<<"auxiliary[1] = normal\n";
	//     normalfile<<196.0<<'\n'<<"Au\n";
	//     normalfile<<AtomPosition[surface_atoms[40]][0]/h[0][0]<<' '<<AtomPosition[surface_atoms[40]][1]/h[1][1]<<' '<<AtomPosition[surface_atoms[40]][2]/h[2][2]<<' '<<1<<' '<<1<<'\n';
	//     normalfile<<normal[0]<<' '<<normal[1]<<' '<<normal[2]<<' '<<1<<' '<<2<<'\n';
	//     normalfile.close();
	// printf("file showing a normal finished\n");
	
	//Test for single point//
	cout<<"Tested Atom : SurfaceAtom 657. Particle Identifier="<<SurfaceAtomInfo[26333][0]<<'\n';
	distance=diameter(SurfaceAtomPosition[26333],rc_body/3.0);
	cout<<"Diameter from the point: "<<distance<<'\n';
	point_list=neighbor_expansion_v2(SurfaceAtomPosition[26333],3,rc_surface, SurfaceAtomPosition, shead,slscl);
	parameters=fitting_surface(point_list);

    x=SurfaceAtomPosition[26333][0];
    y=SurfaceAtomPosition[26333][1];
    z=SurfaceAtomPosition[26333][2];
	
    phix=parameters[3]+2.0*parameters[0]*x+parameters[2]*y;
    phiy=parameters[4]+2.0*parameters[1]*y+parameters[2]*x;
    phiz=-1;
    phixx=2*parameters[0];
    phiyy=2*parameters[1];
    phizz=0;
    phixy=parameters[2];
    phixz=0;
    phiyz=0;
	
	cur_pos[0]=x;
	cur_pos[1]=y;
	cur_pos[2]=z;
	
	cout<<"Postion discussed: "<<cur_pos[0]<<' '<<cur_pos[1]<<' '<<cur_pos[2]<<'\n';
	
	analytical_normal=bivariate_surface_gradient(cur_pos,parameters);
	cout<<"Analytical_normal: "<<analytical_normal[0]<<' '<<analytical_normal[1]<<' '<<analytical_normal[2]<<'\n';

    k_h=(phix*phix*(phiyy+phizz)+phiy*phiy*(phixx+phizz)+phiz*phiz*(phixx+phiyy))/(2.0*pow(phix*phix+phiy*phiy+phiz*phiz,1.5))-(phix*phiy*phixy+phix*phiz*phixz+phiy*phiz*phiyz)/(pow(phix*phix+phiy*phiy+phiz*phiz,1.5));
	k_h= (double)consistency_of_normal(cur_pos,analytical_normal,rc_body/3.0)*k_h;

    k_k=2*(phix*phiy*(phixz*phiyz-phixy*phizz)+phix*phiz*(phixy*phiyz-phixz*phiyy)+phiy*phiz*(phixy*phixz-phiyz*phixx))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0)+(phix*phix*(phiyy*phizz-phiyz*phiyz)+phiy*phiy*(phixx*phizz-phixz*phixz)+phiz*phiz*(phixx*phiyy-phixy*phixy))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0);

 	//printf("mean and gaussian curvatures:%f and %f\n",k_h,k_k);

    k_1=k_h-sqrt(abs(k_h*k_h-k_k));
    k_2=k_h+sqrt(abs(k_h*k_h-k_k));
	
	printf("K1 and K2 are: %f and %f\n",k_1,k_2);
	
	
	//exit(0);
	
	interval=(k1_max_value-k1_min_value)/double(bin);
	cur_data=new int[bin*bin];
	
	for(i=0;i<bin*bin;i++){
		cur_data[i]=0;
	}
	
	for(i=0;i<SurfaceAtomPosition.size();i++){
		cout<<"Processing Atom "<<i<<'\n';
		distance=diameter(SurfaceAtomPosition[i],rc_body/3.0);
		layers=int(distance/5/rc_surface);
		if (layers<3) layers=3;
		point_list=neighbor_expansion_v2(SurfaceAtomPosition[i],layers,rc_surface, SurfaceAtomPosition, shead,slscl);
		parameters=fitting_surface(point_list);
		if (parameters.size()==0) {
			failed++;
			continue;
		}
	    x=SurfaceAtomPosition[i][0];
	    y=SurfaceAtomPosition[i][1];
	    z=SurfaceAtomPosition[i][2];
		
	 	// x=0;
	 	// y=0;
	 	// z=0;

	    phix=parameters[3]+2.0*parameters[0]*x+parameters[2]*y;
	    phiy=parameters[4]+2.0*parameters[1]*y+parameters[2]*x;
	    phiz=-1;
	    phixx=2*parameters[0];
	    phiyy=2*parameters[1];
	    phizz=0;
	    phixy=parameters[2];
	    phixz=0;
	    phiyz=0;
		
		cur_pos[0]=x;
		cur_pos[1]=y;
		cur_pos[2]=z;
		cout<<"Postion discussed: "<<cur_pos[0]<<' '<<cur_pos[1]<<' '<<cur_pos[2]<<'\n';
		analytical_normal=bivariate_surface_gradient(cur_pos,parameters);
		cout<<"Analytical_normal: "<<analytical_normal[0]<<' '<<analytical_normal[1]<<' '<<analytical_normal[2]<<'\n';
		
	    k_h=(phix*phix*(phiyy+phizz)+phiy*phiy*(phixx+phizz)+phiz*phiz*(phixx+phiyy))/(2.0*pow(phix*phix+phiy*phiy+phiz*phiz,1.5))-(phix*phiy*phixy+phix*phiz*phixz+phiy*phiz*phiyz)/(pow(phix*phix+phiy*phiy+phiz*phiz,1.5));
		k_h= (double)consistency_of_normal(cur_pos,analytical_normal,rc_body/3.0)*k_h;

	    k_k=2*(phix*phiy*(phixz*phiyz-phixy*phizz)+phix*phiz*(phixy*phiyz-phixz*phiyy)+phiy*phiz*(phixy*phixz-phiyz*phixx))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0)+(phix*phix*(phiyy*phizz-phiyz*phiyz)+phiy*phiy*(phixx*phizz-phixz*phixz)+phiz*phiz*(phixx*phiyy-phixy*phixy))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0);

	 	//printf("mean and gaussian curvatures:%f and %f\n",k_h,k_k);

	    k_1=k_h-sqrt(abs(k_h*k_h-k_k));
	    k_2=k_h+sqrt(abs(k_h*k_h-k_k));
		
		
		m=int((k_1-k1_min_value)/interval);
		// if (m>=bin) m=bin-1;
		// if (m<0) m=0;
		if (m>=bin) continue;
		if (m<0) continue;
		n=int((k_2-k2_min_value)/interval);
		// if (n>=bin) n=bin-1;
		// if (n<0) n=0;
		if (n>=bin) continue;
		if (n<0) continue;
		cur_data[m*bin+n]++;
		total_count++;
	 	printf("K1 and K2 are: %f and %f\n",k_1,k_2);
	}
	cur_output.open("k1_vs_k2.txt");
	cur_output<<"K1_interval            k2_interval             Count\n";
	cout<<"Failed points: "<<failed<<'\n';
	cout<<"Total points: "<<total_count++;
	
	
	for(i=0;i<bin;i++){
		for(j=0;j<bin;j++){
			cur_output<<i<<"   "<<j<<"   "<<cur_data[i*bin+j]/((double)(large+small)/(double)(bin))/((double)(large+small)/(double)(bin))/total_count<<'\n';
		}
	}
	
	
	
	// distance=diameter(AtomPosition[surface_atoms[40]],6.5/2);
	// cout<<"Approximate Distance: "<<distance<<'\n';
	
    //====Normalize the normal vector========================================
    // double norm=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
//     normal[0]/=norm;
//     normal[1]/=norm;
//     normal[2]/=norm;
    
    //=====Find nearest neighbor list from surface atoms. Two ways====================
		// =========Just nearest neighbors=======
    // neighbors_for_surface(neighbor_list,surface_atoms[40],60,surface_atoms,AtomPosition);
		//=========Nearest neighbors of nearest neighbors======
	// neighbor_list=neighbor_expansion(surface_atoms[40],int(distance/5/rc_surface),6.5,AtomPosition,shead,slscl);
	// cout<<"Total number: "<<neighbor_list.size()<<'\n';
	// for(i=0;i<neighbor_list.size();i++){
	// 	cout<<"ParticleIdentifier=="<<AtomInfo[neighbor_list[i]][0]<<"||\n";
	// }
	//
	// for(i=0;i<neighbor_list.size();i++){
	// 	//point_list.push_back({AtomPosition[neighbor_list[i]][0]-AtomPosition[neighbor_list[0]][0],AtomPosition[neighbor_list[i]][1]-AtomPosition[neighbor_list[0]][1],AtomPosition[neighbor_list[i]][2]-AtomPosition[neighbor_list[0]][2]});
	// 	point_list.push_back({AtomPosition[neighbor_list[i]][0],AtomPosition[neighbor_list[i]][1],AtomPosition[neighbor_list[i]][2]});
	// }

	//=====bivariate fitting and calculate the curvature=====================
    // parameters=fitting_surface(point_list);
//  	printf("Bivariate surface parameters:%f %f %f %f %f %f\n",parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5]);
//     x=AtomPosition[surface_atoms[40]][0];
//     y=AtomPosition[surface_atoms[40]][1];
//     z=AtomPosition[surface_atoms[40]][2];
//  	// x=0;
//  	// y=0;
//  	// z=0;
//
//     phix=parameters[3]+2.0*parameters[0]*x+parameters[2]*y;
//     phiy=parameters[4]+2.0*parameters[1]*y+parameters[2]*x;
//     phiz=-1;
//     phixx=2*parameters[0];
//     phiyy=2*parameters[1];
//     phizz=0;
//     phixy=parameters[2];
//     phixz=0;
//     phiyz=0;
//
//
//     k_h=(phix*phix*(phiyy+phizz)+phiy*phiy*(phixx+phizz)+phiz*phiz*(phixx+phiyy))/(2.0*pow(phix*phix+phiy*phiy+phiz*phiz,1.5))-(phix*phiy*phixy+phix*phiz*phixz+phiy*phiz*phiyz)/(pow(phix*phix+phiy*phiy+phiz*phiz,1.5));
//
//     k_k=2*(phix*phiy*(phixz*phiyz-phixy*phizz)+phix*phiz*(phixy*phiyz-phixz*phiyy)+phiy*phiz*(phixy*phixz-phiyz*phixx))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0)+(phix*phix*(phiyy*phizz-phiyz*phiyz)+phiy*phiy*(phixx*phizz-phixz*phixz)+phiz*phiz*(phixx*phiyy-phixy*phixy))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0);
//
//  	//printf("mean and gaussian curvatures:%f and %f\n",k_h,k_k);
//
//     k_1=k_h-sqrt(abs(k_h*k_h-k_k));
//     k_2=k_h+sqrt(abs(k_h*k_h-k_k));
//
//  	printf("K1 and K2 are: %f and %f\n",k_1,k_2);
//
//     //new_neighbors=new_points(neighbor_list,{0.0,0.0,1.0},normal);
//
	//============quadratic surface fitting===========================//
	// parameters=fitting_qudratic_surface(point_list);
// 	printf("Quadratic surface parameters:%f %f %f %f %f %f %f %f %f %f %f\n",parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6],parameters[7],parameters[8],parameters[9],parameters[10]);
//     // x=AtomPosition[surface_atoms[40]][0];
//     // y=AtomPosition[surface_atoms[40]][1];
//     // z=AtomPosition[surface_atoms[40]][2];
// 	x=point_list[0][0];
// 	y=point_list[0][1];
// 	z=point_list[0][2];
// 	// x=0;
// 	// y=0;
// 	// z=0;
//
//     phix=parameters[6]+2.0*parameters[0]*x+parameters[5]*y+parameters[4]*z;
//     phiy=parameters[7]+2.0*parameters[1]*y+parameters[5]*x+parameters[3]*z;
//     phiz=parameters[8]+2.0*parameters[2]*z+parameters[3]*y+parameters[4]*x;
//     phixx=2*parameters[0];
//     phiyy=2*parameters[1];
//     phizz=2*parameters[2];
//     phixy=parameters[5];
//     phixz=parameters[4];
//     phiyz=parameters[3];
//
//
//     k_h=(phix*phix*(phiyy+phizz)+phiy*phiy*(phixx+phizz)+phiz*phiz*(phixx+phiyy))/(2.0*pow(phix*phix+phiy*phiy+phiz*phiz,1.5))-(phix*phiy*phixy+phix*phiz*phixz+phiy*phiz*phiyz)/(pow(phix*phix+phiy*phiy+phiz*phiz,1.5));
//
//     k_k=2*(phix*phiy*(phixz*phiyz-phixy*phizz)+phix*phiz*(phixy*phiyz-phixz*phiyy)+phiy*phiz*(phixy*phixz-phiyz*phixx))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0)+(phix*phix*(phiyy*phizz-phiyz*phiyz)+phiy*phiy*(phixx*phizz-phixz*phixz)+phiz*phiz*(phixx*phiyy-phixy*phixy))/pow(phix*phix+phiy*phiy+phiz*phiz,2.0);
//
// 	printf("mean and gaussian curvatures:%f and %f\n",k_h,k_k);
//
//     k_1=k_h-sqrt(abs(k_h*k_h-k_k));
//     k_2=k_h+sqrt(abs(k_h*k_h-k_k));
//
// 	printf("K1 and K2 are: %f and %f\n",k_1*parameters[10],k_2*parameters[10]);
    
	return 0;
}


void initialize_vector(vector<double> & v){
    for(int i=0;i<v.size();i++){
        v[i]=0.0;
    }
}

void print_vector(vector<double> & to_be_print){
    for (int i=0;i<to_be_print.size();i++){
        cout<<to_be_print[i]<<'\n';
    }
}

void read_CFG(ifstream &file){
    string line;
    int i,j,c;
    int entry;
    double rc=10;//cut off
    int Lxyz;
    int mc[3];//single atom location
    vector<int> info{2,0};
    long n;
    double x,y,z;
    char *pch;
    char *pEnd;
    char *pEnd_help;
    char *line_in_char;
    cout<<"file opened\n";
    getline(file,line);
    line_in_char=new char [line.length()+1];
    strcpy(line_in_char,line.c_str());
    pch=strtok(line_in_char,"=");
    pch=strtok(NULL," ");
    n=strtol(pch,&pEnd,10);
    cout<<"total atoms: "<<n<<'\n';
    
    getline(file,line);
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            getline(file,line);
            line_in_char=new char [line.length()+1];
            strcpy(line_in_char,line.c_str());
            pch=strtok(line_in_char,"=");
            pch=strtok(NULL," ");
            h[j][i]=strtod(pch,&pEnd);
        }
    }
    
    L[0]=int(h[0][0]/rc);
    L[1]=int(h[1][1]/rc);
    L[2]=int(h[2][2]/rc);
    
    rc_alpha[0]=h[0][0]/L[0];
    rc_alpha[1]=h[1][1]/L[1];
    rc_alpha[2]=h[2][2]/L[2];

    getline(file,line);
    getline(file,line);
    line_in_char=new char [line.length()+1];
    strcpy(line_in_char,line.c_str());
    strcpy(line_in_char,line.c_str());
    pch=strtok(line_in_char,"=");
    pch=strtok(NULL," ");
    entry=atoi(pch);
    cout<<"Total column: "<<entry<<'\n';
    
    for(i=0;i<entry-3;i++){
        getline(file,line);
        cout<<line<<'\n';
    }
    
    while(getline(file,line)){
        //cout<<"get the line: "<<line<<'\n';
        x=strtod(line.c_str(),&pEnd);
        if(!strtod(pEnd,&pEnd_help)&&!strtod(pEnd_help,&pEnd_help)&&!strtod(pEnd_help,&pEnd_help)){
            getline(file,line);
            //cout<<"jumped "<<line<<'\n';
        }
        
        else{
            y=strtod(pEnd, &pEnd);
            z=strtod(pEnd, &pEnd);
//            cout<<"x y z"<<x<<' '<<y<<' '<<z<<'\n';
//            exit(0);
        //    printf("%f %f %f to %f %f %f\n",x,y,z,coordinates[0],coordinates[1],coordinates[2]);
            info[0]=int(strtod(pEnd,&pEnd));
            info[1]=int(strtod(pEnd,&pEnd));
            AtomPosition.push_back({x*h[0][0]+y*h[1][0]+z*h[2][0],x*h[0][1]+y*h[1][1]+z*h[2][1],x*h[0][2]+y*h[1][2]+z*h[2][2]});
            AtomInfo.push_back(info);
            i++;
            //cout<<info[0]<<' '<<info[1]<<'\n';
        }
    }
    //exit(0);
    //================Create Atom Linked list================================
    head=new int[L[0]*L[1]*L[2]];
    lscl=new int[AtomPosition.size()];
    for (i=0;i<L[0]*L[1]*L[2];i++){
        head[i]=-1;//-1 as EMPTY
    }
    for(i=0;i<AtomPosition.size();i++){
        lscl[i]=-1;
    }
    
    for(i=0;i<AtomPosition.size();i++){
        for(j=0;j<3;j++) mc[j]=AtomPosition[i][j]/rc_alpha[j];
        c=mc[0]*L[1]*L[2]+mc[1]*L[2]+mc[2];
        lscl[i]=head[c];
        head[c]=i;
    }
    cout<<"Linked list created\n";
}

void read_surface_cfg(ifstream &file){
    string line;
    int i,j,c,a;
    int entry;
    int mc[3];//single atom location
    vector<int> info{2,0};
    long n;
    double x,y,z;
    char *pch;
    char *pEnd;
    char *pEnd_help;
    char *line_in_char;
    cout<<"Surface file opened\n";
    getline(file,line);
    line_in_char=new char [line.length()+1];
    strcpy(line_in_char,line.c_str());
    pch=strtok(line_in_char,"=");
    pch=strtok(NULL," ");
    n=strtol(pch,&pEnd,10);
    cout<<"total surface atoms: "<<n<<'\n';
	
    
    getline(file,line);
    for(i=0;i<9;i++){
    	getline(file,line);
    }
	
    getline(file,line);
    getline(file,line);
    line_in_char=new char [line.length()+1];
    strcpy(line_in_char,line.c_str());
    strcpy(line_in_char,line.c_str());
    pch=strtok(line_in_char,"=");
    pch=strtok(NULL," ");
    entry=atoi(pch);
    cout<<"Total column in Surface file: "<<entry<<'\n';
	
    for(i=0;i<entry-3;i++){
        getline(file,line);
        cout<<line<<'\n';
    }
	
    while(getline(file,line)){
        //cout<<"get the line: "<<line<<'\n';
        x=strtod(line.c_str(),&pEnd);
        if(!strtod(pEnd,&pEnd_help)&&!strtod(pEnd_help,&pEnd_help)&&!strtod(pEnd_help,&pEnd_help)){
            getline(file,line);
            //cout<<"jumped "<<line<<'\n';
        }
        
        else{
            y=strtod(pEnd, &pEnd);
            z=strtod(pEnd, &pEnd);
//            cout<<"x y z"<<x<<' '<<y<<' '<<z<<'\n';
//            exit(0);
        //    printf("%f %f %f to %f %f %f\n",x,y,z,coordinates[0],coordinates[1],coordinates[2]);
            info[0]=int(strtod(pEnd,&pEnd));
            info[1]=int(strtod(pEnd,&pEnd));
            SurfaceAtomPosition.push_back({x*h[0][0]+y*h[1][0]+z*h[2][0],x*h[0][1]+y*h[1][1]+z*h[2][1],x*h[0][2]+y*h[1][2]+z*h[2][2]});
            SurfaceAtomInfo.push_back(info);
            i++;
            //cout<<info[0]<<' '<<info[1]<<'\n';
        }
    }
	
	//=======create surface atom linked list=========
    shead=new int[L[0]*L[1]*L[2]];
    slscl=new int[AtomPosition.size()];
    for(i=0;i<L[0]*L[1]*L[2];i++) shead[i]=-1;
    for(i=0;i<SurfaceAtomPosition.size();i++) slscl[i]=-1;
    for(i=0;i<SurfaceAtomPosition.size();i++){
        for(a=0;a<3;a++) mc[a]=SurfaceAtomPosition[i][a]/rc_alpha[a];
        c=mc[0]*L[1]*L[2]+mc[1]*L[2]+mc[2];
        slscl[i]=shead[c];
        shead[c]=i;
    }
	cout<<"Surface Atoms Linked List Created\n";
	
}
void read_internal_cfg(ifstream &file){
    string line;
    int i,j,c,a;
    int entry;
    int mc[3];//single atom location
    vector<int> info{2,0};
    long n;
    double x,y,z;
    char *pch;
    char *pEnd;
    char *pEnd_help;
    char *line_in_char;
    cout<<"Internal file opened\n";
    getline(file,line);
    line_in_char=new char [line.length()+1];
    strcpy(line_in_char,line.c_str());
    pch=strtok(line_in_char,"=");
    pch=strtok(NULL," ");
    n=strtol(pch,&pEnd,10);
    cout<<"total internal atoms: "<<n<<'\n';
	
    
    getline(file,line);
    for(i=0;i<9;i++){
    	getline(file,line);
    }
	
    getline(file,line);
    getline(file,line);
    line_in_char=new char [line.length()+1];
    strcpy(line_in_char,line.c_str());
    strcpy(line_in_char,line.c_str());
    pch=strtok(line_in_char,"=");
    pch=strtok(NULL," ");
    entry=atoi(pch);
    cout<<"Total column in Internal file: "<<entry<<'\n';
	
    for(i=0;i<entry-3;i++){
        getline(file,line);
        cout<<line<<'\n';
    }
	
    while(getline(file,line)){
        //cout<<"get the line: "<<line<<'\n';
        x=strtod(line.c_str(),&pEnd);
        if(!strtod(pEnd,&pEnd_help)&&!strtod(pEnd_help,&pEnd_help)&&!strtod(pEnd_help,&pEnd_help)){
            getline(file,line);
            //cout<<"jumped "<<line<<'\n';
        }
        
        else{
            y=strtod(pEnd, &pEnd);
            z=strtod(pEnd, &pEnd);
//            cout<<"x y z"<<x<<' '<<y<<' '<<z<<'\n';
//            exit(0);
        //    printf("%f %f %f to %f %f %f\n",x,y,z,coordinates[0],coordinates[1],coordinates[2]);
            info[0]=int(strtod(pEnd,&pEnd));
            info[1]=int(strtod(pEnd,&pEnd));
            InternalAtomPosition.push_back({x*h[0][0]+y*h[1][0]+z*h[2][0],x*h[0][1]+y*h[1][1]+z*h[2][1],x*h[0][2]+y*h[1][2]+z*h[2][2]});
            InternalAtomInfo.push_back(info);
            i++;
            //cout<<info[0]<<' '<<info[1]<<'\n';
        }
    }
	
	//=======create surface atom linked list=========
    inhead=new int[L[0]*L[1]*L[2]];
    inlscl=new int[InternalAtomPosition.size()];
    for(i=0;i<L[0]*L[1]*L[2];i++) inhead[i]=-1;
    for(i=0;i<InternalAtomPosition.size();i++) inlscl[i]=-1;
    for(i=0;i<InternalAtomPosition.size();i++){
        for(a=0;a<3;a++) mc[a]=InternalAtomPosition[i][a]/rc_alpha[a];
        c=mc[0]*L[1]*L[2]+mc[1]*L[2]+mc[2];
        inlscl[i]=inhead[c];
        inhead[c]=i;
    }
	cout<<"Internal Atoms Linked List Created\n";
	
}

void find_nearest_neighbor(vector<int> &neighbor_list,int num, double r,vector<vector<double> > &AtomPosition,int *head,int *lscl){
    neighbor_list.clear();
    int mc[3];
    int center,c;
    int i;
    int dmx,dmy,dmz,mx,my,mz;
    double x,y,z,distance;

    
    for(i=0;i<3;i++) mc[i]=AtomPosition[num][i]/rc_alpha[i];
    center=mc[0]*L[1]*L[2]+mc[1]*L[2]+mc[2];
    for(dmx=-1;dmx<2;dmx++){
        mx=mc[0]+dmx;
        if(mx<0){mx=L[0]-1;}
        else if(mx==L[0]){
            mx=0;
        }
        for(dmy=-1;dmy<2;dmy++){
            my=mc[1]+dmy;
            if(my<0){my=L[1]-1;}
            else if(mx==L[1]){
                my=0;
            }
            for(dmz=-1;dmz<2;dmz++){
                mz=mc[2]+dmz;
                if(mz<0){mx=L[2]-1;}
                else if(mx==L[2]){
                    mz=0;
                }
                c=mx*L[1]*L[2]+my*L[2]+mz;
                if(head[c]==-1) continue;
                if(c==center){
                    i=head[c];
                    while(lscl[i]!=-1){
                        if(i!=num){
                            x=min(abs(AtomPosition[num][0]-AtomPosition[i][0]),h[0][0]-abs(AtomPosition[num][0]-AtomPosition[i][0]));
                            y=min(abs(AtomPosition[num][1]-AtomPosition[i][1]),h[1][1]-abs(AtomPosition[num][1]-AtomPosition[i][1]));
                            z=min(abs(AtomPosition[num][2]-AtomPosition[i][2]),h[2][2]-abs(AtomPosition[num][2]-AtomPosition[i][2]));
                            distance=sqrt(x*x+y*y+z*z);
                            if (distance<r) neighbor_list.push_back(i);
                        }
                        i=lscl[i];
                    }
                    if(i!=num){
                        x=min(abs(AtomPosition[num][0]-AtomPosition[i][0]),h[0][0]-abs(AtomPosition[num][0]-AtomPosition[i][0]));
                        y=min(abs(AtomPosition[num][1]-AtomPosition[i][1]),h[1][1]-abs(AtomPosition[num][1]-AtomPosition[i][1]));
                        z=min(abs(AtomPosition[num][2]-AtomPosition[i][2]),h[2][2]-abs(AtomPosition[num][2]-AtomPosition[i][2]));
                        distance=sqrt(x*x+y*y+z*z);
                        if (distance<r) neighbor_list.push_back(i);
                    }
                    
                }
                else{
                    i=head[c];
                    while(lscl[i]!=-1){
                        x=min(abs(AtomPosition[num][0]-AtomPosition[i][0]),h[0][0]-abs(AtomPosition[num][0]-AtomPosition[i][0]));
                        y=min(abs(AtomPosition[num][1]-AtomPosition[i][1]),h[1][1]-abs(AtomPosition[num][1]-AtomPosition[i][1]));
                        z=min(abs(AtomPosition[num][2]-AtomPosition[i][2]),h[2][2]-abs(AtomPosition[num][2]-AtomPosition[i][2]));
                        distance=sqrt(x*x+y*y+z*z);
                        if (distance<r) neighbor_list.push_back(i);
                        i=lscl[i];
                    }
                    x=min(abs(AtomPosition[num][0]-AtomPosition[i][0]),h[0][0]-abs(AtomPosition[num][0]-AtomPosition[i][0]));
                    y=min(abs(AtomPosition[num][1]-AtomPosition[i][1]),h[1][1]-abs(AtomPosition[num][1]-AtomPosition[i][1]));
                    z=min(abs(AtomPosition[num][2]-AtomPosition[i][2]),h[2][2]-abs(AtomPosition[num][2]-AtomPosition[i][2]));
                    distance=sqrt(x*x+y*y+z*z);
                    if (distance<r) neighbor_list.push_back(i);
                }
            }
        }
    }
}

void write_surface_cfg(){
    int i,j;
    ofstream myfile;
    myfile.open("surface_atoms.cfg");
    myfile<<"Number of particles = "<<surface_atoms.size()<<'\n';
    myfile<<"A = 1.0 Angstrom (basic length-scale)\n";
    for (i=0;i<3;i++){
        for(j=0;j<3;j++){
            myfile<<"H0("<<j+1<<','<<i+1<<") = "<<h[j][i]<<" A\n";
        }
    }
    myfile<<".NO_VELOCITY.\n";
    myfile<<"entry_count = 5\n"; //Information Just Type and Grain ID
    //What is element type, type 1 for what? Type 2 for what?
    myfile<<"auxiliary[0] = type\n";
    myfile<<"auxiliary[1] = id\n";
    if (AtomInfo[surface_atoms[0]][0]==1){
        myfile<<196.0<<'\n'<<"Au\n";
    }
    else{
        myfile<<190.0<<'\n'<<"Zr\n";
    }
    myfile<<AtomPosition[surface_atoms[0]][0]/h[0][0]<<' '<<AtomPosition[surface_atoms[0]][1]/h[1][1]<<' '<<AtomPosition[surface_atoms[0]][2]/h[2][2]<<' '<<"1 "<<AtomInfo[surface_atoms[0]][0]<<'\n';
    
    for(i=1;i<surface_atoms.size();i++){
        if (AtomInfo[surface_atoms[i]][0]!=AtomInfo[surface_atoms[i-1]][0]){
            if(AtomInfo[surface_atoms[i]][0]==1){
                myfile<<196.0<<'\n'<<"Au\n";
            }
            else{
                myfile<<190.0<<'\n'<<"Zr\n";
            }
        }
    myfile<<AtomPosition[surface_atoms[i]][0]/h[0][0]<<' '<<AtomPosition[surface_atoms[i]][1]/h[1][1]<<' '<<AtomPosition[surface_atoms[i]][2]/h[2][2]<<' '<<"1 "<<AtomInfo[surface_atoms[i]][0]<<'\n';
    }
    myfile.close();
}

vector<vector<double> > new_points(vector<vector<double> > &old_points, double vertical_axe[3], double new_normal[3]){
    //This function convert coordinates of old_points into new coordination system, rotate by certain angle.
    vector<vector<double> > result;
    //Calculate Rotation
    double v[3];
    double r;
    double s[3];
    double c;
    double R[8],V[8];
    double alpha;
    double Rinv[8];
    vector<double> aux(3,0.0);
    MKL_INT n=3,i;
    MKL_INT info;
    MKL_INT ipiv[n];
    v[0]=vertical_axe[1]*new_normal[2]-vertical_axe[2]*new_normal[1];
    v[1]=vertical_axe[2]*new_normal[0]-vertical_axe[0]*new_normal[2];
    v[2]=vertical_axe[0]*new_normal[1]-vertical_axe[1]*new_normal[0];
    
    r=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    c=vertical_axe[0]*new_normal[0]+vertical_axe[1]*new_normal[1]+vertical_axe[2]*new_normal[2];
    if (abs(c+1)<0.0000001){
        for(i=0;i<old_points.size();i++){
            result.push_back({old_points[i][0],old_points[i][1],-(old_points[i][2])});
        }
        return result;
    }
    

    V[0]=0,V[1]=-v[2],V[2]=v[1],V[3]=v[2],V[4]=0,V[5]=-v[0],V[6]=-v[1],V[7]=v[0],V[8]=0;
    cout<<"Calculate Rotation matrix\n";
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1/(1+c),V,n,V,n,0.0,R,n);
    printf("Matrix Multiplication Complete\n");
    for(i=0;i<9;i++){
        if(i==0||i==4||i==8){
            R[i]=R[i]+V[i]+1;
        }
        else{
            R[i]=R[i]+V[i];
        }
    }
    printf("Rotation Matrix Contructed\n");
    print_matrix("Rotation Matrix",n,n,R,n);
    copy(R,R+9,Rinv);
    info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR,n,n,Rinv,n,ipiv);
    
    if(i>0){
        printf("The diagonal element of the triangular factor of A,\n");
        printf("U(%i,%i) is zero, so that A is singular.\n",info,info);
        return old_points;
    }
        
    printf("R is LU factored with permutation if needed\n");
    print_int_vector("Row permutation info",n,ipiv);
    info=LAPACKE_dgetri(LAPACK_ROW_MAJOR,n,Rinv,n,ipiv);
    printf("Invese Rotation Matrix Obtained");
    print_matrix("Inv(R)",n,n,Rinv,n);
    
    for(i=0;i<old_points.size();i++){
        aux[0]=Rinv[0]*old_points[i][0]+Rinv[1]*old_points[i][1]+Rinv[2]*old_points[i][2];
        aux[1]=Rinv[3]*old_points[i][0]+Rinv[4]*old_points[i][1]+Rinv[5]*old_points[i][2];
        aux[2]=Rinv[6]*old_points[i][0]+Rinv[7]*old_points[i][1]+Rinv[8]*old_points[i][2];
        result.push_back(aux);
    }
    
    printf("New coordinates are found\n");
    return result;
        
}

double diameter(vector<double> &position,double movement){
	vector<double> cur_pos;
	double total_movement=0.0;
	vector<double> normal_vector;
	vector<vector<double>> neighbors;
	int num_old;
	int num_new=0;
	int go_through=0;
	normal_vector=find_normal(position,rc_surface);
	//cout<<"Normal is:"<<normal_vector[0]<<' '<<normal_vector[1]<<' '<<normal_vector[2]<<'\n';
	cur_pos.push_back(position[0]);
	cur_pos.push_back(position[1]);
	cur_pos.push_back(position[2]);
	num_old=find_atoms_around(position,rc_surface,SurfaceAtomPosition,shead,slscl).size();
	num_new=num_old;
	//cout<<"initial number of neighbors: "<<num_old<<'\n';
	while(num_old>=num_new){
		cur_pos[0]+=movement*normal_vector[0];
		cur_pos[1]+=movement*normal_vector[1];
		cur_pos[2]+=movement*normal_vector[2];	
		total_movement+=movement;
		if(cur_pos[0]>h[0][0]) {cur_pos[0]=cur_pos[0]-h[0][0]; go_through=1;}
		if(cur_pos[0]<0) {cur_pos[0]=h[0][0]+cur_pos[0];go_through=1;}
		if(cur_pos[1]>h[1][1]) {cur_pos[1]=cur_pos[1]-h[1][1]; go_through=1;}
		if(cur_pos[1]<0) {cur_pos[1]=h[1][1]+cur_pos[1]; go_through=1;}
		if(cur_pos[2]>h[2][2]) {cur_pos[2]=cur_pos[2]-h[2][2]; go_through=1;}
		if(cur_pos[2]<0) {cur_pos[2]=h[2][2]+cur_pos[2]; go_through=1;}
		//cout<<"cur_pos: "<<cur_pos[0]<<' '<<cur_pos[1]<<' '<<cur_pos[2]<<'\n';
		num_old=num_new;
		num_new=find_atoms_around(cur_pos,rc_surface,SurfaceAtomPosition,shead,slscl).size();
		//cout<<"Total Movement: "<<total_movement<<'\n';
		//cout<<"Neighbors: "<<num_new<<'\n';
	}
	neighbors=find_atoms_around_v2(cur_pos,rc_surface,SurfaceAtomPosition,shead,slscl);
	// for(int i=0;i<neighbors.size();i++){
	// 	cout<<"ParticleIdentifier=="<<SurfaceAtomInfo[(int)neighbors[i][0]][0]<<'\n';
	// }
	//cout<<"go though boundary? "<<go_through<<'\n';
	return total_movement;	
}

vector<double> find_normal(vector<double> &position,double rc){
	vector<vector<double>> neighbors;
	vector<double> normal_vector;
	double norm;
	normal_vector.push_back(0.0);
	normal_vector.push_back(0.0);
	normal_vector.push_back(0.0);
	//cout<<"Finding neighbors...\n";
	neighbors=find_atoms_around(position,rc,InternalAtomPosition,inhead,inlscl);
	//cout<<"Internal Atoms Number: "<<neighbors.size()<<'\n';
	//cout<<"Neighbors found\n";
	if(neighbors.size()==0) return {1,0,0};
	for(int i=0;i<neighbors.size();i++){
		normal_vector[0]+=neighbors[i][0];
		normal_vector[1]+=neighbors[i][1];
		normal_vector[2]+=neighbors[i][2];
	}
	norm=normal_vector[0]*normal_vector[0]+normal_vector[1]*normal_vector[1]+normal_vector[2]*normal_vector[2];
	norm=sqrt(norm);
	normal_vector[0]=normal_vector[0]/norm;
	normal_vector[1]=normal_vector[1]/norm;
	normal_vector[2]=normal_vector[2]/norm;
	//cout<<"Normal is: "<<normal_vector[0]<<' '<<normal_vector[1]<<' '<<normal_vector[2]<<"\n";
	return normal_vector;
	
}

vector<vector<double> > find_atoms_around(vector<double> &position, double rc,vector<vector<double> > &AtomPosition,int* head, int* lscl){
    //Function returns a set of relative coordinates;
	vector<vector<double> > neighbors;
	int mc[3];
    int center,c;
    int i;
    int dmx,dmy,dmz,mx,my,mz;
    double x,y,z,distance;
	double box_position[3];
	double relative_position[3];
	
	//cout<<"rc_alpha : "<<rc_alpha[0]<<' '<<rc_alpha[1]<<' '<<rc_alpha[2]<<'\n';
	//cout<<"Position: "<<position[0]<<' '<<position[1]<<' '<<position[2]<<"\n";
    for(i=0;i<3;i++) mc[i]=position[i]/rc_alpha[i];
    center=mc[0]*L[1]*L[2]+mc[1]*L[2]+mc[2];
	//cout<<"Center box Id: "<<center<<'\n';
	box_position[0]=rc_alpha[0]*mc[0];
	box_position[1]=rc_alpha[1]*mc[1];
	box_position[2]=rc_alpha[2]*mc[2];
	relative_position[0]=position[0]-box_position[0];
	relative_position[1]=position[1]-box_position[1];
	relative_position[2]=position[2]-box_position[2];
	for(dmx=-1;dmx<2;dmx++){
        mx=mc[0]+dmx;
        if(mx<0){mx=L[0]-1;}
        else if(mx==L[0]){
            mx=0;
        }
		for(dmy=-1;dmy<2;dmy++){
	        my=mc[1]+dmy;
	        if(my<0){my=L[1]-1;}
	        else if(my==L[1]){
	            my=0;
	        }
			for(dmz=-1;dmz<2;dmz++){
		        mz=mc[2]+dmz;
		        if(mz<0){mz=L[2]-1;}
		        else if(mz==L[2]){
		            mz=0;
		        }
                c=mx*L[1]*L[2]+my*L[2]+mz;
				//cout<<"Box ID: "<<c<<'\n';
                if(head[c]==-1) continue;
				i=head[c];
				box_position[0]=rc_alpha[0]*mx;
				box_position[1]=rc_alpha[1]*my;
				box_position[2]=rc_alpha[2]*mz;
				while(lscl[i]!=-1){
					x=AtomPosition[i][0]-box_position[0]-relative_position[0]+dmx*rc_alpha[0];
					y=AtomPosition[i][1]-box_position[1]-relative_position[1]+dmy*rc_alpha[1];
					z=AtomPosition[i][2]-box_position[2]-relative_position[2]+dmz*rc_alpha[2];
                    distance=sqrt(x*x+y*y+z*z);
                    if (distance<rc) neighbors.push_back({x,y,z});
					i=lscl[i];
				}
				x=AtomPosition[i][0]-box_position[0]-relative_position[0]+dmx*rc_alpha[0];
				y=AtomPosition[i][1]-box_position[1]-relative_position[1]+dmy*rc_alpha[1];
				z=AtomPosition[i][2]-box_position[2]-relative_position[2]+dmz*rc_alpha[2];
                distance=sqrt(x*x+y*y+z*z);
                if (distance<rc) neighbors.push_back({x,y,z});
			}
		}
	}
	return neighbors;
}

vector<vector<double> > find_atoms_around_v2(vector<double> &position, double rc,vector<vector<double> > &AtomPosition,int* head, int* lscl){//return also the number of each atoms, saved in 0 position
	vector<vector<double> > neighbors;
	int mc[3];
    int center,c;
    int i;
    int dmx,dmy,dmz,mx,my,mz;
    double x,y,z,distance;
	double box_position[3];
	double relative_position[3];
    for(i=0;i<3;i++) mc[i]=position[i]/rc_alpha[i];
    center=mc[0]*L[1]*L[2]+mc[1]*L[2]+mc[2];
	//cout<<"Center box Id: "<<center<<'\n';
	box_position[0]=rc_alpha[0]*mc[0];
	box_position[1]=rc_alpha[1]*mc[1];
	box_position[2]=rc_alpha[2]*mc[2];
	relative_position[0]=position[0]-box_position[0];
	relative_position[1]=position[1]-box_position[1];
	relative_position[2]=position[2]-box_position[2];
	for(dmx=-1;dmx<2;dmx++){
        mx=mc[0]+dmx;
        if(mx<0){mx=L[0]-1;}
        else if(mx>=L[0]){
            mx=0;
        }
		for(dmy=-1;dmy<2;dmy++){
	        my=mc[1]+dmy;
	        if(my<0){my=L[1]-1;}
	        else if(my>=L[1]){
	            my=0;
	        }
			for(dmz=-1;dmz<2;dmz++){
		        mz=mc[2]+dmz;
		        if(mz<0){mz=L[2]-1;}
		        else if(mz>=L[2]){
		            mz=0;
		        }
                c=mx*L[1]*L[2]+my*L[2]+mz;
				//cout<<"mx,my,mz:"<<mx<<' '<<my<<' '<<mz<<'\n';
				//cout<<"Box ID: "<<c<<'\n';
                if(head[c]==-1) continue;
				i=head[c];
				box_position[0]=rc_alpha[0]*mx;
				box_position[1]=rc_alpha[1]*my;
				box_position[2]=rc_alpha[2]*mz;
				while(lscl[i]!=-1){
					x=AtomPosition[i][0]-box_position[0]-relative_position[0]+dmx*rc_alpha[0];
					y=AtomPosition[i][1]-box_position[1]-relative_position[1]+dmy*rc_alpha[1];
					z=AtomPosition[i][2]-box_position[2]-relative_position[2]+dmz*rc_alpha[2];
                    distance=sqrt(x*x+y*y+z*z);
                    if (distance<rc) neighbors.push_back({double(i)+0.000001,x+position[0],y+position[1],z+position[2]});
					i=lscl[i];
				}
				x=AtomPosition[i][0]-box_position[0]-relative_position[0]+dmx*rc_alpha[0];
				y=AtomPosition[i][1]-box_position[1]-relative_position[1]+dmy*rc_alpha[1];
				z=AtomPosition[i][2]-box_position[2]-relative_position[2]+dmz*rc_alpha[2];
                distance=sqrt(x*x+y*y+z*z);
                if (distance<rc) neighbors.push_back({double(i)+0.000001,x+position[0],y+position[1],z+position[2]});
			}
		}
	}
	return neighbors;
	
}

void print_int_vector(char *desc,MKL_INT n, MKL_INT *a){
    MKL_INT j;
    printf("\n %s \n",desc);
    for (j=0;j<n;j++) printf(" %6i",a[j]);
    printf( "\n\n");
}

void print_ge_vector(char *name, MKL_INT n,double *b){
    MKL_INT j;
    printf("\n %s \n",name);
    for (j=0;j<n;j++) printf("%10.2f ",b[j]);
    printf("\n\n");
}

void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %12.10f", a[i*lda+j] );
                printf( "\n" );
        }
        printf("\n");
}

vector<double> fitting_surface(vector<vector<double> > &point_list){
    //fitting a bivariate surface with given neighbor points;
    //z=ax*x+by*y+cx*y+dx+ey+f
    //cout<<"Sovle Parameters using lapacke_dgesv\n";
    vector<double> result;
    MKL_INT i;
    MKL_INT length=point_list.size();
    double b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    MKL_INT n=6,nrhs=1,lda=6,ldb=1,info;
    MKL_INT ipiv[n];
    double A[lda*n];
    double B[ldb*n];
	double cost=0,difference;
    for(i=0;i<length;i++){
        b0+=point_list[i][0]*point_list[i][0]*point_list[i][2];//x_square[i]*z[i];
        b1+=point_list[i][1]*point_list[i][1]*point_list[i][2];//y_square[i]*z[i];
        b2+=point_list[i][0]*point_list[i][1]*point_list[i][2];//xy[i]*z[i];
        b3+=point_list[i][0]*point_list[i][2];//x[i]*z[i];
        b4+=point_list[i][1]*point_list[i][2];//y[i]*z[i];
        b5+=point_list[i][2];//z[i];
    }
    B[0]=b0,B[1]=b1,B[2]=b2,B[3]=b3,B[4]=b4,B[5]=b5;
    //cout<<"Vector B prepared and it is \n";
    //print_ge_vector("Vector B",n,B);

    b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    for(i=0;i<length;i++){
        b0+=point_list[i][0]*point_list[i][0]*point_list[i][0]*point_list[i][0];//x_square[i]*x_square[i];
        b1+=point_list[i][0]*point_list[i][0]*point_list[i][1]*point_list[i][1];//x_square[i]*y_square[i];
        b2+=point_list[i][0]*point_list[i][0]*point_list[i][0]*point_list[i][1];//x_square[i]*xy[i];
        b3+=point_list[i][0]*point_list[i][0]*point_list[i][0];//x_square[i]*x[i];
        b4+=point_list[i][0]*point_list[i][0]*point_list[i][1];//x_square[i]*y[i];
        b5+=point_list[i][0]*point_list[i][0];//x_square[i];
    }
    A[0]=b0,A[1]=b1,A[2]=b2;A[3]=b3,A[4]=b4,A[5]=b5;

    b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    for(i=0;i<length;i++){
        b0+=point_list[i][1]*point_list[i][1]*point_list[i][0]*point_list[i][0];
        b1+=point_list[i][1]*point_list[i][1]*point_list[i][1]*point_list[i][1];
        b2+=point_list[i][1]*point_list[i][1]*point_list[i][0]*point_list[i][1];
        b3+=point_list[i][1]*point_list[i][1]*point_list[i][0];
        b4+=point_list[i][1]*point_list[i][1]*point_list[i][1];
        b5+=point_list[i][1]*point_list[i][1];
//        b0+=y_square[i]*x_square[i];
//        b1+=y_square[i]*y_square[i];
//        b2+=y_square[i]*xy[i];
//        b3+=y_square[i]*x[i];
//        b4+=y_square[i]*y[i];
//        b5+=y_square[i];
    }
    A[6]=b0,A[7]=b1,A[8]=b2;A[9]=b3,A[10]=b4,A[11]=b5;

    b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    for(i=0;i<length;i++){
        b0+=point_list[i][0]*point_list[i][1]*point_list[i][0]*point_list[i][0];
        b1+=point_list[i][0]*point_list[i][1]*point_list[i][1]*point_list[i][1];
        b2+=point_list[i][0]*point_list[i][1]*point_list[i][0]*point_list[i][1];
        b3+=point_list[i][0]*point_list[i][1]*point_list[i][0];
        b4+=point_list[i][0]*point_list[i][1]*point_list[i][1];
        b5+=point_list[i][0]*point_list[i][1];
//        b0+=xy[i]*x_square[i];
//        b1+=xy[i]*y_square[i];
//        b2+=xy[i]*xy[i];
//        b3+=xy[i]*x[i];
//        b4+=xy[i]*y[i];
//        b5+=xy[i];
    }
    A[12]=b0,A[13]=b1,A[14]=b2;A[15]=b3,A[16]=b4,A[17]=b5;

    b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    for(i=0;i<length;i++){
        b0+=point_list[i][0]*point_list[i][0]*point_list[i][0];
        b1+=point_list[i][0]*point_list[i][1]*point_list[i][1];
        b2+=point_list[i][0]*point_list[i][0]*point_list[i][1];
        b3+=point_list[i][0]*point_list[i][0];
        b4+=point_list[i][0]*point_list[i][1];
        b5+=point_list[i][0];
//        b0+=x[i]*x_square[i];
//        b1+=x[i]*y_square[i];
//        b2+=x[i]*xy[i];
//        b3+=x[i]*x[i];
//        b4+=x[i]*y[i];
//        b5+=x[i];
    }
    A[18]=b0,A[19]=b1,A[20]=b2;A[21]=b3,A[22]=b4,A[23]=b5;

    b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    for(i=0;i<length;i++){
        b0+=point_list[i][1]*point_list[i][0]*point_list[i][0];
        b1+=point_list[i][1]*point_list[i][1]*point_list[i][1];
        b2+=point_list[i][1]*point_list[i][0]*point_list[i][1];
        b3+=point_list[i][1]*point_list[i][0];
        b4+=point_list[i][1]*point_list[i][1];
        b5+=point_list[i][1];
//        b5+=point_list[i][0];
//        b0+=y[i]*x_square[i];
//        b1+=y[i]*y_square[i];
//        b2+=y[i]*xy[i];
//        b3+=y[i]*x[i];
//        b4+=y[i]*y[i];
//        b5+=y[i];
    }
    A[24]=b0,A[25]=b1,A[26]=b2;A[27]=b3,A[28]=b4,A[29]=b5;

    b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    for(i=0;i<length;i++){
        b0+=point_list[i][0]*point_list[i][0];
        b1+=point_list[i][1]*point_list[i][1];
        b2+=point_list[i][0]*point_list[i][1];
        b3+=point_list[i][0];
        b4+=point_list[i][1];
        b5+=1;
//        b0+=x_square[i];
//        b1+=y_square[i];
//        b2+=xy[i];
//        b3+=x[i];
//        b4+=y[i];
//        b5+=1;
    }
    A[30]=b0,A[31]=b1,A[32]=b2;A[33]=b3,A[34]=b4,A[35]=b5;
    //cout<<"A value assigned and is \n";
    //print_matrix("Initial A",n,n,A,lda);
    //cout<<"LAPACKE_degsv (row-major, high-level) calculation result\n";

    info=LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,nrhs,A,lda,ipiv,B,ldb);

    if(info>0){
        printf("The diagonal element of the triangular factor of A,\n");
        printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
        printf("No solution obtained\n");
        return {};
    }
    //print_matrix("Solution",n,nrhs,B,ldb);
    //print_matrix("Details of LU decomposition",n,n,A,lda);
    //print_int_vector("Pivot indices",n,ipiv);
    for(i=0;i<n;i++){
        result.push_back(B[i]);
    }
	for(i=0;i<length;i++){
		difference=(B[0]*point_list[i][0]*point_list[i][0]+B[1]*point_list[i][1]*point_list[i][1]+B[2]*point_list[i][0]*point_list[i][1]+B[3]*point_list[i][0]+B[4]*point_list[i][1]+B[5]-point_list[i][2]);
		cost+=difference*difference;
	}
	cost=sqrt(cost)/length;
	cout<<"Remaining Cost: "<<cost<<'\n';
	if (cost>fit_threshold) return {};
    return result;
}

vector<double> fitting_qudratic_surface(vector<vector<double> > &point_list){
	double scale;
	double actual_max_x,actual_max_y,actual_max_z;
	double actual_min_x,actual_min_y,actual_min_z;
	double max_range;
	double x_mid,y_mid,z_mid;
	double target_range=4;
    MKL_INT i,j,k=1;
	MKL_INT info;
	MKL_INT nr=10,ldr=10;
	MKL_INT ldp=1,np=1,mp=10;
	MKL_INT ldK_2=6,nK_2=6;
	MKL_INT ldc=6,nc=6;
	MKL_INT ldb=4,nb=6;
	MKL_INT lda=4,na=4;
	MKL_INT ldm_prime=6,nm_prime=6;
	MKL_INT ldm=6,nm=6;
    double R[nr*nr];
	double p[np*mp];
	double C[ldc*nc];
	double B[ldb*nb];
	double A[lda*na];
	double Ainv[lda*na];
	double alpha[lda];
	double beta[nb];
	double beta_prime[nb];
	double K_2[nK_2*ldK_2];
	double M_prime[ldm_prime,nm_prime];
	double M[ldm*nm];
	double Hinv[6*6];
	double AinvBtrans[4*6];
	double MHinv[6*6];
	double d[nm_prime],e[nm_prime-1],tau[nm_prime-1];
	double TM_prime[3*nm_prime];
	double Z[nm_prime*nm_prime];
	double W[nm_prime];
	vector<double> parameters(11,0.0);
	double factor;
	MKL_INT min_tag;

	MKL_INT ipiv[lda];
	
	actual_max_x=point_list[0][0];
	actual_min_x=point_list[0][0];
	actual_max_y=point_list[0][1];
	actual_min_y=point_list[0][1];
	actual_max_z=point_list[0][2];
	actual_min_z=point_list[0][2];
	
	for(i=0;i<point_list.size();i++){
		actual_max_x=point_list[i][0]>actual_max_x?point_list[i][0]:actual_max_x;
		actual_min_x=point_list[i][0]<actual_min_x?point_list[i][0]:actual_min_x;
		actual_max_y=point_list[i][1]>actual_max_y?point_list[i][1]:actual_max_y;
		actual_min_y=point_list[i][1]<actual_min_y?point_list[i][1]:actual_min_y;
		actual_max_z=point_list[i][2]>actual_max_z?point_list[i][2]:actual_max_z;
		actual_min_z=point_list[i][2]<actual_min_z?point_list[i][2]:actual_min_z;
	}
	
	max_range=actual_max_x-actual_min_x;
	max_range=(actual_max_y-actual_min_y)>max_range?(actual_max_y-actual_min_y):max_range;
	max_range=(actual_max_z-actual_min_z)>max_range?(actual_max_z-actual_min_z):max_range;
	
	scale=max_range/(target_range);
	x_mid=(actual_max_x-actual_min_x)/2;
	y_mid=(actual_max_y-actual_min_y)/2;
	z_mid=(actual_max_z-actual_min_z)/2;
	
	for(i=0;i<point_list.size();i++){
		point_list[i][0]=(point_list[i][0]-x_mid)/scale;
		point_list[i][1]=(point_list[i][1]-y_mid)/scale;
		point_list[i][2]=(point_list[i][2]-z_mid)/scale;
	}
	
	
	
	double x,y,z;
	for(i=0;i<nr*nr;i++){
		R[i]=0;
	}
	for(i=0;i<nK_2*ldK_2;i++){
		K_2[i]=0;
	}
	
	for(i=0;i<point_list.size();i++){
		x=point_list[i][0];
		y=point_list[i][1];
		z=point_list[i][2];
		p[0]=x*x,p[1]=y*y,p[2]=z*z,p[3]=y*z,p[4]=z*x,p[5]=x*y,p[6]=x,p[7]=y,p[8]=z,p[9]=1;
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,mp,mp,np,1,p,ldp,p,ldp,1,R,ldr);
	}
	
	print_matrix("R",nr,ldr,R,ldr);
	K_2[0]=1,K_2[7]=1,K_2[14]=1,K_2[21]=0.5,K_2[28]=0.5,K_2[35]=0.5;
	print_matrix("K_2",nK_2,ldK_2,K_2,ldK_2);
	
	
	for(i=0;i<nc;i++){
		for(j=0;j<ldc;j++){
			C[i*ldc+j]=R[i*ldr+j];
		}
	}
	
	for(i=0;i<nb;i++){
		for(j=0;j<ldb;j++){
			B[i*ldb+j]=R[i*ldr+j+ldc];
		}
	}
	
	for(i=0;i<na;i++){
		for(j=0;j<lda;j++){
			A[i*lda+j]=R[(i+nb)*ldr+j+ldc];
		}
	}
	
	//solving for beta;
	//Get inverse of A;
	copy(A,A+lda*na,Ainv);
	info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR,na,na,Ainv,lda,ipiv);
	if(info>0){
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular.\n",info,info);
		printf("consider if the surface is planar\n");
	}
	info=LAPACKE_dgetri(LAPACK_ROW_MAJOR,na,Ainv,na,ipiv);
	print_matrix("Inverse of A",na,lda,Ainv,lda);
	
	
	copy(C,C+nc*ldc,M);
	print_matrix("C before",6,6,C,6);
	print_matrix("M before",6,6,M,6);
	for(i=0;i<4*6;i++){
		AinvBtrans[i]=0;
	}
	
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,na,nb,lda,1,Ainv,lda,B,ldb,0,AinvBtrans,6);
	print_matrix("AinvBtrans",4,6,AinvBtrans,6);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nb,6,4,-1,B,ldb,AinvBtrans,6,1,M,ldm);
	print_matrix("M",6,6,M,6);
	
	for(i=0;i<36;i++){
		Hinv[i]=0;
	}
	
	Hinv[0]=1,Hinv[7]=1,Hinv[14]=1,Hinv[21]=1.414213562373095,Hinv[28]=1.414213562373095,Hinv[35]=1.414213562373095;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,6,6,6,1,M,6,Hinv,6,0,MHinv,6);
	print_matrix("MHinv",6,6,MHinv,6);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,6,6,6,1,Hinv,6,MHinv,6,0,M_prime,6);
	print_matrix("M_prime",6,6,M_prime,6);
	info=LAPACKE_dsytrd(LAPACK_ROW_MAJOR,'U',6,M_prime,6,d,e,tau);
	if(info<0){
		printf("The %d th parameter has illegal value",-info);
		exit(0);
	}
	// print_matrix("Matrix after decomposition",n,n,A,lda);
	// print_ge_vector("Diagonal elements",n,d);
	// print_ge_vector("Superdiagonal elements",n-1,e);
	// print_ge_vector("Elementary reflector",n-1,tau);
	// printf("QTQ decomposition done\n");
	
	
	for(i=0;i<nm_prime;i++){
		for(j=i;j<i+2&&j<nm_prime;j++){
			TM_prime[(1+i-j)*nm_prime+j]=M_prime[i*nm_prime+j];
		}
	}
	print_matrix("Band of M",6,6,M_prime,6);
	
	info=LAPACKE_dsbevd(LAPACK_ROW_MAJOR,'V','U',nm_prime,1,TM_prime,ldm_prime,W,Z,ldm_prime);
	info=LAPACKE_dormtr(LAPACK_ROW_MAJOR,'L','U','N',nm_prime,nm_prime,M_prime,ldm_prime,tau,Z,nm_prime);
	print_ge_vector("Eigenvalues of M_prime", 6,W);
	print_matrix("Eigenvectors of M_prime",6,6,Z,6);
	
	
	min_tag=0;
	for(i=0;i<nm_prime;i++){
		if (abs(W[i])<abs(W[min_tag])){
			min_tag=i;
		}
	}
	
	for(i=0;i<nb;i++){
		beta_prime[i]=Z[i*nm_prime+min_tag];
	}
	print_ge_vector("Chosen Beta_prime",6,beta_prime);
	
	
	for(i=0;i<nb;i++){
		beta[i]=Hinv[i*nm_prime+i]*beta_prime[i];
	}
	print_ge_vector("Beta",6,beta);
	
	
	for(i=0;i<lda;i++){
		alpha[i]=0;
	}
	
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,4,1,6,-1,AinvBtrans,6,beta,1,0,alpha,1);
	print_ge_vector("Alpha",lda,alpha);
	
	parameters[6]=alpha[0],parameters[7]=alpha[1],parameters[8]=alpha[2],parameters[9]=alpha[3];
	parameters[0]=beta[0],parameters[1]=beta[1],parameters[2]=beta[2],parameters[3]=beta[3],parameters[4]=beta[4],parameters[5]=beta[5];
	
	printf("Parameters:%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6],parameters[7],parameters[8],parameters[9]);
	
	//Factor = a'Ka
	factor=parameters[0]*parameters[0]+parameters[1]*parameters[1]+parameters[2]*parameters[2]+0.5*parameters[3]*parameters[3]+0.5*parameters[4]*parameters[4]+0.5*parameters[5]*parameters[5];
	printf("Normalization Factor: %f\n",factor);
	parameters[10]=scale;
	
	return parameters;
	
}

void neighbors_for_surface(vector<int> &neighbor_list,int num,double r){
	int i;
	double x,y,z;
	neighbor_list.clear();
	for (i=0;i<surface_atoms.size();i++){
		//cout<<surface_atoms[i]<<'\n';
		if (surface_atoms[i]==num) continue;
		x=abs(AtomPosition[num][0]-AtomPosition[surface_atoms[i]][0]);
		if (x>r) continue;
		y=abs(AtomPosition[num][1]-AtomPosition[surface_atoms[i]][1]);
		if (y>r) continue;
		z=abs(AtomPosition[num][2]-AtomPosition[surface_atoms[i]][2]);
		if (z>r) continue;
		if(sqrt(x*x+y*y+z*z)<r){
			neighbor_list.push_back(surface_atoms[i]);
		}
	}
}

vector<int> neighbor_expansion(int num,int steps,double r,vector<vector<double> > &AtomPosition, int *head, int *lscl){
	vector<int> result;
	vector<int> neighbor_list;
	vector<int> :: iterator it;
	int i,j,p,start=1,end;
	result.push_back(num);
	for (i=0;i<steps;i++){
		cout<<"steps: "<<i<<'\n';
		if(i==0){
			find_nearest_neighbor(neighbor_list,num,r,AtomPosition,head,lscl);
			//cout<<"Number of nearest neighbor:" << neighbor_list.size()<<'\n';
			for (j=0;j<neighbor_list.size();j++){
				result.push_back(neighbor_list[j]);
			}
			end=result.size();
			cout<<"Start and End: "<<start<<' '<<end<<'\n';
		}
		else{
			for(p=start;p<end;p++){
				find_nearest_neighbor(neighbor_list,result[p],r,AtomPosition,head,lscl);
				for(j=0;j<neighbor_list.size();j++){
					it=find(result.begin(),result.end(),neighbor_list[j]);
					if(it==result.end()) {
						result.push_back(neighbor_list[j]);
					}
				}
			}
			start=end;
			end=result.size();
		}	
	}
	return result;
}

vector<vector<double> > neighbor_expansion_v2(vector<double> &position,int steps,double r, vector<vector<double> > &AtomPosition, int *head, int *lscl){
	vector<vector<double> > result;
	vector<vector<double> > ins_memory;
	vector<int> visited;
	vector<int> :: iterator it;
	vector<double> ins_pos;
	
	ins_pos.push_back(0);
	ins_pos.push_back(0);
	ins_pos.push_back(0);
	int i,j,p,start=0,end;
	for(i=0;i<steps;i++){
		//cout<<"Steps "<<i<<'\n';
		if (i==0){
			ins_memory=find_atoms_around_v2(position, r, AtomPosition, head, lscl);
			for (j=0;j<ins_memory.size();j++){
				visited.push_back(int(ins_memory[j][0]));
				result.push_back({ins_memory[j][1],ins_memory[j][2],ins_memory[j][3]});
			}
			end=visited.size();
		}
		else{
			for(p=start;p<end;p++){
				//cout<<"p="<<p<<"\n";
				ins_pos[0]=result[p][0];
				ins_pos[1]=result[p][1];
				ins_pos[2]=result[p][2];
				ins_memory=find_atoms_around_v2(ins_pos,r,AtomPosition,head,lscl);
				//cout<<"ins_memory saved\n";
				for(j=0;j<ins_memory.size();j++){
					//cout<<"ins_memory[0]: "<<ins_memory[j][0]<<'\n';
					it=find(visited.begin(),visited.end(),int(ins_memory[j][0]));
					if(it==visited.end()){
						//cout<<"Not found\n";
						visited.push_back(int(ins_memory[j][0]));
						result.push_back({ins_memory[j][1],ins_memory[j][2],ins_memory[j][3]});
					}
					//cout<<"existed\n";
				}
				//cout<<"all new neighbors checked\n";
			}
			start=end;
			end=visited.size();
		}
		//cout<<"Current numbers: "<<end<<'\n';
		// for(int i=0;i<visited.size();i++){
		// 	cout<<"ParticleIdentifier=="<<SurfaceAtomInfo[(int)visited[i]][0]<<"||\n";
		// }
	}
	return result;
}

vector<double> bivariate_surface_gradient(vector<double>cur_pos,vector<double> &parameters){
	vector<double> result;
	double norm;
	result.push_back(parameters[3]+2*parameters[0]*cur_pos[0]+parameters[2]*cur_pos[1]);
	result.push_back(parameters[4]+2*parameters[1]*cur_pos[1]+parameters[2]*cur_pos[0]);
	result.push_back(1.0);
	norm=result[0]*result[0]+result[1]*result[1]+result[2]*result[2];
	norm=sqrt(norm);
	result[0]=result[0]/norm;
	result[1]=result[1]/norm;
	result[2]=result[2]/norm;
	
	return result;	
}

int consistency_of_normal(vector<double> &cur_pos, vector<double> &normal_direction, double movement){
	vector<double> follow_normal(3,0.0);
	vector<double> anti_normal(3,0.0);
	cout<<find_atoms_around_v2(cur_pos,rc_surface,InternalAtomPosition,inhead,inlscl).size()<<'\n';
	cout<<movement<<'\n';
	if(find_atoms_around_v2(cur_pos,rc_surface,InternalAtomPosition,inhead,inlscl).size()==0||movement>10*rc_surface){
		return 1;
	}
	follow_normal[0]=(cur_pos[0]+movement*normal_direction[0]);
	if(follow_normal[0]>h[0][0]) follow_normal[0]=follow_normal[0]-h[0][0];
	else if(follow_normal[0]<0) follow_normal[0]=h[0][0]+follow_normal[0];
	
	follow_normal[1]=(cur_pos[1]+movement*normal_direction[1]);
	if(follow_normal[1]>h[1][1]) follow_normal[1]=follow_normal[1]-h[1][1];
	else if(follow_normal[1]<0) follow_normal[1]=h[1][1]+follow_normal[1];
	
	follow_normal[2]=(cur_pos[2]+movement*normal_direction[2]);
	if(follow_normal[2]>h[2][2]) follow_normal[2]=follow_normal[2]-h[2][2];
	else if(follow_normal[2]<0) follow_normal[2]=h[2][2]+follow_normal[2];
	
	anti_normal[0]=(cur_pos[0]-movement*normal_direction[0]);
	if(anti_normal[0]>h[0][0]) anti_normal[0]=anti_normal[0]-h[0][0];
	else if(anti_normal[0]<0) anti_normal[0]=h[0][0]+anti_normal[0];
	
	anti_normal[1]=(cur_pos[1]-movement*normal_direction[1]);
	if(anti_normal[1]>h[1][1]) anti_normal[1]=anti_normal[1]-h[1][1];
	else if(anti_normal[1]<0) anti_normal[1]=h[1][1]+anti_normal[1];
	
	anti_normal[2]=(cur_pos[2]-movement*normal_direction[2]);
	if(anti_normal[2]>h[2][2]) anti_normal[2]=anti_normal[2]-h[2][2];
	else if(anti_normal[2]<0) anti_normal[2]=h[2][2]+anti_normal[2];
	
	
	if(find_atoms_around_v2(follow_normal,rc_surface,InternalAtomPosition,inhead,inlscl).size()>find_atoms_around_v2(anti_normal,rc_surface,InternalAtomPosition,inhead,inlscl).size()) return -1;
	else if(find_atoms_around_v2(follow_normal,rc_surface,InternalAtomPosition,inhead,inlscl).size()<find_atoms_around_v2(anti_normal,rc_surface,InternalAtomPosition,inhead,inlscl).size()) return 1;
	else {
		cout<<"Same amount: "<<find_atoms_around_v2(follow_normal,rc_surface,InternalAtomPosition,inhead,inlscl).size()<<'\n';
	}return(consistency_of_normal(cur_pos,normal_direction,2*movement));
	
}