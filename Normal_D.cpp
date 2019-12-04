#include "Normal_D.h"
using namespace std;
int main(int argc, const char *argv[]){
	int i,j,a,m,n;
	double mm,nn;
	int mc[3];
	int c;
	ifstream cfgfile,surfacecfgfile,internalcfgfile;
    vector<vector<double> > point_list;
    vector<int> neighbor_list;
    double x,y,z;
    vector<double> normal(3,0);
	vector<double> cuml_normal(3.0);
    vector<vector<double> > new_neighbors;
	double distance;
	ofstream xyfile,xzfile,xplanefile;
	int bin=22,bin2=11;
	int xycount=0,xzcount=0;
	double interval;
	int *xydata,*xzdata,*xplanedata;
	
	
	cfgfile.open("Frame0.cfg");
	// surfacecfgfile.open("diameter_100_surface.cfg");
	internalcfgfile.open("Frame0_internal.cfg");
	//cfgfile.open("hemi_cylinder.cfg");
	surfacecfgfile.open("Frame0_surface.cfg");
	//internalcfgfile.open("hemi_internal.cfg");
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
	xyfile.open("x.vs.y.txt");
	xzfile.open("x.vs.z.txt");
	xplanefile.open("x.vs.plane.txt");
	xyfile<<"X_interval    Y_interval     Count\n";
	xzfile<<"X_interval    Z_interval     Count\n";
	interval=2.2/double(bin);
	xydata=new int[bin*bin];
	xzdata=new int[bin*bin];
	xplanedata=new int[bin2*bin2];
	
	for(i=0;i<bin*bin;i++){
		xydata[i]=0;
		xzdata[i]=0;
	}
	for(i=0;i<bin2*bin2;i++){
		xplanedata[i]=0;
	}
	
	
	for(i=0;i<SurfaceAtomPosition.size();i++){
		if (i%500==0) cout<<"Processing Atom "<<SurfaceAtomInfo[i][0]<<'\n';
		normal=find_normal(SurfaceAtomPosition[i],10);
		cuml_normal[0]+=abs(normal[0]);
		cuml_normal[1]+=abs(normal[1]);
		cuml_normal[2]+=abs(normal[2]);
		m=int((normal[0]+1.1)/interval);
		//if(m==bin) m--;
		n=int((normal[1]+1.1)/interval);
		//if(n==bin) n--;
		xydata[m*bin+n]++;
		n=int((normal[2]+1.1)/interval);
		//if(n==bin) n--;
		xzdata[m*bin+n]++;
		mm=abs(normal[0]);
		nn=sqrt(normal[1]*normal[1]+normal[2]*normal[2]);
		//cout<<mm*mm+nn*nn<<'\n';
		m=int(mm/interval);
		n=int(nn/interval);
		xplanedata[m*bin2+n]++;
	}
	
	for(i=0;i<bin;i++){
		for(j=0;j<bin;j++){
			xyfile<<i<<"   "<<j<<"   "<<double(xydata[i*bin+j])/double(SurfaceAtomPosition.size())/0.1/0.1<<'\n';
			xzfile<<i<<"   "<<j<<"   "<<double(xzdata[i*bin+j])/double(SurfaceAtomPosition.size())/0.1/0.1<<'\n';
		}
	}
	
	for(i=0;i<bin2;i++){
		for(j=0;j<bin2;j++){
			xplanefile<<i<<"   "<<j<<"   "<<double(xplanedata[i*bin2+j])/double(SurfaceAtomPosition.size())/0.1/0.1<<'\n';
		}
	}
	
	xplanefile.close();
	xyfile.close();
	xzfile.close();

	cout<<"Normal is "<<cuml_normal[0]<<' '<<cuml_normal[1]<<' '<<cuml_normal[2]<<'\n';
	
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
vector<double> find_normal(vector<double> &position,double rc){
	vector<vector<double>> neighbors;
	vector<double> normal_vector;
	double norm;
	normal_vector.push_back(0.0);
	normal_vector.push_back(0.0);
	normal_vector.push_back(0.0);
	//cout<<"Finding neighbors...\n";
	neighbors=find_atoms_around(position,rc,InternalAtomPosition,inhead,inlscl);
	//cout<<"Neighbors found\n";
	//cout<<"Number of Neighbors: "<<neighbors.size()<<'\n';
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
double diameter(vector<double> &position,double movement,double rc){
	vector<double> cur_pos;
	double total_movement=0.0;
	vector<double> normal_vector;
	int num_old;
	int num_new=0;
	normal_vector=find_normal(position,rc);
	cur_pos.push_back(position[0]);
	cur_pos.push_back(position[1]);
	cur_pos.push_back(position[2]);
	num_old=find_atoms_around(position,rc_surface,SurfaceAtomPosition,shead,slscl).size();
	num_new=num_old;
	cout<<"initial number of neighbors: "<<num_old<<'\n';
	while(num_old>=num_new){
		cur_pos[0]+=movement*normal_vector[0];
		cur_pos[1]+=movement*normal_vector[1];
		cur_pos[2]+=movement*normal_vector[2];	
		total_movement+=movement;
		if(cur_pos[0]>h[0][0]) cur_pos[0]=cur_pos[0]-h[0][0];
		if(cur_pos[0]<0) cur_pos[0]=h[0][0]+cur_pos[0];
		if(cur_pos[1]>h[1][1]) cur_pos[1]=cur_pos[1]-h[1][1];
		if(cur_pos[1]<0) cur_pos[1]=h[1][1]+cur_pos[1];
		if(cur_pos[2]>h[2][2]) cur_pos[2]=cur_pos[2]-h[2][2];
		if(cur_pos[2]<0) cur_pos[2]=h[2][2]+cur_pos[2];
		num_old=num_new;
		num_new=find_atoms_around(cur_pos,rc_surface,SurfaceAtomPosition,shead,slscl).size();
	}
	return total_movement;	
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