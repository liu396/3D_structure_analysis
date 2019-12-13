#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAXSTRING 1024
#define NUMSTRING 11
#define NUMSKIP 15
#define NPOINT 233371
#define NTRI 488442
#define NATOM 2145536

int main()
{
	double pos_x[NPOINT], pos_y[NPOINT], pos_z[NPOINT];
	double atom_x[NATOM], atom_y[NATOM], atom_z[NATOM];
	double xb = 740.0, xdelta = 169.0, yb = 909.0, ydelta = 169.0, zb = 1063.0, zdelta = 287.0;
	double ab_x = 2.4745301520011589e+02, ab_y = 2.4745301520011589e+02, ab_z = 4.2212573181196234e+02;
	double s, t, area[NTRI], s1[NTRI], t1[NTRI], resX, resY, resZ, X, Y, Z;
	double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, n_x, n_y, n_z, cal, pz, cputime;
	char info[NUMSTRING][MAXSTRING], info2[NUMSKIP][MAXSTRING];
	int i, j, geo, tri[NTRI][3], svatom[NATOM] = {0}, arsign[NTRI] = {0};
	int aid[NATOM], atype[NATOM], aft0[NATOM], aft1[NATOM], aft2[NATOM];
	int tri_l, tri_r, acount = 0;
    clock_t time1,time2;
	time1 = clock();
	FILE *fp;
	fp = fopen("rock.ply", "r");
	if (fp != NULL)
	{
		for (i = 0; i < NUMSTRING; i++)
		{
			fgets(info[i], MAXSTRING, fp);
		}
		for (i = 0; i < NPOINT; i++)
		{
			fscanf(fp, "%lf %lf %lf\n", &pos_x[i], &pos_y[i], &pos_z[i]);
		}
		for (i = 0; i < NTRI; i++)
		{
			fscanf(fp, "%d %d %d %d\n", &geo, &tri[i][0], &tri[i][1], &tri[i][2]);
		}
	}
	fclose(fp);
	
	fp = fopen("crystal.data", "r");
	if (fp != NULL)
	{
		for (i = 0; i < NUMSKIP; i++)
		{
			fgets(info2[i], MAXSTRING, fp);
		}
		for (i = 0; i < NATOM; i++)
		{
			fscanf(fp, "%d %d %lf %lf %lf %d %d %d\n", &aid[i], &atype[i], &atom_x[i], &atom_y[i], &atom_z[i], &aft0[i], &aft1[i], &aft2[i]);
		}
	}
	fclose(fp);
	
	// Barycentric coordinate preparations
	// Details see https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
	for (i = 0; i < NTRI; i++)
	{
		// following = 2 * singed area
		area[i] = -(pos_y[tri[i][1]] * pos_x[tri[i][2]]) \
		 	      + pos_y[tri[i][0]] * (pos_x[tri[i][2]] - pos_x[tri[i][1]]) \
			      + pos_x[tri[i][0]] * (pos_y[tri[i][1]] - pos_y[tri[i][2]]) \
			      + pos_x[tri[i][1]] * pos_y[tri[i][2]];	
		if (area[i] > 0) arsign[i] = 1;
		s1[i] = pos_y[tri[i][0]] * pos_x[tri[i][2]] - pos_x[tri[i][0]] * pos_y[tri[i][2]];
		t1[i] = pos_x[tri[i][0]] * pos_y[tri[i][1]] - pos_y[tri[i][0]] * pos_x[tri[i][1]];
	}
	
	// Check with a point 
	resX = xdelta / ab_x;
	resY = ydelta / ab_y;
	resZ = zdelta / ab_z;
	for (j = 0; j < NATOM; j++)
	{
		X = atom_x[j] * resX + xb;
		Y = atom_y[j] * resY + yb;
		Z = atom_z[j] * resZ + zb;
		tri_l = 0;
		tri_r = 0;
		for (i = 0; i < NTRI; i++)
		{
			s = s1[i] + (pos_y[tri[i][2]] - pos_y[tri[i][0]]) * X + (pos_x[tri[i][0]] - pos_x[tri[i][2]]) * Y;
			t = t1[i] + (pos_y[tri[i][0]] - pos_y[tri[i][1]]) * X + (pos_x[tri[i][1]] - pos_x[tri[i][0]]) * Y;
			if (arsign[i] == 1)
			{
				if (s > 0 && t > 0 && s+t < area[i])
				{
					if (pos_z[tri[i][0]] <= Z && pos_z[tri[i][1]] <= Z && pos_z[tri[i][2]] <= Z)
						tri_l ++;
					else if (pos_z[tri[i][0]] > Z && pos_z[tri[i][1]] > Z && pos_z[tri[i][2]] > Z)
						tri_r ++;
					else
					{
						v1_x = pos_x[tri[i][1]] - pos_x[tri[i][0]];
						v1_y = pos_y[tri[i][1]] - pos_y[tri[i][0]];
						v1_z = pos_z[tri[i][1]] - pos_z[tri[i][0]];
						v2_x = pos_x[tri[i][2]] - pos_x[tri[i][0]];
						v2_y = pos_y[tri[i][2]] - pos_y[tri[i][0]];
						v2_z = pos_z[tri[i][2]] - pos_z[tri[i][0]];
						n_x = v1_y * v2_z - v1_z * v2_y;
						n_y = v1_x * v2_z - v1_z * v2_x;
						n_z = v1_x * v2_y - v1_y * v2_x;
						cal = n_x * (X - pos_x[tri[i][0]]) + n_y * (Y - pos_y[tri[i][0]]);
						pz = pos_z[tri[i][0]] - cal / n_z;
						if (pz <= Z)
							tri_l ++;
						else
							tri_r ++;
					}
				}
			}
			else
			{
				if (s < 0 && t < 0 && s+t > area[i])
				{
					if (pos_z[tri[i][0]] <= Z && pos_z[tri[i][1]] <= Z && pos_z[tri[i][2]] <= Z)
						tri_l ++;
					else if (pos_z[tri[i][0]] > Z && pos_z[tri[i][1]] > Z && pos_z[tri[i][2]] > Z)
						tri_r ++;
					else
					{
						v1_x = pos_x[tri[i][1]] - pos_x[tri[i][0]];
						v1_y = pos_y[tri[i][1]] - pos_y[tri[i][0]];
						v1_z = pos_z[tri[i][1]] - pos_z[tri[i][0]];
						v2_x = pos_x[tri[i][2]] - pos_x[tri[i][0]];
						v2_y = pos_y[tri[i][2]] - pos_y[tri[i][0]];
						v2_z = pos_z[tri[i][2]] - pos_z[tri[i][0]];
						n_x = v1_y * v2_z - v1_z * v2_y;
						n_y = v1_x * v2_z - v1_z * v2_x;
						n_z = v1_x * v2_y - v1_y * v2_x;
						cal = n_x * (X - pos_x[tri[i][0]]) + n_y * (Y - pos_y[tri[i][0]]);
						pz = pos_z[tri[i][0]] - cal / n_z;
						if (pz <= Z)
							tri_l ++;
						else
							tri_r ++;
					}
				}
			}						
		}
		if (tri_l % 2)
		{
			acount ++;
			svatom[j] = 1;		
		}
	}
	
	fp = fopen("newcrystal.data", "w");
	if (fp != NULL)
	{
		for (i = 0; i < NUMSKIP; i++)
		{
			fputs(info2[i], fp);
		}
		for (i = 0; i < NATOM; i++)
		{
			if (svatom[i] == 1)
			{
				fprintf(fp, "%d %d %lf %lf %lf %d %d %d\n", aid[i], atype[i], atom_x[i], atom_y[i], atom_z[i], aft0[i], aft1[i], aft2[i]);
			}
		}
		fprintf(fp, "Please change the number of atoms to: %d. \n", acount);
	}
	fclose(fp);
    time2 = clock();
    cputime = (time2 - time1) / (double)CLOCKS_PER_SEC;
	printf("===== Total cpu time = %lf =====\n", cputime);
	
	return 0;
}