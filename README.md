# 3D Atomic Structure Projection, MD Simulation and Analysis
- This project aims to measure the geometric characteristic value, curvatures, normals and genus of atomic structure after MD simulation. 

## **Project Member** ##
Chang Liu (liu396@usc.edu)
Suyue Yuan (suyueyua@usc.edu)

### Problem Description
- Project any 3D structure obtained from experiments, in the form of voxel, or surface mesh file, into structure at atomic scale with desired crystal lattice. 
- Use MD simulation to relax and perform mechanical tests on it.
- Analyze the morphology including curvatures, normals and genus, change after the simulation. 
   
### Methods and Algorithm
- Atomic Structure is seperated into surface atoms and internal atoms.
- Normal direction can be found summing all the vectors pointing from internal atoms around certain surface atoms within a cutoff.  
- Normal direction can be used to evaluate how many atoms needed to measure curvature.  
- Curvatures are found by fitting a bivariate [1] or a quadratic surface [2]. 
- Genus is calculated from faces and edges of surface mesh. 
- Paralled computing can achieved by spatial decomposition.

### Expected Results
1. Atomic structures from experiments.

![Image of a porous rock](https://github.com/liu396/CS653/blob/master/zhaxiong.png)

2. Normal direction distribution of structure. 

![Image of Normal](https://github.com/liu396/CS653/blob/master/xy_frame0.png)

3. Principle curvature distribution of structure.
![Image of Curvature](https://github.com/liu396/CS653/blob/master/45RD_matrix.png)
4. Genus 

### Some results
1.The atomic model of obtained from experimental sample:
![Image of Curvature](https://github.com/liu396/CS653/blob/master/nanoporous_MG.png)

2.The curvature distribution of the above sample: 
![Image of Curvature](https://github.com/liu396/CS653/blob/master/new_md_curvature.png)


### Reference
[1]: Yokoya, N. and Levine, M.D., 1989. Range image segmentation based on differential geometry: A hybrid approach. IEEE Transactions on Pattern Analysis and Machine Intelligence, 11(6), pp.643-649.

[2]: Groshong, B., Bilbro, G. and Snyder, W., 1989. Fitting a Quadratic Surface to Three Dimensional Data.
