# CSCI-653 Project
## 3D Atomic Structure Projection, MD Simulation and Analysis
- This project aims to measure the geometric characteristic value, curvatures, normals and genus of atomic structure after MD simulation. 

## **Project Member** ##
Chang Liu (liu396@usc.edu)
Suyue Yuan (suyueyua@usc.edu)

### Problem Description
- Project any 3D structure obtained from experiments, in the form of voxel, VTK file or others, into structure at atomic scale with desired crystal lattice. 
- Use MD simulation to relax and perform mechanical tests on it.
- Analyze the morphology including curvatures, normals and genus, change after the simulation. 
   
### Methods and Algorithm
- Atomic Structure is seperated into surface atoms and internal atoms.
- Normal direction can be found summing all the vectors pointing from internal atoms around certain surface atoms within a cutoff.  
- Normal direction can be used to evaluate how many atoms needed to measure curvature.  
- Curvatures are found by fitting a bivariate [] or a quadratic surface []. 
- Genus is calculated from faces and edges of surface mesh. 
- Paralled computing can achieved by spatial decomposition.

### Expected Results
1. Atomic structures from experiments.

![Image of a porous rock](https://github.com/liu396/CS653/blob/master/zhazha.jpg | width="400" height="790")
![Image of the projected atomic structure](https://github.com/liu396/CS653/blob/master/zhazha.png){:height="50%" width="50%"}

2. Normal direction distribution of structure. 

![Image of Normal](https://github.com/liu396/CS653/blob/master/xy_frame0.png)

3. Principle curvature distribution of structure.
![Image of Curvature](https://github.com/liu396/CS653/blob/master/45RD_matrix.png)
4. Genus 

### Reference
[1]: Lee, Ta-Chih, Rangasami L. Kashyap, and Chong-Nam Chu. "Building skeleton models via 3-D medial surface axis thinning algorithms." CVGIP: Graphical Models and Image Processing 56.6 (1994): 462-478.

[2]: Stuckner, Joshua, et al. "AQUAMI: An open source Python package and GUI for the automatic quantitative analysis of morphologically complex multiphase materials." Computational Materials Science 139 (2017): 320-329.

[3]: Soyarslan, Celal, et al. "3D stochastic bicontinuous microstructures: Generation, topology and elasticity." Acta materialia 149 (2018): 326-340.

[4]: Liu, Chang, and Paulo S. Branicio. "Efficient generation of non-cubic stochastic periodic bicontinuous nanoporous structures." Computational Materials Science 169 (2019): 109101.
