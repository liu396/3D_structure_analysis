# CSCI-653 Project
## Bicontinuous Structure Based on Leveled Field Method and Geometric Measurement
- This project aims to generate periodic or non periodic bicontinuous atomistic sample in arbitary parallelpiped, and to measure the geometric characteristic value, especially the ligament size. 

## **Project Member** ##
Chang Liu (liu396@usc.edu)
Suyue Yuan (suyueyua@usc.edu)

### Problem Description
There exists several methods to model of bicontinuous structure, such as Phase Field and Kinetic Monte Carlos methods. These simulation based methods usually need considerable computational resourses and time. This project aims to implement a direct method based on Cahn's original idea and the following work [3-4]. The atomistic structure generated with this user friendly code can reside in cubic or other arbitary parallelepipeds with or without periodic boudary conditions.   

The next challenge is to measure the ligament size, which is one of the most important characteristics of these structures. The prevalent method to measure such quantity is to take slices from the 3D structure and use tools for 2D pictures, such as Aquami [2], measuring the ligament size in 2D. Several slices along different directions would be take and the average ligament distribution would indicate the overal ligament size. This project aims to 3D skeletonization [1] and distance map on the target structure. The combination of this two transforms will dictate accurately the ligament size distribution.  

### Methods and Algorithm
- The leveled field method is efficiently realized by using OpenMPI. 
- To determine the ligament size, the multiplication of distance map transformation and skeletonization using 3D digital thinning technique is the key.

### Expected Results
1. Bicontinuous structrue in arbitary parallelepiped with controllable topological and geometrical featrues.
2. Skeleton of 3D bicontinuous structure without change of topological information.
3. Distance Map of 3D bicontinuos structure.
4. Ligament size distribution. 

### Reference
[1]: Lee, Ta-Chih, Rangasami L. Kashyap, and Chong-Nam Chu. "Building skeleton models via 3-D medial surface axis thinning algorithms." CVGIP: Graphical Models and Image Processing 56.6 (1994): 462-478.

[2]: Stuckner, Joshua, et al. "AQUAMI: An open source Python package and GUI for the automatic quantitative analysis of morphologically complex multiphase materials." Computational Materials Science 139 (2017): 320-329.

[3]: Soyarslan, Celal, et al. "3D stochastic bicontinuous microstructures: Generation, topology and elasticity." Acta materialia 149 (2018): 326-340.

[4]: Liu, Chang, and Paulo S. Branicio. "Efficient generation of non-cubic stochastic periodic bicontinuous nanoporous structures." Computational Materials Science 169 (2019): 109101.
