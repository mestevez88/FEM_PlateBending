**FEM_PlateBending** is a Python implementation of the Finite Element Method applied to a Kirchoff plate subjected to forces acting normal to its mid-surface plane. The theoretical background is based on the book "The Implementation of the Finite Element Method" by Viktor Hristovski. **FEM_PlateBending** consists on the following:

The Jupyter Notebook "FEM_PlateBending.ipynb" is a step-by-step guide explaining the theoretical background and its Python implementation. Each of the steps are followed by a code section in which the reader can directly relate the theoretical background with the implementation. The reader can try different examples by modifying the desired geometrical and mechanical properties as well as the structural problem's boundary conditions.

The "main.py" script allows a more flexible use of the implementation. Here, the user can controll the following properties:
1) Plate Geometric Properties
2) Rectangular Mesh
3) Material Properties
4) Load Vector
5) Boundary Conditions

The "LocalStiffness.py" script contains a class called FEM_PlateBending. The class is initialized by specifying the following parameters:
LX: Length of the plate in the X-direction.
LY: Length of the plate in the Y-direction.
EX: Number of rectangular elements in the X-direction.
EY: Number of rectangular elements in the Y-direction.
Gpx: Number of Gauss Points for numerical integration in the X-direction.
Gpy: Number of Gauss Points for numerical integration in the Y-direction.
t: Plate thickness.
materials: a 2 by 2 material matrix specifying the Elasticity Modulus (E) and Poisson's ratio (\nu) of the Solid and Voided materials.

**Installation** 
1) Download this repository to your local drive.
2) Create a virtual environment.
3) Activate your virtual environment.
4) Install rquirements.txt file typing "pip install -r requirements.txt".

**Usage**
1) 
