**FEM_PlateBending** is a Python implementation of the Finite Element Method applied to a Kirchoff plate subjected to forces acting normal to its mid-surface plane. The theoretical background is based on the book "The Implementation of the Finite Element Method" by Viktor Hristovski. **FEM_PlateBending** consists on the following:

The Jupyter Notebook "FEM_PlateBending.ipynb" is a step-by-step guide explaining the theoretical background and its Python implementation. Each of the steps are followed by a code section in which the reader can directly relate the theoretical background with the implementation. The reader can try different examples by modifying the desired geometrical and mechanical properties as well as the structural problem's boundary conditions.

The "main.py" script allows a more flexible use of the implementation. Here, the user can controll the following properties:
1) Plate Geometric Properties
2) Rectangular Mesh
3) Material Properties
4) Load Vector
5) Boundary Conditions

The **LocalStiffness.py** script calculates the local stifness matrices for each rectangular mesh element considering either a solid and a voided material.

The **GlobalStiffness.py** scripts constructs the global stiffness matrix $K$ by allocating the nodal stiffness provided by each adjacent plate element having either solid or voided material properties.

**Installation** 
1) Download this repository to your local drive.
2) Create a new virtual environment and activate it.
3) Install the requirements specified in rquirements.txt by typing "pip install -r requirements.txt".
