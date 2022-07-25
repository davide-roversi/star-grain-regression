# star-grain-regression
A MATLAB 2D simulation of the regression in a star grain solid rocket motor

--------------------------------------------------------------------------------------------------------------------------------
This code was written for a Bachelor's final project in Aerospace Engineering at Politecnico di Torino.
The code simulates the regression of a star shaped grain in a solid rocket motor.

Several assumptions are made:
1. the grain is homogeneous and has constant density;
2. A 0D combustion chamber model is used;
3. Flow is assumed isentropic and chamber pressure and temperature values are considered as total values.

The profile of the star is discretized through nodes that are moved backwards at each iteration.
At each iteration the code calculates the burning surface, the chamber pressure and the regression ratio for the next iteration.

For more information and explanations about what the code does and how to use it, see the .pptx and .pdf files in the repo.
