This code was written for a Bachelor's final project in Aerospace Engineering at Politecnico di Torino.
The code simulates the regression of a star grain section in a solid rocket motor.

Several assumptions are made:
1. the grain is homogeneous and has constant density
2. A 0D combustion chamber model is used 
3. Flow is assumed isentropic and chamber pressure values are considered trotal values (for p and T)

The profile of the star is discretized through nodes that are moved backwards at each iteration.
At each iteration the code calculates the burning surface, the chamber pressure and the regression ratio for the next iteration.

For visual information on the geometry inputs required by the code, see attached .pdf files.
