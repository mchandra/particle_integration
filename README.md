# particle_integration

##Points to remember before attempting to use the code.

* Edit Parameters required in params.py(Set no of particles, number of zones along x and y, Length of domain along x and y etc).

* b0 Spline shape factors for current and charges has been used in the code.

* For a 1 D simulation, set number of zones along required dimension to be 1.

* There is no Poisson solver, Hence the current code will give erroreneous results.

* Run diagnostics.py to test all the modules. Uncomment the file writing code found near the end of the file to write the data required for particles in your working directory

* currentdepositer.py contains the functions for charge and current deposition using a b0 spline.

* fields.py has been modified to include currents due to particles.

* While running the code make sure all the numerical stability criterion are satisfied(described in the document). Keep charges low and velocites low in general.

* currentdepositer.py contains the commented scripts that show the inner workings of the current depositor function.
