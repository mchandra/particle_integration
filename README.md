# particle_integration

1. Edit Parameters required in params.py 

2. Run diagnostics.py to test all the modules. Uncomment the file writing code found near the end of the file to write the data required for particles in your working directory

3. currentdepositer.py contains the functions for charge and current deposition using a b0 spline.

4. fields.py has been modified to include currents due to particles.

5. While running the code make sure all the numerical stability criterion are satisfied(described in the document). Keep charges low and velocites low in general.

6. currentdepositer.py contains the commented scripts that show the inner workings of the current depositor function.
