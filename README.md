# Particle-Based Solver

This package contains of Python scripts that allow for solving a system consisting of n-particles. Currently the framework can be used to study collisions(short-range interactions) between the particles, as well as model the thermodynamic behavior. The work is still in development, and shall be updated with additional modules such as long-range interactions.

## Starting Off

Most of the simulation parameters can be modified in the params.py file. Appropriate commenting has been added in the file to allow for easy understanding of the code. The variables have been named aptly to their served function. All simulation parameters can be changed in this file, except the time parameters for the simulation which must be changed from initialize.py. Once desired parameters are set, all that is left is the execution:

```
python initialize.py
python run.py
```

This will generate the output files in data_files folder. The simulation data may be retrieved from these files. The usual plots which are necessary can be easily generated from these files by making use of the post-processing scripts from the post_processing file.


### Prerequisites

The code is a purely Python based, and makes use of the numpy,scipy,h5py, and pylab modules. The current code also has arrayfire listed as a module. However, the packaged framework does not make use of this currently in this version, and is a feature in development. That being said, the stand-alone files are implemented using arrayfire.

## Authors

* **Shyam Sankaran** - [GitHub Profile](https://github.com/ShyamSS-95)
* **Tejas Mane** - [GitHub Profile](https://github.com/TejasMane)
* **Mani Chandra** - [GitHub Profile](https://github.com/mchandra)


## Current Features and Development

* In 2D, all current features are fully functional
* In 3D, only purely collisionless cases may be run. Collisions in 3D are still in development
* Speedup using ArrayFire are also yet to be implemented, although options are provided