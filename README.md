# MolecularDynamics
A Code to compare between serial and parallel version of Molecular Dynamics.

Running Serial Version (md.cpp) :-
```
g++ -o md -fopenmp md.cpp
./md nd np num_steps
```

Running Serial Version (md.cpp) :-
```
export OMP_NUM_THREADS=NUMBER_OF_THREADS
g++ -o mdParallel -fopenmp mdParallel.cpp
./mdParallel nd np num_steps
```
Video Presentation
```
https://drive.google.com/file/d/1zOTlu3Gy7qIeaYmTMc7NwS_nT6K5blpb/view?usp=sharing
```
