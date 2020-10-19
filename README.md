# ising_hmc
C-Program simulating the Ising model with the HMC

To run the code, download it, e.g. with

```
git clone git@github.com:HISKP-LQCD/ising_hmc.git
```

go into the newly created directory `ising_hmc` and execute:

```
make
./ising_hmc example_input_2d.dat example_input_alltoall.dat
```

You are encouraged to replace the example-files with alternative ones including
different parameters.

Every run generates an output file containing five columns, namely the 
magnetisation, its absolute value, its square, the energy per site and an 
indicator if this trajectory has been accepted (1 means accepted, 0 rejected). 
Every line corresponds to a different measured trajectory.
