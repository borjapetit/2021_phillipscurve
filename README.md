
# **Costain, Nakov and Petit (2021) - Replication codes**

This repostory includes all the codes required to replicate the results in the paper **"Flattening of the Phillips curve with state-dependent prices and wages"**, by with [James Costain](https://sites.google.com/site/jimcostain/) (Banco de España), [Anton Nakov](https://sites.google.com/site/antonnakov/) (ECB and CEPR) and [Borja Petit](https://borjapetit.github.io) (CUNEF).

## Struture of the codes

The model is solved it two steps. First, we use Fortran to solve for the steady state and to compute the Jacobian of the dynamic system. Then, we upload the Jacobian matrix into Matlab to make the linearization using the Klein's methodology, and to compute the Impulse Response Functions.

### Fortran codes

### Matlab codes

### Others


## Compiling and running the codes


To compìle the Fortran codes the command is:

```
cd  [set your own working directory]

gfortran [optional: compilation flags] -J $(pwd)/compiledfiles toolkit.f90 parameters.f90 solution.f90 dynamics.f90 main.f90 -o lpw
```

You tend o use the following compilation flags:

```
-fopenmp
-O3
-ffixed-line-length-150
-fmax-stack-var-size=1000000
```

