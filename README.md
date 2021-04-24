
# **Costain, Nakov and Petit (2021) - Replication codes**

This repostory includes all the codes required to replicate the results in the paper **"Flattening of the Phillips curve with state-dependent prices and wages"**, by with [James Costain](https://sites.google.com/site/jimcostain/) (Banco de España), [Anton Nakov](https://sites.google.com/site/antonnakov/) (ECB and CEPR) and [Borja Petit](https://borjapetit.github.io) (CUNEF).

## Structure of the codes

The model is solved it two steps. First, we use Fortran to solve for the steady state and to linearize the dynamic system. Then, we use the Kelin's complex QZ decomposition solution method to solve the model and compute the Impulse Response FUnction. 

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

The code offers 7 possibilities:
1. Solve only the steady state
2. Solve the steady state and the dynamics
3. Calibrate the parameters of the model
4. Solve the dynamics for the 6 possible combinations of decision-making costs (with baseline inflation rate)
5. Solve the dynamics for an inflation rate of -2%, -1%, 0%, +4% and +8%. The rate of 2% is the baseline inflation rate) and for the sticky version, the flexible -prices case, the flexible-wages case and the flexible prices and wages case.
6.Solve the steady-state for all possible cases (different combinations of inflation rate and adjustment cost parameters)

| Adjustment cost / Inflation rate   | 2%   | -1%  |   0%  |   4%   |  8%  |  -2%   |  1% |
|-----|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
| Baseline                     |   V10  |  V11  |  V12  |  V13  |  V14  |  V15  |  V16 |
| Semi-flex P and sticky W     |   V20  |  V21  |  V22  |  V23  |  V24  |  V25  |  V26 |
| Flexible P and sticky W      |   V30  |  V31  |  V32  |  V33  |  V34  |  V35  |  V36 |
| Sticky P and semi-flex W     |   V40  |  V41  |  V42  |  V43  |  V44  |  V45  |  V46 |
| Sticky P and flexible W      |   V50  |  V51  |  V52  |  V53  |  V54  |  V55  |  V56 |
| Flexible P and flexible W    |   V60  |  V61  |  V62  |  V63  |  V64  |  V65  |  V66 |

To manage the different versions of the model, we label each of them as Vab, where "a" and "b" depend on the different combinations of adjustment cost parameters and steady-state inflation rates. In particular, the adjustment cost parameters can take on 6 different combinations:

| Adjustment cost   | a   |
|-----|:-----:|
| Baseline                              |   1  |
| Semi-flex prices and sticky wages     |   2  |
| Flexible prices and sticky wages      |   3  |
| Sticky prices and semi-flex wages     |   4  |
| Sticky prices and flexible wages      |   5  |
| Flexible prices and flexible wages    |   6  |

The inflation rate can take on 7 values:

| Inflation rate        | b      |
|-----                  |:------:|
| Inflation rate of 2%  |   0  |  
| Inflation rate of -1% |   1  |
| Inflation rate of 0%  |   2  |
| Inflation rate of 4%  |   3  |
| Inflation rate of 8%  |   4  |  
| Inflation rate of -2% |   5  |
| Inflation rate of 1%  |   6  |

Beyond this, we allow the inflation rate to take two different extra values to compute the Phillips curve slope presented in the paper. In particular:

| Inflation rate        | b      |
|-----                  |:------:|
| Inflation rate of 4.63% (US, 1980-2000)  |   7  |
| Inflation rate of 2.01% (US, 2000-2020)  |   8  |
