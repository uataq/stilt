## Getting started 

Familiar with the workflow and ready to start a new STILT project?

1. Ensure the UATAQ R package is up to date

```bash
Rscript -e "install.packages('devtools'); devtools::install_github('benfasoli/uataq')"
```

2. Initialize STILT project using the `uataq::stilt_init()` R function

```bash
Rscript -e "uataq::stilt_init('myproject', branch='hysplit-merge')"
# Cloning into 'myproject'...
# remote: Enumerating objects: 60, done.
# remote: Counting objects: 100% (60/60), done.
# remote: Compressing objects: 100% (56/56), done.
# remote: Total 60 (delta 3), reused 25 (delta 2), pack-reused 0
# Unpacking objects: 100% (60/60), 2.26 MiB | 4.09 MiB/s, done.
# Compiling footprint kernel aggregation subroutine...
#
# STILT installation successful.
#
# Relevant manuscripts:
# 1. Fasoli, B., Lin, J. C., Bowling, D. R., Mitchell, L., and Mendoza, D.: 
#    Simulating atmospheric tracer concentrations for spatially distributed 
#    receptors: updates to the Stochastic Time-Inverted Lagrangian Transport 
#    model's R interface (STILT-R version 2), Geosci. Model Dev., 11, 2813-2824, 
#    [10.5194/gmd-11-2813-2018](https://doi.org/10.5194/gmd-11-2813-2018), 2018.
# 2. Lin, J. C., Gerbig, C., Wofsy, S. C., Andrews, A. E., Daube, B. C., Davis,
#    K. J. and Grainger, C. A.: A near-field tool for simulating the upstream 
#    influence of atmospheric observations: The Stochastic Time-Inverted Lagrangian
#    Transport (STILT) model, J. Geophys. Res., 108(D16), ACH 2-1-ACH 2-17, 
#    [10.1029/2002JD003161](https://doi.org/10.1029/2002JD003161), 2003.
#
# We strongly suggest you subscribe to the critical update notifications at
# https://uataq.github.io/stilt/
# to be notified if important STILT model updates updates.

cd myproject
ls
# Dockerfile  README.md  bin/  exe/  r/  setup  test/
```

3. Edit model configuration in `r/run_stilt.r`

4. Run model with `Rscript r/run_stilt.r`

```bash
Rscript r/run_stilt.r
# Initializing STILT
# Number of receptors: 1
# Number of parallel threads: 1
# Parallelization disabled. Executing simulations sequentially...
#
# Running simulation ID:   2015061822_-111.980323_40.782561_5
```

---

## Next steps

- [Install](install.md) for expanded documentation for installation with required dependencies
- [Tutorial: Stationary simulations](https://github.com/uataq/stilt-tutorials/tree/master/01-wbb)
