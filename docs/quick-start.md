## Getting started

Familiar with the workflow and ready to start a new STILT project?

1. Ensure the UATAQ R package is up to date

```bash
Rscript -e "install.packages('devtools'); devtools::install_github('benfasoli/uataq')"
```

2. Initialize STILT project using the `uataq::stilt_init()` R function

```bash
Rscript -e "uataq::stilt_init('myproject')"
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
- [Tutorial: Stationary simulations](https://github.com/uataq/stilt-tutorials/tree/main/01-wbb)
