# spinchem
My cool new project!

Thanks to @damyanfr for providing the base for this project.

# Build

## Install [fpm](https://github.com/jordansissel/fpm)

The easiest way is by using [uv](https://docs.astral.sh/uv/).

```bash
# check uv is installed
$ uv --version
# install uv (if not already installed)
$ curl -fsSL https://astral.sh/uv/install.sh | sh
```

```bash
$ uv tool install fpm
$ fpm --version
Version:     0.10.1, alpha
Program:     fpm(1)
Description: A Fortran package manager and build system
Home Page:   https://github.com/fortran-lang/fpm
License:     MIT
OS Type:     Linux
```

## Install libomp-dev (for OpenMP support)

```bash
$ sudo apt-get install libomp-dev
```

## Clone the repository 
```bash
$ git clone https://github.com/KenHino/Spin_dynamics.git
$ fpm build
```

## Run
```bash
$ mkdir out
$ build/gfortran_BFCF334AC0CBF7DD/app/spinchem input2.ini
```

## References

- [Trace sampling](https://pubs.aip.org/aip/jcp/article/154/8/084121/76322)
- [Symmetry dynamics](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.220604)
- [Semi-classical](https://pubs.aip.org/aip/jcp/article/139/12/124106/74601)