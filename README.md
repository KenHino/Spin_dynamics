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
$ cd Spin_dynamics
$ git submodule update --init
$ fpm install --profile release --prefix .
```

## Run
```bash
$ mkdir out
$ export PATH=$PWD/bin:$PATH
$ spinchem input2.ini
```

## Build documents
```bash
$ cd Spin_chemistry
$ uv sync # If not already installed
$ uv run mkdocs build
# or
$ uv run mkdocs serve
```

## For Apple Silicon
In Apple Silicon device, gcc is aliased to clang, which does not support openmp.
Thus, one needs to install gcc and set manually.
```bash
brew tap fortran-lang/fortran # will give you v0.12.0 on Apple Silicon
brew install libomp gcc fpm
export PATH="/opt/homebrew/opt/gcc/bin:$PATH"   # let plain “gcc” resolve to gcc-15

# pick compilers
export FPM_FC=gfortran-15
export FPM_CC=gcc-15
export FPM_CXX=g++-15

# add OpenMP once, so it’s on every Fortran *and* C/C++ compile
export FPM_FFLAGS="-fopenmp"
export FPM_CFLAGS="-fopenmp"
export FPM_CXXFLAGS="-fopenmp"

/opt/homebrew/bin/fpm install --profile release --prefix .
```

## References

- [Trace sampling](https://pubs.aip.org/aip/jcp/article/154/8/084121/76322)
- [Symmetry dynamics](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.220604)
- [Semi-classical](https://pubs.aip.org/aip/jcp/article/139/12/124106/74601)
