# QEDprobing.jl

[![Dev](https://img.shields.io/badge/docs-main-blue.svg)](https://qedjl-applications.pages.hzdr.de/QEDprobing.jl)
[![Build Status](https://github.com/qedapplications/QEDprobing.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/qedapplications/QEDprobing.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://codebase.helmholtz.cloud/qedapplications/QEDprobing.jl/badges/main/pipeline.svg)](https://codebase.helmholtz.cloud/qedapplications/QEDprobing.jl/pipelines)
[![Coverage](https://codebase.helmholtz.cloud/qedapplications/QEDprobing.jl/badges/main/coverage.svg)](https://codebase.helmholtz.cloud/qedapplications/QEDprobing.jl/commits/main)

`QEDprobing.jl` is an experimental implementation of the interaction of photons or laser
fields with a cloud of distributed electrons. The code will change rapidly without proper
announcement. For more curated version of parts of the code, see [`QED.jl`](https://github.com/QEDjl-project).

For more background information on concepts and algorithms used in the development, see
the [wiki](https://codebase.helmholtz.cloud/qedjl-applications/QEDprobing.jl/-/wikis/Overview) (only for members).

# Installation

## Building from source

Since `QEDprobing.jl` is not registered yet, one needs to build is from source. First one
needs to clone the repository and go to the source folder.

Then, the package can be build locally by using the `Pkg`:

```bash
julia --project -e "import Pkg; Pkg.build()"
```

Since `QEDprobing.jl` relies on the latest development of `QED.jl`, one needs to add the
respective feature branches by hand. For that, you can use the script `add_dev_packages.jl`:

```bash
julia --project add_dev_packages.jl
```

## Building the docs locally

Building the docs locally involves the following steps:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs add_dev_packages.jl
julia --project=docs --color=yes docs/make.jl
```

The website with the documentation can then be accessed using the browser of your choice

```bash
open docs/build/index.html
```

Here `open` stands for your browser of the open command on macOS.
