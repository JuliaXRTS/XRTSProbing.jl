# XRTSProbing

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaXRTS.github.io/XRTSProbing.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaXRTS.github.io/XRTSProbing.jl/dev)
[![Build Status](https://github.com/JuliaXRTS/XRTSProbing.jl/workflows/Test/badge.svg)](https://github.com/JuliaXRTS/XRTSProbing.jl/actions)
[![Test workflow status](https://github.com/JuliaXRTS/XRTSProbing.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/JuliaXRTS/XRTSProbing.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/JuliaXRTS/XRTSProbing.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/JuliaXRTS/XRTSProbing.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/JuliaXRTS/XRTSProbing.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/JuliaXRTS/XRTSProbing.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaXRTS/XRTSProbing.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaXRTS/XRTSProbing.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/JuliaXRTS/XRTSProbing.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)

`XRTSProbing.jl` is an experimental implementation of the interaction of photons or laser
fields with a cloud of distributed electrons. The code will change rapidly without proper
announcement. For more curated version of parts of the code, see [`JuliaXRTS`](https://github.com/JuliaXRTS).

If you use XRTSProbing.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/JuliaXRTS/XRTSProbing.jl/blob/main/CITATION.cff).

## Installation

### Building from source

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

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

`XRTSProbing.jl` is an experimental implementation of the interaction of photons or laser
fields with a cloud of distributed electrons. The code will change rapidly without proper
announcement. For more curated version of parts of the code, see [`JuliaXRTS`](https://github.com/JuliaXRTS).

For more background information on concepts and algorithms used in the development, see
the [wiki](https://codebase.helmholtz.cloud/qedjl-applications/XRTSProbing.jl/-/wikis/Overview) (only for members).

## Installation

### Building from source

Since `XRTSProbing.jl` is not registered yet, one needs to build is from source. First one
needs to clone the repository and go to the source folder.

Then, the package can be build locally by using the `Pkg`:

```bash
julia --project -e "import Pkg; Pkg.build()"
```

Since `XRTSProbing.jl` relies on the latest development of `QuantumElectrodynamics.jl`, one needs to add the
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

## How to Cite

If you use XRTSProbing.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/JuliaXRTS/XRTSProbing.jl/blob/main/CITATION.cff).

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://JuliaXRTS.github.io/XRTSProbing.jl/dev/90-contributing/)

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
