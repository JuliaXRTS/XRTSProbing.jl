# XRTSProbing

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaXRTS.github.io/XRTSProbing.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaXRTS.github.io/XRTSProbing.jl/dev)
[![Build Status](https://github.com/JuliaXRTS/XRTSProbing.jl/workflows/Test/badge.svg)](https://github.com/JuliaXRTS/XRTSProbing.jl/actions)

XRTSProbing.jl is an experimental Julia package to simulate and analyze the interaction of
photons (or laser fields) with a cloud of distributed electrons, with an eye toward X-ray
Thomson Scattering (XRTS) applications.

## Installation

### Building from source

Since `XRTSProbing.jl` is not registered yet, one needs to build is from source. First one
needs to clone the repository and go to the source folder.

Then, the package can be build locally by using the `Pkg`:

```bash
julia --project -e "import Pkg; Pkg.build()"
```

## Building the docs locally

Building the docs locally involves the following steps:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs --color=yes docs/make.jl
```

The website with the documentation can then be accessed using the browser of your choice

```bash
open docs/build/index.html
```

Here `open` stands for your browser of the open command on macOS.

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://JuliaXRTS.github.io/XRTSProbing.jl/dev/90-contributing/)

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
