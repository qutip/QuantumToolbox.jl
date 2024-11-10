# How to build documentation and start Vitepress site locally ?

## Working Directory
All the commands should be run under the root folder of the package: `/path/to/QuantumToolbox.jl/`

The document pages will be generated in the directory: `/path/to/QuantumToolbox.jl/docs/build/` (which is ignored by git).

## Method 1: Run with `make` command
Run the following command to instantiate and build the documentation:
```shell
make docs
```

Run the following command to start Vitepress site of documentation:
> [!NOTE]
> You need to install `Node.js` and `npm` first.
```shell
make vitepress
```

## Method 2: Run commands manually

### Build Pkg
Run the following command:
```shell
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```
> **_NOTE:_** `Pkg.develop(PackageSpec(path=pwd()))` adds the local version of `QuantumToolbox` as dev-dependency instead of pulling from the registered version.

### Build Documentation
Run the following command:
```shell
julia --project=docs docs/make.jl
```

### Start a local Vitepress site
> [!NOTE]
> You need to install `Node.js` and `npm` first.

Install `npm` dependencies:
```shell
npm --prefix docs i
```

Run the following command:
```shell
npm --prefix docs run docs:dev
```