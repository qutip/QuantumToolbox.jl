# How to build documentation locally ?

## Working Directory
All the commands should be run under the root folder of the package: `/path/to/QuantumToolbox.jl/`

The document pages will be generated in the directory: `/path/to/QuantumToolbox.jl/docs/build/` (which is ignored by git).

## Method 1: Run with `make` command
Run the following command:
```shell
make docs
```

## Method 2: Run `julia` command manually

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

# How to start a local Vitepress site ?

> [!NOTE]
> You need to install `Node.js` and `npm` first. 

Enter `docs` directory first:
```shell
cd /path/to/QuantumToolbox.jl/docs
```

Install `npm` dependencies:
```shell
npm i
```

Run the following command:
```shell
npm run docs:dev
```