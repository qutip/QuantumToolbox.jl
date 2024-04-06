# QuantumToolbox

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://albertomercurio.github.io/QuantumToolbox.jl/dev)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://albertomercurio.github.io/QuantumToolbox.jl/stable)
[![Build Status](https://github.com/albertomercurio/QuantumToolbox.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/albertomercurio/QuantumToolbox.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Coverage](https://codecov.io/gh/albertomercurio/QuantumToolbox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/albertomercurio/QuantumToolbox.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10822817.svg)](https://doi.org/10.5281/zenodo.10822817)

## Introduction
[QuantumToolbox.jl](https://github.com/albertomercurio/QuantumToolbox.jl) is a cutting-edge Julia package designed for quantum physics simulations, closely emulating the popular Python [QuTiP](https://github.com/qutip/qutip) package. It uniquely combines the simplicity and power of Julia with advanced features like GPU acceleration and distributed computing, making simulation of quantum systems more accessible and efficient.

## Features
QuantumToolbox.jl is equipped with a robust set of features:

- **Quantum State and Operator Manipulation:** Easily handle quantum states and operators with a rich set of tools.
- **Dynamical Evolution:** Advanced solvers for time evolution of quantum systems.
- **Measurement and Statistics:** Comprehensive quantum measurement simulation and analysis.
- **GPU and Distributed Computing:** Leverage GPU and distributed resources for high-performance computing.

## Installation
```julia
using Pkg
Pkg.add("QuantumToolbox")
