<div align="center">
  <img src="./docs/src/assets/logo.png" alt="QuantumToolbox.jl logo" width="150">
</div>
# QuantumToolbox.jl

[![Release](https://img.shields.io/github/release/albertomercurio/QuantumToolbox.jl.svg)](https://github.com/albertomercurio/QuantumToolbox.jl/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11096277.svg)](https://doi.org/10.5281/zenodo.11096277)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FQuantumToolbox&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/QuantumToolbox)

[![Build Status](https://github.com/albertomercurio/QuantumToolbox.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/albertomercurio/QuantumToolbox.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Coverage](https://codecov.io/gh/albertomercurio/QuantumToolbox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/albertomercurio/QuantumToolbox.jl)

[![Doc-Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://albertomercurio.github.io/QuantumToolbox.jl/stable)
[![Doc-Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://albertomercurio.github.io/QuantumToolbox.jl/dev)

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
```
