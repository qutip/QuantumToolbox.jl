# [Intensive parallelization on a Cluster](@id doc:Intensive-parallelization-on-a-Cluster)

## Introduction

In this example, we will demonstrate how to seamlessly perform intensive parallelization on a cluster using the **QuantumToolbox.jl** package. Indeed, thanks to the [**Distributed.jl**](https://docs.julialang.org/en/v1/manual/distributed-computing/) and [**ClusterManagers.jl**](https://github.com/JuliaParallel/ClusterManagers.jl) packages, it is possible to parallelize on a cluster with minimal effort. The following examples are applied to a cluster with the [SLURM](https://slurm.schedmd.com/documentation.html) workload manager, but the same principles can be applied to other workload managers, as the [**ClusterManagers.jl**](https://github.com/JuliaParallel/ClusterManagers.jl) package is very versatile.

## SLURM batch script

To submit a batch script to [SLURM](https://slurm.schedmd.com/documentation.html), we start by creating a file named `run.batch` with the following content:

```bash
#!/bin/bash
#SBATCH --job-name=example
#SBATCH --output=output.out
#SBATCH --account=your_account
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=128GB
#SBATCH --time=0:10:00
#SBATCH --qos=parallel

# Set PATH to include the directory of your custom Julia installation
export PATH=/home/username/.juliaup/bin:$PATH

# Now run Julia
julia --project script.jl
```

where we have to replace `your_account` with the name of your account. This script will be used to submit the job to the cluster by using the following command in terminal:

```shell
sbatch run.batch
```

Here, we are requesting `10` nodes with `72` threads each (`720` parallel jobs). The `--time` flag specifies the maximum time that the job can run. To see all the available options, you can check the [SLURM documentation](https://slurm.schedmd.com/documentation.html). We also export the path to the custom Julia installation, which is necessary to run the script (replace `username` with your username). Finally, we run the script `script.jl` with the command `julia --project script.jl`.

In the following, we will consider two examples:

1. **Parallelization of a Monte Carlo quantum trajectories**
2. **Parallelization of a Master Equation by sweeping over parameters**

## Monte Carlo Quantum Trajectories

Let's consider a `2`-dimensional transverse field Ising model with `4x3` spins. The Hamiltonian is given by

```math
\hat{H} = \frac{J_z}{2} \sum_{\langle i,j \rangle} \hat{\sigma}_i^z \hat{\sigma}_j^z + h_x \sum_i \hat{\sigma}_i^x \, ,
```

where the sums are over nearest neighbors, and the collapse operators are given by 

```math
\hat{c}_i = \sqrt{\gamma} \hat{\sigma}_i^- \, .
```

In this case, the `script.jl` contains the following content:

```julia
using Distributed
using ClusterManagers

const SLURM_NUM_TASKS = parse(Int, ENV["SLURM_NTASKS"])
const SLURM_CPUS_PER_TASK = parse(Int, ENV["SLURM_CPUS_PER_TASK"])

exeflags = ["--project=.", "-t $SLURM_CPUS_PER_TASK"]
addprocs(SlurmManager(SLURM_NUM_TASKS); exeflags=exeflags, topology=:master_worker)


println("################")
println("Hello! You have $(nworkers()) workers with $(remotecall_fetch(Threads.nthreads, 2)) threads each.")

println("----------------")


println("################")

flush(stdout)

@everywhere begin
    using QuantumToolbox
    using OrdinaryDiffEq

    BLAS.set_num_threads(1)
end

# Define lattice

Nx = 4
Ny = 3
latt = Lattice(Nx = Nx, Ny = Ny)

# Define Hamiltonian and collapse operators
Jx = 0.0
Jy = 0.0
Jz = 1.0
hx = 0.2
hy = 0.0
hz = 0.0
γ = 1

Sx = mapreduce(i -> multisite_operator(latt, i=>sigmax()), +, 1:latt.N)
Sy = mapreduce(i -> multisite_operator(latt, i=>sigmay()), +, 1:latt.N)
Sz = mapreduce(i -> multisite_operator(latt, i=>sigmaz()), +, 1:latt.N)

H, c_ops = DissipativeIsing(Jx, Jy, Jz, hx, hy, hz, γ, latt; boundary_condition = Val(:periodic_bc), order = 1)
e_ops = [Sx, Sy, Sz]

# Time Evolution

ψ0 = fock(2^latt.N, 0, dims = ntuple(i->2, Val(latt.N)))

tlist = range(0, 10, 100)

sol_mc = mcsolve(H, ψ0, tlist, c_ops, e_ops=e_ops, ntraj=5000, ensemblealg=EnsembleSplitThreads())

##

println("FINISH!")

rmprocs(workers())
```

In this script, we first load the necessary packages for distributed computing on the cluster (`Distributed.jl` and `ClusterManagers.jl`). Thanks to the environment variables (previously defined in the SLURM script), we can define the number of tasks and the number of CPUs per task. Then, we initialize the distributed network with the `addprocs(SlurmManager(SLURM_NUM_TASKS); exeflags=exeflags, topology=:master_worker)` command. We then import the packages with the `@everywhere` macro, meaning to load them in all the workers. Moreover, in order to avoid conflicts between the multithreading of the BLAS library and the native Julia multithreading, we set the number of threads of the BLAS library to 1 with the `BLAS.set_num_threads(1)` command. More information about this can be found [here](https://docs.julialang.org/en/v1/manual/performance-tips/#man-multithreading-linear-algebra).

With the

```julia
println("Hello! You have $(nworkers()) workers with $(remotecall_fetch(Threads.nthreads, 2)) threads each.")
```

command, we test that the distributed network is correctly initialized. The `remotecall_fetch(Threads.nthreads, 2)` command returns the number of threads of the worker with ID `2`.

We then write the main part of the script, where we define the lattice through the [`Lattice`](@ref) function. We set the parameters and define the Hamiltonian and collapse operators with the [`DissipativeIsing`](@ref) function. We also define the expectation operators `e_ops` and the initial state `ψ0`. Finally, we perform the Monte Carlo quantum trajectories with the [`mcsolve`](@ref) function. The `ensemblealg=EnsembleSplitThreads()` argument is used to parallelize the Monte Carlo quantum trajectories, by splitting the ensemble of trajectories among the workers. For a more detailed explanation of the different ensemble methods, you can check the [official documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/) of the [**DifferentialEquations.jl**](https://github.com/SciML/DifferentialEquations.jl/) package. Finally, the `rmprocs(workers())` command is used to remove the workers after the computation is finished.

The output of the script will be printed in the `output.out` file, which contains an output similar to the following:

```
################
Hello! You have 10 workers with 72 threads each.
----------------
################

Progress: [==============================] 100.0% --- Elapsed Time: 0h 00m 21s (ETA: 0h 00m 00s)

FINISH!
```

where we can see that the computation **lasted only 21 seconds**.

## Master Equation by Sweeping Over Parameters

In this example, we will consider a driven Jaynes-Cummings model, describing a two-level atom interacting with a driven cavity mode. The Hamiltonian is given by

```math
\hat{H} = \omega_c \hat{a}^\dagger \hat{a} + \frac{\omega_q}{2} \hat{\sigma}_z + g (\hat{a} \hat{\sigma}_+ + \hat{a}^\dagger \hat{\sigma}_-) + F \cos(\omega_d t) (\hat{a} + \hat{a}^\dagger) \, ,
```

and the collapse operators are given by

```math
\hat{c}_1 = \sqrt{\gamma} \hat{a} \, , \quad \hat{c}_2 = \sqrt{\gamma} \hat{\sigma}_- \, .
```

The SLURM batch script file is the same as before, but the `script.jl` file now contains the following content:

```julia
using Distributed
using ClusterManagers

const SLURM_NUM_TASKS = parse(Int, ENV["SLURM_NTASKS"])
const SLURM_CPUS_PER_TASK = parse(Int, ENV["SLURM_CPUS_PER_TASK"])

exeflags = ["--project=.", "-t $SLURM_CPUS_PER_TASK"]
addprocs(SlurmManager(SLURM_NUM_TASKS); exeflags=exeflags, topology=:master_worker)


println("################")
println("Hello! You have $(nworkers()) workers with $(remotecall_fetch(Threads.nthreads, 2)) threads each.")

println("----------------")


println("################")

flush(stdout)

@everywhere begin
    using QuantumToolbox
    using OrdinaryDiffEq

    BLAS.set_num_threads(1)
end

@everywhere begin
    const Nc = 20
    const ωc = 1.0
    const g = 0.05
    const γ = 0.01
    const F = 0.01

    const a = tensor(destroy(Nc), qeye(2))

    const σm = tensor(qeye(Nc), sigmam())
    const σp = tensor(qeye(Nc), sigmap())

    H(ωq) = ωc * a' * a + ωq * tensor(num(Nc), qeye(2)) + g * (a' * σm + a * σp)

    coef(p, t) = p.F * cos(p.ωd * t) # coefficient for the time-dependent term

    const c_ops = [sqrt(γ) * a, sqrt(γ) * σm]
    const e_ops = [a' * a]
end

# Define the ODE problem and the EnsembleProblem generation function

@everywhere begin
    ωq_list = range(ωc - 3*g, ωc + 3*g, 100)
    ωd_list = range(ωc - 3*g, ωc + 3*g, 100)

    const iter = collect(Iterators.product(ωq_list, ωd_list))

    function my_prob_func(prob, i, repeat, channel)
        ωq, ωd = iter[i]
        H_i = H(ωq)
        H_d_i = H_i + QobjEvo(a + a', coef) # Hamiltonian with a driving term

        L = liouvillian(H_d_i, c_ops).data # Make the Liouvillian

        put!(channel, true) # Update the progress bar channel

        remake(prob, f=L, p=(F = F, ωd = ωd))
    end
end

ωq, ωd = iter[1]
H0 = H(ωq) + QobjEvo(a + a', coef)
ψ0 = tensor(fock(Nc, 0), basis(2, 1)) # Ground State
tlist = range(0, 20 / γ, 1000)

prob = mesolveProblem(H0, ψ0, tlist, c_ops, e_ops=e_ops, progress_bar=Val(false), params=(F = F, ωd = ωd))

### Just to print the progress bar
progr = ProgressBar(length(iter))
progr_channel::RemoteChannel{Channel{Bool}} = RemoteChannel(() -> Channel{Bool}(1))
###
ens_prob = EnsembleProblem(prob.prob, prob_func=(prob, i, repeat) -> my_prob_func(prob, i, repeat, progr_channel))


@sync begin
    @async while take!(progr_channel)
        next!(progr)
    end

    @async begin
        sol = solve(ens_prob, Tsit5(), EnsembleSplitThreads(), trajectories = length(iter))
        put!(progr_channel, false)
    end
end
```

We are using the [`mesolveProblem`](@ref) function to define the master equation problem. We added some code to manage the progress bar, which is updated through a `RemoteChannel`. The `prob_func` argument of the `EnsembleProblem` function is used to define the function that generates the problem for each iteration. The `iter` variable contains the product of the `ωq_list` and `ωd_list` lists, which are used to sweep over the parameters. The `sol = solve(ens_prob, Tsit5(), EnsembleDistributed(), trajectories=length(iter))` command is used to solve the problem with the distributed ensemble method. The output of the script will be printed in the `output.out` file, which contains an output similar to the previous example.


