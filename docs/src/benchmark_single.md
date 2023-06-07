```@meta
CurrentModule = QuPhys
```

```@eval
using Pkg
Pkg.add("PkgBenchmark")
Pkg.add("BenchmarkTools")

using PkgBenchmark
using Markdown

nthreads = div(Sys.CPU_THREADS, 2)

config1 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => 1))
config4 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => nthreads))

result1 = benchmarkpkg("QuPhys", config1)
judgement_14 = judge("QuPhys", config4, config1)

mdtext = sprint(io -> PkgBenchmark.export_markdown(io, result1))
Markdown.parse(mdtext)
```
