```@meta
CurrentModule = QuPhys
```

```@eval
using PkgBenchmark
using Markdown

nthreads = div(Sys.CPU_THREADS, 2)

config1 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => 1))
config4 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => nthreads))

judgement_14 = judge("QuPhys", config4, config1)

mdtext = sprint(io -> PkgBenchmark.export_markdown(io, judgement_14))
Markdown.parse(mdtext)
```