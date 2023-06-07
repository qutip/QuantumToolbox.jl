```@meta
CurrentModule = QuPhys
```

```@eval
using PkgBenchmark
using Markdown

nthreads = Sys.CPU_THREADS

config1 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => 1))
config2 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => nthreads))

judgement_12 = judge("QuPhys", config2, config1)

mdtext = sprint(io -> PkgBenchmark.export_markdown(io, judgement_12))
Markdown.parse(mdtext)
```