```@meta
CurrentModule = QuPhys
```

```@eval
using PkgBenchmark
using Markdown


config1 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => 1))

result1 = benchmarkpkg("QuPhys", config1)

mdtext = sprint(io -> PkgBenchmark.export_markdown(io, result1))
Markdown.parse(mdtext)
```
