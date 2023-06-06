# Benchmarks

```@example
using Pkg
Pkg.add("PkgBenchmark")
Pkg.add("BenchmarkTools")

using PkgBenchmark

result = benchmarkpkg("QuPhys")

mdtext = PkgBenchmark.export_markdown(result)

println(mdtext)
```