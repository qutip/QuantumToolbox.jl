window.BENCHMARK_DATA = {
  "lastUpdate": 1715021220106,
  "repoUrl": "https://github.com/albertomercurio/QuantumToolbox.jl",
  "entries": {
    "Benchmark results": [
      {
        "commit": {
          "author": {
            "email": "61953577+albertomercurio@users.noreply.github.com",
            "name": "Alberto Mercurio",
            "username": "albertomercurio"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c79aae033cf1d92e3ded58411b3ceb740118a10c",
          "message": "Merge pull request #96 from albertomercurio/dev/patch-3\n\nAdd benchmark tracking",
          "timestamp": "2024-05-05T11:56:44+02:00",
          "tree_id": "3f7caf2742be57bc22b286a04ae0d2138fbe04a3",
          "url": "https://github.com/albertomercurio/QuantumToolbox.jl/commit/c79aae033cf1d92e3ded58411b3ceb740118a10c"
        },
        "date": 1714903359775,
        "tool": "julia",
        "benches": [
          {
            "name": "steadystate/driven-dissipative harmonic oscillator",
            "value": 7402800.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8686048\nallocs=250\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "61953577+albertomercurio@users.noreply.github.com",
            "name": "Alberto Mercurio",
            "username": "albertomercurio"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bb7fb7ca4374cd49918145e75d7e6bf8ddc94922",
          "message": "Merge pull request #97 from albertomercurio/dev/patch-3\n\nChange README Zenodo Badge",
          "timestamp": "2024-05-05T12:25:33+02:00",
          "tree_id": "875d7271688776efb5302232bfca8dbf34f4c403",
          "url": "https://github.com/albertomercurio/QuantumToolbox.jl/commit/bb7fb7ca4374cd49918145e75d7e6bf8ddc94922"
        },
        "date": 1714904791570,
        "tool": "julia",
        "benches": [
          {
            "name": "steadystate/driven-dissipative harmonic oscillator",
            "value": 7271233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8686048\nallocs=250\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ],
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "61953577+albertomercurio@users.noreply.github.com",
            "name": "Alberto Mercurio",
            "username": "albertomercurio"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "82caa759ade1210c600d0d586cd9a09121d6ff3a",
          "message": "Merge pull request #98 from albertomercurio/dev/patch-3\n\nAdd more benchmarks",
          "timestamp": "2024-05-05T23:02:56+02:00",
          "tree_id": "2b8f18c49f36966eb3978c0658a47231ebb21070",
          "url": "https://github.com/albertomercurio/QuantumToolbox.jl/commit/82caa759ade1210c600d0d586cd9a09121d6ff3a"
        },
        "date": 1714943483430,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6781976,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11510696\nallocs=344\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602917388.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 42230481.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6809344\nallocs=226\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 12316949,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1360888\nallocs=968\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Exponential Series",
            "value": 43159591,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5167600\nallocs=8732\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1415559050,
            "unit": "ns",
            "extra": "gctime=26973874\nmemory=3019323376\nallocs=435121\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4479039,
            "unit": "ns",
            "extra": "gctime=0\nmemory=102064\nallocs=114\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 24914752,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3695568\nallocs=468\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 184610666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3429840\nallocs=9583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 93746425,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3431024\nallocs=9594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "61953577+albertomercurio@users.noreply.github.com",
            "name": "Alberto Mercurio",
            "username": "albertomercurio"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d04d1da461982ec42324bdf5a16001928eab49b5",
          "message": "Merge pull request #100 from albertomercurio/dependabot/github_actions/julia-actions/cache-2\n\nBump julia-actions/cache from 1 to 2",
          "timestamp": "2024-05-06T20:38:29+02:00",
          "tree_id": "afbfd01c85e2b4eb9da5720bd5052d836f24663f",
          "url": "https://github.com/albertomercurio/QuantumToolbox.jl/commit/d04d1da461982ec42324bdf5a16001928eab49b5"
        },
        "date": 1715021219613,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6801200,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11510696\nallocs=344\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593703698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 48227138.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6809344\nallocs=226\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 12223917.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1360888\nallocs=968\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Exponential Series",
            "value": 43128621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5167600\nallocs=8732\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1372021563,
            "unit": "ns",
            "extra": "gctime=18395318\nmemory=3019323376\nallocs=435121\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4535977.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=102064\nallocs=114\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 24932344.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3695568\nallocs=468\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 183919863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3429840\nallocs=9583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 94075418,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3431024\nallocs=9594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}