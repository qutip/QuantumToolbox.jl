window.BENCHMARK_DATA = {
  "lastUpdate": 1759045295077,
  "repoUrl": "https://github.com/qutip/QuantumToolbox.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5b08c046bff9a758ceac9d4f0ced77d98afffbdf",
          "message": "Improve `show` for `Qobj`/`QobjEvo` lists (#365)\n\n* improve `show` for `Qobj` lists\r\n\r\n* fix doctests",
          "timestamp": "2025-01-07T11:53:27+01:00",
          "tree_id": "3595f9102fd607dff850fafdaba21552c93a4a79",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5b08c046bff9a758ceac9d4f0ced77d98afffbdf"
        },
        "date": 1736247667316,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6902396.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 607241840,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67420246.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 34595901.5,
            "unit": "ns",
            "extra": "gctime=11251472.5\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42962141,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14128636.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1733012409.5,
            "unit": "ns",
            "extra": "gctime=387348642\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4392683,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26912489,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221149562,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113007043.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3564304\nallocs=12474\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4f81980df99e0f89edb7139aeb6c72fa71cee239",
          "message": "Introduce `Space`, `Dimensions`, and `GeneralDimensions` (#360)",
          "timestamp": "2025-01-13T10:30:16+01:00",
          "tree_id": "69ca9ab0fafb160acb2319707ac8343f03bd6a32",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/4f81980df99e0f89edb7139aeb6c72fa71cee239"
        },
        "date": 1736780983345,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6817421,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606599208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67842324.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31041822,
            "unit": "ns",
            "extra": "gctime=8311574\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44037111.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13898284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863496\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1787204018,
            "unit": "ns",
            "extra": "gctime=365657391\nmemory=3066309904\nallocs=303537\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4350116,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106208\nallocs=105\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27514798,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104272\nallocs=606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222322539,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563792\nallocs=12483\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112564014,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565088\nallocs=12495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2db0455da2960a78b2d5a2d791cff8ada4e58cf5",
          "message": "CompatHelper: bump compat for CairoMakie in [weakdeps] to 0.13, (keep existing compat) (#369)\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>",
          "timestamp": "2025-01-15T14:40:53+01:00",
          "tree_id": "e3e2724a4418a50de858f6baf09af9f285b42b2b",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2db0455da2960a78b2d5a2d791cff8ada4e58cf5"
        },
        "date": 1736948850667,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6874906,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602618002,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 81002589.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30119506,
            "unit": "ns",
            "extra": "gctime=7862591.5\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43715024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13950394,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863496\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1740049284,
            "unit": "ns",
            "extra": "gctime=343367536\nmemory=3066309904\nallocs=303537\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4390009.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106208\nallocs=105\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27732087,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104272\nallocs=606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220643964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563792\nallocs=12483\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112350923.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565088\nallocs=12495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "88bf11b14f47ff94c5cc246d9b28ca6c03b97e51",
          "message": "Improve lazy tensor warning for `SciMLOperators` (#370)",
          "timestamp": "2025-01-18T14:36:41+09:00",
          "tree_id": "b66716b819af255c1f6215182fd5e157bd19f3fe",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/88bf11b14f47ff94c5cc246d9b28ca6c03b97e51"
        },
        "date": 1737179037226,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6759204,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593360144,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78965886,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28046372,
            "unit": "ns",
            "extra": "gctime=6969463\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43551347,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13938629,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863496\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1673086304,
            "unit": "ns",
            "extra": "gctime=333698686\nmemory=3066309904\nallocs=303537\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4380146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106208\nallocs=105\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27517426.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104272\nallocs=606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219901730,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563792\nallocs=12483\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112135367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565088\nallocs=12495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b69f2b29a143d9d109e83eebbc4251ddd84d64eb",
          "message": "Change order of `AbstractQuantumObject` data type (#371)",
          "timestamp": "2025-01-20T11:43:40+09:00",
          "tree_id": "40724186e6acc065df8bd4ec2963609180c88daa",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/b69f2b29a143d9d109e83eebbc4251ddd84d64eb"
        },
        "date": 1737341459364,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6850953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 605059035,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 71908137.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30722480.5,
            "unit": "ns",
            "extra": "gctime=8417763\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43893746,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14229577,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863496\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1779853843,
            "unit": "ns",
            "extra": "gctime=341341005\nmemory=3066309904\nallocs=303537\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4431837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106208\nallocs=105\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27709359.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104272\nallocs=606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 223139830,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563792\nallocs=12483\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113168361,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565088\nallocs=12495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d0f0b11620c709cdcb423275425d890279aa3b55",
          "message": "bump version to `v0.25.0` (#372)",
          "timestamp": "2025-01-20T11:52:14+09:00",
          "tree_id": "9b4e548c815574c4f944f8c2e6f469f533cbd145",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d0f0b11620c709cdcb423275425d890279aa3b55"
        },
        "date": 1737341770072,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6978183,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 607033759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80511824.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29233662,
            "unit": "ns",
            "extra": "gctime=8529001.5\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43167712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13808393.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863496\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1691150043,
            "unit": "ns",
            "extra": "gctime=357143073\nmemory=3066309904\nallocs=303537\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4325850.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106208\nallocs=105\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26968512,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104272\nallocs=606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221424056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563792\nallocs=12483\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112343181,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565088\nallocs=12495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "9cc02248eb7fefb8442d4db9bc545090fbd88842",
          "message": "Fix dynamical fock dimension saving (#375)",
          "timestamp": "2025-01-26T09:21:13+01:00",
          "tree_id": "32d3ddf4c137d9c1410db8af28bbc819762d8ec4",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9cc02248eb7fefb8442d4db9bc545090fbd88842"
        },
        "date": 1737880245031,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6923272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 607116206,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 69427647.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 36240949,
            "unit": "ns",
            "extra": "gctime=12507349\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43029804,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13848113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863496\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1776101519,
            "unit": "ns",
            "extra": "gctime=402709034\nmemory=3066311728\nallocs=303594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4344759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106208\nallocs=105\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27076977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104272\nallocs=606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220777566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563792\nallocs=12483\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112450095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565088\nallocs=12495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2f0ce8efd16d2e091a121993af8a1d1964f50d27",
          "message": "Support a list of observables for `expect` (#376)",
          "timestamp": "2025-01-29T11:45:41+01:00",
          "tree_id": "ba9f8b3c848fb07750dd199bee6b3a24b6eec1d1",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2f0ce8efd16d2e091a121993af8a1d1964f50d27"
        },
        "date": 1738148087698,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6891925,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596543171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67179878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31129524,
            "unit": "ns",
            "extra": "gctime=9883046\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43746439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14186129,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863496\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1736022329,
            "unit": "ns",
            "extra": "gctime=344771940\nmemory=3066311728\nallocs=303594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4397247,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106208\nallocs=105\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27712780,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104272\nallocs=606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 223737812,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563792\nallocs=12483\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113499960.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565088\nallocs=12495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "de31edb0d738fa39a49f368792ce7f80bd99a933",
          "message": "Check `tlist` properties (#378)",
          "timestamp": "2025-01-29T12:10:37+01:00",
          "tree_id": "f201ceb1bd724d969b04b769a9ab8270be200016",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/de31edb0d738fa39a49f368792ce7f80bd99a933"
        },
        "date": 1738149272782,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6883661,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597407846,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80269042.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30916159.5,
            "unit": "ns",
            "extra": "gctime=9965774\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43960208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14214751,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863560\nallocs=715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1778668580,
            "unit": "ns",
            "extra": "gctime=379289928\nmemory=3066311792\nallocs=303595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4407680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27132751,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104336\nallocs=607\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222396520,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112803999,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "9ca0815ad0fd8b913d3c0b479fbffe45b6b114a4",
          "message": "Bump version to v0.25.1 (#379)",
          "timestamp": "2025-01-29T12:17:00+01:00",
          "tree_id": "d77f1ac790c31bc5576571bed9a299862cd999c5",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9ca0815ad0fd8b913d3c0b479fbffe45b6b114a4"
        },
        "date": 1738149653497,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6872685,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598479128,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 76994118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30196395,
            "unit": "ns",
            "extra": "gctime=9974678\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43824716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13935342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863560\nallocs=715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1777486052,
            "unit": "ns",
            "extra": "gctime=357131325\nmemory=3066311792\nallocs=303595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4358653.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27671143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104336\nallocs=607\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221885770,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112999346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bf7145d7b27a6dcb1c205fe6d6b8d97ab30fb48f",
          "message": "Move code quality dependencies to separate environment (#380)",
          "timestamp": "2025-01-30T11:15:42+09:00",
          "tree_id": "00d701f33f4b1fcee35b8f43830bb8491ce81df6",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/bf7145d7b27a6dcb1c205fe6d6b8d97ab30fb48f"
        },
        "date": 1738203752674,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6876493.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597266323,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 81137636.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30442698.5,
            "unit": "ns",
            "extra": "gctime=9529053\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43784282.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14178030.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863560\nallocs=715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1731331745,
            "unit": "ns",
            "extra": "gctime=345895856\nmemory=3066311792\nallocs=303595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4414297,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27620065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104336\nallocs=607\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221968264,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113116425,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "8cbd5499f0473b638f8c58eb15c5532789473a3d",
          "message": "Add state normalization during evolution in ssesolve (#383)",
          "timestamp": "2025-02-02T13:35:26+09:00",
          "tree_id": "3baf0d7f8b0fe8c1f6bdb7faea22ea0f2f00954b",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/8cbd5499f0473b638f8c58eb15c5532789473a3d"
        },
        "date": 1738471360758,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6740104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598672497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80515386,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30727944,
            "unit": "ns",
            "extra": "gctime=9711174\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43786163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13864349,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863560\nallocs=715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1773708310,
            "unit": "ns",
            "extra": "gctime=365263886\nmemory=3066311792\nallocs=303595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4317628.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27404420.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104336\nallocs=607\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221182059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112125017,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "504ac18fc145e355a5237244b098d71f30e07c06",
          "message": "Bump version to v0.25.2 (#384)",
          "timestamp": "2025-02-02T13:51:08+09:00",
          "tree_id": "c57c2a7eec9adbaaae8ecfbb55bc4308e4a1f785",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/504ac18fc145e355a5237244b098d71f30e07c06"
        },
        "date": 1738472098562,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6742426,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602894639,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79859770,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29981110,
            "unit": "ns",
            "extra": "gctime=9221631\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43724931,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13903907,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863560\nallocs=715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1708693424,
            "unit": "ns",
            "extra": "gctime=355713460\nmemory=3066311792\nallocs=303595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4350758,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27230173,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104336\nallocs=607\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221313892,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112074807,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "04080696f00bd6228dec1b076efa4814ed3010ba",
          "message": "Fix CUDA sparse_to_dense (#386)",
          "timestamp": "2025-02-08T16:14:48+09:00",
          "tree_id": "a0f1764412772978c6f27af390a891fb0d7c8846",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/04080696f00bd6228dec1b076efa4814ed3010ba"
        },
        "date": 1738999306628,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6705479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594585619,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 68734680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29483623,
            "unit": "ns",
            "extra": "gctime=8536860\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43544687.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13970998,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863560\nallocs=715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1724744366,
            "unit": "ns",
            "extra": "gctime=339877682\nmemory=3066311792\nallocs=303595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4328496.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27353646,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104336\nallocs=607\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219917439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 111660886,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "c5ee41780da5d775071274305f40212c33328c87",
          "message": "Improve pseudo inverse spectrum (#388)\n\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>",
          "timestamp": "2025-02-08T16:23:45+09:00",
          "tree_id": "1a0d8bf92ec0ff1231d796df97825ae2827e3e28",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/c5ee41780da5d775071274305f40212c33328c87"
        },
        "date": 1738999660878,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6851314,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284712\nallocs=469\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598276947,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79690995.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 34729581.5,
            "unit": "ns",
            "extra": "gctime=11697393\nmemory=121473128\nallocs=104744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43542036.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791888\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13872473,
            "unit": "ns",
            "extra": "gctime=0\nmemory=863560\nallocs=715\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1802566420,
            "unit": "ns",
            "extra": "gctime=369466121\nmemory=3066311792\nallocs=303595\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4325502,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27309901.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104336\nallocs=607\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221026481,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113204735.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "1886ccaa8bd08d18c8d3149620a1ab59784bcfb3",
          "message": "Add smesolve solver (#389)\n\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>",
          "timestamp": "2025-02-09T12:57:10+09:00",
          "tree_id": "9db5687fd3cd6de2cbc0231745930951b687986a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/1886ccaa8bd08d18c8d3149620a1ab59784bcfb3"
        },
        "date": 1739073791720,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6857963,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601002439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67860102,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30080978,
            "unit": "ns",
            "extra": "gctime=9282489\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43716291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13841563,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1759954126,
            "unit": "ns",
            "extra": "gctime=344889331\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4327564.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27304904,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221594699,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112232737,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "8c728f08107c8ef5dcbe006b8d044eed907fdd70",
          "message": "Bump to v0.26.0 (#390)",
          "timestamp": "2025-02-09T13:04:52+09:00",
          "tree_id": "075fc4f052ef0f911bd556982ecd29d41f46693c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/8c728f08107c8ef5dcbe006b8d044eed907fdd70"
        },
        "date": 1739074124044,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6805072,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602025733,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80334308,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 33484409.5,
            "unit": "ns",
            "extra": "gctime=11401549.5\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43718351,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14245719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1721632973,
            "unit": "ns",
            "extra": "gctime=347641225\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4411059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27738216,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221315239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112692155,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "162c76be2c6f9db79d653a077aa19727b78e4c63",
          "message": "Rename function `sparse_to_dense` as `to_dense` and `dense_to_sparse` as `to_sparse` (#392)",
          "timestamp": "2025-02-10T14:41:36+09:00",
          "tree_id": "7a2357ca6736240770f54d39fbbdbcc2fe681c6c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/162c76be2c6f9db79d653a077aa19727b78e4c63"
        },
        "date": 1739166521655,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6810687,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606627522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80600050,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 32244914.5,
            "unit": "ns",
            "extra": "gctime=10484146.5\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44051548,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13961362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1751442562,
            "unit": "ns",
            "extra": "gctime=350929068\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4376753,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27690116,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222551849,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112798500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1a2479e5eb4dc966f0bff1ecbfc3a540549f6e69",
          "message": "Minor changes in documentation and docstrings (#391)",
          "timestamp": "2025-02-10T14:43:02+09:00",
          "tree_id": "82a6b7242f6f49802131e6771cea4f95d2981f12",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/1a2479e5eb4dc966f0bff1ecbfc3a540549f6e69"
        },
        "date": 1739166605793,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6862120,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598414402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67991653,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31222097,
            "unit": "ns",
            "extra": "gctime=9734524\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43810103.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13889840,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1783952729,
            "unit": "ns",
            "extra": "gctime=354892462\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4340722,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27508823,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222558226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112578647.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2b6e0d67cb2dda093935be465a0541832548e790",
          "message": "Fix of erroneous definition of stochastic term in  `smesolve` (#393)",
          "timestamp": "2025-02-11T10:15:06+09:00",
          "tree_id": "6881832a99022cb79209f54ec6cdf35ad37df7f5",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2b6e0d67cb2dda093935be465a0541832548e790"
        },
        "date": 1739237050402,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6892255,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601594709,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80804148,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30620490,
            "unit": "ns",
            "extra": "gctime=9409963\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43803440,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13951573,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1683403277,
            "unit": "ns",
            "extra": "gctime=334288469\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4370896,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 28327076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222950014,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113006637,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "c984f9536814e4bd61f419bca4f0283164ffe716",
          "message": "Change MultiSiteOperator function name to multisite_operator (#394)\n\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>",
          "timestamp": "2025-02-11T14:45:29+09:00",
          "tree_id": "ed94bed84d2c301764849867102a18c9644c7239",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/c984f9536814e4bd61f419bca4f0283164ffe716"
        },
        "date": 1739253074289,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6856302.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597794347,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80941383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29216003,
            "unit": "ns",
            "extra": "gctime=8891496\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43586660,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13963475,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1733213671,
            "unit": "ns",
            "extra": "gctime=336137321\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4368613,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26941367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220668852,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113435718.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2fbace49430ad480386e3b979e39f199ba1ae2b0",
          "message": "Fix stochastic solvers and add documentation page (#395)",
          "timestamp": "2025-02-11T19:09:38+09:00",
          "tree_id": "1b951f75af39152c3521692970b8190d909fd92a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2fbace49430ad480386e3b979e39f199ba1ae2b0"
        },
        "date": 1739268813638,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6955019.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600510495,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79409963,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31517944,
            "unit": "ns",
            "extra": "gctime=9895908\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43781303.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14060275,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1716233919,
            "unit": "ns",
            "extra": "gctime=341754625\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4380332.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27839670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219747543,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 111792766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "aaaf0d333107e2d24513c57bc03dabf42b3a2367",
          "message": "Fix time evolution output when using `saveat` (#398)",
          "timestamp": "2025-02-13T09:27:34+09:00",
          "tree_id": "07990e96e1f4f28ee5e9f256f34a363d64e5b200",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/aaaf0d333107e2d24513c57bc03dabf42b3a2367"
        },
        "date": 1739407003837,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6826307,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601952668,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80272504,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30944146,
            "unit": "ns",
            "extra": "gctime=9590087\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43697402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13990552,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1766694674,
            "unit": "ns",
            "extra": "gctime=362549879\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4355299,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27427689,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222266137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112397283,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "28fe4a902cb986209f6d2bb1c0f7816406bc271e",
          "message": "CompatHelper: bump compat for LinearSolve to 3, (keep existing compat) (#387)",
          "timestamp": "2025-02-13T09:53:18+09:00",
          "tree_id": "8ea237ad25180b49e20f1253b0104333d4ed3932",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/28fe4a902cb986209f6d2bb1c0f7816406bc271e"
        },
        "date": 1739408229778,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6964646.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595966726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80241267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30613674.5,
            "unit": "ns",
            "extra": "gctime=9855160\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43978617.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13941044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1772883788,
            "unit": "ns",
            "extra": "gctime=347603615\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4362220,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27364027.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220955752,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112340743,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "0ce0b452128c569dc01291b4f94e6bba9bc01b07",
          "message": "Add more generic solver for steadystate_floquet (#396)",
          "timestamp": "2025-02-13T19:23:55+09:00",
          "tree_id": "a3237f894701b08cd6ec4b11a6d57019ace17dcc",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/0ce0b452128c569dc01291b4f94e6bba9bc01b07"
        },
        "date": 1739442472406,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6818899,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596735502,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80230743,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30100394.5,
            "unit": "ns",
            "extra": "gctime=9568983\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44392363,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14199982.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1680229825,
            "unit": "ns",
            "extra": "gctime=336733682\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4426159,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27724231,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221588777,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112914201,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "dff6aa17cb6a5e9effd228d086dcbdf67001dcde",
          "message": "Improve ensemble generation of ssesolve and change parameters on stochastic processes (#403)",
          "timestamp": "2025-02-13T23:27:56+09:00",
          "tree_id": "1b63a9289d44658c9035c410f55d13a96c339ecf",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/dff6aa17cb6a5e9effd228d086dcbdf67001dcde"
        },
        "date": 1739457311958,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6888923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600982401,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 66779266,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 34850278.5,
            "unit": "ns",
            "extra": "gctime=12033364\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42954039,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13893777.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1766435152,
            "unit": "ns",
            "extra": "gctime=378326856\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4340230,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27106390,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219297499,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112892920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "ed5b0da12f273282e11d939196e3aba386aeb51d",
          "message": "Align some attributes of `mcsolve`, `ssesolve` and `smesolve` results with QuTiP (#402)",
          "timestamp": "2025-02-14T00:12:31+09:00",
          "tree_id": "3ecbe9e375fd40c9a338eb7683588702f6a20097",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/ed5b0da12f273282e11d939196e3aba386aeb51d"
        },
        "date": 1739459941371,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6846988,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598213841,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67672640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 37168819.5,
            "unit": "ns",
            "extra": "gctime=13353248\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44073579,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13950125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1670853117,
            "unit": "ns",
            "extra": "gctime=360277105\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4336149,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27449478,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 224340380,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112236454,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "81724addb8acf774911452066ebddbc4107d3e42",
          "message": "Set default trajectories to 500 and rename the keyword argument `ensemble_method` to `ensemblealg` (#405)",
          "timestamp": "2025-02-14T15:37:16+09:00",
          "tree_id": "664ff2d9a56034718990e1231d708af7b2d8a430",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/81724addb8acf774911452066ebddbc4107d3e42"
        },
        "date": 1739516378524,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6712795.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596443049,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78519195.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29312525,
            "unit": "ns",
            "extra": "gctime=9019861\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44075987,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13892447.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1774352998,
            "unit": "ns",
            "extra": "gctime=364918384\nmemory=3212054576\nallocs=304450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4361036.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27508381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220254454,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112153386.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "49e9ffc109c4f3ec86dcaf58efaceeb3f75cdee9",
          "message": "Introduce measurement on `ssesolve` and `smesolve` (#404)",
          "timestamp": "2025-02-14T17:47:12+09:00",
          "tree_id": "3aa0ff1e663a777619db965be16fd4c62e77d4f8",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/49e9ffc109c4f3ec86dcaf58efaceeb3f75cdee9"
        },
        "date": 1739523064811,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6713795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601222491,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 69125916.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30747695,
            "unit": "ns",
            "extra": "gctime=9562548\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43727926,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13950805,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1734953500,
            "unit": "ns",
            "extra": "gctime=349301112\nmemory=3212054688\nallocs=304452\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4368963,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27715347.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221014570,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112882305.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "3029a54a22313b48c692c67d5a55fb5bfce47d80",
          "message": "Bump version to v0.27.0 (#406)",
          "timestamp": "2025-02-14T18:08:32+09:00",
          "tree_id": "ea9ffd25a0448e188645b9e7046f0c585a22e438",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/3029a54a22313b48c692c67d5a55fb5bfce47d80"
        },
        "date": 1739524341617,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6920385,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596822327,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80419239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30234789.5,
            "unit": "ns",
            "extra": "gctime=9460070.5\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43944294.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13980287,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1716913669,
            "unit": "ns",
            "extra": "gctime=350558417\nmemory=3212054688\nallocs=304452\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4325340,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27545489,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219583334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112276775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "9e3ddf36fe8060415e4395baf93e455286dcea16",
          "message": "Support for single `sc_ops` for faster specific method in `ssesolve` and `smesolve` (#408)",
          "timestamp": "2025-02-17T10:07:55+01:00",
          "tree_id": "bac794af9984810a13aa0a8e4616f32942682d56",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9e3ddf36fe8060415e4395baf93e455286dcea16"
        },
        "date": 1739784075281,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7080188,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603712553,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 77521669.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 33597822,
            "unit": "ns",
            "extra": "gctime=10814494\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44086885,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14296764,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1745582002,
            "unit": "ns",
            "extra": "gctime=351357158\nmemory=3212054688\nallocs=304452\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4405393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27751107,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222945898,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113429383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "d48371cc48d2138ca3d65bee396950a9e53f43bf",
          "message": "Align `eigenstates` and `eigenenergies` to QuTiP (#411)",
          "timestamp": "2025-02-18T08:56:03+09:00",
          "tree_id": "a46403c35a2ceb4abb284a01678240ab3cabcc3b",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d48371cc48d2138ca3d65bee396950a9e53f43bf"
        },
        "date": 1739836952997,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6859860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606848330,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 72587226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31706574,
            "unit": "ns",
            "extra": "gctime=9908464\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44085953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13878717,
            "unit": "ns",
            "extra": "gctime=0\nmemory=882856\nallocs=728\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1760647848,
            "unit": "ns",
            "extra": "gctime=343666486\nmemory=3212054688\nallocs=304452\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4352861,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106272\nallocs=106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27640099,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3247760\nallocs=622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220689079,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563856\nallocs=12484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112370711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3565152\nallocs=12496\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "300fd5de584c86eaca47e89a685e693e9a75c55f",
          "message": "Change save callbacks from `PresetTimeCallback` to `FunctionCallingCallback` (#410)\n\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>",
          "timestamp": "2025-02-18T11:50:24+09:00",
          "tree_id": "35fcf3356e072bb12a06991fa746507b4d32ff68",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/300fd5de584c86eaca47e89a685e693e9a75c55f"
        },
        "date": 1739847258344,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7076538,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597076856,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67728222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 37961272,
            "unit": "ns",
            "extra": "gctime=13900099\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42728315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7841779.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1443894236,
            "unit": "ns",
            "extra": "gctime=360099631\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1685210,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17788954,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152635570,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77437273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7d58e4690efbe53580553e9c4659ebc82d0d1eed",
          "message": "Introduce `operator_to_vector` and `vector_to_operator` (#413)",
          "timestamp": "2025-02-18T08:30:42+01:00",
          "tree_id": "52110a81a4febfb7f221a2a73981c4f5cb581bdf",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/7d58e4690efbe53580553e9c4659ebc82d0d1eed"
        },
        "date": 1739864073599,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6787981,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596508446,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80183277,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30725135,
            "unit": "ns",
            "extra": "gctime=9854693.5\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44279399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7814820,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1448933940,
            "unit": "ns",
            "extra": "gctime=340549246\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1656445,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17895211,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150794197.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77381630,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "46b98de36105478d87b81acd5de02e9d04837c5e",
          "message": "Introduce some entropy related functions (#416)",
          "timestamp": "2025-02-20T10:41:54+09:00",
          "tree_id": "c898a1d463d7736a3f251686a22ee4e06634399a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/46b98de36105478d87b81acd5de02e9d04837c5e"
        },
        "date": 1740016148014,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6804716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598441749,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79704489,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29640886,
            "unit": "ns",
            "extra": "gctime=9281480\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43072286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7698190,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1378714578,
            "unit": "ns",
            "extra": "gctime=332995474.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1686175,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17777326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150013493,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76943295,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "71b8d8cbaaa190e2a0c2abf97fd663d51bade668",
          "message": "Fix `entanglement` and introduce `concurrence` (#419)",
          "timestamp": "2025-02-20T21:58:34+01:00",
          "tree_id": "8aecfa5aece96b9173185b0a18e5838f48185334",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/71b8d8cbaaa190e2a0c2abf97fd663d51bade668"
        },
        "date": 1740098760660,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6852068.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598981326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80013397,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28305156,
            "unit": "ns",
            "extra": "gctime=8748710\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43698032,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7918793,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1433350329.5,
            "unit": "ns",
            "extra": "gctime=331135682\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1675948,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17805021,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151340833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77143085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6e14238233b8db30a3d080faad920739aeadf74b",
          "message": "Introduce some metric functions (#420)",
          "timestamp": "2025-02-21T10:04:25+01:00",
          "tree_id": "4936197cc897f190ca209c993896a017d90eaa9b",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/6e14238233b8db30a3d080faad920739aeadf74b"
        },
        "date": 1740129089647,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7030069,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600451364,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80231368.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 39386742,
            "unit": "ns",
            "extra": "gctime=14514486\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43034957,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7719638,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1477139199.5,
            "unit": "ns",
            "extra": "gctime=364097839.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1642504,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17870367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151189384.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77740143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "2d352f9012c94bb96f82fa3804f3ef3b39b4a266",
          "message": "Align `steadystate` ODE solver and improve GPU support (#421)",
          "timestamp": "2025-02-22T10:59:50+09:00",
          "tree_id": "d5bfdbe2662334da32fd2103ca43fe9ceae53f39",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2d352f9012c94bb96f82fa3804f3ef3b39b4a266"
        },
        "date": 1740189822060,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6849852.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593677636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79480127.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 27466071,
            "unit": "ns",
            "extra": "gctime=7820776.5\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43529873,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7717178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1415336215,
            "unit": "ns",
            "extra": "gctime=341072813\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1651628,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17778889,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151632838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77538181,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "35e9c6693151724c31bd9ae8a71fa7044b46dc18",
          "message": "bump version to `v0.28.0` (#422)",
          "timestamp": "2025-02-22T11:08:42+09:00",
          "tree_id": "74df98a2bf3fc3cea7cacaf0d23786d7d023cb8c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/35e9c6693151724c31bd9ae8a71fa7044b46dc18"
        },
        "date": 1740190350958,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6898298,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597192135,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80012947,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30570171,
            "unit": "ns",
            "extra": "gctime=10148443\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43070146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7741294,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1428427095.5,
            "unit": "ns",
            "extra": "gctime=334685338\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1650493.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17924603,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150870557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77353622,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "5b65086439eed6dfc283df1168227dba333b3166",
          "message": "Add support for `OperatorKet` state input for `mesolve` and `smesolve` (#423)",
          "timestamp": "2025-02-23T14:17:58+09:00",
          "tree_id": "a636a01355c924017958401d820a4c6bb6cc2fca",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5b65086439eed6dfc283df1168227dba333b3166"
        },
        "date": 1740288288320,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6956465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601745446,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80129602,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30307036,
            "unit": "ns",
            "extra": "gctime=9664758\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42713164.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7834478.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1381737318.5,
            "unit": "ns",
            "extra": "gctime=325545066\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1662654,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17798281.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152380657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77504126,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "157601901+TendonFFF@users.noreply.github.com",
            "name": "Li-Xun Cai",
            "username": "TendonFFF"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b5b770f3368dab13beff4cdf4cc6cb577a6c488f",
          "message": "Introduce `plot_fock_distribution` (#428)",
          "timestamp": "2025-03-07T13:17:17+08:00",
          "tree_id": "0732cb8753ba270684dc52650cf5a5d05435abe3",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/b5b770f3368dab13beff4cdf4cc6cb577a6c488f"
        },
        "date": 1741325137042,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6944076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602904720,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79994039,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30816705,
            "unit": "ns",
            "extra": "gctime=9897352.5\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43847324.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7802396.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1410289514.5,
            "unit": "ns",
            "extra": "gctime=315863011.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1659063,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17918706,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150617072.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76802023,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e76fcbeeebc8143f5128e9be18af769721cb973f",
          "message": "Bump version to `v0.29.0` (#429)",
          "timestamp": "2025-03-07T13:22:51+08:00",
          "tree_id": "0a94c8e2160dd6969a84b3e868f7faad1a153031",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e76fcbeeebc8143f5128e9be18af769721cb973f"
        },
        "date": 1741325572941,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6847271.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596375996,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80081850,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30907415,
            "unit": "ns",
            "extra": "gctime=9696814\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43024291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7666892,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1433495119,
            "unit": "ns",
            "extra": "gctime=363328509\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1652677,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17737749.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150967977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76846914,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "157601901+TendonFFF@users.noreply.github.com",
            "name": "Li-Xun Cai",
            "username": "TendonFFF"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e22908f29ecc3e601e371e5eabc714e1dbadeb0b",
          "message": "move and rename eltype conversion (#430)",
          "timestamp": "2025-03-07T16:05:48+08:00",
          "tree_id": "b6f7c84a20736e71b060921a4083a8b8d6195cd4",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e22908f29ecc3e601e371e5eabc714e1dbadeb0b"
        },
        "date": 1741334974954,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6830614,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595080420,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79686242.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28477563,
            "unit": "ns",
            "extra": "gctime=8453970\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42523883,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7670757,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1410305755.5,
            "unit": "ns",
            "extra": "gctime=324162699.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1646331,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17619745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150709006,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76770672,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1c2688850740cef0fc7bf445de5d5b2a0f305924",
          "message": "Bump version to `v0.29.1` (#431)",
          "timestamp": "2025-03-07T16:11:17+08:00",
          "tree_id": "ed2788219311c5f8e3ad5ec0b3bc72c54b488523",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/1c2688850740cef0fc7bf445de5d5b2a0f305924"
        },
        "date": 1741335301327,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6875639,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 608259276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67417027,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 35276453,
            "unit": "ns",
            "extra": "gctime=12868202\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43156828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7704308,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1395470697.5,
            "unit": "ns",
            "extra": "gctime=336470017.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1652597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17906300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150151656.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77319855,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6d4e7d26f3e1d93b406cc45a8ebfcfdc8dedb9cf",
          "message": "Add reference for `concurrence` (#433)",
          "timestamp": "2025-03-11T17:15:45+09:00",
          "tree_id": "ec07b609fb38f0972b770e91c8cd70860b45a962",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/6d4e7d26f3e1d93b406cc45a8ebfcfdc8dedb9cf"
        },
        "date": 1741681435930,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6907606,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600782831,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79508908,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29046247,
            "unit": "ns",
            "extra": "gctime=8718284\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42625371,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7671381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1392260882,
            "unit": "ns",
            "extra": "gctime=313611886\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1648467.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 18352877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150745953.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76874939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "d8f83ccb7cfd7b2084eeec644e3e057a9b6ad0ef",
          "message": "Make CUDA conversion more general (#437)",
          "timestamp": "2025-04-05T08:37:55+09:00",
          "tree_id": "90101715999b3fa511b026dee4ea7164b934c365",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d8f83ccb7cfd7b2084eeec644e3e057a9b6ad0ef"
        },
        "date": 1743810369131,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6970607,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596132643,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 71153566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28628157,
            "unit": "ns",
            "extra": "gctime=8308253\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42815336,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7685650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1388349351,
            "unit": "ns",
            "extra": "gctime=323501315\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1651452.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17721428.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150711077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78007714,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "5fbce6d6e1143e547d623f93e79d9b22bd63c00d",
          "message": "Make `fock` non-mutating (#438)",
          "timestamp": "2025-04-05T17:05:55+02:00",
          "tree_id": "94a21ddd2bb7f27f649328b5105f6fa3a01fc7c9",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5fbce6d6e1143e547d623f93e79d9b22bd63c00d"
        },
        "date": 1743865852906,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6848866,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593662946,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78976653.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28287063.5,
            "unit": "ns",
            "extra": "gctime=8089167\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42508339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7674599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1382652014.5,
            "unit": "ns",
            "extra": "gctime=318905604.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1644572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17717753,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150281245.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76937035,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "7fd9e80a3460f547812e1019b244656ca7d253d6",
          "message": "Format the code according to the new version of JuliaFormatter.jl (#439)",
          "timestamp": "2025-04-06T11:13:40+02:00",
          "tree_id": "30670f6b3bfb103e87ad838bff714f623a0f13ed",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/7fd9e80a3460f547812e1019b244656ca7d253d6"
        },
        "date": 1743931050603,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6846176,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 591553430,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 66924913,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 26555213,
            "unit": "ns",
            "extra": "gctime=7239852\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42384481.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7632191.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1290892651,
            "unit": "ns",
            "extra": "gctime=308195834.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1639742,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17578639.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151211734.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77247720,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "54145bf845f96403d07413958566190bf277e989",
          "message": "Remove Re-export (#443)",
          "timestamp": "2025-04-10T11:55:57+02:00",
          "tree_id": "8fb35447dfe7ddab15720c237e97bf9618ba20a1",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/54145bf845f96403d07413958566190bf277e989"
        },
        "date": 1744279360040,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7025849,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603078457,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 68931697,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30504847,
            "unit": "ns",
            "extra": "gctime=9116500\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42792949,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7700614.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1423913697.5,
            "unit": "ns",
            "extra": "gctime=339313035\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1657870,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17735612,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151654690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76783433,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "a07b9918ca0f965510af60552688ba4dce7c7d15",
          "message": "Add support to automatic differentiation for `sesolve` and `mesolve` (#440)",
          "timestamp": "2025-04-12T14:49:03+09:00",
          "tree_id": "12a0429be8db7835715a6bdfe9b77223aba955f0",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/a07b9918ca0f965510af60552688ba4dce7c7d15"
        },
        "date": 1744437369548,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6982292,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 604044957,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80477143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30262443,
            "unit": "ns",
            "extra": "gctime=9041724\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43006483,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7697976,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1430082545.5,
            "unit": "ns",
            "extra": "gctime=326842617.5\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1645739,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17912685.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150071588.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77385965,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "6bc6e383e414cd37f567e36533bbe3848e4656ac",
          "message": "Bump to v0.30.0 (#446)",
          "timestamp": "2025-04-12T08:29:20+02:00",
          "tree_id": "957888419c5db9dc731b69354c536e3d0077394e",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/6bc6e383e414cd37f567e36533bbe3848e4656ac"
        },
        "date": 1744439589984,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6919208,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500456\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596864884,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80227453.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30301930,
            "unit": "ns",
            "extra": "gctime=9594979\nmemory=121492424\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42925074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4811184\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7818618.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852008\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1382380882.5,
            "unit": "ns",
            "extra": "gctime=326425760\nmemory=2719158672\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1675398,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75424\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17744127,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3244704\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152964259,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3327216\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77707721.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328512\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "ad8b0b88993a566537f409048907ce93d519e9e0",
          "message": "Support different length for `to` and `from` on GeneralDImensions (#448)",
          "timestamp": "2025-04-16T09:43:16+02:00",
          "tree_id": "7a1055d6511f3f10b09c207aea5505694775a48d",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/ad8b0b88993a566537f409048907ce93d519e9e0"
        },
        "date": 1744789908206,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7036822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606288005,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79590186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22951720,
            "unit": "ns",
            "extra": "gctime=5816012\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42743323,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7654986,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1367598209,
            "unit": "ns",
            "extra": "gctime=327071916\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1661197,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17731043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151968620,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77168376.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "130eafee74fe94ce9ed9ab487db3ef86cb571df8",
          "message": "Make Makie extension more general (#450)",
          "timestamp": "2025-04-19T16:22:31+02:00",
          "tree_id": "30f0534755f10daedb73a0f372225890a39d334c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/130eafee74fe94ce9ed9ab487db3ef86cb571df8"
        },
        "date": 1745072972932,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6782759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500384\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593954051,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79907648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22065941.5,
            "unit": "ns",
            "extra": "gctime=5517494\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42488446.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7655287,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1351859713,
            "unit": "ns",
            "extra": "gctime=317231481\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1640314,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17725002,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150340688,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76949723,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "68ed7c3fb810a80e54a4ee0053fa50ad1d975bc9",
          "message": "Fix definition of noise derivative in stochastic solvers (#453)",
          "timestamp": "2025-04-24T23:56:01+02:00",
          "tree_id": "082634b8d5f269b9684274810f157736c9ec1519",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/68ed7c3fb810a80e54a4ee0053fa50ad1d975bc9"
        },
        "date": 1745532275402,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7015218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603281006,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80056910.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22620200,
            "unit": "ns",
            "extra": "gctime=5675067\nmemory=121708456\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42661094,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7790338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1405633944,
            "unit": "ns",
            "extra": "gctime=339749677.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1649796,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17707420,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151852183,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78190656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "2b5e86efb5ca00a65277a4b02bda07531cd529cf",
          "message": "Bump to v0.30.1 (#454)",
          "timestamp": "2025-04-24T23:58:33+02:00",
          "tree_id": "ddfa1e3317b352d6a1a9360e9ae27531c130b41f",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2b5e86efb5ca00a65277a4b02bda07531cd529cf"
        },
        "date": 1745532430667,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7040358,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596679883,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80021013,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22444693,
            "unit": "ns",
            "extra": "gctime=5604232.5\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42733770,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7816350.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1373500362.5,
            "unit": "ns",
            "extra": "gctime=321625711.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1654627,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17815510.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152386415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78262206.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "85d23efaa18525754b8516ebf9c327f9c8d9b667",
          "message": "Return `sesolve` when `mesolve` allows it (#455)",
          "timestamp": "2025-04-30T11:36:59+09:00",
          "tree_id": "5c376665e27bf09d2f3b407ec03bd57707d71351",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/85d23efaa18525754b8516ebf9c327f9c8d9b667"
        },
        "date": 1745981031153,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6841057,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593860188,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67734599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 21510758.5,
            "unit": "ns",
            "extra": "gctime=5395468.5\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42606607.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7612554.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1379464179,
            "unit": "ns",
            "extra": "gctime=317984325\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1656903,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17633860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151211012,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77316326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "201e066b57116075923b9781e1cad3c9ecc670b3",
          "message": "Unify structure of `QuantumObjectType` (#456)",
          "timestamp": "2025-05-03T18:32:12+09:00",
          "tree_id": "3c9d1bbabbe4f34ab1cea967662e5d2a3f2ed6f5",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/201e066b57116075923b9781e1cad3c9ecc670b3"
        },
        "date": 1746265243599,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6891800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594852215,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 68238339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22250564,
            "unit": "ns",
            "extra": "gctime=5570470.5\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42595403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7645738,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1322671701.5,
            "unit": "ns",
            "extra": "gctime=314669531.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1649446,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17778444.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151088012,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77109400.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "243bff734b1e13aed78847ba32a778b924c9dc8d",
          "message": "Bump to v0.31.0 (#459)",
          "timestamp": "2025-05-03T22:02:21+02:00",
          "tree_id": "94774624c56798c512e0fc9c7c5340641f9f0646",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/243bff734b1e13aed78847ba32a778b924c9dc8d"
        },
        "date": 1746302769474,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6895283.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500400\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606194843,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 74228607,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22173310.5,
            "unit": "ns",
            "extra": "gctime=5559624\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42952582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7808010,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1413191639.5,
            "unit": "ns",
            "extra": "gctime=352918766.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1666048,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17820333.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151894239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77332384,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e2b34f1d9f464c5f592a3cd4b237cdf7511109b2",
          "message": "Introduce `QuantumToolbox.settings` and `auto_tidyup` (#460)",
          "timestamp": "2025-05-05T11:28:22+02:00",
          "tree_id": "cd789c2181ff73352e33eade66b73d1f4bb43cfe",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e2b34f1d9f464c5f592a3cd4b237cdf7511109b2"
        },
        "date": 1746437826418,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7062323.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 608753403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80322523,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 26675578,
            "unit": "ns",
            "extra": "gctime=8331494\nmemory=121708264\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42992995.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7689375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1454507820.5,
            "unit": "ns",
            "extra": "gctime=365552201\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1644058,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17711042,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 153362492,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77552587,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e4ff0f7dd88f7aa911f0796383a86d326c85ceea",
          "message": "CompatHelper: bump compat for SciMLOperators to 0.4, (keep existing compat) (#466)",
          "timestamp": "2025-05-10T18:24:34+09:00",
          "tree_id": "d2eccc22b96474e8c0002f172a7e867b56f9b61b",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e4ff0f7dd88f7aa911f0796383a86d326c85ceea"
        },
        "date": 1746869515086,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6872203,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595116370,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80013707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22307872,
            "unit": "ns",
            "extra": "gctime=5719316\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42741946.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7640403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1417764649.5,
            "unit": "ns",
            "extra": "gctime=356984834.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1645389,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17675482.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150983467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77242801,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a1d1627acbba053653808e0f5002b3191e80a25a",
          "message": "Bump to `v0.31.1` (#469)",
          "timestamp": "2025-05-16T19:45:35+08:00",
          "tree_id": "fec890dd81b7e3f81e67950ecf21043eeb2ccf28",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/a1d1627acbba053653808e0f5002b3191e80a25a"
        },
        "date": 1747396466436,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6778522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600232583,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 71153618.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3661352\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22916900.5,
            "unit": "ns",
            "extra": "gctime=5748059\nmemory=121708248\nallocs=104757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42496074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7626058,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1371432912.5,
            "unit": "ns",
            "extra": "gctime=330783192.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1636113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17648899,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152720004,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76939388,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "9fff0169a3c44fea82219f6ec62c026075b97ac8",
          "message": "Simplify `QobEvo` creation and time evolution (#477)",
          "timestamp": "2025-05-29T08:56:17+08:00",
          "tree_id": "80f58a80872419e2da29f29d0ae2fe32af899b03",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9fff0169a3c44fea82219f6ec62c026075b97ac8"
        },
        "date": 1748480696488,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6868924,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500400\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596548138,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80243306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22258319.5,
            "unit": "ns",
            "extra": "gctime=5586394\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42790927,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7667507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1387180556.5,
            "unit": "ns",
            "extra": "gctime=335057831.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1658148,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17815593,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151420590.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76749424,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a965621b17b90a10920afa88deccb932391d12e6",
          "message": "Switch `Core` tests to \"Test Item Framework\" (#475)",
          "timestamp": "2025-06-02T11:24:29+02:00",
          "tree_id": "eb6849b187ff42f3c1dd21fef69b7bc8dd30a456",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/a965621b17b90a10920afa88deccb932391d12e6"
        },
        "date": 1748856693123,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6853913.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601562848,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79916339.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22652315,
            "unit": "ns",
            "extra": "gctime=5633364\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42610148,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7808914.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1390611455,
            "unit": "ns",
            "extra": "gctime=333745558\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1652690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17762629,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151988095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77849753.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "secli.matteo@gmail.com",
            "name": "Matteo Seclì",
            "username": "matteosecli"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6b6f674d073b5fa581195b5381d1bbdf8b08e5ab",
          "message": "Lanczos-based `spectrum` method (#476)",
          "timestamp": "2025-06-04T10:32:36+02:00",
          "tree_id": "367b693cdff47dbe4586e881a369fcdfdd06661d",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/6b6f674d073b5fa581195b5381d1bbdf8b08e5ab"
        },
        "date": 1749026491389,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6913556,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594960902,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80775532,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 430329.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22703277,
            "unit": "ns",
            "extra": "gctime=5604165\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42581157,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7634778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1397980205,
            "unit": "ns",
            "extra": "gctime=343301478.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1648327,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17811679,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150520539.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77481250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "157601901+TendonFFF@users.noreply.github.com",
            "name": "Li-Xun Cai",
            "username": "TendonFFF"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fbe80da8ead53d3f5a46c45c6eccc0710d22318c",
          "message": "Bloch Redfield master equation implementation (#473)",
          "timestamp": "2025-06-04T14:14:34+02:00",
          "tree_id": "afd1c6d44eb83fd71da9cf6faf6238cb1812bc82",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/fbe80da8ead53d3f5a46c45c6eccc0710d22318c"
        },
        "date": 1749039528339,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6984460,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 608276072,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67882682,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 440413,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23277525,
            "unit": "ns",
            "extra": "gctime=6130137\nmemory=121724280\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43159185,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7802219,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1490310820,
            "unit": "ns",
            "extra": "gctime=376940901\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1657107,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17926390,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 153928665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78559123.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "feroz.ahmed.mian@gmail.com",
            "name": "Feroz Ahmed Mian",
            "username": "Fe-r-oz"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "47deb2c2d702c3b3ca83ad3b9f6d7add5b24bffb",
          "message": "Implement Bloch Sphere rendering (#472)\n\nCo-authored-by: Alberto Mercurio <61953577+albertomercurio@users.noreply.github.com>\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>",
          "timestamp": "2025-06-05T09:27:38+08:00",
          "tree_id": "960a3d7e493c4aa084287bcd734f8415ae5dd154",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/47deb2c2d702c3b3ca83ad3b9f6d7add5b24bffb"
        },
        "date": 1749087332701,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7059240,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603014790,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 68690993.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 449433.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 25912523,
            "unit": "ns",
            "extra": "gctime=7429308\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43027081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7808881,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1421569237.5,
            "unit": "ns",
            "extra": "gctime=345936428.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1653807,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17924478,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151747220,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77911188.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "64ec15dea9f27fb8bae7c1558b6edcbda6c474bd",
          "message": "Improve `Bloch` sphere code structure and docs (#480)\n\nCo-authored-by: Alberto Mercurio <61953577+albertomercurio@users.noreply.github.com>",
          "timestamp": "2025-06-06T01:26:23+08:00",
          "tree_id": "9b1395fb92bf2be98a114dc6082279c0e2981d3e",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/64ec15dea9f27fb8bae7c1558b6edcbda6c474bd"
        },
        "date": 1749144650500,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7061877.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595628426,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79558867,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 436483,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 24644022.5,
            "unit": "ns",
            "extra": "gctime=6563070\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42843341.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7809630.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1424080638.5,
            "unit": "ns",
            "extra": "gctime=338567665.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1662786.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17728548,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152012864,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77377653.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "826a03d71a24df641d5c6a2b3f16e44b09adf36d",
          "message": "Implement `Base.copy` for `AbstractQuantumObject` (#486)",
          "timestamp": "2025-06-09T07:54:56+08:00",
          "tree_id": "2e96e53621d43669bb5249f2895446e25a7a8199",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/826a03d71a24df641d5c6a2b3f16e44b09adf36d"
        },
        "date": 1749427456471,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7165173,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603486046,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 75440656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 456171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 25706893,
            "unit": "ns",
            "extra": "gctime=7449855\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43044316,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7664895,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1449461818.5,
            "unit": "ns",
            "extra": "gctime=366869439.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1649896.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17721757.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151303535.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76683599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "22efb110b3cd3ea0a04e95b1966285cb80dfcef0",
          "message": "Use `LScene` instead of `Axis3` for `Bloch` Sphere (#485)\n\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>",
          "timestamp": "2025-06-09T11:37:02+08:00",
          "tree_id": "3f18591435d657f8b579936344f46fc647d7c746",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/22efb110b3cd3ea0a04e95b1966285cb80dfcef0"
        },
        "date": 1749440469300,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6758587,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594993152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79921384.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 429327,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22475648.5,
            "unit": "ns",
            "extra": "gctime=5609938\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42583042,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7603323,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1369314238,
            "unit": "ns",
            "extra": "gctime=345116656\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1638077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17629343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150357584.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 76558768,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "dc4ecde44c69ea77969d832839332f4c555de769",
          "message": "Fix `Bloch` docstrings and improve visualization of the sphere (#487)",
          "timestamp": "2025-06-09T22:40:10+02:00",
          "tree_id": "e3419fd7c53582b83fac33bdd3774474e6b920d6",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/dc4ecde44c69ea77969d832839332f4c555de769"
        },
        "date": 1749502048820,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6973065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603590507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 70325205,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 440196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23283950.5,
            "unit": "ns",
            "extra": "gctime=5756835.5\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42607464,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7681159.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1412827960.5,
            "unit": "ns",
            "extra": "gctime=361647082\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1654496,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17813751,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 154092477,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77891604,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a4620a940ffcca1471f9922912fa5cf3df66a88a",
          "message": "Improve `Bloch` alignment with qutip (#489)",
          "timestamp": "2025-06-12T00:04:10+02:00",
          "tree_id": "65e3566993e87ea7d36fddc595fbce9a9ed55852",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/a4620a940ffcca1471f9922912fa5cf3df66a88a"
        },
        "date": 1749679886057,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6849471,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595537136,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80527658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 431218.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22876742,
            "unit": "ns",
            "extra": "gctime=5566547\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42658202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7629914,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1365650303,
            "unit": "ns",
            "extra": "gctime=328446134\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1649111.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17765945,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152033775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77293683,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b5e8af1cca91557df525d5a267c43f6580748b33",
          "message": "CompatHelper: bump compat for Makie in [weakdeps] to 0.23, (keep existing compat) (#490)\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>",
          "timestamp": "2025-06-13T14:41:26+08:00",
          "tree_id": "55782d6a664c0b2dc543df22cd178d19cb245b2e",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/b5e8af1cca91557df525d5a267c43f6580748b33"
        },
        "date": 1749797287546,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6881328,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595539914,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78466693,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 428463.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 21809996,
            "unit": "ns",
            "extra": "gctime=5369216.5\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42445056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7779597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1385520609,
            "unit": "ns",
            "extra": "gctime=345216021.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1653962,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17660630,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152810597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77706787.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "157601901+TendonFFF@users.noreply.github.com",
            "name": "Li-Xun Cai",
            "username": "TendonFFF"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "29a3785ebe236860db9bff3a3f6d2eed48a1fab6",
          "message": "fix typos in doc strings (#495)\n\nCo-authored-by: cailixun <cailixun23@gmail.com>",
          "timestamp": "2025-06-20T09:36:52+08:00",
          "tree_id": "f4d12fcbbbf2493eef8bf5b93f2b9053f46e19cc",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/29a3785ebe236860db9bff3a3f6d2eed48a1fab6"
        },
        "date": 1750383770883,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6989752,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603443212,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 70313522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 438627,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22613867,
            "unit": "ns",
            "extra": "gctime=5934230\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42622769,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7632787,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1387816926.5,
            "unit": "ns",
            "extra": "gctime=337702659\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1651793,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17922245.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151952173,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77577128,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ee92a07610f62622a847200bf17f2015a96394f9",
          "message": "CompatHelper: bump compat for Makie in [weakdeps] to 0.24, (keep existing compat) (#496)",
          "timestamp": "2025-06-23T09:52:16+08:00",
          "tree_id": "244bc5770574096a3ec816132c96faffa477ac94",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/ee92a07610f62622a847200bf17f2015a96394f9"
        },
        "date": 1750643981457,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6876419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600332475,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67273383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 430972.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22306699.5,
            "unit": "ns",
            "extra": "gctime=5619230.5\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42551077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7781834.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1398583141,
            "unit": "ns",
            "extra": "gctime=339729069.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1672621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17764820,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151856251,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77818432.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3028fe0a4af5b9ec3952fbefa737e6c0e0c369b3",
          "message": "Bump version to v0.32.1 (#498)",
          "timestamp": "2025-06-24T10:13:11+08:00",
          "tree_id": "df044c7ca959c77ebcfbe287dfe86e3ea6210291",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/3028fe0a4af5b9ec3952fbefa737e6c0e0c369b3"
        },
        "date": 1750731585132,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6885834.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 592356491,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80196972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 426931,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 21542640,
            "unit": "ns",
            "extra": "gctime=5453915\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42313903,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7784729,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1353919591.5,
            "unit": "ns",
            "extra": "gctime=324283283.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1668880,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17759134,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152245135,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77500040,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "secli.matteo@gmail.com",
            "name": "Matteo Seclì",
            "username": "matteosecli"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0066cdbb468119a25cd79dc051013234ce2bd9d7",
          "message": "Check for orthogonality breakdown in `Lanczos` (#502)",
          "timestamp": "2025-07-05T11:44:09+02:00",
          "tree_id": "620c0085ddcb268ea5234d1a42006091752ecace",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/0066cdbb468119a25cd79dc051013234ce2bd9d7"
        },
        "date": 1751709197499,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7013266,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500384\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597147832,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 70104481,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 454262,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23626736.5,
            "unit": "ns",
            "extra": "gctime=6316071.5\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42876613.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7645966,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1377481936,
            "unit": "ns",
            "extra": "gctime=333852240.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1639885,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17655564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150855324.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77357243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "157601901+TendonFFF@users.noreply.github.com",
            "name": "Li-Xun Cai",
            "username": "TendonFFF"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a2ade17784e425bc72ed81eb3d6655a993baefc5",
          "message": "Excitation number restricted (`ENR`) state space implementation  (#500)",
          "timestamp": "2025-07-07T11:11:25+02:00",
          "tree_id": "0d934cf725bc293f2508f2540e6269413a2f5de7",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/a2ade17784e425bc72ed81eb3d6655a993baefc5"
        },
        "date": 1751879737686,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6808085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594755537,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79856153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 439250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22911909,
            "unit": "ns",
            "extra": "gctime=5679996\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42575971,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7639119,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1462501614.5,
            "unit": "ns",
            "extra": "gctime=356648630.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1653365,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17841029,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151100922,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328888\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77003487,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330184\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3057e3e002096e30b356fe88228c2edb58d888f8",
          "message": "Fix errors in `Julia v1.12` (#507)",
          "timestamp": "2025-07-16T16:49:40-04:00",
          "tree_id": "9e6247c34a0859c9a2ebb8af2419ec1ae318620d",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/3057e3e002096e30b356fe88228c2edb58d888f8"
        },
        "date": 1752699517483,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7074041,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500400\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 599113589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67308783,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 458587,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23437725,
            "unit": "ns",
            "extra": "gctime=5989940\nmemory=121724312\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42898359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7786662,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1376512371.5,
            "unit": "ns",
            "extra": "gctime=345830872.5\nmemory=2719211208\nallocs=232097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1664888,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17759907,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152370500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328904\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77628031,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330200\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5e1707df94a157d1bc80175e58303102991dae83",
          "message": "Store both `times` and `times_states` in time evolution solutions (#506)",
          "timestamp": "2025-07-22T00:15:24-04:00",
          "tree_id": "adab95ee8d5bde462d5d1a5b0f0ec63eb3e5f00a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5e1707df94a157d1bc80175e58303102991dae83"
        },
        "date": 1753158167732,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6890848,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598116283,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79528935.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 449852,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23697088.5,
            "unit": "ns",
            "extra": "gctime=6255950\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42868031,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7682248,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1425059016,
            "unit": "ns",
            "extra": "gctime=345920844\nmemory=2719211480\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1648634,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17790917,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152049442,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328904\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77263616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330200\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "bc76e9f90b190fe6bd1a45f403ccbbe45539abf9",
          "message": "Bump to v0.33.0 (#508)",
          "timestamp": "2025-07-22T00:24:53-04:00",
          "tree_id": "11ea14eed43d889a3bca81cf5b972ebe43eeee37",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/bc76e9f90b190fe6bd1a45f403ccbbe45539abf9"
        },
        "date": 1753158729505,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7117232,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601110659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67472336,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 466983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 27974105,
            "unit": "ns",
            "extra": "gctime=8126348\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42956153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7686585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1401600155,
            "unit": "ns",
            "extra": "gctime=377156686.5\nmemory=2719211480\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1648724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17938519.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151256766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328904\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78139984.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330200\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "7a6cc5034c47956f95668fecb7430cc920bf5092",
          "message": "Improve efficiency of generation of Bloch-Redfield tensor and fix documentation typos (#509)",
          "timestamp": "2025-07-25T21:44:30-04:00",
          "tree_id": "06622715445df77e455d8ea94a483b2c00e9422a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/7a6cc5034c47956f95668fecb7430cc920bf5092"
        },
        "date": 1753494808700,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7015536,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602852903,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80101296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 440795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22888090.5,
            "unit": "ns",
            "extra": "gctime=5675337.5\nmemory=121724360\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42673093,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7655576,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1404749119,
            "unit": "ns",
            "extra": "gctime=347600146\nmemory=2719211480\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1641878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17631400,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151098719.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328904\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77322441,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330200\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e2c6d3e99435f41b8f38286f45aff2724fd20d36",
          "message": "CompatHelper: bump compat for SciMLOperators to 1, (keep existing compat) (#470)\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>\nCo-authored-by: Yi-Te Huang <y.t.d.huang@phys.ncku.edu.tw>\nCo-authored-by: Alberto Mercurio <alberto.mercurio96@gmail.com>",
          "timestamp": "2025-07-25T21:54:00-04:00",
          "tree_id": "99bf3d2854821da4b1b464ede0268216099da3cb",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e2c6d3e99435f41b8f38286f45aff2724fd20d36"
        },
        "date": 1753495741987,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7192305.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595101428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80281612,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 444793,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22718397,
            "unit": "ns",
            "extra": "gctime=5710036\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42767801.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7804909.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1384206522.5,
            "unit": "ns",
            "extra": "gctime=351002544.5\nmemory=2719211480\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1671581,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17791811,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151181480,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328904\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77515017,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330200\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c6283feebeacbb8b0446f16e9fa39dc1d8a67fff",
          "message": "CompatHelper: bump compat for Makie in [weakdeps] to 0.24 (#513)\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>\nCo-authored-by: Alberto Mercurio <alberto.mercurio96@gmail.com>",
          "timestamp": "2025-07-27T17:50:21+08:00",
          "tree_id": "5e0880ef15125d0ad43eb7537f99e2ddb0d912ea",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/c6283feebeacbb8b0446f16e9fa39dc1d8a67fff"
        },
        "date": 1753610220098,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6812574,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500400\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600769342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 73887242,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 436263.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22126793,
            "unit": "ns",
            "extra": "gctime=5541569\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42625074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7759088,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1372253871.5,
            "unit": "ns",
            "extra": "gctime=330765850\nmemory=2719211480\nallocs=232098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1659655.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17793156,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152170403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328904\nallocs=11787\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77358554,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3330200\nallocs=11799\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0bc7621facad7d691298d7ee056c417567845677",
          "message": "Add `keep_runs_results` option for multi-trajectory solvers to align with qutip (#512)\n\nCo-authored-by: Alberto Mercurio <alberto.mercurio96@gmail.com>",
          "timestamp": "2025-07-27T09:46:14-04:00",
          "tree_id": "0a356b99bd0619b86951f42c01756593353394c8",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/0bc7621facad7d691298d7ee056c417567845677"
        },
        "date": 1753624408557,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6983030,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600860732,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78823863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 443647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22743856,
            "unit": "ns",
            "extra": "gctime=5631274.5\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42684886.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7640854,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1424476605,
            "unit": "ns",
            "extra": "gctime=349721920\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1648167.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17836567,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152520353,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78251810.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "5ba859ac5231a58a44193dba6475c2821f8ae9c5",
          "message": "Add support to ForwarDiff differentiation (#515)",
          "timestamp": "2025-07-27T11:31:23-04:00",
          "tree_id": "795e406f6a3e7d10831759b9ca46200525592e0c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5ba859ac5231a58a44193dba6475c2821f8ae9c5"
        },
        "date": 1753630531113,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6908706,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601042168,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79680383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 450474.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 24020537,
            "unit": "ns",
            "extra": "gctime=6231227\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42383275.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7616237,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1374378570.5,
            "unit": "ns",
            "extra": "gctime=366320735\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1648593,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17642327,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151581625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78333391,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8fb7406e55faba207c233ba72816ff6b1731b2a3",
          "message": "Bump compat of `SciMLBase` to fix piracy error (#516)",
          "timestamp": "2025-07-27T13:19:24-04:00",
          "tree_id": "d8af150a2224289710c9600398f6077a09fd4839",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/8fb7406e55faba207c233ba72816ff6b1731b2a3"
        },
        "date": 1753637198378,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7021644.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600543941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79696161,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 445182,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22732302,
            "unit": "ns",
            "extra": "gctime=5872502\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42537667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7630211.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1388137268,
            "unit": "ns",
            "extra": "gctime=341357983\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1662708,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17687582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152249737,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78914155,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "13709a17a1269f45dba6742843d75bdd76c8653f",
          "message": "Bump to v0.34.0 (#518)",
          "timestamp": "2025-07-29T16:21:00+02:00",
          "tree_id": "5b074898a55f461549bd0c33eba77bd36ecae351",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/13709a17a1269f45dba6742843d75bdd76c8653f"
        },
        "date": 1753799312007,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7102514,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598806932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80322211,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 456358.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 25213553,
            "unit": "ns",
            "extra": "gctime=7143272\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42742828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7784693.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1414536600.5,
            "unit": "ns",
            "extra": "gctime=347829316\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1658583.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17790660,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 154477908,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78293614,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "157601901+TendonFFF@users.noreply.github.com",
            "name": "Li-Xun Cai",
            "username": "TendonFFF"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5110d93ce0389d116d30c88984c6267d03112023",
          "message": "Improve Bloch sphere rendering for animation (#520)",
          "timestamp": "2025-08-08T11:42:09+02:00",
          "tree_id": "ca5a10872a53d7d586e12dbeb9f2b0dd29dc2512",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5110d93ce0389d116d30c88984c6267d03112023"
        },
        "date": 1754646673430,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6851863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597170657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79946696.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 443606,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22739848,
            "unit": "ns",
            "extra": "gctime=5748095\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42502510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7643095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1417974408.5,
            "unit": "ns",
            "extra": "gctime=350312254\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1661534,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17836805.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 153385725,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78358403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "0b9a9f042e28922f244e17af3f67d34e8009dceb",
          "message": "Fix Zygote autodiff errors (#530)",
          "timestamp": "2025-08-20T15:13:31+02:00",
          "tree_id": "1c0e44468aaf3e1c641f592c8fa84e9ac5fdb8ca",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/0b9a9f042e28922f244e17af3f67d34e8009dceb"
        },
        "date": 1755696113307,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6858707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600279465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80196303,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662776\nallocs=446\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 434117,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23704217,
            "unit": "ns",
            "extra": "gctime=5969087.5\nmemory=121884376\nallocs=112761\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42797764.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7643693.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1389980030.5,
            "unit": "ns",
            "extra": "gctime=339279967.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1641234,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17645972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152279432,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78005163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fab15f71c2469f51c254f07f8b4c10c8c58a45f0",
          "message": "Bump actions/checkout from 4 to 5 (#526)\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2025-08-20T16:16:03+02:00",
          "tree_id": "735338aca81ad5070ada92925acef8ac53d9d23c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/fab15f71c2469f51c254f07f8b4c10c8c58a45f0"
        },
        "date": 1755699621810,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6858459.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597232272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67447773,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662776\nallocs=446\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 448286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 24077701,
            "unit": "ns",
            "extra": "gctime=6226065\nmemory=121884376\nallocs=112761\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42733947,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7643763.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1448369828.5,
            "unit": "ns",
            "extra": "gctime=383330639\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1645228.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17866063,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152489526,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78562976,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "d95bf63db2628fdbbb37597ad605f24c5d40ba31",
          "message": "Add support for Enzyme autodiff (#531)",
          "timestamp": "2025-08-21T10:31:40+02:00",
          "tree_id": "2ff1d9e378a99dddaa936ce32f39547e3f5d5cc0",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d95bf63db2628fdbbb37597ad605f24c5d40ba31"
        },
        "date": 1755765471155,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6915252.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597438588,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80401803,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662776\nallocs=446\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 440697,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22434370,
            "unit": "ns",
            "extra": "gctime=5557368\nmemory=121884376\nallocs=112761\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42839983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7613265.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1382888173,
            "unit": "ns",
            "extra": "gctime=336372780.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1642175,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17655362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 150856681.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78552300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "6edcca11d39fd8875cef6b3689199e0010303174",
          "message": "Bump version to v0.34.1 (#534)",
          "timestamp": "2025-08-23T15:22:22+02:00",
          "tree_id": "cfd981d2a0814b1b694115d74f7e556af4c4fe5c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/6edcca11d39fd8875cef6b3689199e0010303174"
        },
        "date": 1755955781138,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7120265,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606321542,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79542794,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 449624.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23703710.5,
            "unit": "ns",
            "extra": "gctime=5822801.5\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42896568,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7611741,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1423098726,
            "unit": "ns",
            "extra": "gctime=386323391\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1638541,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17727247,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152044120,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78512766.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "157601901+TendonFFF@users.noreply.github.com",
            "name": "Li-Xun Cai",
            "username": "TendonFFF"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "15954c9bf89e87b4fe94f785aeef1247337a9f00",
          "message": "Add `QobjEvo` support for `SteadyStateODESolver` (#536)\n\nCo-authored-by: cailixun <cailixun23@gmail.com>",
          "timestamp": "2025-08-29T13:40:22+08:00",
          "tree_id": "9e94a4c39241d891df5abbdcf3e043f8bc6a765e",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/15954c9bf89e87b4fe94f785aeef1247337a9f00"
        },
        "date": 1756446941199,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7158527,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602960473,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67934097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 463097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 27309116,
            "unit": "ns",
            "extra": "gctime=7786438.5\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43075651.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7665757,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1366450659,
            "unit": "ns",
            "extra": "gctime=370598232.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1640180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17676173.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152236920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77927306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "75f33c468cabba15f020e42777c1471776071312",
          "message": "Fix keyword argument handling for `SteadyStateODESolver` (#537)",
          "timestamp": "2025-09-01T09:04:01+02:00",
          "tree_id": "1df2f97ebdc7debfbd3e69e563d1f1adb07a8f98",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/75f33c468cabba15f020e42777c1471776071312"
        },
        "date": 1756710775887,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6873374.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596930020,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80352062,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 440519,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22653257.5,
            "unit": "ns",
            "extra": "gctime=5839274.5\nmemory=121724328\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42669995,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7595605,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1357242264,
            "unit": "ns",
            "extra": "gctime=320927337\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1645119,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17640572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152846349,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 77750650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "326c946b2edc9377350ff29ce36a3c7a31316523",
          "message": "Fix incorrect `negativity` and `partial_transpose` for arbitrary subsystem dimension (#539)",
          "timestamp": "2025-09-03T12:38:34+02:00",
          "tree_id": "d472815ec27efc99836b2d8ea2bd2eba4fce5e97",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/326c946b2edc9377350ff29ce36a3c7a31316523"
        },
        "date": 1756896448155,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6961213,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606590244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80154899,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 443511,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22881533,
            "unit": "ns",
            "extra": "gctime=5677748\nmemory=121724280\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42628321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7660851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1372564014.5,
            "unit": "ns",
            "extra": "gctime=333365985.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1659500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17830457.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152524500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78448650.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "be147e59740c7f4bea34b1a66c1aa8eb9f8035a1",
          "message": "Bump to v0.35.0 (#540)",
          "timestamp": "2025-09-03T12:46:01+02:00",
          "tree_id": "238a3dba38c2eac6db72b88de3d36d273d79fcff",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/be147e59740c7f4bea34b1a66c1aa8eb9f8035a1"
        },
        "date": 1756896910964,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7075522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 605040118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 75268352.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662392\nallocs=434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 451982,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 23651215,
            "unit": "ns",
            "extra": "gctime=5934065\nmemory=121724248\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42946042.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7651971.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1421182887.5,
            "unit": "ns",
            "extra": "gctime=354532612.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1647981,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 18015586,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 151788545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78034565,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "252e3621f4f3459dee5a5b00dbef70f36c5ebb14",
          "message": "Add `sortby` and `rev` arguments to eigensolver functions (#546)",
          "timestamp": "2025-09-27T13:07:22+02:00",
          "tree_id": "5ac91689706e1346999162d8491988dcb26fc6a5",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/252e3621f4f3459dee5a5b00dbef70f36c5ebb14"
        },
        "date": 1758971770512,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7077189,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 599362725,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 69671761,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662552\nallocs=442\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 455836,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 24129259,
            "unit": "ns",
            "extra": "gctime=6545799\nmemory=121740264\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42971273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7780763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1448445023.5,
            "unit": "ns",
            "extra": "gctime=374331423.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1659597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17942402.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 153157161,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78229718.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3f6d8ed8dc8087d501817bbff97091d8de0956de",
          "message": "Fix `Makie` warning message (#528)",
          "timestamp": "2025-09-27T13:55:52+02:00",
          "tree_id": "23653d15888afc188c0c3c978ab7382cd7cd1567",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/3f6d8ed8dc8087d501817bbff97091d8de0956de"
        },
        "date": 1758974401978,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7021310,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597500709,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79826103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662552\nallocs=442\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 461106,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 26078266.5,
            "unit": "ns",
            "extra": "gctime=7904825\nmemory=121740312\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42826086.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7655193,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1436515764,
            "unit": "ns",
            "extra": "gctime=388359711.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1645116,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17826679.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152818687,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78483185,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "64a727c4b1d8756746410268f70ce9f38cc20bf9",
          "message": "Improve `QuantumObjectEvolution` generation function and fix `smesolve` type instabilities (#548)",
          "timestamp": "2025-09-28T10:01:13+08:00",
          "tree_id": "07df2f175a97d9af374883e75d923d3f32c496b8",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/64a727c4b1d8756746410268f70ce9f38cc20bf9"
        },
        "date": 1759025128224,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6922586.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595988932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80302945,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662552\nallocs=442\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 505099,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 24818174,
            "unit": "ns",
            "extra": "gctime=7102750\nmemory=121740264\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43063748.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7823153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1441165678,
            "unit": "ns",
            "extra": "gctime=381312395\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1680182,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17984432,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152395061,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78655586.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "44385685+ytdHuang@users.noreply.github.com",
            "name": "Yi-Te Huang",
            "username": "ytdHuang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b2c90b509c951b7550090f3f344588a87e91c906",
          "message": "Run code-quality CI also on `Julia 1` (#401)\n\nCo-authored-by: Alberto Mercurio <alberto.mercurio96@gmail.com>",
          "timestamp": "2025-09-28T09:37:00+02:00",
          "tree_id": "6d23d1961ae395054a15fcd1bea114fcf4c314ae",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/b2c90b509c951b7550090f3f344588a87e91c906"
        },
        "date": 1759045270073,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6776129,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10500368\nallocs=484\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598423923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091528\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80639996,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3662552\nallocs=442\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Lanczos",
            "value": 441736,
            "unit": "ns",
            "extra": "gctime=0\nmemory=955600\nallocs=580\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 22676592,
            "unit": "ns",
            "extra": "gctime=5980251.5\nmemory=121740264\nallocs=105757\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42686142.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813312\nallocs=1036\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 7591145,
            "unit": "ns",
            "extra": "gctime=0\nmemory=852160\nallocs=714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1396343008.5,
            "unit": "ns",
            "extra": "gctime=338374753.5\nmemory=2719211352\nallocs=232096\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 1639145,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75464\nallocs=92\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 17732375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3245504\nallocs=613\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 152343352,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13677208\nallocs=14195\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 78514083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13678504\nallocs=14207\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}