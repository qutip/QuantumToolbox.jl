window.BENCHMARK_DATA = {
  "lastUpdate": 1739524351227,
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
          "id": "5af51a267ec83daaaec6c0a81cffe5768916c369",
          "message": "Fix type conversion of `tlist` in time evolution (#229)",
          "timestamp": "2024-09-20T18:23:21+08:00",
          "tree_id": "28847ca5c5ff4c191566e21c2b25d77b48f52b73",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5af51a267ec83daaaec6c0a81cffe5768916c369"
        },
        "date": 1726828276481,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7176103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601754380,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 81553370.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13849581,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447624\nallocs=3384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43112583,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1377301111,
            "unit": "ns",
            "extra": "gctime=31125832\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4724203,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26693507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700336\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 277389699,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8328296\nallocs=101391\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 141232300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8222176\nallocs=99897\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "7605b69e2093386eea0f634a20cb617093233b7f",
          "message": "fix `rand_unitary` (#230)",
          "timestamp": "2024-09-22T18:22:42+02:00",
          "tree_id": "dbaf50a388ab56e35e7eec6521650f68ff272f61",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/7605b69e2093386eea0f634a20cb617093233b7f"
        },
        "date": 1727022550585,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7158065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601307851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80352519,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13931871,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447624\nallocs=3384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43495323.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1375342645,
            "unit": "ns",
            "extra": "gctime=33079065\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4747916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26673694.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700336\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 275531905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8307384\nallocs=101296\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 141226245,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8211240\nallocs=99971\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "75fefe8b85eb6ea12d319ed6f7438e13fc474269",
          "message": "Automatically convert TimeDependentOperatorSum to `liouvillian` (#235)",
          "timestamp": "2024-09-27T13:22:26+08:00",
          "tree_id": "40fb0be2257b6d5e73e5ef161b053f40ca5fb814",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/75fefe8b85eb6ea12d319ed6f7438e13fc474269"
        },
        "date": 1727414944326,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7088910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606386045,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 76221774.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13931158,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447624\nallocs=3384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43265774,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1387381348,
            "unit": "ns",
            "extra": "gctime=34047647\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4702379.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26788826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700336\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 278845763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8352808\nallocs=102014\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140932010,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8274864\nallocs=100708\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "fd098b3e0e9a439d245fa118a76eedcf1760ae4b",
          "message": "Change version to v0.14.1",
          "timestamp": "2024-09-27T10:44:12+02:00",
          "tree_id": "d29c89dd4ee66009d2241fc69c5750242ec1ad70",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/fd098b3e0e9a439d245fa118a76eedcf1760ae4b"
        },
        "date": 1727426871045,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7167371,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600473767.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79789233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13899998,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447624\nallocs=3384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43092805,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1362714675,
            "unit": "ns",
            "extra": "gctime=32456617\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4715703,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26612903,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700336\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 273721997,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8310000\nallocs=101409\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140881587,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8192560\nallocs=99864\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "9172e829f1f35f5007373e20b470a90ca7282849",
          "message": "Implement Stochastic Schrodinger equation (#238)\n\n* Implement Stochastic Schrodinger equation\r\n\r\n* Minor changes\r\n\r\n* Minor changes on the import\r\n\r\n* Fix wrong link\r\n\r\n* Fix instabilities\r\n\r\n* Avoid safety_copy",
          "timestamp": "2024-09-27T19:18:34+02:00",
          "tree_id": "bc0c69573ecfabb7958f24927ecff4a064f5aeca",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9172e829f1f35f5007373e20b470a90ca7282849"
        },
        "date": 1727457936650,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7265814,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 604620453.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 69947797.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13966852.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447624\nallocs=3384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 45184742,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1451754190,
            "unit": "ns",
            "extra": "gctime=37443575\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4734089,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 28009972.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700336\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 277218439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8263960\nallocs=100689\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 141431380,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8229192\nallocs=100285\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "d6342342964fe58e1014abf693ce9c8be2473168",
          "message": "Change to version v0.15.0",
          "timestamp": "2024-09-27T21:11:41+02:00",
          "tree_id": "e162ffc0bafc9053eb05033ac2c6522f5f502a60",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d6342342964fe58e1014abf693ce9c8be2473168"
        },
        "date": 1727464517367,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7107863.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601655219.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80036781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13906958,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447624\nallocs=3384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43223566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1399290150,
            "unit": "ns",
            "extra": "gctime=35615831\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4701226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26877282,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700336\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 274390838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8299928\nallocs=101116\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 143598915,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8223408\nallocs=100277\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "c4df4976b2b018c6b780fbda23441081fcda0a7d",
          "message": "CompatHelper: bump compat for DiffEqCallbacks to 4, (keep existing compat) (#239)\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>",
          "timestamp": "2024-09-28T10:51:56+02:00",
          "tree_id": "21b5f2d74d891dccf6f608b25653f84a3346a277",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/c4df4976b2b018c6b780fbda23441081fcda0a7d"
        },
        "date": 1727513908002,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7085924,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 604565737.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67861987.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13892069,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447624\nallocs=3384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43346530.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1397567328,
            "unit": "ns",
            "extra": "gctime=33962593\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4690859,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26760784.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700336\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 274784839,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8061432\nallocs=98145\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 141636644,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8217120\nallocs=100086\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "bdac5f4ac254830eaf9b902d3273b87dd1d7dc0f",
          "message": "[Docs] Update docstring and documentation for `steadystate` (#240)\n\n* add missing docstrings for `steadystate` solvers\r\n\r\n* update documentation for `steadystate`\r\n\r\n* format files\r\n\r\n* improve `steadystate` documentation",
          "timestamp": "2024-09-28T11:58:20+02:00",
          "tree_id": "29a8b00ee3e12401a059766441ae042ccecdaf82",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/bdac5f4ac254830eaf9b902d3273b87dd1d7dc0f"
        },
        "date": 1727517881891,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7169319,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 599655314,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80572780,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13779589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447432\nallocs=3372\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43298460,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1382966960,
            "unit": "ns",
            "extra": "gctime=35184720\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4764281,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26473987,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700320\nallocs=770\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 275569005,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8391976\nallocs=102442\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 141674559,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8170688\nallocs=99608\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "1913fd9e6c6cc0e09dfa4e0be9757964b3d9f83d",
          "message": "Replace `n_traj` with `ntraj` for QuTiP compatibility (#243)",
          "timestamp": "2024-09-29T07:26:03+08:00",
          "tree_id": "98fa2982261c735492172a826de7d84ea51971b1",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/1913fd9e6c6cc0e09dfa4e0be9757964b3d9f83d"
        },
        "date": 1727566366448,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7054963,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597205950.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80040505.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13761233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447432\nallocs=3372\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43170125.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1404467421,
            "unit": "ns",
            "extra": "gctime=37470666\nmemory=3468571712\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4732486,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237376\nallocs=3098\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27285797,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700320\nallocs=770\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 273560594,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8170568\nallocs=99714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 142102016,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8249600\nallocs=100482\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "905e658a18d80a4ab61c5cc6993caeff5b6d4384",
          "message": "fix incorrect `times` in time evolution solutions (#244)",
          "timestamp": "2024-09-29T10:48:21+02:00",
          "tree_id": "ae049ffcb55370d8d6bbd28efdb4bc415d1e92a2",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/905e658a18d80a4ab61c5cc6993caeff5b6d4384"
        },
        "date": 1727600021758,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7090311,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 605557013.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 71892277,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13750448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447944\nallocs=3374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44711195.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1377935680,
            "unit": "ns",
            "extra": "gctime=27315010\nmemory=3468572864\nallocs=163544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4719473,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27579465.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 273461296.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8289880\nallocs=99853\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140346620,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8382952\nallocs=100898\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "88519acc33f93406a0200896a3111a8f20e5ce14",
          "message": "Bump version to v0.16.0",
          "timestamp": "2024-09-29T10:52:13+02:00",
          "tree_id": "7620fed66b935246d7b0ecb692f773954d38df50",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/88519acc33f93406a0200896a3111a8f20e5ce14"
        },
        "date": 1727600253386,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7097242,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 605160970.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79924971,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13694692,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447944\nallocs=3374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43211763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1368076509,
            "unit": "ns",
            "extra": "gctime=26679452\nmemory=3468572864\nallocs=163544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4737937,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26484884.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 273871143.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8428360\nallocs=101741\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 141464377,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8197368\nallocs=98957\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "d21041fe403693092e5977e53b3276c3f7eafbf4",
          "message": "Extend the support to `Tuple` in time_evolution (#248)\n\n* Extend the support to `Tuple` in time_evolution\r\n\r\n* Add missing cases",
          "timestamp": "2024-10-02T14:13:13+02:00",
          "tree_id": "b1f869e7629597a886ad526d72fc5ad2c72c2ef2",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d21041fe403693092e5977e53b3276c3f7eafbf4"
        },
        "date": 1727871592393,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7036705,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606935083.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80063438,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13691649,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447944\nallocs=3374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43253677,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1370896331,
            "unit": "ns",
            "extra": "gctime=30016186\nmemory=3468572864\nallocs=163544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4770605,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26515686,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 275063348,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8361888\nallocs=100869\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140845584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8259952\nallocs=99384\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "a3068b7a6af6caa70ddc6f9595612eb2f7d066dc",
          "message": "Replace `n_th` with `n_thermal` for QuTiP compatibility (#249)\n\n* make `n_th` to `n_thermal` and align with qutip\r\n\r\n* fix fontsize in docs\r\n\r\n* improve type stability",
          "timestamp": "2024-10-03T09:39:40+02:00",
          "tree_id": "e50310714f3d6f13c13b1f89e6ded7bad5a3fc26",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/a3068b7a6af6caa70ddc6f9595612eb2f7d066dc"
        },
        "date": 1727941404750,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7067665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602771693.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67618277,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13782633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447944\nallocs=3374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43147564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1357943015,
            "unit": "ns",
            "extra": "gctime=23805988\nmemory=3468572864\nallocs=163544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4728479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26529152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 271670957.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8355792\nallocs=100791\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 139958343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8335960\nallocs=100549\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "036c765b715c7e5d9152406f027727fb61e1bef1",
          "message": "Remove Luca from list of maintainers (#253)",
          "timestamp": "2024-10-04T02:10:56+08:00",
          "tree_id": "0bd0fc0aa2c0a993d4833ef34c945aa061b0bee7",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/036c765b715c7e5d9152406f027727fb61e1bef1"
        },
        "date": 1727979400658,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7124381.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600082139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67231721,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13803002.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447944\nallocs=3374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44558509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1386924601,
            "unit": "ns",
            "extra": "gctime=33370518\nmemory=3468572864\nallocs=163544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4682601.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27577536,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 273611281.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8408864\nallocs=101382\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140275844,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8189464\nallocs=98564\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "bde91428875a70c4ee5c68c1ae90f0891e3b4de7",
          "message": "Fix return type precision for `n_thermal` (#251)",
          "timestamp": "2024-10-04T16:39:05+08:00",
          "tree_id": "9c95280428cd714bb24a41676cf4e6c645452c4a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/bde91428875a70c4ee5c68c1ae90f0891e3b4de7"
        },
        "date": 1728031504488,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7433619,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 609858777,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79848851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13869070.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1447944\nallocs=3374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43401987,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1376933909,
            "unit": "ns",
            "extra": "gctime=25487659\nmemory=3468572864\nallocs=163544\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4727196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26600663.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 276245733.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8369560\nallocs=100793\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 142070013,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8241816\nallocs=98705\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "48f8c41840c623d791789bedd4e7272c0522139b",
          "message": "Add progress_bar in mcsolve, ssesolve and dsf_mcsolve (#254)",
          "timestamp": "2024-10-04T10:42:14+02:00",
          "tree_id": "c7c4f4338889f2acbd958bc003ba7e811862145a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/48f8c41840c623d791789bedd4e7272c0522139b"
        },
        "date": 1728031696581,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7310152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603554758,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79783516,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13874417,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1446136\nallocs=3346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43618412.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1413386972,
            "unit": "ns",
            "extra": "gctime=29335162\nmemory=3468572752\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4676098.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26683758,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 275084816.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8318792\nallocs=100086\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 190480199.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8352040\nallocs=100283\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "51684ab99544e57b1ff8404e5dffc84f2e8acc8b",
          "message": "Introduce `convert_unit` (#255)\n\n* introduce `convert_unit`\r\n\r\n* add tests for utility functions\r\n\r\n* fix value of physical constants",
          "timestamp": "2024-10-04T16:17:16+02:00",
          "tree_id": "033e64d3358d2c1dfbd4dcb889605892751e2a1f",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/51684ab99544e57b1ff8404e5dffc84f2e8acc8b"
        },
        "date": 1728051669598,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7257198,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600154158,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80188160,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13807640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1446136\nallocs=3346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43597351.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1388960936,
            "unit": "ns",
            "extra": "gctime=29533105\nmemory=3468572752\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4728797,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26600211,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 275397987.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8389088\nallocs=100995\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 188659868,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8386224\nallocs=100753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "c7a34832a65a7a99a62ffa588a61830b96db8b06",
          "message": "Bump version to v0.17.0",
          "timestamp": "2024-10-04T16:17:52+02:00",
          "tree_id": "72b862449d077339ced60044b8574fbeecc98429",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/c7a34832a65a7a99a62ffa588a61830b96db8b06"
        },
        "date": 1728051696526,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7115283,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597974875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 81461912.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13717018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1446136\nallocs=3346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43254920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1374313243,
            "unit": "ns",
            "extra": "gctime=30993387\nmemory=3468572752\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4723344,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26699229,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 275800186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8376720\nallocs=100323\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 190621351.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8345840\nallocs=100148\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "2f92e5bde0d2847c0e4dd75cedcceb18e44d6ad2",
          "message": "introduce `PhysicalConstants` and new energy units (#257)",
          "timestamp": "2024-10-05T09:39:07+02:00",
          "tree_id": "8cd8072dcfc53a062500867216b3c76240dff94d",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2f92e5bde0d2847c0e4dd75cedcceb18e44d6ad2"
        },
        "date": 1728114263097,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7049381.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10733272\nallocs=350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 599619107.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20092336\nallocs=25\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 71139818.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3666032\nallocs=187\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13699594.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1446136\nallocs=3346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43422646,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4822592\nallocs=465\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1358764246,
            "unit": "ns",
            "extra": "gctime=29980164\nmemory=3468572752\nallocs=163541\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4713865,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237568\nallocs=3099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26383309,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3700512\nallocs=771\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 278270629.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8318616\nallocs=99772\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 187295763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8341016\nallocs=100099\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "87ce64e13b82c10f3886e1dcda67e9c1edfb1afe",
          "message": "Dispatch progress bar method in EnsembleProblems (#259)",
          "timestamp": "2024-10-08T13:58:47+08:00",
          "tree_id": "93fae91c400dfa6f9b5c7eab76ca865b3182bed1",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/87ce64e13b82c10f3886e1dcda67e9c1edfb1afe"
        },
        "date": 1728367695673,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6885455,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 592395330.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79009714,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13818746,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1428248\nallocs=3753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42640064,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1990841790.5,
            "unit": "ns",
            "extra": "gctime=632101101.5\nmemory=3504189664\nallocs=303327\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4496676.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=228032\nallocs=3135\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26495943,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3695120\nallocs=984\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 268643571,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8602720\nallocs=114216\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 139026165,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8490928\nallocs=112404\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "bdedbda99469df8014bbaf4920fbb61328ae5ed2",
          "message": "Fix default value of keyword parameter `saveat` (#260)\n\n* fix default value of keyword parameter `saveat`\r\n\r\n* fix typos in docstrings",
          "timestamp": "2024-10-08T11:01:24+02:00",
          "tree_id": "82b9f61c1c9a66cec3b1a1766db95ed15b45fcaa",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/bdedbda99469df8014bbaf4920fbb61328ae5ed2"
        },
        "date": 1728378658624,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7033206.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601501716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80592688,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14033772,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1428248\nallocs=3753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42606916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2039508129.5,
            "unit": "ns",
            "extra": "gctime=670855463\nmemory=3504189664\nallocs=303327\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4522832,
            "unit": "ns",
            "extra": "gctime=0\nmemory=228032\nallocs=3135\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27010364.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3695120\nallocs=984\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 271310412.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8623152\nallocs=114717\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140586724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8275552\nallocs=108934\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "ddc98fcb962dd8d7bc7598e2a0a0757c53dd358c",
          "message": "Improve random number generation on mcsolve and ssesolve (#263)",
          "timestamp": "2024-10-10T01:48:02+02:00",
          "tree_id": "c5e01c8f514df50c9e373bf48d0700a0ab2eea48",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/ddc98fcb962dd8d7bc7598e2a0a0757c53dd358c"
        },
        "date": 1728518107531,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6815775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601376402.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 77223460,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14012226.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1428248\nallocs=3753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42471923.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1921048656,
            "unit": "ns",
            "extra": "gctime=635530264.5\nmemory=3504189664\nallocs=303327\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4532766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=228032\nallocs=3135\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26710220,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3695120\nallocs=984\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 270799800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8582496\nallocs=113068\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140017352,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8509440\nallocs=111737\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "110c2e5715646c3f7dcfe73d9efd073131e261ff",
          "message": "Bump version to v0.18.0",
          "timestamp": "2024-10-10T01:48:51+02:00",
          "tree_id": "86e2edc4cba84a21280d359092332486ba749f32",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/110c2e5715646c3f7dcfe73d9efd073131e261ff"
        },
        "date": 1728518163030,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6861450.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602697687.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78717411,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14000776,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1428248\nallocs=3753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42538970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2026697339,
            "unit": "ns",
            "extra": "gctime=658902665\nmemory=3504189664\nallocs=303327\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4494780.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=228032\nallocs=3135\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26796075.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3695120\nallocs=984\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 270318835,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8614344\nallocs=113846\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 140057209.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8603192\nallocs=113280\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "5e94453a40808b18ac3c1a7d2e39892fd57c7966",
          "message": "Improve mcsolve performance (#265)",
          "timestamp": "2024-10-13T03:06:34+02:00",
          "tree_id": "ae125c7e0caf8560d78052c018142caf2dfc044a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5e94453a40808b18ac3c1a7d2e39892fd57c7966"
        },
        "date": 1728782016911,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6898592,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595683266,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 70106938,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13955121,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1428248\nallocs=3753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42572618,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1914946002,
            "unit": "ns",
            "extra": "gctime=622207401\nmemory=3504189664\nallocs=303327\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4555812,
            "unit": "ns",
            "extra": "gctime=0\nmemory=228032\nallocs=3135\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27070665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3695120\nallocs=984\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 271249799,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8652176\nallocs=113067\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 138418373,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8670800\nallocs=113974\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "d3f131354c3030eb906bb0f4bcfc434c46d11f0a",
          "message": "Make ProgressBar mutable (#268)",
          "timestamp": "2024-10-15T22:29:16+09:00",
          "tree_id": "502604f39fbd13deff7a04a14e5cba09e9b44160",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d3f131354c3030eb906bb0f4bcfc434c46d11f0a"
        },
        "date": 1728999393998,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7015455,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603172117.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80189992.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13979589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1347128\nallocs=2753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42628931,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2054901418.5,
            "unit": "ns",
            "extra": "gctime=676248757\nmemory=3504179968\nallocs=303227\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4472658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147552\nallocs=2135\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26857354,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686640\nallocs=884\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 263431403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7343880\nallocs=102772\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 134791215,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7379608\nallocs=103390\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "33577062136b89dde85b9aefbbc93a9b4889b93d",
          "message": "Introduce `QobjEvo` and use `SciMLOperators` for time evolution (#266)\n\n* First working implementation\r\n\r\n* Minor changes\r\n\r\n* First working case of sesolve\r\n\r\n* Minor changes\r\n\r\n* Rebase commits\r\n\r\n* Apply Yi-Te comments\r\n\r\n* Working mesolve\r\n\r\n* Make dsf_mesolve and dfd_mesolve work\r\n\r\n* Minor changes\r\n\r\n* Working version of ssesolve\r\n\r\n* Remove sleep in runtests\r\n\r\n* Remove OperatorSum and TimeDependentOperatorSum\r\n\r\n* Remove alg as argument from all problem definitions\r\n\r\n* First working  runtests\r\n\r\n* Add tests for QObjEvo and time evolution\r\n\r\n* Add docstrings to documentation\r\n\r\n* Reduce tolerance for stochastic runtests\r\n\r\n* Make runtests multithreaded and reduce time for timeevolution tests\r\n\r\n* Add QobjEvo tests and minor changes\r\n\r\n* Update docstrings\r\n\r\n* Add comment on the use of MatrixOperator",
          "timestamp": "2024-10-22T10:39:16+02:00",
          "tree_id": "68be237c99ee29b09fabc25d5908931dc041b81e",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/33577062136b89dde85b9aefbbc93a9b4889b93d"
        },
        "date": 1729586916636,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6992229,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595600215.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79635126.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13937868,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42319039,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2022806646,
            "unit": "ns",
            "extra": "gctime=664562230.5\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4415825,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26914644,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 248924176,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4939360\nallocs=66892\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 126430577,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4945424\nallocs=67208\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "f19741e478a23be9961e10de97542fdcc05b855f",
          "message": "add deprecate file (#274)",
          "timestamp": "2024-10-25T13:53:35+02:00",
          "tree_id": "7e7b88a0747ac05d0ce8451527371335584449d8",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/f19741e478a23be9961e10de97542fdcc05b855f"
        },
        "date": 1729857760179,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7134419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 590969008.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79688312,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14341680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42735011,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2133673457,
            "unit": "ns",
            "extra": "gctime=781984091\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4494467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27005417,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 245504144.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4941728\nallocs=67040\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 125210846,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4948688\nallocs=67412\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "f5b87fdb92de9119c8e6ce1025d694e332417c5d",
          "message": "fix error message in `mesolve` (#273)",
          "timestamp": "2024-10-25T14:00:36+02:00",
          "tree_id": "a137a9a593f41dd95ef389b3c34b896348784065",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/f5b87fdb92de9119c8e6ce1025d694e332417c5d"
        },
        "date": 1729858171649,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7131233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600879551.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 70783621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13996363,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42829652,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2198616985.5,
            "unit": "ns",
            "extra": "gctime=793865889\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4463572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27149948.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 251162536,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4934624\nallocs=66596\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 125459069,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4943232\nallocs=67071\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "df3b4f72f31ac0070b142a2a2e86323c80547980",
          "message": "Add `show` method for `QobjEvo` (#272)",
          "timestamp": "2024-10-25T21:27:46+09:00",
          "tree_id": "e61883eb77a44fbbd1a1e55b039c3e9eaedb589a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/df3b4f72f31ac0070b142a2a2e86323c80547980"
        },
        "date": 1729859495092,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7564324,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598547692.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79429773,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14346540,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42178235,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2110305427.5,
            "unit": "ns",
            "extra": "gctime=727225675.5\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4471303,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27040224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 245359703.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4931104\nallocs=66376\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 124711853,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4942448\nallocs=67022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "b5ce558695260b2b38030ee826fa7cbd1541ccac",
          "message": "bump version to `0.19.0`",
          "timestamp": "2024-10-25T21:28:36+09:00",
          "tree_id": "0549731c76f77f3e3153a6c5c9636491b961190e",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/b5ce558695260b2b38030ee826fa7cbd1541ccac"
        },
        "date": 1729859628540,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7132906.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 592899654.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79707346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13933318.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42481449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2121443800,
            "unit": "ns",
            "extra": "gctime=761136416\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4489385.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26985998.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 247526415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4942112\nallocs=67064\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 125708363.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4939872\nallocs=66861\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "5820791493c66813774edacc52339859976dc1ed",
          "message": "add `SciMLOperators` to `versioninfo` (#275)",
          "timestamp": "2024-10-27T15:04:10+01:00",
          "tree_id": "fd8736a6639e491e56daaebdabdc0ff5f35faaa6",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5820791493c66813774edacc52339859976dc1ed"
        },
        "date": 1730038315780,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6986048,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602403017,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67468353,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13926738,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42385285.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2191374309.5,
            "unit": "ns",
            "extra": "gctime=763311858.5\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4470694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26894682,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 248197686.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4931040\nallocs=66372\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 125237725,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4935712\nallocs=66601\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "e579ada202d1583abca2fe3161c6d1199ebeca14",
          "message": "Fix issue in mesolve for nothing c_ops (#276)",
          "timestamp": "2024-10-29T09:51:30+09:00",
          "tree_id": "fcafb202dca3eff4fccff9bf26b9f958dbcf6ab0",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e579ada202d1583abca2fe3161c6d1199ebeca14"
        },
        "date": 1730163560316,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7304050,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 592426222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79060029,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14000240,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42289899,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2078506082.5,
            "unit": "ns",
            "extra": "gctime=702194960.5\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4466041,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26722243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 248286449.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4946512\nallocs=67339\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 125714790,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4937408\nallocs=66707\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "3207e188de07cdd1e4d8f31b2d596d89a2e505dd",
          "message": "bump version to `0.19.1`",
          "timestamp": "2024-10-29T09:52:09+09:00",
          "tree_id": "45b9421dce704f5357cc9242447483561795f7db",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/3207e188de07cdd1e4d8f31b2d596d89a2e505dd"
        },
        "date": 1730163592631,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6911471,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10730792\nallocs=477\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596909097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78672987,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13990099.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1346568\nallocs=2760\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42437437.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4809216\nallocs=1023\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 2031231604,
            "unit": "ns",
            "extra": "gctime=702443720\nmemory=3504117648\nallocs=302381\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4464025.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27361236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3686176\nallocs=879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 254883300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4938144\nallocs=66816\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 128893626.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4942096\nallocs=67000\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "25bee9fed9c106ccb8e5ffa515e1eb9b1598a055",
          "message": "Support time-dependent `liouvillian` (#277)\n\n* support time-dependent `liouvillian`\r\n\r\n* improve codecov\r\n\r\n* format files\r\n\r\n* minor changes",
          "timestamp": "2024-10-30T14:44:40+01:00",
          "tree_id": "67fd422f662615c25ac78882c20057bf20069904",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/25bee9fed9c106ccb8e5ffa515e1eb9b1598a055"
        },
        "date": 1730296496460,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6980143.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284488\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603132625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 68571666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13878648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293480\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42731935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791520\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1829130853,
            "unit": "ns",
            "extra": "gctime=553024246\nmemory=3066240208\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4421059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26750520,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109120\nallocs=858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 252931477,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4945056\nallocs=67248\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 125574431,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4939472\nallocs=66836\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "dea0bd7d3d55e3da7ae9dbb2ef3c3ed54bfb5ae2",
          "message": "Clean low-rank mesolve and add Documentation bibliography (#269)\n\n* change Lattice structure\r\n\r\n* Change lr_mesolve and add documentation\r\n\r\n* Fix Documentation errors",
          "timestamp": "2024-10-31T11:15:39+01:00",
          "tree_id": "7c7a32d7f77a82e760a1f6050bad63466be50071",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/dea0bd7d3d55e3da7ae9dbb2ef3c3ed54bfb5ae2"
        },
        "date": 1730370222539,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7430376.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284488\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594686640.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78504598.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14046907.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293480\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42372795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791520\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1799561083,
            "unit": "ns",
            "extra": "gctime=547093898\nmemory=3066240208\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4442013,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27042207,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109120\nallocs=858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 248607303.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4943744\nallocs=67166\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 125120044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4948496\nallocs=67400\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "892c496803c9e097e9461f5426ada4438234bb70",
          "message": "Make clenshaw method as default for wigner (#279)\n\n* Make clenshaw method as default for wigner\r\n\r\n* Minor changes",
          "timestamp": "2024-10-31T11:15:56+01:00",
          "tree_id": "967ce303a7a86d1b7077a4cb10bf08ac6313b450",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/892c496803c9e097e9461f5426ada4438234bb70"
        },
        "date": 1730370590440,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7072112,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284488\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595690533,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80103754,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13950590.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293480\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42672209,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791520\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1872882356,
            "unit": "ns",
            "extra": "gctime=568178874\nmemory=3066240208\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4425230,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26610741,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109120\nallocs=858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 248323492,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4944608\nallocs=67220\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 126499581,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4945328\nallocs=67202\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "33944ac66fceba16584f2acaddbb54d15c148717",
          "message": "Bump to v0.20.0",
          "timestamp": "2024-10-31T11:16:39+01:00",
          "tree_id": "b8c5a8c22ce24e95448b313306be5025c2726b02",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/33944ac66fceba16584f2acaddbb54d15c148717"
        },
        "date": 1730370633890,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7020610.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284488\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 592494024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 66823668.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13891035,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293480\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42344430,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791520\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1962334289.5,
            "unit": "ns",
            "extra": "gctime=613669189\nmemory=3066240208\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4450354,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147408\nallocs=2136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26866836.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109120\nallocs=858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 248148036.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4950000\nallocs=67557\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 126284257,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4933008\nallocs=66432\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "380d017eaffd589da418245018e5e2e34ecac003",
          "message": "Support time-dependent `lindblad_dissipator` (#280)",
          "timestamp": "2024-11-02T09:10:52+09:00",
          "tree_id": "0c21a36227c237bb82df0acee49e03c68d8478e3",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/380d017eaffd589da418245018e5e2e34ecac003"
        },
        "date": 1730506881288,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7087162,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596069218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79407125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14119133,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42614593,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1970584550,
            "unit": "ns",
            "extra": "gctime=627970138\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4449376,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27112646,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 254065142,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4937312\nallocs=66764\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 128776301.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4931696\nallocs=66350\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "54ded38dc069dd7e3c4e30533003636a4a185ccd",
          "message": "A lazy method to construct `QobjEvo` (#282)\n\n* add a lazy method to construct `QobjEvo`\r\n\r\n* format files\r\n\r\n* improve tests to make sure time-independent problems are `MatrixOperator`\r\n\r\n* modify docstrings\r\n\r\n* fix typo\r\n\r\n* fix extra `ScalarOperator` in `mcsolve`\r\n\r\n* rebase to `ScaledOperator` in time-independent `mcsolve`\r\n\r\n* add missing import in test",
          "timestamp": "2024-11-02T18:16:23+01:00",
          "tree_id": "238c3138a35774a7a87300a28d478c96724ce215",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/54ded38dc069dd7e3c4e30533003636a4a185ccd"
        },
        "date": 1730568007459,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7154031.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 587998097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79085044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13892479.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42095602.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1791917625.5,
            "unit": "ns",
            "extra": "gctime=530117898.5\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4493217,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26896204,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 246643256.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4945744\nallocs=67291\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 124562257,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4931920\nallocs=66364\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "70209fe55f7f1c316ed9c390a82cdf4529196b75",
          "message": "Improve SciMLOperator handling of sesolveProblem (#283)",
          "timestamp": "2024-11-03T11:35:36+09:00",
          "tree_id": "eac11161b6303dd27f1341c16faa047675aefbfd",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/70209fe55f7f1c316ed9c390a82cdf4529196b75"
        },
        "date": 1730601790944,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6918559,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600051113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79564608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13971026.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42255864,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1796810581,
            "unit": "ns",
            "extra": "gctime=532126995.5\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4457545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26782084,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 233421916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4827232\nallocs=66468\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 117792258.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4821824\nallocs=66067\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "c63e627568f4349baa783d7efa32ff4f49f74f39",
          "message": "Bump version to `0.21.0`",
          "timestamp": "2024-11-04T11:54:27+09:00",
          "tree_id": "42e8181b5d453593958c7783413680a7f3e02555",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/c63e627568f4349baa783d7efa32ff4f49f74f39"
        },
        "date": 1730689094324,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6985397,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601251371.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67793469,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14002972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42535174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1922591015,
            "unit": "ns",
            "extra": "gctime=582522486.5\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4471162,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27145550,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 231292466,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4813376\nallocs=65602\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 117828376.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4824896\nallocs=66259\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "2adcd8a94f4d4376ae3be945170a6e12cc356acc",
          "message": "Optimize precompilation during runtests (#284)",
          "timestamp": "2024-11-05T09:23:57+09:00",
          "tree_id": "c6b7f1416ae91fa867a01c5451f15bafe5649b8b",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/2adcd8a94f4d4376ae3be945170a6e12cc356acc"
        },
        "date": 1730766855328,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7602417,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594238617.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 78540075,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13989658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42516183,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1853018015,
            "unit": "ns",
            "extra": "gctime=559036350.5\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4459825.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26619174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 231145535,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4830880\nallocs=66696\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 118016231.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4821360\nallocs=66038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "9b23c6f10a71bbf6c90574ad5b5dab8261e6f61b",
          "message": "Add normalize_states option in mcsolve (#285)",
          "timestamp": "2024-11-05T08:32:49+01:00",
          "tree_id": "8169672696a06fbf36a1fda10bc22ef40a1b6ef4",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9b23c6f10a71bbf6c90574ad5b5dab8261e6f61b"
        },
        "date": 1730792198424,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6884477,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598176721.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 72091589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13967296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42339590,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1827993734.5,
            "unit": "ns",
            "extra": "gctime=543716718\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4464604,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27132002,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 232769662,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4827024\nallocs=66455\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 117363728,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4817248\nallocs=65781\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "54e85c9fd2e26870395e4e50712c12a6635568d0",
          "message": "Add spell check CI and fix typos (#286)\n\n* add spell check CI\r\n\r\n* fix version number\r\n\r\n* fix version number\r\n\r\n* add typo ignore\r\n\r\n* fix typos\r\n\r\n* minor changes",
          "timestamp": "2024-11-05T13:55:46+01:00",
          "tree_id": "13b1bf6862f479b381998ba6a27fe7e5122a5851",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/54e85c9fd2e26870395e4e50712c12a6635568d0"
        },
        "date": 1730811804787,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6899676.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596349675,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79265270,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13956660,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42288333.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1829785526.5,
            "unit": "ns",
            "extra": "gctime=552750345.5\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4464497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26814059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 231713487,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4833520\nallocs=66861\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 117416245,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4829776\nallocs=66564\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "c8c9479fa868b0ac4c2eb687a8eba15b20174bd5",
          "message": "Bump version to v0.21.1",
          "timestamp": "2024-11-05T13:56:30+01:00",
          "tree_id": "7aef8e33b22c838cc6ce44e74a40becafad8bbfd",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/c8c9479fa868b0ac4c2eb687a8eba15b20174bd5"
        },
        "date": 1730811849377,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6949424,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 603748065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79685368,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14033448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1293960\nallocs=2730\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42704384,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1830064707,
            "unit": "ns",
            "extra": "gctime=556340609.5\nmemory=3066234224\nallocs=301583\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4448682,
            "unit": "ns",
            "extra": "gctime=0\nmemory=147056\nallocs=2129\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26942299.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3109280\nallocs=856\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 232510702,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4824656\nallocs=66307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 117903342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4830176\nallocs=66589\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "ca872f65c9631799a086efcafab22b6b5d3b3e3e",
          "message": "Add doctest in all docstrings (#294)",
          "timestamp": "2024-11-11T18:14:17+09:00",
          "tree_id": "2032289e80664ddb0baa235829045b7bf3605e0a",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/ca872f65c9631799a086efcafab22b6b5d3b3e3e"
        },
        "date": 1731317072002,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7369056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 604147020,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79741562,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13829655,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262056\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42465066,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1831649592.5,
            "unit": "ns",
            "extra": "gctime=543093271\nmemory=3066228928\nallocs=301318\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4368973,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106640\nallocs=113\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27123227,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104240\nallocs=602\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 228949357,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3903952\nallocs=14763\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 116328369,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3905152\nallocs=14775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "d11841dff0b3cc551e77563a27597910316e4817",
          "message": "Bump to v0.21.2",
          "timestamp": "2024-11-11T11:25:57+01:00",
          "tree_id": "7ae0bbd40f740ca17fac623a94218ec5d088ab2f",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/d11841dff0b3cc551e77563a27597910316e4817"
        },
        "date": 1731320977834,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7683661.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596996780.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79665185,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14250255.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262056\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42187830.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1852019938,
            "unit": "ns",
            "extra": "gctime=552706637\nmemory=3066228928\nallocs=301318\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4388149,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106640\nallocs=113\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27084732,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104240\nallocs=602\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 228612991,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3903952\nallocs=14763\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 116708785,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3905152\nallocs=14775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "eaddad0f2a504dbf0eb33b3bcd9e182a1f4bcd1a",
          "message": "Bump version to `0.21.3`",
          "timestamp": "2024-11-11T22:39:42+09:00",
          "tree_id": "573c013196618f2068982c4373ab20cf8f0d0314",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/eaddad0f2a504dbf0eb33b3bcd9e182a1f4bcd1a"
        },
        "date": 1731332727412,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6965824,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593019995,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80393078.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14003112,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262056\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42392082.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1777946114.5,
            "unit": "ns",
            "extra": "gctime=517869483.5\nmemory=3066228928\nallocs=301318\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4360640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106640\nallocs=113\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27359718,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104240\nallocs=602\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 228744254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3903952\nallocs=14763\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 116599201,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3905152\nallocs=14775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "43535b01d625ce171c83c0c1eccf2720a97a384e",
          "message": "Bump version to v0.21.4",
          "timestamp": "2024-11-13T14:52:58+01:00",
          "tree_id": "e5fabfc7e7d7cbfd2d96132b4e28c353d4552e33",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/43535b01d625ce171c83c0c1eccf2720a97a384e"
        },
        "date": 1731506579289,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6956255,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593405093,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 73497948.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14279537,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262056\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42500702,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1806628669,
            "unit": "ns",
            "extra": "gctime=532303876.5\nmemory=3066228928\nallocs=301318\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4412423,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106640\nallocs=113\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27294163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104240\nallocs=602\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 229857635,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3903952\nallocs=14763\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 117167609,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3905152\nallocs=14775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "afded9ad739b13ad3e35e26bf6ccff239f361c46",
          "message": "Bump version to `v0.21.5` (#309)\n\n* bump version to `v0.21.5`\r\n\r\n* bump version to `v0.21.5`",
          "timestamp": "2024-11-14T16:39:46+01:00",
          "tree_id": "b9ae2b42fab6235b67098bc78fb74b211c5dbdca",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/afded9ad739b13ad3e35e26bf6ccff239f361c46"
        },
        "date": 1731599241616,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7120779.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598768279.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79510065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13956614,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262056\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42911079,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1879591085,
            "unit": "ns",
            "extra": "gctime=563693644\nmemory=3066228928\nallocs=301318\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4351607.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106640\nallocs=113\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27023018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3104240\nallocs=602\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 230123157,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3903952\nallocs=14763\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 116437426.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3905152\nallocs=14775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "283669648d286d20f24568fc1ed7143dd3a7d60d",
          "message": "Make time evolution solvers compatible with automatic differentiation (#311)\n\n* Working sesolve\r\n\r\n* add `inplace` keywork argument\r\n\r\n* add SciMLStructures and relax params type\r\n\r\n* Working mcsolve (no type-stability)\r\n\r\n* Fix type-instabilities for mcsolve\r\n\r\n* Add SciMLStructures.jl methods\r\n\r\n* Add callbacks helpers\r\n\r\n* Fix dsf_mcsolve\r\n\r\n* Remove ProgressBar from ODE parameters\r\n\r\n* Fix abstol and reltol extraction\r\n\r\n* Use Base allequal function\r\n\r\n* Remove expvals from TimeEvolutionParameters\r\n\r\n* Make NullParameters as default for params\r\n\r\n* Remove custom PresetTimeCallback\r\n\r\n* Update description of `inplace` argument\r\n\r\n* Working mesolve\r\n\r\n* Fix dfd_mesolve and dsf_mesolve\r\n\r\n* Remove TimeEvolutionParameters (type-unstable)\r\n\r\n* Fix type instabilities\r\n\r\n* Fix type instabilities on Julia v1.10\r\n\r\n* Format document",
          "timestamp": "2024-11-18T10:07:28+01:00",
          "tree_id": "ea2e73ba1849bdc87bfc3589470132d3ab00f1c3",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/283669648d286d20f24568fc1ed7143dd3a7d60d"
        },
        "date": 1731921329317,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6955854,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 597999838.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 66830792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13896581,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42341561.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1746646608,
            "unit": "ns",
            "extra": "gctime=480969286\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4340992,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27238196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222837067,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112558950.5,
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
          "id": "9810db4732d5f9bcefbc73ea2dfd72121b8b5e14",
          "message": "Fix type instability for `liouvillian` (#318)\n\n* add type stability tests for liouvillian\r\n\r\n* fix type instability and reduce memory alloc in `liouvillian`\r\n\r\n* fix type instability\r\n\r\n* minor changes\r\n\r\n* update changelog",
          "timestamp": "2024-11-20T09:23:34+01:00",
          "tree_id": "959093f4e6ad0da2b5d6d49f02bb34e9ecd25fa7",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9810db4732d5f9bcefbc73ea2dfd72121b8b5e14"
        },
        "date": 1732091446259,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6816159,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 592208409,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 66435478,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13932501,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42408045.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1793671902.5,
            "unit": "ns",
            "extra": "gctime=480018488\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4334263,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27068465.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220853910.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 111705642,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3564304\nallocs=12474\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "9d1cb6b183751df92f34d0ac032f5beb43b12ca8",
          "message": "Bump to v0.22.0 (#319)\n\n* Bump to v0.22.0\r\n\r\n* Add Unreleased section",
          "timestamp": "2024-11-20T09:54:50+01:00",
          "tree_id": "6160ca9f48c14305edae0cba44092efbbbf9ec5f",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9d1cb6b183751df92f34d0ac032f5beb43b12ca8"
        },
        "date": 1732093101140,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7546725.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 593339272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79755175,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13976563,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42612583,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1725245513,
            "unit": "ns",
            "extra": "gctime=455259190.5\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4338599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27314011,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220308601.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112162872,
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
          "id": "ab66af58a6be14debf8c58f56bb1ebafc59127de",
          "message": "extend methods for `_FType` and `_CType` (#320)",
          "timestamp": "2024-11-22T11:56:22+01:00",
          "tree_id": "f7a36a683c3ce0fb4a1580cad4de5fd05fdd397c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/ab66af58a6be14debf8c58f56bb1ebafc59127de"
        },
        "date": 1732273410065,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6997342.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601392120.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 66299775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14144653,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42421759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1747913778,
            "unit": "ns",
            "extra": "gctime=472559841\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4390956,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26788767,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 223022429.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112247250,
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
          "id": "aee73e6065aa043f6cb81999e01d82b476fe3e3c",
          "message": "fix latex in docstrings (#323)",
          "timestamp": "2024-11-29T16:04:56+01:00",
          "tree_id": "1da2a1f733f1cf7be164d1a1f9dc11bcbcdf9a55",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/aee73e6065aa043f6cb81999e01d82b476fe3e3c"
        },
        "date": 1732893126425,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7163140,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 605614267.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 73422587,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14090971.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43129778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1815882135,
            "unit": "ns",
            "extra": "gctime=515140010\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4394016,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27605436,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220326527,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112317613,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3564304\nallocs=12474\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "b314479ce4e49b52bc25a4816f0d0ba67054ce83",
          "message": "Change SingleSiteOperator with MultiSiteOperator (#324)\n\n* Change SingleSiteOperator with MultiSiteOperator\r\n\r\n* Update changelog",
          "timestamp": "2024-11-29T16:05:29+01:00",
          "tree_id": "27dfd8c1b07e6e0c8ab31813fbdb54dbd86dda60",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/b314479ce4e49b52bc25a4816f0d0ba67054ce83"
        },
        "date": 1732893148555,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6826780.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 592702690.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80549198,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13902741.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42405007.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1762149931,
            "unit": "ns",
            "extra": "gctime=475248826\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4383414,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27078305,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220282250.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112497559,
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
          "id": "9a9a80d7593c46d60607de3daa3cfe1c31ed0168",
          "message": "fix typo in `spectrum` docstrings (#325)",
          "timestamp": "2024-11-29T16:52:05+01:00",
          "tree_id": "27d79e4a4867708427305040436aa06d449d1c7c",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/9a9a80d7593c46d60607de3daa3cfe1c31ed0168"
        },
        "date": 1732895748139,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6798182,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 598669639.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 68927662.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13962619,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42751585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1749275087,
            "unit": "ns",
            "extra": "gctime=483759719\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4357331,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26965089.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219801398.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112532290,
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
          "id": "e56afac54b56e6878244d01da5f3eedd939b7889",
          "message": "dummy change for benchmark CI",
          "timestamp": "2024-11-30T14:19:32+09:00",
          "tree_id": "d9ed161b5a6ddb550af779a5ab84f6976b6dd372",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e56afac54b56e6878244d01da5f3eedd939b7889"
        },
        "date": 1732944401000,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 7035316.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602870004,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 68036713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13831745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1262456\nallocs=1246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 0,
            "unit": "ns",
            "extra": ""
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42361084,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1924207724,
            "unit": "ns",
            "extra": "gctime=604543176\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4337312,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26761101.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219377635,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112067082,
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
          "id": "cfe79439bac2327f5ae85be0db08da48c5b20eaf",
          "message": "Make `spectrum` and `correlation` align with QuTiP (#330)\n\n* remove `MySolver` and don't return frequency list from `spectrum`\r\n\r\n* remove `FFTCorrelation` and introduce `spectrum_correlation_fft`\r\n\r\n* format files\r\n\r\n* introduce `PseudoInverse<:SpectrumSolver`\r\n\r\n* update changelog\r\n\r\n* fix benchmark\r\n\r\n* introduce `correlation_3op_1t`\r\n\r\n* format files\r\n\r\n* fix missing links in docstring of `correlation` functions\r\n\r\n* update documentation page for two-time correlation function\r\n\r\n* update doc page of two-time correlation functions\r\n\r\n* make the output correlation `Matrix` align with `QuTiP`\r\n\r\n* minor changes\r\n\r\n* add deprecated methods to `deprecated.jl`\r\n\r\n* update changelog\r\n\r\n* improve runtests for `spectrum_correlation_fft`",
          "timestamp": "2024-12-04T11:08:55+01:00",
          "tree_id": "3e69e0beed20a6d2957a6dad1be22389175c0f64",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/cfe79439bac2327f5ae85be0db08da48c5b20eaf"
        },
        "date": 1733307494812,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6995544,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 605266682.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80185281.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29902466,
            "unit": "ns",
            "extra": "gctime=7615494\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44055402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14134923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1743696774.5,
            "unit": "ns",
            "extra": "gctime=361447582\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4381472,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27496141,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220581329.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112850562,
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
          "id": "e2e22492efb2a3ff61c125a69a1456fd891fbb42",
          "message": "bump version to `v0.23.0` (#331)",
          "timestamp": "2024-12-04T20:09:19+09:00",
          "tree_id": "5d70a9c53a0242fbb07edc924d37ded783e410ba",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e2e22492efb2a3ff61c125a69a1456fd891fbb42"
        },
        "date": 1733310788922,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6707483.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594116651,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79302055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 27716352,
            "unit": "ns",
            "extra": "gctime=6904427\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43778654,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14095415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1680260354.5,
            "unit": "ns",
            "extra": "gctime=336497135.5\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4408151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27446744,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220894674.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112498184,
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
          "id": "0f2c4a1978d4697abcbf2ec61b98470943f735dd",
          "message": "Fix `[compat]` of `DiffEqCallbacks.jl` (#335)\n\n* fix `[compat]` of `DiffEqCallbacks`\r\n\r\n* fix changelog",
          "timestamp": "2024-12-06T09:55:38+01:00",
          "tree_id": "b58452dbeb977a0ab3dc07d19d8f0b2b97de49d2",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/0f2c4a1978d4697abcbf2ec61b98470943f735dd"
        },
        "date": 1733475795515,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6894988,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 599541149,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 67677757,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 29491399,
            "unit": "ns",
            "extra": "gctime=8353526\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42616579,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13732447,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1736128843.5,
            "unit": "ns",
            "extra": "gctime=379546725.5\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4341170.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26729549.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 219995521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112440453.5,
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
          "id": "5e108e19a52715915e236f37f7a300069d342ef3",
          "message": "Improve the construction of `QobjEvo` (#339)",
          "timestamp": "2024-12-10T09:18:21+01:00",
          "tree_id": "3e110eae40d2c0ab56633d76f36ddc1a8170282e",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/5e108e19a52715915e236f37f7a300069d342ef3"
        },
        "date": 1733819129377,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6716573.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 600316893,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80105244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28285589,
            "unit": "ns",
            "extra": "gctime=6754705\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43691885,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13823003,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1652145126,
            "unit": "ns",
            "extra": "gctime=343272210\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4357683,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27643024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220222325.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112698532,
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
          "id": "e754148481e5fe4d69d276cf00068303b8c20cbd",
          "message": "Support `zero` and `one` for `AbstractQuantumObject` (#346)\n\n* support `zero` and `one` for `AbstractQuantumObject`\r\n\r\n* update changelog\r\n\r\n* fix changelog",
          "timestamp": "2024-12-11T09:08:32+01:00",
          "tree_id": "0f42b5b6f37539f7e02b6adc60a3fa33383545c5",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/e754148481e5fe4d69d276cf00068303b8c20cbd"
        },
        "date": 1733904748694,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6729962.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 599661863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 66490216,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31562550.5,
            "unit": "ns",
            "extra": "gctime=7921962\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43628898,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13958328,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1731855377.5,
            "unit": "ns",
            "extra": "gctime=358889985.5\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4359105.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27659864,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222023978.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112279351,
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
          "id": "38a473c91e4854144449e89ccb3079ea8bc3624c",
          "message": "Introduce visualization (#347)\n\n* Add Wigner plotting via CairoMakie\r\n\r\n* CairoMakie default lib\r\n\r\n* docs\r\n\r\n* formatting\r\n\r\n* Fix typo\r\n\r\n* Add error handling for unavailable plotting libraries in plot_wigner function\r\n\r\n* Docs\r\n\r\n* Update CairoMakie dependency and add CairoMakie extension tests\r\n\r\n* Update documentation to reflect changes in Wigner function plotting and computation\r\n\r\n* formatting\r\n\r\n* Update changelog\r\n\r\n* Improve docs\r\n\r\n* Update CairoMakie imports and change axis labels\r\n\r\n* fix CHANGELOG\r\n\r\n* setup CI pipeline config for `CairoMakie` extension\r\n\r\n* rebase accidently removed lines\r\n\r\n* fix typo\r\n\r\n* add `CairoMakie` extension doc page\r\n\r\n* fix typo\r\n\r\n* minor change in docs\r\n\r\n* minor change in docs\r\n\r\n---------\r\n\r\nCo-authored-by: Lorenzo Fioroni <lorenzo.fioroni@epfl.ch>\r\nCo-authored-by: Alberto Mercurio <alberto.mercurio96@gmail.com>",
          "timestamp": "2024-12-13T10:08:01+01:00",
          "tree_id": "583bf4f76d8345807b1fa066f0e4313b7aec1bf1",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/38a473c91e4854144449e89ccb3079ea8bc3624c"
        },
        "date": 1734081440236,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6784381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595270504,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79876401,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 26810910,
            "unit": "ns",
            "extra": "gctime=7332686.5\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42462115,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13698073,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1675977917,
            "unit": "ns",
            "extra": "gctime=347583102.5\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4344897.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26546053,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222315173,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 111631355,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3564304\nallocs=12474\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "91cf21140bc5f67e61efaa5a85695b81b75dcb19",
          "message": "Bump version to v0.24.0 (#348)",
          "timestamp": "2024-12-13T10:12:49+01:00",
          "tree_id": "fc37dc4ea0e3fbe273d41bfcb64c529aa8e80372",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/91cf21140bc5f67e61efaa5a85695b81b75dcb19"
        },
        "date": 1734081731271,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6854770,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 602958967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80039096,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31569756,
            "unit": "ns",
            "extra": "gctime=9004550\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42949096,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14102218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1697182439.5,
            "unit": "ns",
            "extra": "gctime=367453953\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4383453,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26833151.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220782556,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113398837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3564304\nallocs=12474\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "3ffdd3debfdb31e7b21794626467e4f93c5f1586",
          "message": "Change block diagonal form structure (#349)\n\n* Change block diagonal form structure\r\n\r\n* Add changelog",
          "timestamp": "2024-12-15T10:19:41+01:00",
          "tree_id": "11744d51984a7447470134554e4c2802effc3eb6",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/3ffdd3debfdb31e7b21794626467e4f93c5f1586"
        },
        "date": 1734254943172,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6897935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 606142434.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79284961,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 31518004.5,
            "unit": "ns",
            "extra": "gctime=9337617.5\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 42775718,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13824144,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1745321341.5,
            "unit": "ns",
            "extra": "gctime=375423469\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4351347,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 26780515,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 222312689,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112013142,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3564304\nallocs=12474\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "1f857be3d80aff90b5d643820726db1f7790e04e",
          "message": "Add `ptrace` support for any `GPUArray` (#350)\n\n* Add `ptrace` support for any `GPUArray`\r\n\r\n* Add changelog\r\n\r\n---------\r\n\r\nCo-authored-by: Yi-Te Huang <44385685+ytdHuang@users.noreply.github.com>",
          "timestamp": "2024-12-15T13:53:49+01:00",
          "tree_id": "38b18034f77257fc319a07ba17efe84a613331b4",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/1f857be3d80aff90b5d643820726db1f7790e04e"
        },
        "date": 1734267465362,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6803867,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 595464391,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 79102849,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28991918,
            "unit": "ns",
            "extra": "gctime=7410263\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43953029,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13909178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1686960121,
            "unit": "ns",
            "extra": "gctime=364759837\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4354668,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27761594,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 221674016,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112415762,
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
          "id": "ade8b45b24dc0b030c35b3b27035f5f9f0f63417",
          "message": "CompatHelper: bump compat for GPUArrays in [weakdeps] to 11, (keep existing compat) (#351)",
          "timestamp": "2024-12-16T11:41:25+09:00",
          "tree_id": "252a6113a8ffa2435dec4e3b404bfedca1e6dad4",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/ade8b45b24dc0b030c35b3b27035f5f9f0f63417"
        },
        "date": 1734317277607,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6943284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 596830778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 77135458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 34296217,
            "unit": "ns",
            "extra": "gctime=9059951\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 44222946,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 14246733.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1717367406.5,
            "unit": "ns",
            "extra": "gctime=369319841\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4409043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27698364,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 224210133.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 113362507,
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
          "id": "270d3f69f79a171b9c5feec9cfb12191c8553753",
          "message": "Simplify synonyms (#355)",
          "timestamp": "2024-12-20T15:40:22+01:00",
          "tree_id": "47d026cd853cad3db215c1c4bbfcfc9df3f2daa2",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/270d3f69f79a171b9c5feec9cfb12191c8553753"
        },
        "date": 1734706177580,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6829185.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 594484185,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80101207.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 30972807,
            "unit": "ns",
            "extra": "gctime=8042264\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43752400,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13906517,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1724159307,
            "unit": "ns",
            "extra": "gctime=368893524.5\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4365680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27491260,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 220646537.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112221142,
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
          "id": "1b48458065ad4e6addd72623ac4eebc3efa57ac0",
          "message": "Remove `CairoMakie` extension in local tests (#362)\n\n* remove `CairoMakie` extension in local tests\r\n\r\n* format files\r\n\r\n* remove `CairoMakie` from `[targets]` in Project.toml`",
          "timestamp": "2025-01-02T10:32:40+01:00",
          "tree_id": "37a49ebca8ba93f777cf237e6dbd97e29675d9c3",
          "url": "https://github.com/qutip/QuantumToolbox.jl/commit/1b48458065ad4e6addd72623ac4eebc3efa57ac0"
        },
        "date": 1735810925342,
        "tool": "julia",
        "benches": [
          {
            "name": "Steadystate/Direct",
            "value": 6818341.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10284584\nallocs=463\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/dense",
            "value": 601490110.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20091440\nallocs=37\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Eigenvalues/eigenstates/sparse",
            "value": 80931704,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3660592\nallocs=369\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Pseudo Inverse",
            "value": 28816496.5,
            "unit": "ns",
            "extra": "gctime=7455877.5\nmemory=121472856\nallocs=104731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/Spectrum/Exponential Series",
            "value": 43643851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4791680\nallocs=1013\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Correlations and Spectrum/FFT Correlation",
            "value": 13885562.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=862984\nallocs=690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/Dynamical Fock Dimension",
            "value": 1679029119.5,
            "unit": "ns",
            "extra": "gctime=351177919\nmemory=3066226656\nallocs=301307\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/sesolve",
            "value": 4375712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106096\nallocs=102\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mesolve",
            "value": 27515486.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3103712\nallocs=591\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Serial",
            "value": 223611224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3563008\nallocs=12462\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Time Evolution/time-independent/mcsolve/Multithreaded",
            "value": 112515115,
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
      }
    ]
  }
}