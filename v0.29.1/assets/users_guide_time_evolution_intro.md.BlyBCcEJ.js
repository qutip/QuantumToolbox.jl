import{_ as o,c as t,o as l,ai as i}from"./chunks/framework.DCbwEeOm.js";const v=JSON.parse('{"title":"Time Evolution and Quantum System Dynamics","description":"","frontmatter":{},"headers":[],"relativePath":"users_guide/time_evolution/intro.md","filePath":"users_guide/time_evolution/intro.md","lastUpdated":null}'),a={name:"users_guide/time_evolution/intro.md"};function n(s,e,u,r,m,d){return l(),t("div",null,e[0]||(e[0]=[i('<h1 id="doc:Time-Evolution-and-Quantum-System-Dynamics" tabindex="-1">Time Evolution and Quantum System Dynamics <a class="header-anchor" href="#doc:Time-Evolution-and-Quantum-System-Dynamics" aria-label="Permalink to &quot;Time Evolution and Quantum System Dynamics {#doc:Time-Evolution-and-Quantum-System-Dynamics}&quot;">​</a></h1><p><strong>Table of contents</strong></p><ul><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/intro#doc-TE:Introduction">Introduction</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/solution#doc-TE:Time-Evolution-Solutions">Time Evolution Solutions</a></p><ul><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/solution#doc-TE:Solution">Solution</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/solution#doc-TE:Accessing-data-in-solutions">Accessing data in solutions</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/solution#doc-TE:Multiple-trajectories-solution">Multiple trajectories solution</a></p></li></ul></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/sesolve#doc-TE:Schrödinger-Equation-Solver">Schrödinger Equation Solver</a></p><ul><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/sesolve#doc-TE:Unitary-evolution">Unitary evolution</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/sesolve#doc-TE:Example:Spin-dynamics">Example: Spin dynamics</a></p></li></ul></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mesolve#doc-TE:Lindblad-Master-Equation-Solver">Lindblad Master Equation Solver</a></p><ul><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mesolve#doc-TE:Von-Neumann-equation">Von Neumann equation</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mesolve#doc-TE:The-Lindblad-master-equation">The Lindblad master equation</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mesolve#doc-TE:Example:Dissipative-Spin-dynamics">Example: Dissipative Spin dynamics</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mesolve#doc-TE:Example:Harmonic-oscillator-in-thermal-bath">Example: Harmonic oscillator in thermal bath</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mesolve#doc-TE:Example:Two-level-atom-coupled-to-dissipative-single-mode-cavity">Example: Two-level atom coupled to dissipative single-mode cavity</a></p></li></ul></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mcsolve#doc-TE:Monte-Carlo-Solver">Monte Carlo Solver</a></p><ul><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mcsolve#doc-TE:Monte-Carlo-wave-function">Monte Carlo wave-function</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mcsolve#doc-TE:Example:Two-level-atom-coupled-to-dissipative-single-mode-cavity-(MC)">Example: Two-level atom coupled to dissipative single-mode cavity (MC)</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/mcsolve#doc-TE:Running-trajectories-in-parallel">Running trajectories in parallel</a></p></li></ul></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/stochastic#doc-TE:Stochastic-Solver">Stochastic Solver</a></p><ul><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/stochastic#doc-TE:Stochastic-Schrödinger-equation">Stochastic Schrödinger equation</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/stochastic#doc-TE:Stochastic-master-equation">Stochastic master equation</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/stochastic#doc-TE:Example:Homodyne-detection">Example: Homodyne detection</a></p></li></ul></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/time_dependent#doc-TE:Solving-Problems-with-Time-dependent-Hamiltonians">Solving Problems with Time-dependent Hamiltonians</a></p><ul><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/time_dependent#doc-TE:Generate-QobjEvo">Generate QobjEvo</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/time_dependent#doc-TE:QobjEvo-fields-(attributes)">QobjEvo fields (attributes)</a></p></li><li><p><a href="/QuantumToolbox.jl/v0.29.1/users_guide/time_evolution/time_dependent#doc-TE:Using-parameters">Using parameters</a></p></li></ul></li></ul><h1 id="doc-TE:Introduction" tabindex="-1">Introduction <a class="header-anchor" href="#doc-TE:Introduction" aria-label="Permalink to &quot;Introduction {#doc-TE:Introduction}&quot;">​</a></h1><p>Although in some cases, we want to find the stationary states of a quantum system, often we are interested in the dynamics: how the state of a system or an ensemble of systems evolves with time. <code>QuantumToolbox</code> provides many ways to model dynamics.</p><p>There are two kinds of quantum systems: open systems that interact with a larger environment and closed systems that do not. In a closed system, the state can be described by a state vector. When we are modeling an open system, or an ensemble of systems, the use of the density matrix is mandatory.</p><p>The following table lists the solvers provided by <code>QuantumToolbox</code> for dynamic quantum systems and the corresponding type of solution returned by the solver:</p><table tabindex="0"><thead><tr><th style="text-align:left;"><strong>Equation</strong></th><th style="text-align:left;"><strong>Function Call</strong></th><th style="text-align:left;"><strong>Problem</strong></th><th style="text-align:left;"><strong>Returned Solution</strong></th></tr></thead><tbody><tr><td style="text-align:left;">Unitary evolution, Schrödinger equation</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.sesolve"><code>sesolve</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.sesolveProblem"><code>sesolveProblem</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.TimeEvolutionSol"><code>TimeEvolutionSol</code></a></td></tr><tr><td style="text-align:left;">Lindblad master eqn. or Von Neuman eqn.</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.mesolve"><code>mesolve</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.mesolveProblem"><code>mesolveProblem</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.TimeEvolutionSol"><code>TimeEvolutionSol</code></a></td></tr><tr><td style="text-align:left;">Monte Carlo evolution</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.mcsolve"><code>mcsolve</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.mcsolveProblem"><code>mcsolveProblem</code></a> <a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.mcsolveEnsembleProblem"><code>mcsolveEnsembleProblem</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.TimeEvolutionMCSol"><code>TimeEvolutionMCSol</code></a></td></tr><tr><td style="text-align:left;">Stochastic Schrödinger equation</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.ssesolve"><code>ssesolve</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.ssesolveProblem"><code>ssesolveProblem</code></a> <a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.ssesolveEnsembleProblem"><code>ssesolveEnsembleProblem</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.TimeEvolutionStochasticSol"><code>TimeEvolutionStochasticSol</code></a></td></tr><tr><td style="text-align:left;">Stochastic master equation</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.smesolve"><code>smesolve</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.smesolveProblem"><code>smesolveProblem</code></a> <a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.smesolveEnsembleProblem"><code>smesolveEnsembleProblem</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.29.1/resources/api#QuantumToolbox.TimeEvolutionStochasticSol"><code>TimeEvolutionStochasticSol</code></a></td></tr></tbody></table><div class="tip custom-block"><p class="custom-block-title">Solving dynamics with pre-defined problems</p><p><code>QuantumToolbox</code> provides two different methods to solve the dynamics. One can use the function calls listed above by either taking all the operators (like Hamiltonian and collapse operators, etc.) as inputs directly, or generating the <code>prob</code>lems by yourself and take it as an input of the function call, e.g., <code>sesolve(prob)</code>.</p></div>',9)]))}const p=o(a,[["render",n]]);export{v as __pageData,p as default};
