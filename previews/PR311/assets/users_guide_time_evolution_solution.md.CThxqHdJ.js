import{_ as i,o as a,c as t,a5 as e}from"./chunks/framework.8NlStuqe.js";const n="/QuantumToolbox.jl/previews/PR311/assets/qwhtrdd.Y3DCNuer.svg",g=JSON.parse('{"title":"Time Evolution Solutions","description":"","frontmatter":{},"headers":[],"relativePath":"users_guide/time_evolution/solution.md","filePath":"users_guide/time_evolution/solution.md","lastUpdated":null}'),l={name:"users_guide/time_evolution/solution.md"};function h(p,s,k,o,d,r){return a(),t("div",null,s[0]||(s[0]=[e(`<h1 id="doc-TE:Time-Evolution-Solutions" tabindex="-1">Time Evolution Solutions <a class="header-anchor" href="#doc-TE:Time-Evolution-Solutions" aria-label="Permalink to &quot;Time Evolution Solutions {#doc-TE:Time-Evolution-Solutions}&quot;">​</a></h1><h2 id="doc-TE:Solution" tabindex="-1">Solution <a class="header-anchor" href="#doc-TE:Solution" aria-label="Permalink to &quot;Solution {#doc-TE:Solution}&quot;">​</a></h2><p><code>QuantumToolbox</code> utilizes the powerful <a href="https://docs.sciml.ai/DiffEqDocs/stable/" target="_blank" rel="noreferrer"><code>DifferentialEquation.jl</code></a> to simulate different kinds of quantum system dynamics. Thus, we will first look at the data structure used for returning the solution (<code>sol</code>) from <a href="https://docs.sciml.ai/DiffEqDocs/stable/" target="_blank" rel="noreferrer"><code>DifferentialEquation.jl</code></a>. The solution stores all the crucial data needed for analyzing and plotting the results of a simulation. A generic structure <a href="/QuantumToolbox.jl/previews/PR311/resources/api#QuantumToolbox.TimeEvolutionSol"><code>TimeEvolutionSol</code></a> contains the following properties for storing simulation data:</p><table tabindex="0"><thead><tr><th style="text-align:left;"><strong>Fields (Attributes)</strong></th><th style="text-align:left;"><strong>Description</strong></th></tr></thead><tbody><tr><td style="text-align:left;"><code>sol.times</code></td><td style="text-align:left;">The time list of the evolution.</td></tr><tr><td style="text-align:left;"><code>sol.states</code></td><td style="text-align:left;">The list of result states.</td></tr><tr><td style="text-align:left;"><code>sol.expect</code></td><td style="text-align:left;">The expectation values corresponding to each time point in <code>sol.times</code>.</td></tr><tr><td style="text-align:left;"><code>sol.alg</code></td><td style="text-align:left;">The algorithm which is used during the solving process.</td></tr><tr><td style="text-align:left;"><code>sol.abstol</code></td><td style="text-align:left;">The absolute tolerance which is used during the solving process.</td></tr><tr><td style="text-align:left;"><code>sol.reltol</code></td><td style="text-align:left;">The relative tolerance which is used during the solving process.</td></tr><tr><td style="text-align:left;"><code>sol.retcode</code> (or <code>sol.converged</code>)</td><td style="text-align:left;">The returned status from the solver.</td></tr></tbody></table><h2 id="doc-TE:Accessing-data-in-solutions" tabindex="-1">Accessing data in solutions <a class="header-anchor" href="#doc-TE:Accessing-data-in-solutions" aria-label="Permalink to &quot;Accessing data in solutions {#doc-TE:Accessing-data-in-solutions}&quot;">​</a></h2><p>To understand how to access the data in solution, we will use an example as a guide, although we do not worry about the simulation details at this stage. The Schrödinger equation solver (<a href="/QuantumToolbox.jl/previews/PR311/resources/api#QuantumToolbox.sesolve"><code>sesolve</code></a>) used in this example returns <a href="/QuantumToolbox.jl/previews/PR311/resources/api#QuantumToolbox.TimeEvolutionSol"><code>TimeEvolutionSol</code></a>:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">H </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.5</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> *</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sigmay</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ψ0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> basis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">e_ops </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    proj</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">basis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)),</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    proj</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">basis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)),</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    basis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> basis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&#39;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">tlist </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> LinRange</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sol </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sesolve</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(H, ψ0, tlist, e_ops </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> e_ops, progress_bar </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Val</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><p>To see what is contained inside the solution, we can use the <code>print</code> function:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sol)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Solution of time evolution</span></span>
<span class="line"><span>(return code: Success)</span></span>
<span class="line"><span>--------------------------</span></span>
<span class="line"><span>num_states = 1</span></span>
<span class="line"><span>num_expect = 3</span></span>
<span class="line"><span>ODE alg.: OrdinaryDiffEqTsit5.Tsit5{typeof(OrdinaryDiffEqCore.trivial_limiter!), typeof(OrdinaryDiffEqCore.trivial_limiter!), Static.False}(OrdinaryDiffEqCore.trivial_limiter!, OrdinaryDiffEqCore.trivial_limiter!, static(false))</span></span>
<span class="line"><span>abstol = 1.0e-8</span></span>
<span class="line"><span>reltol = 1.0e-6</span></span></code></pre></div><p>It tells us the number of expectation values are computed and the number of states are stored. Now we have all the information needed to analyze the simulation results. To access the data for the three expectation values, one can do:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">expt1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> real</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sol</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">expect[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,:])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">expt2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> real</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sol</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">expect[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,:])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">expt3 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> real</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sol</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">expect[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,:])</span></span></code></pre></div><p>Recall that <code>Julia</code> uses <code>Fortran</code>-style indexing that begins with one (i.e., <code>[1,:]</code> represents the 1-st observable, where <code>:</code> represents all values corresponding to <code>tlist</code>).</p><p>Together with the array of times at which these expectation values are calculated:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">times </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> sol</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">times</span></span></code></pre></div><p>we can plot the resulting expectation values:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CairoMakie</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CairoMakie</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">enable_only_mime!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">MIME</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;image/svg+xml&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">500</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">350</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ax </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], xlabel </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;t&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, times, expt1, label </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\l</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">angle 0 | </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\r</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ho(t) | 0 </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\r</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">angle&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, times, expt2, label </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\l</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">angle 1 | </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\r</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ho(t) | 1 </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\r</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">angle&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">lines!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, times, expt3, label </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> L</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\l</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">angle 0 | </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\r</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ho(t) | 1 </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\r</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">angle&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ylims!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, (</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">axislegend</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ax, position </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :lb</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig</span></span></code></pre></div><p><img src="`+n+`" alt=""></p><p>State vectors, or density matrices, are accessed in a similar manner:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sol</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">states</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>1-element Vector{QuantumObject{Vector{ComplexF64}, KetQuantumObject, 1}}:</span></span>
<span class="line"><span> Quantum Object:   type=Ket   dims=[2]   size=(2,)</span></span>
<span class="line"><span>2-element Vector{ComplexF64}:</span></span>
<span class="line"><span> 0.28366218546240907 + 0.0im</span></span>
<span class="line"><span> -0.9589242745990754 + 0.0im</span></span></code></pre></div><p>Here, the solution contains only one (final) state. Because the <code>states</code> will be saved depend on the keyword argument <code>saveat</code> in <code>kwargs</code>. If <code>e_ops</code> is empty, the default value of <code>saveat=tlist</code> (saving the states corresponding to <code>tlist</code>), otherwise, <code>saveat=[tlist[end]]</code> (only save the final state). One can also specify <code>e_ops</code> and <code>saveat</code> separately.</p><p>Some other solvers can have other output.</p><h2 id="doc-TE:Multiple-trajectories-solution" tabindex="-1">Multiple trajectories solution <a class="header-anchor" href="#doc-TE:Multiple-trajectories-solution" aria-label="Permalink to &quot;Multiple trajectories solution {#doc-TE:Multiple-trajectories-solution}&quot;">​</a></h2><p>This part is still under construction, please visit <a href="/QuantumToolbox.jl/previews/PR311/resources/api#doc-API">API</a> first.</p>`,25)]))}const c=i(l,[["render",h]]);export{g as __pageData,c as default};
