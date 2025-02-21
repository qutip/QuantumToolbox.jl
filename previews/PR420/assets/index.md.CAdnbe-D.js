import{_ as e,c as a,o as i,ai as s}from"./chunks/framework.CmV4vDCh.js";const c=JSON.parse('{"title":"Introduction","description":"","frontmatter":{"layout":"home","hero":{"name":"QuantumToolbox.jl","tagline":"A pure Julia framework designed for High-performance quantum physics simulations","image":{"src":"/logo.png","alt":"QuantumToolbox"},"actions":[{"theme":"brand","text":"Getting Started","link":"/getting_started/brief_example"},{"theme":"alt","text":"Users Guide","link":"/users_guide/QuantumObject/QuantumObject"},{"theme":"alt","text":"Tutorials","link":"https://qutip.org/qutip-julia-tutorials/"},{"theme":"alt","text":"API","link":"/resources/api"},{"theme":"alt","text":"View on Github","link":"https://github.com/qutip/QuantumToolbox.jl"},{"theme":"alt","text":"Visit QuTiP.org","link":"https://qutip.org/"}]},"features":[{"icon":"<img width=\\"64\\" height=\\"64\\" src=\\"https://docs.sciml.ai/DiffEqDocs/stable/assets/logo.png\\" alt=\\"markdown\\"/>","title":"Dynamical Evolution","details":"Advanced solvers for time evolution of quantum systems, thanks to the powerful DifferentialEquations.jl package.","link":"/users_guide/time_evolution/intro"},{"icon":"<img width=\\"64\\" height=\\"64\\" src=\\"https://cuda.juliagpu.org/stable/assets/logo.png\\" />","title":"GPU Computing","details":"Leverage GPU resources for high-performance computing. Simulate quantum dynamics directly on the GPU with the same syntax as the CPU case.","link":"/users_guide/extensions/cuda"},{"icon":"<img width=\\"64\\" height=\\"64\\" src=\\"https://img.icons8.com/?size=100&id=1W4Bkj363ov0&format=png&color=000000\\" />","title":"Distributed Computing","details":"Distribute the computation over multiple nodes (e.g., a cluster). Simulate hundreds of quantum trajectories in parallel on a cluster, with, again, the same syntax as the simple case.","link":"/users_guide/cluster"}]},"headers":[],"relativePath":"index.md","filePath":"index.md","lastUpdated":null}'),n={name:"index.md"};function o(l,t,r,u,h,p){return i(),a("div",null,t[0]||(t[0]=[s(`<h1 id="doc:Introduction" tabindex="-1">Introduction <a class="header-anchor" href="#doc:Introduction" aria-label="Permalink to &quot;Introduction {#doc:Introduction}&quot;">​</a></h1><p><a href="https://github.com/qutip/QuantumToolbox.jl" target="_blank" rel="noreferrer"><code>QuantumToolbox.jl</code></a> is a cutting-edge <a href="https://julialang.org/" target="_blank" rel="noreferrer"><code>Julia</code></a> package designed for quantum physics simulations, closely emulating the popular <a href="https://github.com/qutip/qutip" target="_blank" rel="noreferrer"><code>Python QuTiP</code></a> package. It uniquely combines the simplicity and power of Julia with advanced features like GPU acceleration and distributed computing, making simulation of quantum systems more accessible and efficient. Taking advantage of the <a href="https://julialang.org/" target="_blank" rel="noreferrer"><code>Julia</code></a> language features (like multiple dispatch and metaprogramming), <a href="https://github.com/qutip/QuantumToolbox.jl" target="_blank" rel="noreferrer"><code>QuantumToolbox.jl</code></a> is designed to be easily extendable, allowing users to build upon the existing functionalities.</p><p><em><strong>With this package, moving from Python to Julia for quantum physics simulations has never been easier</strong></em>, due to the similar syntax and functionalities.</p><h1 id="doc:Installation" tabindex="-1">Installation <a class="header-anchor" href="#doc:Installation" aria-label="Permalink to &quot;Installation {#doc:Installation}&quot;">​</a></h1><div class="tip custom-block"><p class="custom-block-title">Requirements</p><p><code>QuantumToolbox.jl</code> requires <code>Julia 1.10+</code>.</p></div><p>To install <code>QuantumToolbox.jl</code>, run the following commands inside Julia&#39;s interactive session (also known as REPL):</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;QuantumToolbox&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Alternatively, this can also be done in <code>Julia</code>&#39;s <a href="https://julialang.github.io/Pkg.jl/v1/getting-started/" target="_blank" rel="noreferrer">Pkg REPL</a> by pressing the key <code>]</code> in the REPL to use the package mode, and then type the following command:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> add QuantumToolbox</span></span></code></pre></div><p>More information about <code>Julia</code>&#39;s package manager can be found at <a href="https://julialang.github.io/Pkg.jl/v1/" target="_blank" rel="noreferrer"><code>Pkg.jl</code></a>.</p><p>To load the package and check the version information, use either <a href="/QuantumToolbox.jl/previews/PR420/resources/api#QuantumToolbox.versioninfo"><code>QuantumToolbox.versioninfo()</code></a> or <a href="/QuantumToolbox.jl/previews/PR420/resources/api#QuantumToolbox.about"><code>QuantumToolbox.about()</code></a>, namely</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> QuantumToolbox</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">QuantumToolbox</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">versioninfo</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">QuantumToolbox</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">about</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span></code></pre></div>`,12)]))}const g=e(n,[["render",o]]);export{c as __pageData,g as default};
