import{_ as a,c as e,a5 as t,o as n}from"./chunks/framework.auZU70jD.js";const h=JSON.parse('{"title":"Functions operating on Qobj","description":"","frontmatter":{},"headers":[],"relativePath":"users_guide/QuantumObject/QuantumObject_functions.md","filePath":"users_guide/QuantumObject/QuantumObject_functions.md","lastUpdated":null}'),i={name:"users_guide/QuantumObject/QuantumObject_functions.md"};function l(p,s,o,d,r,c){return n(),e("div",null,s[0]||(s[0]=[t(`<h1 id="doc:Functions-operating-on-Qobj" tabindex="-1">Functions operating on Qobj <a class="header-anchor" href="#doc:Functions-operating-on-Qobj" aria-label="Permalink to &quot;Functions operating on Qobj {#doc:Functions-operating-on-Qobj}&quot;">​</a></h1><p><code>QuantumToolbox</code> also provide functions (methods) that operates on <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a>.</p><p>You can click the function links and see the corresponding docstring for more information.</p><h2 id="Linear-algebra-and-attributes" tabindex="-1">Linear algebra and attributes <a class="header-anchor" href="#Linear-algebra-and-attributes" aria-label="Permalink to &quot;Linear algebra and attributes {#Linear-algebra-and-attributes}&quot;">​</a></h2><p>Here is a table that summarizes all the supported linear algebra functions and attribute functions operating on a given <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> <code>Q</code>:</p><table tabindex="0"><thead><tr><th style="text-align:left;"><strong>Description</strong></th><th style="text-align:left;"><strong>Function call</strong></th><th style="text-align:left;"><strong>Synonyms</strong></th></tr></thead><tbody><tr><td style="text-align:left;">zero-like array</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.zero"><code>zero(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">identity-like matrix</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.one"><code>one(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">conjugate</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.conj"><code>conj(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">transpose</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.transpose"><code>transpose(Q)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.trans"><code>trans(Q)</code></a></td></tr><tr><td style="text-align:left;">conjugate transposition</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.adjoint"><code>adjoint(Q)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.adjoint"><code>Q&#39;</code></a>, <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.dag"><code>dag(Q)</code></a></td></tr><tr><td style="text-align:left;">partial transpose</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.partial_transpose"><code>partial_transpose(Q, mask)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">dot product</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.dot"><code>dot(Q1, Q2)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">generalized dot product</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.dot"><code>dot(Q1, Q2, Q3)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.matrix_element"><code>matrix_element(Q1, Q2, Q3)</code></a></td></tr><tr><td style="text-align:left;">trace</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.tr"><code>tr(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">partial trace</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.ptrace"><code>ptrace(Q, sel)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">singular values</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.svdvals"><code>svdvals(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">standard vector <code>p</code>-norm or <a href="https://en.wikipedia.org/wiki/Schatten_norm" target="_blank" rel="noreferrer">Schatten</a> <code>p</code>-norm</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.norm"><code>norm(Q, p)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">normalization</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.normalize"><code>normalize(Q, p)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.unit"><code>unit(Q, p)</code></a></td></tr><tr><td style="text-align:left;">normalization (in-place)</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.normalize!"><code>normalize!(Q, p)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">matrix inverse</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.inv"><code>inv(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">matrix square root</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.sqrt"><code>sqrt(Q)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.sqrt"><code>√(Q)</code></a>, <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.sqrtm"><code>sqrtm(Q)</code></a></td></tr><tr><td style="text-align:left;">matrix logarithm</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.log"><code>log(Q)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.logm"><code>logm(Q)</code></a></td></tr><tr><td style="text-align:left;">matrix exponential</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.exp"><code>exp(Q)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.expm"><code>expm(Q)</code></a></td></tr><tr><td style="text-align:left;">matrix sine</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.sin"><code>sin(Q)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.sinm"><code>sinm(Q)</code></a></td></tr><tr><td style="text-align:left;">matrix cosine</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#Base.cos"><code>cos(Q)</code></a></td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.cosm"><code>cosm(Q)</code></a></td></tr><tr><td style="text-align:left;">diagonal elements</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.diag"><code>diag(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">projector</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.proj"><code>proj(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">purity</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.purity"><code>purity(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">permute</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.permute"><code>permute(Q, order)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">remove small elements</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.tidyup"><code>tidyup(Q, tol)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">remove small elements (in-place)</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.tidyup!"><code>tidyup!(Q, tol)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">get data</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.get_data"><code>get_data(Q)</code></a></td><td style="text-align:left;">-</td></tr><tr><td style="text-align:left;">get coherence</td><td style="text-align:left;"><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.get_coherence"><code>get_coherence(Q)</code></a></td><td style="text-align:left;">-</td></tr></tbody></table><h2 id="Eigenvalue-decomposition" tabindex="-1">Eigenvalue decomposition <a class="header-anchor" href="#Eigenvalue-decomposition" aria-label="Permalink to &quot;Eigenvalue decomposition {#Eigenvalue-decomposition}&quot;">​</a></h2><ul><li><p><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.eigenenergies"><code>eigenenergies</code></a>: return eigenenergies (eigenvalues)</p></li><li><p><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.eigenstates"><code>eigenstates</code></a>: return <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.EigsolveResult"><code>EigsolveResult</code></a> (contains eigenvalues and eigenvectors)</p></li><li><p><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.eigvals"><code>eigvals</code></a>: return eigenvalues</p></li><li><p><a href="/QuantumToolbox.jl/previews/PR364/resources/api#LinearAlgebra.eigen"><code>eigen</code></a>: using dense eigen solver and return <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.EigsolveResult"><code>EigsolveResult</code></a> (contains eigenvalues and eigenvectors)</p></li><li><p><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.eigsolve"><code>eigsolve</code></a>: using sparse eigen solver and return <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.EigsolveResult"><code>EigsolveResult</code></a> (contains eigenvalues and eigenvectors)</p></li><li><p><a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.eigsolve_al"><code>eigsolve_al</code></a>: using the Arnoldi-Lindblad eigen solver and return <a href="/QuantumToolbox.jl/previews/PR364/resources/api#QuantumToolbox.EigsolveResult"><code>EigsolveResult</code></a> (contains eigenvalues and eigenvectors)</p></li></ul><h2 id="examples" tabindex="-1">Examples <a class="header-anchor" href="#examples" aria-label="Permalink to &quot;Examples&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ψ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> normalize</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">basis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> basis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Ket   dims=[4]   size=(4,)</span></span>
<span class="line"><span>4-element Vector{ComplexF64}:</span></span>
<span class="line"><span>                0.0 + 0.0im</span></span>
<span class="line"><span> 0.7071067811865475 + 0.0im</span></span>
<span class="line"><span> 0.7071067811865475 + 0.0im</span></span>
<span class="line"><span>                0.0 + 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ψ</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&#39;</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Bra   dims=[4]   size=(1, 4)</span></span>
<span class="line"><span>1×4 adjoint(::Vector{ComplexF64}) with eltype ComplexF64:</span></span>
<span class="line"><span> 0.0-0.0im  0.707107-0.0im  0.707107-0.0im  0.0-0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ρ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> coherent_dm</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[5]   size=(5, 5)   ishermitian=true</span></span>
<span class="line"><span>5×5 Matrix{ComplexF64}:</span></span>
<span class="line"><span> 0.367911+0.0im   0.367744+0.0im  …  0.146207+0.0im   0.088267+0.0im</span></span>
<span class="line"><span> 0.367744+0.0im   0.367577+0.0im      0.14614+0.0im  0.0882269+0.0im</span></span>
<span class="line"><span> 0.261054+0.0im   0.260936+0.0im     0.103742+0.0im  0.0626306+0.0im</span></span>
<span class="line"><span> 0.146207+0.0im    0.14614+0.0im     0.058102+0.0im   0.035077+0.0im</span></span>
<span class="line"><span> 0.088267+0.0im  0.0882269+0.0im     0.035077+0.0im  0.0211765+0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">diag</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ρ)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span>   0.3679111729923387 + 0.0im</span></span>
<span class="line"><span>   0.3675770456232403 + 0.0im</span></span>
<span class="line"><span>  0.18523331233838003 + 0.0im</span></span>
<span class="line"><span>  0.05810197190350208 + 0.0im</span></span>
<span class="line"><span> 0.021176497142538265 + 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_data</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ρ)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>5×5 Matrix{ComplexF64}:</span></span>
<span class="line"><span> 0.367911+0.0im   0.367744+0.0im  …  0.146207+0.0im   0.088267+0.0im</span></span>
<span class="line"><span> 0.367744+0.0im   0.367577+0.0im      0.14614+0.0im  0.0882269+0.0im</span></span>
<span class="line"><span> 0.261054+0.0im   0.260936+0.0im     0.103742+0.0im  0.0626306+0.0im</span></span>
<span class="line"><span> 0.146207+0.0im    0.14614+0.0im     0.058102+0.0im   0.035077+0.0im</span></span>
<span class="line"><span> 0.088267+0.0im  0.0882269+0.0im     0.035077+0.0im  0.0211765+0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">norm</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ρ)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>0.9999999999999996</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sqrtm</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ρ)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[5]   size=(5, 5)   ishermitian=true</span></span>
<span class="line"><span>5×5 Matrix{ComplexF64}:</span></span>
<span class="line"><span> 0.367911+0.0im   0.367744+0.0im  …  0.146207+0.0im   0.088267+0.0im</span></span>
<span class="line"><span> 0.367744-0.0im   0.367577+0.0im      0.14614+0.0im  0.0882269+0.0im</span></span>
<span class="line"><span> 0.261054-0.0im   0.260936-0.0im     0.103742+0.0im  0.0626306+0.0im</span></span>
<span class="line"><span> 0.146207-0.0im    0.14614-0.0im     0.058102+0.0im   0.035077+0.0im</span></span>
<span class="line"><span> 0.088267-0.0im  0.0882269-0.0im     0.035077-0.0im  0.0211765+0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">tr</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ρ)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>0.9999999999999993 + 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">eigenenergies</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ρ)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>5-element Vector{Float64}:</span></span>
<span class="line"><span> -3.188562758433179e-17</span></span>
<span class="line"><span> -7.70881124130806e-18</span></span>
<span class="line"><span>  2.0684351320241445e-17</span></span>
<span class="line"><span>  1.0235138990670176e-16</span></span>
<span class="line"><span>  0.9999999999999996</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> eigenstates</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ρ)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>EigsolveResult:   type=Operator   dims=[5]</span></span>
<span class="line"><span>values:</span></span>
<span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span> -2.8177392874225097e-17 + 0.0im</span></span>
<span class="line"><span>   3.416070845000482e-17 + 0.0im</span></span>
<span class="line"><span>   7.728341127110703e-17 + 0.0im</span></span>
<span class="line"><span>   8.881784197001252e-16 + 0.0im</span></span>
<span class="line"><span>      0.9999999999999994 + 0.0im</span></span>
<span class="line"><span>vectors:</span></span>
<span class="line"><span>5×5 Matrix{ComplexF64}:</span></span>
<span class="line"><span>  0.349081+0.0im   0.699533+0.0im  …  -0.0892167+0.0im  -0.606557+0.0im</span></span>
<span class="line"><span> -0.698381+0.0im  -0.132457+0.0im     -0.0891762+0.0im  -0.606281+0.0im</span></span>
<span class="line"><span>  0.596321+0.0im  -0.648729+0.0im     -0.0633045+0.0im  -0.430387+0.0im</span></span>
<span class="line"><span> -0.186568+0.0im  -0.268811+0.0im     -0.0354544+0.0im  -0.241044+0.0im</span></span>
<span class="line"><span>       0.0+0.0im        0.0+0.0im       0.989355+0.0im  -0.145521-0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">λ, ψ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> result</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">λ </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># eigenvalues</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span> -2.8177392874225097e-17 + 0.0im</span></span>
<span class="line"><span>   3.416070845000482e-17 + 0.0im</span></span>
<span class="line"><span>   7.728341127110703e-17 + 0.0im</span></span>
<span class="line"><span>   8.881784197001252e-16 + 0.0im</span></span>
<span class="line"><span>      0.9999999999999994 + 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ψ </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># eigenvectors</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>5-element Vector{QuantumObject{Vector{ComplexF64}, KetQuantumObject, 1}}:</span></span>
<span class="line"><span> Quantum Object:   type=Ket   dims=[5]   size=(5,)</span></span>
<span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span>  0.34908090726365426 + 0.0im</span></span>
<span class="line"><span>  -0.6983812275352181 + 0.0im</span></span>
<span class="line"><span>   0.5963209517322131 + 0.0im</span></span>
<span class="line"><span> -0.18656769210014237 + 0.0im</span></span>
<span class="line"><span>                  0.0 + 0.0im</span></span>
<span class="line"><span> Quantum Object:   type=Ket   dims=[5]   size=(5,)</span></span>
<span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span>  0.6995328227558701 + 0.0im</span></span>
<span class="line"><span> -0.1324572834510804 + 0.0im</span></span>
<span class="line"><span> -0.6487290623177562 + 0.0im</span></span>
<span class="line"><span> -0.2688112751584079 + 0.0im</span></span>
<span class="line"><span>                 0.0 + 0.0im</span></span>
<span class="line"><span> Quantum Object:   type=Ket   dims=[5]   size=(5,)</span></span>
<span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span>  0.11369058957774669 + 0.0im</span></span>
<span class="line"><span> -0.34523801249281383 + 0.0im</span></span>
<span class="line"><span> -0.18523271117272722 + 0.0im</span></span>
<span class="line"><span>   0.9130027422100534 + 0.0im</span></span>
<span class="line"><span>                  0.0 + 0.0im</span></span>
<span class="line"><span> Quantum Object:   type=Ket   dims=[5]   size=(5,)</span></span>
<span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span>  -0.08921674125585327 + 0.0im</span></span>
<span class="line"><span>   -0.0891762198903574 + 0.0im</span></span>
<span class="line"><span>  -0.06330447537681073 + 0.0im</span></span>
<span class="line"><span> -0.035454413343867425 + 0.0im</span></span>
<span class="line"><span>    0.9893550944213416 + 0.0im</span></span>
<span class="line"><span> Quantum Object:   type=Ket   dims=[5]   size=(5,)</span></span>
<span class="line"><span>5-element Vector{ComplexF64}:</span></span>
<span class="line"><span>  -0.6065568176126116 + 0.0im</span></span>
<span class="line"><span>   -0.606281325477901 + 0.0im</span></span>
<span class="line"><span> -0.43038739797812414 + 0.0im</span></span>
<span class="line"><span> -0.24104350624628368 + 0.0im</span></span>
<span class="line"><span> -0.14552146626026788 - 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">λ, ψ, T </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> result</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">T </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># transformation matrix</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>5×5 Matrix{ComplexF64}:</span></span>
<span class="line"><span>  0.349081+0.0im   0.699533+0.0im  …  -0.0892167+0.0im  -0.606557+0.0im</span></span>
<span class="line"><span> -0.698381+0.0im  -0.132457+0.0im     -0.0891762+0.0im  -0.606281+0.0im</span></span>
<span class="line"><span>  0.596321+0.0im  -0.648729+0.0im     -0.0633045+0.0im  -0.430387+0.0im</span></span>
<span class="line"><span> -0.186568+0.0im  -0.268811+0.0im     -0.0354544+0.0im  -0.241044+0.0im</span></span>
<span class="line"><span>       0.0+0.0im        0.0+0.0im       0.989355+0.0im  -0.145521-0.0im</span></span></code></pre></div>`,35)]))}const m=a(i,[["render",l]]);export{h as __pageData,m as default};
