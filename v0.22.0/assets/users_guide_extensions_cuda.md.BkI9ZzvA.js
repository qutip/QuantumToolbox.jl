import{_ as a,c as e,a5 as i,o as n}from"./chunks/framework.BKgRIcLl.js";const u=JSON.parse('{"title":"Extension for CUDA.jl","description":"","frontmatter":{},"headers":[],"relativePath":"users_guide/extensions/cuda.md","filePath":"users_guide/extensions/cuda.md","lastUpdated":null}'),t={name:"users_guide/extensions/cuda.md"};function p(l,s,o,d,c,h){return n(),e("div",null,s[0]||(s[0]=[i(`<h1 id="doc:CUDA" tabindex="-1">Extension for CUDA.jl <a class="header-anchor" href="#doc:CUDA" aria-label="Permalink to &quot;Extension for CUDA.jl {#doc:CUDA}&quot;">​</a></h1><h2 id="introduction" tabindex="-1">Introduction <a class="header-anchor" href="#introduction" aria-label="Permalink to &quot;Introduction&quot;">​</a></h2><p>This is an extension to support <code>QuantumObject.data</code> conversion from standard dense and sparse CPU arrays to GPU (<a href="https://github.com/JuliaGPU/CUDA.jl" target="_blank" rel="noreferrer"><code>CUDA.jl</code></a>) arrays.</p><p>This extension will be automatically loaded if user imports both <code>QuantumToolbox.jl</code> and <a href="https://github.com/JuliaGPU/CUDA.jl" target="_blank" rel="noreferrer"><code>CUDA.jl</code></a>:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> QuantumToolbox</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CUDA</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CUDA</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CUSPARSE</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CUDA</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">allowscalar</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Avoid unexpected scalar indexing</span></span></code></pre></div><p>We wrapped several functions in <code>CUDA</code> and <code>CUDA.CUSPARSE</code> in order to not only converting <code>QuantumObject.data</code> into GPU arrays, but also changing the element type and word size (<code>32</code> and <code>64</code>) since some of the GPUs perform better in <code>32</code>-bit. The functions are listed as follows (where input <code>A</code> is a <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a>):</p><ul><li><p><code>cu(A; word_size=64)</code>: return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA</code> arrays and specified <code>word_size</code>.</p></li><li><p><code>CuArray(A)</code>: If <code>A.data</code> is a dense array, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CuArray</code>.</p></li><li><p><code>CuArray{T}(A)</code>: If <code>A.data</code> is a dense array, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CuArray</code> under element type <code>T</code>.</p></li><li><p><code>CuSparseVector(A)</code>: If <code>A.data</code> is a sparse vector, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CUSPARSE.CuSparseVector</code>.</p></li><li><p><code>CuSparseVector{T}(A)</code>: If <code>A.data</code> is a sparse vector, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CUSPARSE.CuSparseVector</code> under element type <code>T</code>.</p></li><li><p><code>CuSparseMatrixCSC(A)</code>: If <code>A.data</code> is a sparse matrix, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CUSPARSE.CuSparseMatrixCSC</code>.</p></li><li><p><code>CuSparseMatrixCSC{T}(A)</code>: If <code>A.data</code> is a sparse matrix, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CUSPARSE.CuSparseMatrixCSC</code> under element type <code>T</code>.</p></li><li><p><code>CuSparseMatrixCSR(A)</code>: If <code>A.data</code> is a sparse matrix, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CUSPARSE.CuSparseMatrixCSR</code>.</p></li><li><p><code>CuSparseMatrixCSR{T}(A)</code>: If <code>A.data</code> is a sparse matrix, return a new <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a> with <code>CUDA.CUSPARSE.CuSparseMatrixCSR</code> under element type <code>T</code>.</p></li></ul><p>We suggest to convert the arrays from CPU to GPU memory by using the function <code>cu</code> because it allows different <code>data</code>-types of input <a href="/QuantumToolbox.jl/v0.22.0/resources/api#QuantumToolbox.QuantumObject"><code>QuantumObject</code></a>.</p><p>Here are some examples:</p><h2 id="Converting-dense-arrays" tabindex="-1">Converting dense arrays <a class="header-anchor" href="#Converting-dense-arrays" aria-label="Permalink to &quot;Converting dense arrays {#Converting-dense-arrays}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">V </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> fock</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># CPU dense vector</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Ket   dims=[2]   size=(2,)</span></span>
<span class="line"><span>2-element Vector{ComplexF64}:</span></span>
<span class="line"><span> 1.0 + 0.0im</span></span>
<span class="line"><span> 0.0 + 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(V)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Ket   dims=[2]   size=(2,)</span></span>
<span class="line"><span>2-element CuArray{ComplexF64, 1, CUDA.DeviceMemory}:</span></span>
<span class="line"><span> 1.0 + 0.0im</span></span>
<span class="line"><span> 0.0 + 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(V; word_size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 32</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Ket   dims=[2]   size=(2,)</span></span>
<span class="line"><span>2-element CuArray{ComplexF32, 1, CUDA.DeviceMemory}:</span></span>
<span class="line"><span> 1.0 + 0.0im</span></span>
<span class="line"><span> 0.0 + 0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">M </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Qobj</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># CPU dense matrix</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=false</span></span>
<span class="line"><span>2×2 Matrix{Int64}:</span></span>
<span class="line"><span> 1  2</span></span>
<span class="line"><span> 3  4</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(M)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=false</span></span>
<span class="line"><span>2×2 CuArray{Int64, 2, CUDA.DeviceMemory}:</span></span>
<span class="line"><span> 1  2</span></span>
<span class="line"><span> 3  4</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(M; word_size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 32</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=false</span></span>
<span class="line"><span>2×2 CuArray{Int32, 2, CUDA.DeviceMemory}:</span></span>
<span class="line"><span> 1  2</span></span>
<span class="line"><span> 3  4</span></span></code></pre></div><h2 id="Converting-sparse-arrays" tabindex="-1">Converting sparse arrays <a class="header-anchor" href="#Converting-sparse-arrays" aria-label="Permalink to &quot;Converting sparse arrays {#Converting-sparse-arrays}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">V </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> fock</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; sparse</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># CPU sparse vector</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Ket   dims=[2]   size=(2,)</span></span>
<span class="line"><span>2-element SparseVector{ComplexF64, Int64} with 1 stored entry:</span></span>
<span class="line"><span>  [1]  =  1.0+0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(V)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Ket   dims=[2]   size=(2,)</span></span>
<span class="line"><span>2-element CuSparseVector{ComplexF64, Int32} with 1 stored entry:</span></span>
<span class="line"><span>  [1]  =  1.0+0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(V; word_size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 32</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Ket   dims=[2]   size=(2,)</span></span>
<span class="line"><span>2-element CuSparseVector{ComplexF32, Int32} with 1 stored entry:</span></span>
<span class="line"><span>  [1]  =  1.0+0.0im</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">M </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sigmax</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># CPU sparse matrix</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true</span></span>
<span class="line"><span>2×2 SparseMatrixCSC{ComplexF64, Int64} with 2 stored entries:</span></span>
<span class="line"><span>     ⋅      1.0+0.0im</span></span>
<span class="line"><span> 1.0+0.0im      ⋅</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(M)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true</span></span>
<span class="line"><span>2×2 CuSparseMatrixCSC{ComplexF64, Int32} with 2 stored entries:</span></span>
<span class="line"><span>     ⋅      1.0+0.0im</span></span>
<span class="line"><span> 1.0+0.0im      ⋅</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cu</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(M; word_size </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 32</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true</span></span>
<span class="line"><span>2×2 CuSparseMatrixCSC{ComplexF32, Int32} with 2 stored entries:</span></span>
<span class="line"><span>     ⋅      1.0+0.0im</span></span>
<span class="line"><span> 1.0+0.0im      ⋅</span></span></code></pre></div>`,35)]))}const k=a(t,[["render",p]]);export{u as __pageData,k as default};
