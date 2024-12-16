import{_ as t,c as a,a5 as i,o}from"./chunks/framework.CdqZ1XjE.js";const h=JSON.parse('{"title":"Extension for CairoMakie.jl","description":"","frontmatter":{},"headers":[],"relativePath":"users_guide/extensions/cairomakie.md","filePath":"users_guide/extensions/cairomakie.md","lastUpdated":null}'),r={name:"users_guide/extensions/cairomakie.md"};function s(n,e,l,d,p,c){return o(),a("div",null,e[0]||(e[0]=[i(`<h1 id="doc:CairoMakie" tabindex="-1">Extension for CairoMakie.jl <a class="header-anchor" href="#doc:CairoMakie" aria-label="Permalink to &quot;Extension for CairoMakie.jl {#doc:CairoMakie}&quot;">​</a></h1><p>This is an extension to support visualization (plotting functions) using <a href="https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie" target="_blank" rel="noreferrer"><code>CairoMakie.jl</code></a> library.</p><p>This extension will be automatically loaded if user imports both <code>QuantumToolbox.jl</code> and <a href="https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie" target="_blank" rel="noreferrer"><code>CairoMakie.jl</code></a>:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> QuantumToolbox</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CairoMakie</span></span></code></pre></div><p>To plot with <a href="https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie" target="_blank" rel="noreferrer"><code>CairoMakie.jl</code></a> library, specify the keyword argument <code>library = Val(:CairoMakie)</code> for the plotting functions.</p><div class="warning custom-block"><p class="custom-block-title">Beware of type-stability!</p><p>If you want to keep type stability, it is recommended to use <code>Val(:CairoMakie)</code> instead of <code>:CairoMakie</code>. See <a href="https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type" target="_blank" rel="noreferrer">this link</a> and the <a href="/QuantumToolbox.jl/v0.24.0/getting_started/type_stability#doc:Type-Stability">related Section</a> about type stability for more details.</p></div><p>The supported plotting functions are listed as follows:</p><table tabindex="0"><thead><tr><th style="text-align:left;"><strong>Plotting Function</strong></th><th style="text-align:left;"><strong>Description</strong></th></tr></thead><tbody><tr><td style="text-align:left;"><a href="/QuantumToolbox.jl/v0.24.0/resources/api#QuantumToolbox.plot_wigner"><code>plot_wigner</code></a></td><td style="text-align:left;"><a href="https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution" target="_blank" rel="noreferrer">Wigner quasipropability distribution</a></td></tr></tbody></table>`,8)]))}const k=t(r,[["render",s]]);export{h as __pageData,k as default};
