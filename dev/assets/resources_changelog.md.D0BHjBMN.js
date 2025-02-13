import{_ as t,c as o,o as a,ai as r}from"./chunks/framework.CufjCsJt.js";const b=JSON.parse('{"title":"ChangeLog","description":"","frontmatter":{},"headers":[],"relativePath":"resources/changelog.md","filePath":"resources/changelog.md","lastUpdated":null}'),l={name:"resources/changelog.md"};function i(s,e,u,n,h,c){return a(),o("div",null,e[0]||(e[0]=[r('<h1 id="changelog" tabindex="-1">ChangeLog <a class="header-anchor" href="#changelog" aria-label="Permalink to &quot;ChangeLog&quot;">​</a></h1><p>All notable changes to <a href="https://github.com/qutip/QuantumToolbox.jl" target="_blank" rel="noreferrer"><code>QuantumToolbox.jl</code></a> will be documented in this file.</p><p>The format is based on <a href="https://keepachangelog.com/en/1.0.0/" target="_blank" rel="noreferrer">Keep a Changelog</a>, and this project adheres to <a href="https://semver.org/spec/v2.0.0.html" target="_blank" rel="noreferrer">Semantic Versioning</a>.</p><h2 id="unreleased" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/tree/main" target="_blank" rel="noreferrer">Unreleased</a> <a class="header-anchor" href="#unreleased" aria-label="Permalink to &quot;[Unreleased](https://github.com/qutip/QuantumToolbox.jl/tree/main)&quot;">​</a></h2><ul><li><p>Rename <code>sparse_to_dense</code> as <code>to_dense</code> and <code>dense_to_sparse</code> as <code>to_sparse</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/392" target="_blank" rel="noreferrer">#392</a>)</p></li><li><p>Fix erroneous definition of the stochastic term in <code>smesolve</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/393" target="_blank" rel="noreferrer">#393</a>)</p></li><li><p>Change name of <code>MultiSiteOperator</code> to <code>multisite_operator</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/394" target="_blank" rel="noreferrer">#394</a>)</p></li><li><p>Fix <code>smesolve</code> for specifying initial state as density matrix. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/395" target="_blank" rel="noreferrer">#395</a>)</p></li><li><p>Add more generic solver for <code>steadystate_floquet</code> to allow more linear solvers. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/396" target="_blank" rel="noreferrer">#396</a>)</p></li><li><p>Fix time evolution output when using <code>saveat</code> keyword argument. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/398" target="_blank" rel="noreferrer">#398</a>)</p></li></ul><h2 id="v0-26-0" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.26.0" target="_blank" rel="noreferrer">v0.26.0</a> <a class="header-anchor" href="#v0-26-0" aria-label="Permalink to &quot;[v0.26.0](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.26.0)&quot;">​</a></h2><p>Release date: 2025-02-09</p><ul><li><p>Fix CUDA <code>sparse_to_dense</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/386" target="_blank" rel="noreferrer">#386</a>)</p></li><li><p>Improve pseudo inverse spectrum solver. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/388" target="_blank" rel="noreferrer">#388</a>)</p></li><li><p>Add <code>smesolve</code> function for stochastic master equation. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/389" target="_blank" rel="noreferrer">#389</a>)</p></li></ul><h2 id="v0-25-2" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.2" target="_blank" rel="noreferrer">v0.25.2</a> <a class="header-anchor" href="#v0-25-2" aria-label="Permalink to &quot;[v0.25.2](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.2)&quot;">​</a></h2><p>Release date: 2025-02-02</p><ul><li><p>Move code quality dependencies to separate environment. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/380" target="_blank" rel="noreferrer">#380</a>)</p></li><li><p>Add additional normalization of the state during time evolution of <code>ssesolve</code>. This improves the numerical stability of the solver. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/383" target="_blank" rel="noreferrer">#383</a>)</p></li></ul><h2 id="v0-25-1" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.1" target="_blank" rel="noreferrer">v0.25.1</a> <a class="header-anchor" href="#v0-25-1" aria-label="Permalink to &quot;[v0.25.1](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.1)&quot;">​</a></h2><p>Release date: 2025-01-29</p><ul><li><p>Fix Dynamical Fock Dimension states saving due to wrong saving of dimensions. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/375" target="_blank" rel="noreferrer">#375</a>)</p></li><li><p>Support a list of observables for <code>expect</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/374" target="_blank" rel="noreferrer">#374</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/376" target="_blank" rel="noreferrer">#376</a>)</p></li><li><p>Add checks for <code>tlist</code> in time evolution solvers. The checks are to ensure that <code>tlist</code> is not empty, the elements are in increasing order, and the elements are unique. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/378" target="_blank" rel="noreferrer">#378</a>)</p></li></ul><h2 id="v0-25-0" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.0" target="_blank" rel="noreferrer">v0.25.0</a> <a class="header-anchor" href="#v0-25-0" aria-label="Permalink to &quot;[v0.25.0](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.0)&quot;">​</a></h2><p>Release date: 2025-01-20</p><ul><li><p>Change the structure of block diagonalization functions, using <code>BlockDiagonalForm</code> struct and changing the function name from <code>bdf</code> to <code>block_diagonal_form</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/349" target="_blank" rel="noreferrer">#349</a>)</p></li><li><p>Add <strong>GPUArrays</strong> compatibility for <code>ptrace</code> function, by using <strong>KernelAbstractions.jl</strong>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/350" target="_blank" rel="noreferrer">#350</a>)</p></li><li><p>Introduce <code>Space</code>, <code>Dimensions</code>, <code>GeneralDimensions</code> structures to support wider definitions and operations of <code>Qobj/QobjEvo</code>, and potential functionalities in the future. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/271" target="_blank" rel="noreferrer">#271</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/353" target="_blank" rel="noreferrer">#353</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/360" target="_blank" rel="noreferrer">#360</a>)</p></li><li><p>Improve lazy tensor warning for <code>SciMLOperators</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/370" target="_blank" rel="noreferrer">#370</a>)</p></li><li><p>Change order of <code>AbstractQuantumObject</code> data type. For example, from <code>QuantumObject{DataType,ObjType,DimsType}</code> to <code>QuantumObject{ObjType,DimsType,DataType}</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/371" target="_blank" rel="noreferrer">#371</a>)</p></li></ul><h2 id="v0-24-0" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.24.0" target="_blank" rel="noreferrer">v0.24.0</a> <a class="header-anchor" href="#v0-24-0" aria-label="Permalink to &quot;[v0.24.0](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.24.0)&quot;">​</a></h2><p>Release date: 2024-12-13</p><ul><li><p>Improve the construction of <code>QobjEvo</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/338" target="_blank" rel="noreferrer">#338</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/339" target="_blank" rel="noreferrer">#339</a>)</p></li><li><p>Support <code>Base.zero</code> and <code>Base.one</code> for <code>AbstractQuantumObject</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/342" target="_blank" rel="noreferrer">#342</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/346" target="_blank" rel="noreferrer">#346</a>)</p></li><li><p>Introduce visualization and function <code>plot_wigner</code> for easy plotting of Wigner functions. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/86" target="_blank" rel="noreferrer">#86</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/292" target="_blank" rel="noreferrer">#292</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/347" target="_blank" rel="noreferrer">#347</a>)</p></li></ul><h2 id="v0-23-1" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.23.1" target="_blank" rel="noreferrer">v0.23.1</a> <a class="header-anchor" href="#v0-23-1" aria-label="Permalink to &quot;[v0.23.1](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.23.1)&quot;">​</a></h2><p>Release date: 2024-12-06</p><ul><li>Update <code>[compat]</code> to fix the incompatibility between <code>QuantumToolbox v0.22.0+</code> and <code>DiffEqCallbacks &lt; v4.2.1</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/335" target="_blank" rel="noreferrer">#335</a>)</li></ul><h2 id="v0-23-0" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.23.0" target="_blank" rel="noreferrer">v0.23.0</a> <a class="header-anchor" href="#v0-23-0" aria-label="Permalink to &quot;[v0.23.0](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.23.0)&quot;">​</a></h2><p>Release date: 2024-12-04</p><ul><li><p>Change <code>SingleSiteOperator</code> with the more general <code>MultiSiteOperator</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/324" target="_blank" rel="noreferrer">#324</a>)</p></li><li><p>Make <code>spectrum</code> and <code>correlation</code> functions align with <code>Python QuTiP</code>, introduce spectrum solver <code>PseudoInverse</code>, remove spectrum solver <code>FFTCorrelation</code>, and introduce <code>spectrum_correlation_fft</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/330" target="_blank" rel="noreferrer">#330</a>)</p></li></ul><h2 id="v0-22-0" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.22.0" target="_blank" rel="noreferrer">v0.22.0</a> <a class="header-anchor" href="#v0-22-0" aria-label="Permalink to &quot;[v0.22.0](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.22.0)&quot;">​</a></h2><p>Release date: 2024-11-20</p><ul><li><p>Change the parameters structure of <code>sesolve</code>, <code>mesolve</code> and <code>mcsolve</code> functions to possibly support automatic differentiation. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/311" target="_blank" rel="noreferrer">#311</a>)</p></li><li><p>Fix type instability and reduce extra memory allocation in <code>liouvillian</code>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/315" target="_blank" rel="noreferrer">#315</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/318" target="_blank" rel="noreferrer">#318</a>)</p></li></ul><h2 id="v0-21-5" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.21.5" target="_blank" rel="noreferrer">v0.21.5</a> <a class="header-anchor" href="#v0-21-5" aria-label="Permalink to &quot;[v0.21.5](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.21.5)&quot;">​</a></h2><p>Release date: 2024-11-15</p><ul><li>This is a demonstration of how to bump version number and also modify <code>CHANGELOG.md</code> before new release. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/309" target="_blank" rel="noreferrer">#309</a>)</li></ul><h2 id="v0-21-4" tabindex="-1"><a href="https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.21.4" target="_blank" rel="noreferrer">v0.21.4</a> <a class="header-anchor" href="#v0-21-4" aria-label="Permalink to &quot;[v0.21.4](https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.21.4)&quot;">​</a></h2><p>Release date: 2024-11-13</p><ul><li>This is just a demonstration about <a href="https://github.com/JuliaDocs/Changelog.jl" target="_blank" rel="noreferrer"><code>Changelog.jl</code></a>. (<a href="https://github.com/qutip/QuantumToolbox.jl/issues/139" target="_blank" rel="noreferrer">#139</a>, <a href="https://github.com/qutip/QuantumToolbox.jl/issues/306" target="_blank" rel="noreferrer">#306</a>)</li></ul>',35)]))}const d=t(l,[["render",i]]);export{b as __pageData,d as default};
