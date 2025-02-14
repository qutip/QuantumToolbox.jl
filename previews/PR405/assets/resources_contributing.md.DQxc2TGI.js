import{_ as t,c as o,o as i,ai as a}from"./chunks/framework.COX_rLcS.js";const p=JSON.parse('{"title":"Contributing to Quantum Toolbox in Julia","description":"","frontmatter":{},"headers":[],"relativePath":"resources/contributing.md","filePath":"resources/contributing.md","lastUpdated":null}'),s={name:"resources/contributing.md"};function r(n,e,l,d,c,h){return i(),o("div",null,e[0]||(e[0]=[a('<h1 id="doc-Contribute" tabindex="-1">Contributing to Quantum Toolbox in Julia <a class="header-anchor" href="#doc-Contribute" aria-label="Permalink to &quot;Contributing to Quantum Toolbox in Julia {#doc-Contribute}&quot;">​</a></h1><h2 id="doc-Contribute:Quick-Start" tabindex="-1">Quick Start <a class="header-anchor" href="#doc-Contribute:Quick-Start" aria-label="Permalink to &quot;Quick Start {#doc-Contribute:Quick-Start}&quot;">​</a></h2><p><code>QuantumToolbox.jl</code> is developed using the <a href="https://git-scm.com/" target="_blank" rel="noreferrer"><code>git</code></a> version-control system, with the <a href="https://github.com/qutip/QuantumToolbox.jl" target="_blank" rel="noreferrer">main repository</a> hosted in the <a href="https://github.com/qutip" target="_blank" rel="noreferrer">qutip organisation on GitHub</a>. You will need to be familiar with <a href="https://git-scm.com/" target="_blank" rel="noreferrer"><code>git</code></a> as a tool, and the <a href="https://docs.github.com/en/get-started/quickstart/github-flow" target="_blank" rel="noreferrer">GitHub Flow</a> workflow for branching and making pull requests. The exact details of environment set-up, build process, and runtests vary by repository are discussed below. In overview, the steps to contribute are:</p><ul><li><p>Consider creating an issue on the GitHub page of the relevant repository, describing the change you think should be made and why, so we can discuss details with you and make sure it is appropriate.</p></li><li><p><em><strong>If this is your first contribution</strong></em>, make a fork of the relevant repository on GitHub (which will be called as <code>origin</code>) and clone it to your local computer. Also add our copy as a remote (let&#39;s call it <code>qutip</code> here): <code>git remote add qutip https://github.com/qutip/&lt;repo-name&gt;</code>.</p></li><li><p>Start from the <code>main</code> branch in your local computer (<code>git checkout main</code>), and pull all the changes from the remote (<code>qutip</code>) repository (on GitHub) to make sure you have an up-to-date copy: <code>git pull qutip main</code>.</p></li><li><p>Switch to a new <code>git</code> branch in your local computer: <code>git checkout -b &lt;branch-name&gt;</code>.</p></li><li><p>Make the changes you want.</p></li><li><p>Go through the following build processes (if the changes you made relates to any of them) locally in your computer to build the final result so you can check your changes work sensibly:</p><ul><li><p>Write and make sure all the runtests pass. See <a href="/QuantumToolbox.jl/previews/PR405/resources/contributing#doc-Contribute:Runtests">here</a> for more details.</p></li><li><p>Make sure all the changes match the <code>Julia</code> code format (style). See <a href="/QuantumToolbox.jl/previews/PR405/resources/contributing#doc-Contribute:Julia-Code-Format">here</a> for more details.</p></li><li><p>Improve and make sure the documentation can be built successfully. See <a href="/QuantumToolbox.jl/previews/PR405/resources/contributing#doc-Contribute:Documentation">here</a> for more details.</p></li><li><p>Update changelog. See <a href="/QuantumToolbox.jl/previews/PR405/resources/contributing#doc-Contribute:Update-ChangeLog">here</a> for more details.</p></li></ul></li><li><p>Add the changed files to the <code>git</code> staging area <code>git add &lt;file1&gt; &lt;file2&gt; ...</code>, and then create some commits with short-but-descriptive names: <code>git commit</code>.</p></li><li><p>Push the changes to your fork (<code>origin</code>): <code>git push -u origin &lt;branch-name&gt;</code>. You won’t be able to push to the remote (<code>qutip</code>) repositories directly.</p></li><li><p>Go to the GitHub website for the repository you are contributing to, click on the “Pull Requests” tab, click the “New Pull Request” button, and follow the instructions there.</p></li></ul><p>Once the pull request (PR) is created, some members of the QuTiP admin team will review the code to make sure it is suitable for inclusion in the library, to check the programming, and to ensure everything meets our standards. For some repositories, several automated CI pipelines will run whenever you create or modify a PR. In general, these will basically be the same ones which you can run locally, and all CI pipelines are required to pass online before your changes are merged to the remote <code>main</code> branch. There might be some feedbacks and requested changes. You can add more commits to address these, and push them to the branch (<code>&lt;branch-name&gt;</code>) of your fork (<code>origin</code>) to update the PR.</p><p>The rest of this document covers programming standards.</p><h2 id="doc-Contribute:Runtests" tabindex="-1">Runtests <a class="header-anchor" href="#doc-Contribute:Runtests" aria-label="Permalink to &quot;Runtests {#doc-Contribute:Runtests}&quot;">​</a></h2><p>All the test scripts should be located in the folder <code>test</code> in the repository. To run the test, use the following command under the <em><strong>root directory of the repository</strong></em> you are working on:</p><div class="language-shell vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">shell</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">make</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> test</span></span></code></pre></div><p>This command will automatically rebuild <code>Julia</code> and run the script located in <code>test/runtests.jl</code> (should cover both the original tests and the new test(s) you add).</p><p>The tests are divided into several test groups, where the group names are defined in the file <code>test/runtests.jl</code> with a variable <code>GROUP</code>. One can also run the test scripts just for a certain test group by adding an argument <code>GROUP=&lt;test-group-name&gt;</code> to the <code>make test</code> command. For example, to run the tests for group <code>Core</code>, one can use the following command:</p><div class="language-shell vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">shell</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">make</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> GROUP=Core</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> test</span></span></code></pre></div><h2 id="doc-Contribute:Julia-Code-Format" tabindex="-1">Julia Code Format <a class="header-anchor" href="#doc-Contribute:Julia-Code-Format" aria-label="Permalink to &quot;Julia Code Format {#doc-Contribute:Julia-Code-Format}&quot;">​</a></h2><p>We use <a href="https://github.com/domluna/JuliaFormatter.jl" target="_blank" rel="noreferrer"><code>JuliaFormatter.jl</code></a> to format all the source codes. The code style and extra formatting options is defined in the file <code>.JuliaFormatter.toml</code> in the repository.</p><p>To format the changed codes, use the following command under the <em><strong>root directory of the repository</strong></em> you are working on:</p><div class="tip custom-block"><p class="custom-block-title">Requirements</p><p>If this is your first time running <code>make</code> command in the local repository you are working on or you just had reinstalled <code>Julia</code>, you should run <code>make setup</code> first.</p></div><div class="language-shell vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">shell</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">make</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> format</span></span></code></pre></div><h2 id="doc-Contribute:Documentation" tabindex="-1">Documentation <a class="header-anchor" href="#doc-Contribute:Documentation" aria-label="Permalink to &quot;Documentation {#doc-Contribute:Documentation}&quot;">​</a></h2><p>All the documentation source files [in markdown (<code>.md</code>) format] and build scripts should be located in the folder <code>docs</code> in the repository.</p><p>The document pages will be generated in the folder <code>docs/build</code> (which is ignored by <code>git</code>) in the repository.</p><p>To instantiate and build the documentation, run the following command under the <em><strong>root directory of the repository</strong></em> you are working on:</p><div class="language-shell vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">shell</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">make</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> docs</span></span></code></pre></div><p>This command will automatically rebuild <code>Julia</code> and run the script located in <code>docs/make.jl</code> (should be able to build the necessary files for the documentation).</p><p>To read the documentation in a browser, you can run the following command:</p><div class="tip custom-block"><p class="custom-block-title">Requirements</p><p>You need to install <code>Node.js</code> and <code>npm</code> first.</p></div><div class="language-shell vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">shell</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">make</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> vitepress</span></span></code></pre></div><p>This will start a local Vitepress site of documentation at <code>http://localhost:5173/QuantumToolbox.jl/</code> in your computer.</p><h2 id="doc-Contribute:Update-ChangeLog" tabindex="-1">Update ChangeLog <a class="header-anchor" href="#doc-Contribute:Update-ChangeLog" aria-label="Permalink to &quot;Update ChangeLog {#doc-Contribute:Update-ChangeLog}&quot;">​</a></h2><p>The changelog is written in the file <code>CHANGELOG.md</code> in the repository. If you add some changes to the repository and made a PR, you should also add some messages or release notes together with the related PRs/issues entries to <code>CHANGELOG.md</code>. For example, add a new line in <code>CHANGELOG.md</code>:</p><div class="language-markdown vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">markdown</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#E36209;--shiki-dark:#FFAB70;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> some messages to describe the changes. ([</span><span style="--shiki-light:#032F62;--shiki-light-text-decoration:underline;--shiki-dark:#DBEDFF;--shiki-dark-text-decoration:underline;">#issue-ID</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], [</span><span style="--shiki-light:#032F62;--shiki-light-text-decoration:underline;--shiki-dark:#DBEDFF;--shiki-dark-text-decoration:underline;">#PR-ID</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><p>See also the <a href="/QuantumToolbox.jl/previews/PR405/resources/changelog#ChangeLog">ChangeLog page</a> for more examples.</p><p>After that, you can run the following command under the <em><strong>root directory of the repository</strong></em> you are working on:</p><div class="tip custom-block"><p class="custom-block-title">Requirements</p><p>If this is your first time running <code>make</code> command in the local repository you are working on or you just had reinstalled <code>Julia</code>, you should run <code>make setup</code> first.</p></div><div class="language-shell vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">shell</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">make</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> changelog</span></span></code></pre></div><p>This will automatically generate the full URLs for the references to PRs/issues by utilizing <a href="https://github.com/JuliaDocs/Changelog.jl" target="_blank" rel="noreferrer"><code>Changelog.jl</code></a>.</p><p>If the changes you made are not necessary to be recorded in <code>CHANGELOG.md</code>, you can add the label <code>[Skip ChangeLog]</code> to the PR you made in the GitHub repository.</p>',36)]))}const m=t(s,[["render",r]]);export{p as __pageData,m as default};
