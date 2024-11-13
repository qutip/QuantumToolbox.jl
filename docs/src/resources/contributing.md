# [Contributing to QuantumToolbox.jl](@id doc-Contribute)

## [Quick Start](@id doc-Contribute:Quick-Start)

`QuantumToolbox.jl` is developed using the [`git`](https://git-scm.com/) version-control system, with the [main repository](https://github.com/qutip/QuantumToolbox.jl) hosted in the [qutip organisation on GitHub](https://github.com/qutip). You will need to be familiar with [`git`](https://git-scm.com/) as a tool, and the [GitHub Flow](https://docs.github.com/en/get-started/quickstart/github-flow) workflow for branching and making pull requests. The exact details of environment set-up, build process, and runtests vary by repository are discussed below. In overview, the steps to contribute are:

1. Consider creating an issue on the GitHub page of the relevant repository, describing the change you think should be made and why, so we can discuss details with you and make sure it is appropriate.
2. *__If this is your first contribution__*, make a fork of the relevant repository on GitHub and clone it to your local computer. Also add our copy as a remote called `qutip`: `git remote add qutip https://github.com/qutip/<repo-name>`.
3. Start from the `main` branch in your local computer (`git checkout main`), and pull all the changes from the remote repository (on GitHub) to make sure you have an up-to-date copy: `git pull qutip main`.
4. Switch to a new `git` branch: `git checkout -b <branch-name>`.
5. Make the changes you want
6. Go through the following build processes (if the changes you made relates to any of them) for the repository to build the final result so you can check your changes work sensibly:
    - Write and make sure all the runtests pass. See [here](@ref doc-Contribute:Runtests) for more details.
    - Make sure all the changes match the `Julia` code format (style). See [here](@ref doc-Contribute:Julia-Code-Format) for more details.
    - Improve and make sure the documentation can be built successfully. See [here](@ref doc-Contribute:Documentation) for more details.
    - Update changelog. See [here](@ref doc-Contribute:Update-ChangeLog) for more details.
7. Add the changed files to the `git` staging area `git add <file1> <file2> ...`, and then create some commits with short-but-descriptive names `git commit`.
8. Push the changes to your fork (`origin`) `git push -u origin <branch-name>`. You won’t be able to push to the remote `qutip` repositories directly.
9. Go to the GitHub website for the repository you are contributing to, click on the “Pull Requests” tab, click the “New Pull Request” button, and follow the instructions there.

Once the pull request (PR) is created, some members of the QuTiP admin team will review the code to make sure it is suitable for inclusion in the library, to check the programming, and to ensure everything meets our standards. For some repositories, several automated CI pipelines will run whenever you create or modify a PR. In general, these will basically be the same ones which you can run locally, and all CI pipelines are required to pass online before your changes are merged to the remote `main` branch. There might be some feedbacks and requested changes. You can add more commits to address these, and push them to the branch (`<branch-name>`) of your fork (`origin`) to update the PR.

The rest of this document covers programming standards.

## [Runtests](@id doc-Contribute:Runtests)

All the test scripts should be located in the folder `test` in the repository. To run the test, use the following command under the *__root directory of the repository__* you are working on:

```shell
make test
```

This command will automatically rebuild `Julia` and run the script located in `test/runtests.jl` (should cover both the original tests and the new test(s) you add).

## [Julia Code Format](@id doc-Contribute:Julia-Code-Format)

We use [`JuliaFormatter.jl`](https://github.com/domluna/JuliaFormatter.jl) to format all the source codes. The code style and extra formatting options is defined in the file `.JuliaFormatter.toml` in the repository.

To format the changed codes, use the following command under the *__root directory of the repository__* you are working on:

!!! note "Requirements"
    If this is your first time running `make` command in the local repository you are working on or you just had reinstalled `Julia`, you should run `make setup` first.

```shell
make format
```

## [Documentation](@id doc-Contribute:Documentation)

All the documentation source files [in markdown (`.md`) format] and build scripts should be located in the folder `docs` in the repository.

The document pages will be generated in the folder `docs/build` (which is ignored by `git`) in the repository.

To instantiate and build the documentation, run the following command under the *__root directory of the repository__* you are working on:

```shell
make docs
```

This command will automatically rebuild `Julia` and run the script located in `docs/make.jl` (should be able to build the necessary files for the documentation).

To read the documentation in a browser, you can run the following command:

!!! note "Requirements"
    You need to install `Node.js` and `npm` first.

```shell
make vitepress
```

This will start a local Vitepress site of documentation at `http://localhost:5173/QuantumToolbox.jl/`

## [Update ChangeLog](@id doc-Contribute:Update-ChangeLog)

The changelog is written in the file `CHANGELOG.md` in the repository. If you add some changes to the package and made a PR, you should also add some messages or release notes together with the related PRs or issues entries to `CHANGELOG.md`, for example:

```markdown
- some messages to describe the changes. ([#issue-ID], [#PR-ID])
```

See also the [Change Log page](@ref ChangeLog) for more examples.

After that, you can run the following command under the *__root directory of the repository__* you are working on:

!!! note "Requirements"
    If this is your first time running `make` command in the local repository you are working on or you just had reinstalled `Julia`, you should run `make setup` first.

```shell
make changelog
```

This will automatically generate the full URLs for the references to PRs or issues by utilizing [`Changelog.jl`](https://github.com/JuliaDocs/Changelog.jl).

If the changes you made are not necessary to be recorded in `CHANGELOG.md`, you can add the label `[Skip ChangeLog]` to the PR you made in the GitHub repository.
