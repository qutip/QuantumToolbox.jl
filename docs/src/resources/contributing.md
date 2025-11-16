# [Contributing to Quantum Toolbox in Julia](@id doc-Contribute)

[![Zulip](https://img.shields.io/badge/chat-on%20Zulip-6f73af.svg)](https://quantumtoolbox-jl.zulipchat.com)

## [Quick Start](@id doc-Contribute:Quick-Start)

### [Git Tutorial](@id doc-Contribute:Git-Tutorial)

`QuantumToolbox.jl` is developed using the [`git`](https://git-scm.com/) version-control system, with the [main repository](https://github.com/qutip/QuantumToolbox.jl) hosted in the [qutip organisation on GitHub](https://github.com/qutip). You will need to be familiar with [`git`](https://git-scm.com/) as a tool, and the [GitHub Flow](https://docs.github.com/en/get-started/quickstart/github-flow) workflow for branching and making pull requests. The exact details of environment set-up, build process, and runtests vary by repository are discussed below. In overview, the steps to contribute are:

- Consider creating an issue on the GitHub page of the relevant repository, describing the change you think should be made and why, so we can discuss details with you and make sure it is appropriate.
- *__If this is your first contribution__*, make a fork of the relevant repository on GitHub (which will be called as `origin`) and clone it to your local computer. Also add our copy as a remote (let's call it `qutip` here): `git remote add qutip https://github.com/qutip/<repo-name>`.
- Start from the `main` branch in your local computer (`git checkout main`), and pull all the changes from the remote (`qutip`) repository (on GitHub) to make sure you have an up-to-date copy: `git pull qutip main`.
- Switch to a new `git` branch in your local computer: `git checkout -b <branch-name>`.
- Make the changes you want.
- Go through the following build processes (if the changes you made relates to any of them) locally in your computer to build the final result so you can check your changes work sensibly:
    - Write and make sure all the runtests pass. See [here](@ref doc-Contribute:Runtests) for more details.
    - Make sure all the changes match the `Julia` code format (style). See [here](@ref doc-Contribute:Julia-Code-Format) for more details.
    - Improve and make sure the documentation can be built successfully. See [here](@ref doc-Contribute:Documentation) for more details.
    - Update changelog. See [here](@ref doc-Contribute:ChangeLog) for more details.
- Add the changed files to the `git` staging area `git add <file1> <file2> ...`, and then create some commits with short-but-descriptive names: `git commit`.
- Push the changes to your fork (`origin`): `git push -u origin <branch-name>`. You won’t be able to push to the remote (`qutip`) repositories directly.
- Go to the GitHub website for the repository you are contributing to, click on the “Pull Requests” tab, click the “New Pull Request” button, and follow the instructions there.

Once the pull request (PR) is created, some members of the QuTiP admin team will review the code to make sure it is suitable for inclusion in the library, to check the programming, and to ensure everything meets our standards. For some repositories, several automated CI pipelines will run whenever you create or modify a PR. In general, these will basically be the same ones which you can run locally, and all CI pipelines are required to pass online before your changes are merged to the remote `main` branch. There might be some feedbacks and requested changes. You can add more commits to address these, and push them to the branch (`<branch-name>`) of your fork (`origin`) to update the PR.

### [Make Command](@id doc-Contribute:Make-Command)

!!! note "Install `make` in Windows"
    The `make` command is not available by default on Windows. You will need to install a compatible version of `make` (for example via Git-for-Windows, MSYS2, or similar tools). A helpful reference for setting this up can be found here: https://stackoverflow.com/questions/66525016/how-to-run-make-command-in-gitbash-in-windows

To simplify common developer tasks, `QuantumToolbox.jl` provides a set of predefined `make` commands. These commands are defined in the repository’s file called `Makefile` and offer a convenient interface for running documentation builds, tests, formatting, and other maintenance routines.

If you are contributing for the first time, or have installed a fresh `julia` version, you should begin with:

```shell
make setup
```

to install necessary dependencies. Here are several other common commands:

```shell
make test      # run the tests
make format    # format codes with JuliaFormatter
make docs      # instantiate and build the documentation
make changelog # generate changelog
make help      # view all commands
```

Each command may accept additional options or environment variables to customize its behavior. The following sections describe these commands in detail and explain the available configuration options.

## [Programming Standards](@id doc-Contribute:Programming-Standards)

The rest of this document covers programming standards for `QuantumToolbox.jl`.

### [Runtests](@id doc-Contribute:Runtests)

All the test scripts should be located in the folder `test` in the repository. To run the test, use the following command under the *__root directory of the repository__* you are working on:

```shell
make test
```

This command will automatically rebuild `Julia` and run the script located in `test/runtests.jl` (should cover both the original tests and the new test(s) you add).

The tests are divided into several test groups, where the group names are defined in the file `test/runtests.jl` with a variable `GROUP`. One can also run the test scripts just for a certain test group by adding an argument `GROUP=<test-group-name>` to the `make test` command. For example, to run the tests for group `Core`, one can use the following command:

```shell
make GROUP=Core test
```

#### [Test Item Framework for Core tests](@id doc-Contribute:Test-Item-Framework-for-Core-tests)

The tests in `GROUP=Core` are provided using the [Test Item Framework](https://www.julia-vscode.org/docs/stable/userguide/testitems/), which structures the test codes into `@testitems` and makes it easier to run individually.

The [VS Code](https://code.visualstudio.com/) and its [Julia extension](https://www.julia-vscode.org/) provides us with options to run individual `@testitems`. It is much easier to find the specific core test that failed since the [Julia extension](https://www.julia-vscode.org/) in [VS Code](https://code.visualstudio.com/) will collect all core test failures and then display them in a structured way, directly at the place in the code where a specific core test failed. See [here](https://www.julia-vscode.org/docs/stable/userguide/testitems/) for more details.

### [Julia Code Format](@id doc-Contribute:Julia-Code-Format)

We use [`JuliaFormatter.jl`](https://github.com/domluna/JuliaFormatter.jl) to format all the source codes. The code style and extra formatting options is defined in the file `.JuliaFormatter.toml` in the repository.

To format the changed codes, use the following command under the *__root directory of the repository__* you are working on:

!!! note "Requirements"
    If this is your first time running `make` command in the local repository you are working on or you just had reinstalled `Julia`, you should run `make setup` first.

```shell
make format
```

### [Documentation](@id doc-Contribute:Documentation)

All the documentation source files [in markdown (`.md`) format] and build scripts should be located in the folder `docs` in the repository.

The document pages will be generated in the folder `docs/build/1/` (which is ignored by `git`) in the repository.

To instantiate and build the documentation, run the following command under the *__root directory of the repository__* you are working on:

!!! note "Requirements"
    You need to install `Node.js` and `npm` first.

```shell
make docs
```

This command will automatically rebuild `Julia` and run the script located in `docs/make.jl` (should be able to build the necessary files for the documentation).

One can also skip several processes in documentation build by adding arguments to the `make docs` command. For example,

```shell
make docs DOCTEST=false # skip doc tests
make docs DRAFT=true    # disable cell evaluation
```

To read the documentation in a browser, you can run the following command:

```shell
make vitepress
```

This will start a local Vitepress site of documentation at `http://localhost:5173` in your computer.

### [ChangeLog](@id doc-Contribute:ChangeLog)

The changelog is written in the file `CHANGELOG.md` in the repository. If you add some changes to the repository and made a PR, you should also add some messages or release notes together with the related PRs/issues entries to `CHANGELOG.md`. For example, add a new line in `CHANGELOG.md`:

```markdown
- some messages to describe the changes. ([#issue-ID], [#PR-ID])
```

See also the [ChangeLog page](@ref ChangeLog) for more examples.

After that, you can run the following command under the *__root directory of the repository__* you are working on:

!!! note "Requirements"
    If this is your first time running `make` command in the local repository you are working on or you just had reinstalled `Julia`, you should run `make setup` first.

```shell
make changelog
```

This will automatically generate the full URLs for the references to PRs/issues by utilizing [`Changelog.jl`](https://github.com/JuliaDocs/Changelog.jl).

If the changes you made are not necessary to be recorded in `CHANGELOG.md`, you can add the label `[Skip ChangeLog]` to the PR you made in the GitHub repository.
