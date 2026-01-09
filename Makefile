JULIA:=julia

default: help

setup:
    ${JULIA} -e 'import Pkg; Pkg.add(["Changelog"])'
    ${JULIA} --project=@runic --startup-file=no -e 'using Pkg; Pkg.add("Runic")'

format:
    ${JULIA} --project=@runic --startup-file=no -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace .

changelog:
	${JULIA} -e 'using Changelog; Changelog.generate(Changelog.CommonMark(), "CHANGELOG.md"; repo = "qutip/QuantumToolbox.jl")'

test:
	${JULIA} --project -e 'using Pkg; Pkg.update(); Pkg.test()'

docs:
	${JULIA} --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.update()'
	${JULIA} --project=docs docs/make.jl

vitepress:
	npm --prefix docs run docs:dev

all: setup format changelog test docs vitepress

help:
	@echo "The following make commands are available:"
	@echo " - make setup: install the dependencies for make command"
	@echo " - make format: format codes with Runic"
	@echo " - make changelog: generate changelog"
	@echo " - make test: run the tests"
	@echo " - make docs: instantiate and build the documentation"
	@echo " - make vitepress: start Vitepress site of documentation"
	@echo " - make all: run every commands in the above order"

.PHONY: default setup format changelog test docs vitepress all help
