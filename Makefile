JULIA:=julia

default: help

docs:
	${JULIA} --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=docs docs/make.jl

format:
	${JULIA} -e 'using JuliaFormatter; format(".")'

test:
	${JULIA} --project -e 'using Pkg; Pkg.resolve(); Pkg.test()'

all: format test docs

help:
	@echo "The following make commands are available:"
	@echo " - make docs: instantiate and build the documentation"
	@echo " - make format: format codes with JuliaFormatter"
	@echo " - make test: run the tests"
	@echo " - make all: run every commands above"

.PHONY: default docs format test all help