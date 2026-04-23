@testitem "Cite" begin

    # citation bibtex
    io_buffer = IOBuffer()
    QuantumToolbox.cite(io_buffer)
    captured_output = String(take!(io_buffer))
    @test captured_output ==
        """@article{QuantumToolbox.jl2025,\n""" *
        """  title = {Quantum{T}oolbox.jl: {A}n efficient {J}ulia framework for simulating open quantum systems},\n""" *
        """  author = {Mercurio, Alberto and Huang, Yi-Te and Cai, Li-Xun and Chen, Yueh-Nan and Savona, Vincenzo and Nori, Franco},\n""" *
        """  journal = {{Quantum}},\n""" *
        """  issn = {2521-327X},\n""" *
        """  publisher = {{Verein zur F{\\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},\n""" *
        """  volume = {9},\n""" *
        """  pages = {1866},\n""" *
        """  month = sep,\n""" *
        """  year = {2025},\n""" *
        """  doi = {10.22331/q-2025-09-29-1866},\n""" *
        """  url = {https://doi.org/10.22331/q-2025-09-29-1866}\n""" *
        """}\n"""
end
