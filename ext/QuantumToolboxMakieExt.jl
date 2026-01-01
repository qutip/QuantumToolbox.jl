module QuantumToolboxMakieExt

using QuantumToolbox
import QuantumToolbox: _state_to_bloch

import LinearAlgebra: cross, deg2rad, normalize, size
import Makie:
    Axis,
    Axis3,
    LScene,
    Colorbar,
    Figure,
    GridLayout,
    heatmap!,
    surface!,
    barplot!,
    GridPosition,
    @L_str,
    Reverse,
    ylims!,
    RGBAf,
    Sphere,
    lines!,
    scatter!,
    arrows3d!,
    text!,
    mesh!,
    RGBf,
    Point3f,
    NoShading,
    cameracontrols,
    update_cam!,
    cam3d!

include("MakieExt/utils.jl")

include("MakieExt/bloch_sphere.jl")
include("MakieExt/fock_distribution.jl")
include("MakieExt/wigner.jl")

end
