module QuantumToolboxMakieExt

using QuantumToolbox
import QuantumToolbox:
    makeVal, getVal, _state_to_bloch, _handle_matrix_plot_data, _gen_default_ket_labels, _gen_default_bra_labels

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
    meshscatter!,
    RGBf,
    Point3f,
    Rect3f,
    Vec3f,
    NoShading,
    cameracontrols,
    update_cam!,
    cam3d!,
    translate!

include("MakieExt/utils.jl")

include("MakieExt/bloch_sphere.jl")
include("MakieExt/fock_distribution.jl")
include("MakieExt/matrix.jl")
include("MakieExt/wigner.jl")

end
