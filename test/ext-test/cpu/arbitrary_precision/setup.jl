# Setup problem for arbitrary precision tests
module TESetup
setprecision(128) # Instead of 256. This speeds up the tests.

using QuantumToolbox
QuantumToolbox.settings.default_eltype_warning = false # these tests create higher precision objects on purpose

N = 20
Δ = 1.0
U = 0.1
F = 0.1
κ = 0.1

tlist = range(0, 40, 100)

# The model does not specify any element type, so it uses `QuantumToolbox.settings.default_eltype`
function model()
    a = destroy(N)
    ψ0 = fock(N, 1)
    H = Δ * a' * a + U / 2 * a' * a' * a * a + F * (a + a')
    c_ops = [sqrt(real(eltype(a))(κ)) * a]
    return a, ψ0, H, c_ops
end

a, ψ0, H, c_ops = model()

QuantumToolbox.settings.default_eltype = Complex{BigFloat}
try
    global a_big, ψ0_big, H_big, c_ops_big = model()
finally
    QuantumToolbox.settings.default_eltype = ComplexF64 # do not affect the other tests running in the same process
end
end
