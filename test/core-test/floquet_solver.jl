using QuantumToolbox
using Test 
using Random

function _matrix_element(c_op_fmmesolve, vp,ep)
    return vp' * c_op_fmmesolve * ep
end

function _convert_c_ops(c_op_fmmesolve, noise_spectrum, vp, ep)
    c_op_mesolve = []
    N = length(vp)
    for i in 1:N
        for j in 1:N
            if i != j
                #calculate the rate 

                gamma = 2 * np.pi * _matrix_element(c_op_fmmesolve, vp[j], vp[i]) * 
                                _matrix_element(c_op_fmmesolve, vp[i], vp[j]) * 
                                noise_spectrum(ep[j] - ep[i])
                
                # add c_op for mesolve
                push!(c_op_mesolve, sqrt(gamma) * (vp[i] * vp[j]'))
                return c_op_mesolve
            end
        end
    end
end



@testitem "Floquet Solver" begin
    N = 10
    a = destroy(N)
    a_d = a'
    coef(p,t) = cos(t)
    H0 = num(N)
    H1 = a + a_d
    H_tuple = (H0,(H1,coef))
    T = 2 * π
    tlist = range(0.0, 3T, length=101)
    floquet_basis = FloquetBasis(H_tuple, T, tlist)
    psi0 = rand_ket(N)
    floquet_psi0 = to_floquet_basis(floquet_basis, psi0)
    sol = sesolve(H_tuple, psi0, tlist, e_ops = [], saveat = tlist)
    states = sol.states
    fse = floquet_sesolve(H_tuple, psi0, tlist, T=T)
    states_fse = fse.states
    for (t,state) in zip(tlist,states)
        from_floquet = from_floquet_basis(floquet_basis, floquet_psi0, t)
        @test overlap(state, from_floquet) ≈ 1.0 atol=8e-5  
    end
    for (state_se,state_fse) in zip(states, states_fse)
        @test overlap(state_se, state_fse) ≈ 1.0 atol=5e-5  
    end    
end

@testitem "Floquet Master equation" begin
    delta = 1.0 * 2π
    eps0 = 1.0 * 2π
    A = 0.5 * 2π
    omega = sqrt(delta^2 + eps0^2)
    T = 2π / omega
    tlist = range(0.0, 2 * T, 101)
    psi0 = rand_ket(2)
    H0 = - eps0 / 2 * sigmaz() - delta / 2 * sigmax()
    H1 = A/2.0 * sigmax()
    H = (H0, (H1, t -> sin(omega * t)))
    e_ops = [num(2)]
    gamma1 = 0

    #collapse operator for Floquet-Markov master equation
    c_op_fmmesolve = sigmax()

    #collapse operator for Lindblad master equation
    @inline spectrum(omega::Real) = ifelse(ω>0, ω * 0.5 * gamma1 / (2π), zero(ω))

    (ep, vp) = eigenstates(H0)
    op0 = vp[1] * vp[1]'
    op1 = vp[2] * vp[2]'

    c_op_mesolve = _convert_c_ops(c_op_fmmesolve, spectrum, vp, ep)

    #Solve the floquet markov master equation
    p_ex = fmmesolve()    
end