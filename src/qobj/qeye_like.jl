export qeye_like

function eye(x :: T) where T<:Union{AbstractQuantumObject{Operator}, AbstractQuantumObject{SuperOperator}}
    qeye_like(x)
end
function qeye_like(x :: T) where T<:QuantumObject{Operator}
    y = 0*x
    for i in 1:size(x)[1]
        y[i,i] = 1
    end
    return y
end

function qeye_like(x :: T) where T<:QobjEvo
    y = 0*x(0)
    for i in 1:size(x)[1]
        y[i,i] = 1
    end
    return y
end

function qeye_like(x :: T) where T<:QuantumObject{Ket}
    qeye_like(x*x')
end
function qeye_like(x :: T) where T<:QuantumObject{Bra}
    qeye_like(x'*x)
end

function qeye_like(x::T) where T<:QuantumObject{SuperOperator}
    y = 0*x
    for i in 1:size(x)[1]
        y[i,i] = 1
    end
    return y
end