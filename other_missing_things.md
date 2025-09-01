# Things that are missing from the QuantumToolbox.jl package.

1) Powers are not supported for `QobjEvo` objects.
   - There is a bug with the one I did, as the following returns a stackoverflow error.
    ```julia
        import Base.^
        import Base.^
        function ^(x::QobjEvo, y::Integer)
            if y< 0
                error("Exponent must be non-negative integer")
            end
            if y == 0
                return QobjEvo(qeye(size(x)[1], dims = x.dims))
            end
            to_return = x
            for i in 2:y
                to_return = to_return * x
            end
            return to_return
        end

        1*(sum([(QobjEvo((qeye(2), (p,t) -> 1)))^n*1 for n in 0:1:2]))
        1*(sum([1*(QobjEvo((qeye(2), (p,t) -> 1)))^n for n in 0:1:2]))
   ```
   - However, the following do not
   ```julia
    1*(sum([(QobjEvo((qeye(2), (p,t) -> 1)))^n for n in 0:1:2]))
    1*(sum([(1*QobjEvo((qeye(2), (p,t) -> 1)))^n for n in 0:1:2]))
    1*(sum([(QobjEvo((qeye(2), (p,t) -> 1))*1)^n for n in 0:1:2]))
    ```

    The error comes from expressions like 
    ```julia
    a = QobjEvo((qeye(2), (p,t) -> 1))
    1*(a+a*a*1)
    ```

    The error is:
    StackOverflowError:
    ```julia
        StackOverflowError:
        
        Stacktrace:
        [1] *(α::ScalarOperator{Int64, SciMLOperators.FilterKwargs{typeof(SciMLOperators.DEFAULT_UPDATE_FUNC), Tuple{}}}, x::UniformScaling{Int64})
        @ SciMLOperators ~/.julia/packages/SciMLOperators/DKrpP/src/scalar.jl:366
        [2] *(α::ScalarOperator{Int64, SciMLOperators.FilterKwargs{typeof(SciMLOperators.DEFAULT_UPDATE_FUNC), Tuple{}}}, x::UniformScaling{Int64})
        @ SciMLOperators ~/.julia/packages/SciMLOperators/DKrpP/src/scalar.jl:369    
    ```