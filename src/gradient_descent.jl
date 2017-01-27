function cumulate_grad!(∇m, m, neigs)
    p = 1.
    @inbounds for si in neigs
        J, i = sign(si), abs(si)
        p *= (1-J*m[i]) / 2
    end
    @inbounds for si in neigs
        J, i = sign(si), abs(si)
        ∇m[i] += J * p/(1-J*m[i]) / (1-p)
        @assert isfinite(∇m[i]) "$(∇m[i]) $(m[i]) p=$p "
    end
end

clip(h::Vector{Float64}) = Int.(sign.(h))

function solve_gd(cnf::CNF;
                    maxiters = 5000, ϵ::Float64 = 1e-4,
                    γ::Float64 = 0.1, γstep::Float64=0.01, γmax = 0.5,
                    η = 0.01, # gd learning rate
                    seed::Int = -1,
                    infotime = 100)

    seed > 0 && srand(seed)
    N, M = cnf.N, cnf.M
    h = randn(N)
    Δ = 0.


    σ = clip(h)
    E = energy(cnf, σ)
    @printf("it=0 γ=%.3f\t  E=%d   \tΔ=%f\n", γ, E, Δ)
    for it=1:maxiters
        m = tanh.(γ .* h)
        ∇m = zeros(N)
        for μ=1:M
            cumulate_grad!(∇m, m, cnf.clauses[μ])
        end

        h .+= γ .* η .* (1 .- m.^2) .* ∇m
        h ./= norm(h)/√N

        σ = clip(h)
        E = energy(cnf, σ)
        # it%infotime == 0 && @printf("\rit=%d γ=%.3f\t  E=%d   \tΔ=%f", it, γ, E, Δ)
        if it % infotime == 0
            @printf("it=%d γ=%.3f\t  E=%d \n", it, γ, E)
            γ += γstep
        end

        E == 0 && break
    end
    @printf("γ=%.3f\t  E=%d\n", γ, E)

    return clip(h)
end
