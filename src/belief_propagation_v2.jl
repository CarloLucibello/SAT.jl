@inline Base.getindex(p::Ptr) = unsafe_load(p)
@inline Base.setindex!(p::Ptr{T}, x::T) where T = unsafe_store!(p, x)

# ζ̂(a→i) = P(σ_i != J_ai)
# ζ(i→a) = P(σ_i != J_ai)

const MessU = Float64  # = ζ̂(a→i) = P(σ_i != J_ai)
const MessH = Float64  # = ζ(i→a) = P(σ_i != J_ai)
@inline getref(v::Vector, i::Integer) = pointer(v, i)

const PU = Ptr{MessU}
const PH = Ptr{MessH}

const VU = Vector{MessU}
const VH = Vector{MessH}
const VRU = Vector{PU}
const VRH = Vector{PH}

mutable struct Fact
    ζlist::Vector{MessH} # incoming messages
    ζ̂list::VRU           # outgoing messages
end
Fact() = Fact(VH(), VRU())

# 
mutable struct Var
    pinned::Int
    ζ̂listp::Vector{MessU}
    ζ̂listm::Vector{MessU}
    ζlistp::VRH
    ζlistm::VRH

    #used only in BP+reinforcement
    ζreinfp::MessU
    ζreinfm::MessU
end

Var() = Var(0, VU(),VU(), VRH(), VRH(), 1., 1.)

abstract type FactorGraph end

struct FactorGraphKSAT <: FactorGraph
    N::Int
    M::Int
    fnodes::Vector{Fact}
    vnodes::Vector{Var}
    cnf::CNF

    function FactorGraphKSAT(cnf::CNF)
        @extract cnf: M N clauses
        println("# read CNF formula")
        println("# N=$N M=$M α=$(M/N)")
        fnodes = [Fact() for i=1:M]
        vnodes = [Var() for i=1:N]
        kf = map(length, clauses)
        kvp = zeros(Int, N)
        kvm = zeros(Int, N)
        for clause in clauses
            for id in clause
                if id > 0
                    kvm[abs(id)] += 1
                else
                    kvp[abs(id)] += 1
                end
            end
        end

        ## Reserve memory in order to avoid invalidation of Refs
        for (a,f) in enumerate(fnodes)
            sizehint!(f.ζlist, kf[a])
            sizehint!(f.ζ̂list, kf[a])
        end
        for (i,v) in enumerate(vnodes)
            sizehint!(v.ζ̂listm, kvm[i])
            sizehint!(v.ζ̂listp, kvp[i])
            sizehint!(v.ζlistm, kvm[i])
            sizehint!(v.ζlistp, kvp[i])
        end

        for (a, clause) in enumerate(clauses)
            for id in clause
                i = abs(id)
                @assert id != 0
                f = fnodes[a]
                v = vnodes[i]
                if id > 0
                    push!(v.ζ̂listm, MessU(0.))
                    push!(f.ζ̂list, getref(v.ζ̂listm,length(v.ζ̂listm)))

                    push!(f.ζlist, MessH(0.))
                    push!(v.ζlistm, getref(f.ζlist,length(f.ζlist)))
                else
                    push!(v.ζ̂listp, MessU(0.))
                    push!(f.ζ̂list, getref(v.ζ̂listp,length(v.ζ̂listp)))

                    push!(f.ζlist, MessH(0.))
                    push!(v.ζlistp, getref(f.ζlist,length(f.ζlist)))
                end
            end
        end
        new(N, M, fnodes, vnodes, cnf)
    end
end

mutable struct ReinfParams
    r::Float64
    rstep::Float64
    γ::Float64
    γstep::Float64
    tγ::Float64
    wait_count::Int

    ReinfParams(r=0.,rstep=0.,γ=0.,γstep=0.) = new(r, rstep, γ, γstep, tanh(γ))
end

deg(f::Fact) = length(f.ζ̂list)
degp(v::Var) = length(v.ζ̂listp)
degm(v::Var) = length(v.ζ̂listm)

function initrand!(g::FactorGraphKSAT)
    for f in g.fnodes
        for k=1:deg(f)
            f.ζlist[k] = rand()
        end
    end
    for v in g.vnodes
        for k=1:degp(v)
            r = 0.5*rand()
            v.ζ̂listp[k] = 1 - (1-2r)/(1-r)
        end
        for k=1:degm(v)
            r = 0.5*rand()
            v.ζ̂listm[k] = 1 - (1-2r)/(1-r)
        end
        v.ζreinfm = 1
        v.ζreinfp = 1
    end
end

function update!(f::Fact, dump=0.5)
    @extract f ζ̂list ζlist
    ζ̂tot = 1.
    eps = 1e-15
    nzeros = 0
    @inbounds for i=1:deg(f)
        if ζlist[i] > eps
            ζ̂tot *= ζlist[i]
        else
            nzeros += 1
        end
    end
    @inbounds for i=1:deg(f)
        if nzeros == 0
            ζ̂i = ζ̂tot / ζlist[i]
        elseif nzeros == 1 && ζlist[i] < eps
            ζ̂i = ζ̂tot
        else
            ζ̂i = 0.
        end
        ζ̂list[i][] = dump * ζ̂list[i][] +(1-dump) * (1 - ζ̂i)/(2-ζ̂i)
    end
end

setfree!(v::Var) = v.pinned = 0

function setpinned!(v::Var, σ::Int)
    #TODO check del denominatore=0
    @extract v  ζlistp ζlistm
    v.pinned = σ
    ### compute cavity fields
    for i=1:degp(v)
        ζlistp[i][] = σ > 0 ? 1. : 0.
    end

    for i=1:degm(v)
        ζlistm[i][] = σ > 0 ? 0. : 1.
    end
end

ispinned(v::Var) = v.pinned != 0
numpinned(g::FactorGraphKSAT) = sum(ispinned, g.vnodes)

# r = fraction of N to assign
function pin_most_biased!(g::FactorGraphKSAT, r::Float64 = 0.02)
    mlist = Vector{Tuple{Int,Float64}}()
    npin = numpinned(g)
    sizehint!(mlist, g.N - npin)
    for (i,v) in enumerate(g.vnodes)
        if !ispinned(v)
            push!(mlist, (i, mag(v)))
        end
    end

    ntopin = min(ceil(Int, r*g.N), length(mlist))
    println("# Pinning $ntopin Variables")
    sort!(mlist, lt = (x,y)->abs(x[2]) > abs(y[2]))
    for k=1:ntopin
        i, m = mlist[k]
        setpinned!(g.vnodes[i], 1-2signbit(m))
    end
end

# r = fraction of N to free
function free_most_frustated!(g::FactorGraphKSAT, r::Float64 = 0.01)
    mlist = Vector{Tuple{Int,Float64}}()
    npin = numpinned(g)
    sizehint!(mlist, npin)
    for (i,v) in enumerate(g.vnodes)
        if ispinned(v)
            σ = v.pinned
            v.pinned = 0. # == setfree!(v)
            push!(mlist, (i, σ*mag(v)))
            v.pinned = σ
        end
    end

    ntofree = min(ceil(Int, r*g.N), length(mlist))
    println("# Freeing $ntofree Variables")
    sort!(mlist, lt = (x,y)->x[2] < y[2])
    for k=1:ntofree
        i, m = mlist[k]
        setfree!(g.vnodes[i])
    end
end

function update!(v::Var, r::Float64 = 0., tγ::Float64 = 0.)
    #TODO check del denominatore=0
    ispinned(v) && return 0.
    @extract v ζ̂listp ζ̂listm ζlistp ζlistm
    Δ = 0.
    ### compute total fields
    πp, πm, fπp, fπm = πpm(v)

    @assert πm >= 0 && πp >= 0

    ### compute cavity fields
    @inbounds for i=1:degp(v)
        πpi = (πp / ζ̂listp[i])*fπm
        πpi /= (πpi*fπm + πm*fπp/(1-ζ̂listp[i]))
        old = ζlistp[i][]
        ζlistp[i][] = πpi
        Δ = max(Δ, abs(πpi- old))
    end

    @inbounds for i=1:degm(v)
        πmi = (πm / ζ̂listm[i])*fπp
        πmi /= (πp*fπm/(1-ζ̂listm[i]) + πmi*fπp)
        old = ζlistm[i][]
        ζlistm[i][] = πmi
        Δ = max(Δ, abs(πmi- old))
    end

    ###############

    if tγ == 0.
        #### reinforcement ######
        #TODO togliere l'if simmetrizzando
        # mi sembra che moltiplicare tutto per πm^r o πp^r  sia lecito, e quindi resterebbe
        # v.ζreinfp = πp^r
        # v.ζreinfm = πm^r
        #
        # In ogni caso da rincontrollare.
        # stesso discorso vale per lo pseudo-reinforcement più sotto.
        # Perché l'ho messo così? problemi numerici?
        if πp < πm
            v.ζreinfp = (πp/πm)^r
            v.ζreinfm = 1
        else
            v.ζreinfm = (πm/πp)^r
            v.ζreinfp = 1
        end
    else
        #### pseudo-reinforcement ######
        πpR = πp / v.ζreinfp
        πmR = πm / v.ζreinfm
        mγ = (πpR-πmR) / (πpR+πmR) * tγ
        pp = (1+mγ)^r
        mm = (1-mγ)^r
        mR = tγ * (pp-mm) / (pp+mm)
        if πp < πm
            v.ζreinfp = (1 + mR) / (1 - mR)
            v.ζreinfm = 1
        else
            v.ζreinfm = (1 - mR) / (1 + mR)
            v.ζreinfp = 1
        end
    end
    Δ
end

function oneBPiter!(g::FactorGraphKSAT, r::Float64=0., tγ::Float64=0.)
    Δ = 0.

    for a=randperm(g.M)
        update!(g.fnodes[a])
    end

    for i=randperm(g.N)
        d = update!(g.vnodes[i], r, tγ)
        Δ = max(Δ, d)
    end

    Δ
end

function update_reinforcement!(reinfpar::ReinfParams)
    # @show reinfpar

    if reinfpar.wait_count < 10
        reinfpar.wait_count += 1
    else
        if reinfpar.γ == 0.
            reinfpar.r = 1 - (1-reinfpar.r) * (1-reinfpar.rstep)
        else
            reinfpar.r *= 1 + reinfpar.rstep
            reinfpar.γ += reinfpar.γstep
            reinfpar.tγ = tanh(reinfpar.γ)
            # @show reinfpar
        end
    end
    # @show reinfpar
end

getσ(mags::Vector) = Int[1-2signbit(m) for m in mags]

function converge!(g::FactorGraph; maxiters::Int = 100, ϵ::Float64=1e-5
        , reinfpar=ReinfParams(), altsolv=false, altconv=true)
    iters=0
    infotime = 10
    while iters<maxiters
        iters += 1
        iters % infotime == 0 && print("it=$iters ... ")
        Δ = oneBPiter!(g, reinfpar.r, reinfpar.tγ)
        E = energy(g)
        fp = numpinned(g) / g.N
        iters % infotime == 0 && @printf("r=%.3f γ=%.3f ρ_pin=%f\t  E=%d   \tΔ=%f \n",reinfpar.r,  reinfpar.γ, fp, E, Δ)
        # σ_noreinf = getσ(mags_noreinf(g))
        # E_noreinf = energy(g.cnf, σ_noreinf)
        # @printf("r=%.3f γ=%.3f \t  E=%d   \tE_noreinf=%d   Δ=%f \n",reinfpar.r, reinfpar.γ, E, E_noreinf, Δ)
        # println(mags(g)[1:10])
        update_reinforcement!(reinfpar)
        if altsolv && E == 0
            println("Found Solution!")
            break
        end
        if altconv && iters > 5 &&  Δ < ϵ
            println("Converged!")
            break
        end
    end
    return iters
end

energy(g::FactorGraphKSAT) = energy(g.cnf, getσ(mags(g)))

function πpm(v::Var)
    @extract v ζ̂listp ζ̂listm
    πp = 1.
    for ζ̂ in ζ̂listp
        πp *= ζ̂
    end
    πm = 1.
    for ζ̂ in ζ̂listm
        πm *= ζ̂
    end
    fπp = 1.
    for ζ̂ in ζ̂listp
        fπp *= 1 - ζ̂                                                       
    end
    fπm = 1.
    for ζ̂ in ζ̂listm
        fπm *= 1 - ζ̂
    end
    # (nzp > 0 && nzm > 0) && exit("contradiction")
    πp *= v.ζ̂reinfp
    πm *= v.ζ̂reinfm

    return πp, πm, fπp, fπm
end

function mag(v::Var)
    ispinned(v) && return float(v.pinned)
    πp, πm, fπp, fπm = πpm(v)
    m = (πp - πm) / (πm + πp)
    @assert isfinite(m)
    return m
end

function mag_noreinf(v::Var)
    ispinned(v) && return float(v.pinned)
    πp, πm, fπp, fπm = πpm(v)
    πp /= v.ζreinfp
    πm /= v.ζreinfm
    m = (πp - πm) / (πm + πp)
    # @assert isfinite(m)
    return m
end

mags(g::FactorGraph) = Float64[mag(v) for v in g.vnodes]
mags_noreinf(g::FactorGraphKSAT) = Float64[mag_noreinf(v) for v in g.vnodes]

"""
    solve(cnf::CNF; maxiters = 5000, ϵ::Float64 = 1e-4,
                method = :reinforcement, #[:reinforcement, :decimation]
                r::Float64 = 0., rstep::Float64= 0.001,
                γ::Float64 = 0., γstep::Float64=0.,
                γmax = 0.5,
                altsolv::Bool = true,
                seed::Int = -1)

Try to solve an instance of boolean satisfiability problem.
"""
function solve(cnf::CNF; maxiters = 5000, ϵ::Float64 = 1e-4,
                method = :reinforcement, #[:reinforcement, :decimation]
                r::Float64 = 0., rstep::Float64= 0.001,
                γ::Float64 = 0., γstep::Float64=0.,
                γmax = 0.5,
                altsolv::Bool = true,
                seed::Int = -1)
    seed > 0 && Random.seed!(seed)
    g = FactorGraphKSAT(cnf)
    initrand!(g)
    E = -1
    if method == :reinforcement
        reinfpar = ReinfParams(r, rstep, γ, γstep)
        converge!(g, maxiters=maxiters, ϵ=ϵ, reinfpar=reinfpar, altsolv=altsolv)
    elseif method == :decimation
        converge!(g, maxiters=maxiters, ϵ=ϵ, altsolv=altsolv)
        while true
            pin_most_biased!(g, r) # r=frac fixed , γ=frac freed
            free_most_frustated!(g, γ) # r=frac fixed , γ=frac freed
            converge!(g, maxiters=maxiters, ϵ=ϵ, altsolv=altsolv)
            E = energy(g)
            numdec = numpinned(g)
            (E == 0 || numdec == g.N) && break
        end
    elseif method == :converge
        # converge than update reinforcement
        reinf = ReinfParams(r, rstep, γ, γstep)
        reinf.wait_count = 100
        while reinf.γ < γmax
            iters = converge!(g, maxiters=maxiters, ϵ=ϵ,
                reinfpar=ReinfParams(reinf.r, 0., reinf.γ, 0.),
                altconv=true, altsolv=true)

            update_reinforcement!(reinf)
            # iters == maxiters && energy(g) != 0 && break
            energy(g) == 0 && break
        end
    else
        error("invalid method ", method)
    end

    return getσ(mags(g))
end
