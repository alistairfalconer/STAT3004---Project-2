using StatsBase, Distributions, Random, LinearAlgebra, Plots, Statistics
gr()

function doobGillespie(T,Q,initProb)
    n = size(Q)[1]
    Pjump  = (Q-diagm(0 => diag(Q)))./-diag(Q)
    lamVec = abs.(-diag(Q))
    #lamVec[1] = 0.0
    state  = sample(1:n,weights(initProb))
    sojournTime = rand(Exponential(1/lamVec[state]))
    t = 0.0
    while t + sojournTime < T
        t += sojournTime
        state = sample(1:n,weights(Pjump[state,:]))
        sojournTime = rand(Exponential(1/lamVec[state]))
    end
    return state
end

function makeQSI(X0::Int,Y0::Int,rSI::Float64)
    N = X0 + Y0
    Q = zeros(X0+1,X0+1)
    for i in 2:X0+1 #first row is absorbing
        Q[i,i-1] = (i-1)*(N-i+1)*rSI
        Q[i,i] = -(i-1)*(N-i+1)*rSI
    end
    InitState = zeros(X0+1,1)
    InitState[X0+1] = 1
    return Q,InitState
end

function Q3b(X0::Int,Y0::Int,rSI::Float64,tSave::Array{Float64,1};M=1e6)
    Q,InitState = makeQSI(X0,Y0,rSI)
    pnSave = []
    for t in tSave
        TrajList = []
        for _ in 1:M
            push!(TrajList,doobGillespie(t,Q,InitState)-1)
        end
        pn = [sum(TrajList .== i)/M for i in 0:X0]
        push!(pnSave,pn)
    end
    p = plot()
    for i in 1:X0+1
        p = scatter!(tSave,[pnSave[t][i] for t in 1:length(tSave)],label="p$(i-1)")
    end
    xlabel!("Time")
    ylabel!("Probability")
    return p
end


function Q3(X0::Int,Y0::Int,rSI::Float64,tSave::Array{Float64,1};M=1e6)
    N = X0 + Y0#population
    Q,InitState = makeQSI(X0,Y0,rSI)
    pnSaveMatrix = []
    pnSaveMonteCarlo = []
    for t in tSave
        TrajList = []
        for _ in 1:M
            push!(TrajList,doobGillespie(t,Q,InitState)-1)
        end
        pn1 = [sum(TrajList .== i)/M for i in 0:X0]
        push!(pnSaveMonteCarlo,pn1)
        P = exp(Q*t)
        pn2 = [(InitState'*P)[i] for i in 1:X0+1]
        push!(pnSaveMatrix,pn2)
    end
    p = plot(legend=:outertopleft)
    col = distinguishable_colors(X0+1,RGB(255/255,127/255,80/255))
    for i in 1:X0+1
        p = scatter!(tSave,[pnSaveMonteCarlo[t][i] for t in 1:length(tSave)],label="p$(i-1)",color=col[i])
        p = plot!(tSave,[pnSaveMatrix[t][i] for t in 1:length(tSave)],#=label="p$(i-1)",=#label=:none,color=col[i])
    end
    xlabel!("Time")
    ylabel!("Probability")
    return p
end


function Q3a(X0::Int,Y0::Int,rSI::Float64,tSave::Array{Float64,1})
    N = X0 + Y0#population
    Q,InitState = makeQSI(X0,Y0,rSI)
    PSave = []
    pnSave = []
    for t in tSave
        P = exp(Q*t)
        push!(PSave,P)
        pn = [(InitState'*P)[i] for i in 1:X0+1]
        push!(pnSave,pn)
    end
    p = plot()
    for i in 1:X0+1
        plot!(tSave,[pnSave[t][i] for t in 1:length(tSave)],label="p$(i-1)")
    end
    xlabel!("Time")
    ylabel!("Probability")
    return p
end
