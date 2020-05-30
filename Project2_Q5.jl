using StatsBase, Distributions, Random, LinearAlgebra, Plots, Statistics, LaTeXStrings
pgfplotsx()

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
#=
function Q5b(X0::Int,Y0::Int,rSI::Float64,rIR::Float64;T=5.0::Float64,a=0.0::Float64,b=0.0::Float64)
    tVals = range(0,stop=T,length=1000)
    dt = tVals[2] - tVals[1]
    pVals = zeros(X0+1,X0+Y0+1,length(tVals))
    pVals[X0+1,Y0+1,1] = 1.0
    for t in 2:length(tVals)
        pVals[X0+1,Y0+1,t] = pVals[X0+1,Y0+1,t] + dt*-Y0*(rSI*X0+rIR)*pVals[X0+1,Y0+1,t-1]
        for i in 1:X0+1
            for j in 1:X0+Y0+1
                if (i+j > X0+Y0+2) || (i,j)==(X0+1,Y0+1)
                    continue
                end
                pVals[i,j,t] = pVals[i,j,t-1]
                if i<X0+1 && j>1
                    pVals[i,j,t] += dt*rSI*(i+1)*(j-1)*pVals[i+1,j-1,t-1]
                end
                pVals[i,j,t] += dt*-j*(rSI*i+rIR)*pVals[i,j,t-1]
                if j<X0+Y0+1
                    pVals[i,j,t] += dt*rIR*(j+1)*pVals[i,j+1,t-1]
                end
            end
        end
    end
    #display(pVals[:,:,end])
    EIT = 0.0
    for i in 1:X0+1
        for j in 1:X0+Y0+1
            EIT += (j-1)*pVals[i,j,end]
        end
    end
    return EIT
end=#
function Q5b(X0::Int,Y0::Int,rSI::Float64,rIR::Float64;T=5.0::Float64,a=0.0::Float64,b=0.0::Float64)
    tVals = range(0,stop=T,length=10000)
    dt = tVals[2] - tVals[1]
    pVals = zeros(X0+1,X0+Y0+1,length(tVals))
    pVals[X0+1,Y0+1,1] = 1.0
    for t in 2:length(tVals)
        pVals[X0+1,Y0+1,t] = pVals[X0+1,Y0+1,t-1] - dt*Y0*(rSI*X0+rIR)*pVals[X0+1,Y0+1,t-1]
        for i in 1:X0+1
            for j in 1:X0+Y0+1
                if (i+j > X0+Y0+2) || (i,j)==(X0+1,Y0+1)
                    continue
                end
                pVals[i,j,t] = pVals[i,j,t-1]
                if i<X0+1 && j>1
                    pVals[i,j,t] += dt*rSI*(i)*(j-2)*pVals[i+1,j-1,t-1]
                end
                pVals[i,j,t] -= dt*(j-1)*(rSI*(i-1)+rIR)*pVals[i,j,t-1]
                if j<X0+Y0+1
                    pVals[i,j,t] += dt*rIR*(j)*pVals[i,j+1,t-1]
                end
            end
        end
    end
    #display(pVals[:,:,end])
    println("SUM:",sum([sum(pVals[:,j,end]) for j in 1:X0+Y0+1]))
    EIT = 0.0
    for i in 1:X0+1
        for j in 1:X0+Y0+1
            EIT += (j-1)*pVals[i,j,end]
        end
    end
    return EIT
end

function Q5bPlot(X0::Int,Y0::Int,rSI::Float64,rIR::Float64,TList::Array)
    p = scatter(TList,[Q5b(X0,Y0,rSI,rIR,T=T) for T in TList],legend=:none)
    xlabel!("Time")
    ylabel!(L" E[I_t]")
    return p
end

function makeQSIR(X0::Int,Y0::Int,rSI::Float64,rIR::Float64)
    N = X0 + Y0
    States = Dict()
    count = 1
    for i in 0:X0
        for j in 0:N-i
            States[count] = (i,j)
            count += 1
        end
    end
    Setats = Dict(States[k]=>k for k in keys(States))#reverse of States
    NStates = length(keys(States))
    Q = zeros(NStates,NStates)
    for from in 1:NStates
        ifrom,jfrom = States[from]
        #println(ifrom,jfrom)
        if ifrom > 0
            Q[from,Setats[ifrom-1,jfrom+1]] = rSI*ifrom*jfrom
            Q[from,from] -= rSI*ifrom*jfrom
        end
        if jfrom > 0
            Q[from,Setats[ifrom,jfrom-1]] = rIR*jfrom
            Q[from,from] -= rIR*jfrom
        end
    end
    InitState = zeros(NStates)
    InitState[Setats[X0,Y0]] = 1
    return Q,InitState
end

function Q5cTime(X0::Int,Y0::Int,rSI::Float64,rIR::Float64;T=5.0::Float64,a=0.0::Float64,b=0.0::Float64,M=1e2)
    N = X0 + Y0
    States = Dict()
    count = 1
    for i in 0:X0
        for j in 0:N-i
            States[count] = (i,j)
            count += 1
        end
    end
    Setats = Dict(States[k]=>k for k in keys(States))#reverse of States
    NStates = length(keys(States))
    Q,InitState = makeQSIR(X0,Y0,rSI,rIR)
    ITs = []
    for _ in 1:M
        FinalState = States[doobGillespie(T,Q,InitState)]
        push!(ITs,FinalState[2])
    end
    return Statistics.mean(ITs)
end

function Q5c(X0::Int,Y0::Int,rSI::Float64,rIR::Float64,TList::Array;a=0.0::Float64,b=0.0::Float64,M=1e4)
    EIt = []
    for T in TList
        push!(EIt,Q5cTime(X0,Y0,rSI,rIR,T=T,a=a,b=b,M=M))
    end
    plot(TList,[Q5b(X0,Y0,rSI,rIR;T=T) for T in TList],label="Forward Equations")
    scatter!(TList,EIt,label="Monte Carlo",legend=:outertopleft)
    xlabel!("Time")
    ylabel!("Expected Number of Infectives")
end

function Q5d(X0::Int,Y0::Int,rSI::Float64,rIR::Float64;T=500.0::Float64,M=1e2)
    N = X0 + Y0
    States = Dict()
    count = 1
    for i in 0:X0
        for j in 0:N-i
            States[count] = (i,j)
            count += 1
        end
    end
    Setats = Dict(States[k]=>k for k in keys(States))#reverse of States
    NStates = length(keys(States))
    Q,InitState = makeQSIR(X0,Y0,rSI,rIR)
    FinalInfections = []
    for _ in 1:M
        FinalState = States[doobGillespie(T,Q,InitState)]
        push!(FinalInfections,X0-FinalState[1])
    end
    return FinalInfections
end
