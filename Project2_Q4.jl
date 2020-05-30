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

function makeQSIS(X0::Int,Y0::Int,rSI::Float64,rIS::Float64)
    N = X0 + Y0
    Q = zeros(N+1,N+1)
    Q[1,1] = -rIS*N
    Q[1,2] = rIS*N
    for i in 2:N #first row is absorbing
        Q[i,i-1] = (i-1)*(N-i+1)*rSI
        Q[i,i] = -(i-1)*(N-i+1)*rSI-(N-i+1)*rIS
        Q[i,i+1] = (N-i+1)*rIS
    end
    InitState = zeros(X0+1,1)
    InitState[X0+1] = 1
    return Q,InitState
end

function Q4aSample(X0::Int,Y0::Int,rSI::Float64,rIS::Float64;M=1e1,Loud=false)
    Q,InitState = makeQSIS(X0,Y0,rSI,rIS)
    N = X0 + Y0
    stateList = []
    for m in 1:M
        push!(stateList,doobGillespie(15.0,Q,InitState)-1)
    end
    prob = sum(stateList .== N)/M
    if Loud
        println("There is a probability $prob that the infection has gone extinct at time 15")
    end
    return prob
end

function Q4a(X0::Int,Y0::Int,rSIList::Array{Float64,1},rIS::Float64;M=1e4)
    problist = []
    for rSI in rSIList
        push!(problist,Q4aSample(X0,Y0,rSI,rIS,M=M))
    end
    #scatter(rSIList,problist)
    #xlabel!("rSI")
    #ylabel!("Probability of extinction by time 15")
    #hline!([0.9],color=:black)
    ind = 0
    while problist[ind+1] > 0.9
        ind += 1
        if ind == length(rSIList)
            println("Out of indices")
            return
        end
    end
    newproblist = []
    newrSIList = collect(range(rSIList[ind]-0.02,stop=rSIList[ind+1]+0.02,length=10))
    for rSI in newrSIList
        push!(newproblist,Q4aSample(X0,Y0,rSI,rIS,M=M))
    end
    push!(problist,newproblist...)
    push!(rSIList,newrSIList...)
    for i in 1:length(problist)
        if problist[i]>0.9 && problist[i+1]<=0.9
            println("The final value of r_{SI} with a greater than 0.9 probability of having the extinction go extinct by time 15 is r_SI}=$(rSIList[i])")
        end
    end
    scatter(rSIList,problist,legend=:none)
    xlabel!(L"r_{SI}")
    ylabel!("Probability of extinction by time 15")
    hline!([0.9],color=:black)
end

function Q4c(X0::Int,Y0::Int,rIS::Float64,rSI::Float64,a::Float64)
    N = X0 + Y0
    pi = [1.0]
    for i in 0:N-1
        push!(pi,pi[end]*rIS*(N-i)/(rSI*(i+1)*(N-i-1)+a))
    end
    return pi./sum(pi)
end

function Q4d(X0::Int,Y0::Int,rIS::Float64,rSI::Float64,a::Float64)
    N = X0 + Y0
    pi = Q4c(X0,Y0,rIS,rSI,a)
    return [(N-i) for i in 0.0:N]'*pi
end
