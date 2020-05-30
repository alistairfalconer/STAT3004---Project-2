using StatsBase, Distributions, Random, LinearAlgebra, Plots, Statistics
gr()

#Given a Q-matrix and and initial proability distribution, will generate a sample trajectory up to time T and return the final state. The state is returned by indexing a basis of the Q matrix
function doobGillespie(T,Q,initProb)
    n = size(Q)[1]
    Pjump  = (Q-diagm(0 => diag(Q)))./-diag(Q) #remonve the diagonal from Q
    lamVec = abs.(-diag(Q)) #parameter for finding time until jump
    state  = sample(1:n,weights(initProb)) #extract current state from initial probabilities
    sojournTime = rand(Exponential(1/lamVec[state])) #find time until next transition
    t = 0.0
    while t + sojournTime < T #Find transitions before time T
        t += sojournTime
        state = sample(1:n,weights(Pjump[state,:])) #find next state
        sojournTime = rand(Exponential(1/lamVec[state])) #find time until next jump
    end
    return state
end

#generate the Q matrix for the SI model given parameters X0,Y0,rSI
function makeQSI(X0::Int,Y0::Int,rSI::Float64)
    N = X0 + Y0#total population
    Q = zeros(X0+1,X0+1)#initialise matrix
    for i in 2:X0+1 #first row is absorbing
        Q[i,i-1] = (i-1)*(N-i+1)*rSI#down transitions
        Q[i,i] = -(i-1)*(N-i+1)*rSI #negative component on diagonal
    end
    InitState = zeros(X0+1,1)
    InitState[X0+1] = 1#find initial state (known distribution)
    return Q,InitState
end

#assume one initial infective, find the expected time at which the i'th susceptible becomes infected
function FindEtj(rSI::Float64,N::Int;M=1e3,Error=false)
    X0 = N-1
    Y0 = 1
    iVals = X0:-1:2
    tVals = []
    tValsCalc = []
    for i in iVals
        itVals = []
        for _ in 1:M
            time = 0.0
            for idash in X0:-1:i
                time += rand(Exponential(1/(rSI*idash*(N-idash)))) #EM - 3.1.3
            end
            push!(itVals,time)

        end
        push!(tValsCalc,1/(rSI*N)*log( (X0*(N+1-i))/(Y0*(i-1)) ))#deterministic approximation
        push!(tVals,sum(itVals)/M) #MC simulation
    end
    if Error #Offset between calculations
        p = scatter(iVals,tVals-tValsCalc,label=:none)
        xlabel!("Infection Number")
        ylabel!("Offset")
        return p
    end
    scatter(iVals,tVals,label="Monte Carlo")
    plot!(iVals,tValsCalc, label="Deterministic Approximation")
    xlabel!("Infection Number")
    ylabel!("Expected Time")
end

function Q2bCollect(rSIList::Array,NList::Array;Error=false,M=1e3)
     pltList = []
     for rSi in rSIList
         for N in NList
             push!(pltList,FindEtj(rSi,N,Error=Error,M=M))
             plot!(legend=:none,title="N=$N, rSI=$rSi")
         end
     end
     plot(pltList...)
end

function Q2cTrace(rSI::Float64,N::Int;Y0=1::Int,ctup=1.0::Float64,ctdown=0.1::Float64,Loud=false)
    X0 = N - Y0
    iVals = X0-1:-1:0
    tVals = []
    time = 0.0
    LastSwitchTime = 0.0
    Distancing = false #switches to indicate if distancing is currently in place
    for i in iVals
        if i%100 == 99
            Distancing = true
            LastSwitchTime = time
        end
        if Distancing
            dtime = rand(Exponential(1/(ctdown*rSI*(i+1)*(N-(i+1)))))#time to transition from i+1 -> i susceptibles
            if time + dtime < LastSwitchTime + 20.0
                time += dtime
            else
                Distancing = false
                time = LastSwitchTime + 20.0 + rand(Exponential(1/(ctup*rSI*(i+1)*(N-(i+1)))))#use memoryless property of exponential
            end
        else
            time += rand(Exponential(1/(ctup*rSI*(i+1)*(N-(i+1)))))
        end
        push!(tVals,time)
    end
    if !Loud
        return tVals,iVals
    end
    scatter(tVals,iVals,label=:none)
    xlabel!("Infection Number")
    ylabel!("Expected Time")

end

function Q2c(M::Int,rSI::Float64,N::Int;Y0=1::Int,ctup=1.0::Float64,ctdown=0.1::Float64)
    p = plot()
    for _ in 1:M
        plot!(Q2cTrace(rSI,N,Y0=Y0,ctup=ctup,ctdown=ctdown,Loud=false)...,legend=:none)
    end
    xlabel!("Time")
    ylabel!("Number of Susceptibles")
end

#=
function Q2c(M::Int,rSI::Float64,N::Int;Y0=1::Int,ct=1.0::Float64,ctdown=0.1::Float64)
    X0 = N - Y0
    iVals = X0-1:-1:0
    tVals = [[] for _ in 1:M]
    for m in 1:M
        time = 0.0
        for i in iVals
            time += rand(Exponential(1/(rSI*(i+1)*(N-(i+1)))))#time to transition from i+1 -> i susceptibles
            push!(tVals[m],time)
        end
    end
    p = plot()
    for m in 1:M
        plot!(tVals[m],iVals,label=:none)
    end
    xlabel!("Time")
    ylabel!("Number of Susceptibles")
    return p
end
=#
