using StatsBase, Distributions, Random, LinearAlgebra, Plots, Statistics, LaTeXStrings
gr()

function Q6b(X0::Int,Y0::Int,rSI::Float64,rIR::Float64,XS::Int,ctdown::Float64;T=1000.0::Float64)
    tVals = range(0,stop=T,length=50000)
    dt = tVals[2] - tVals[1]
    pVals = zeros(X0+1,X0+Y0+1,length(tVals))
    pVals[X0+1,Y0+1,1] = 1.0
    for t in 2:length(tVals)
        pVals[X0+1,Y0+1,t] = pVals[X0+1,Y0+1,t-1] - dt*Y0*(1.0*rSI*X0+rIR)*pVals[X0+1,Y0+1,t-1]
        for i in 1:X0+1
            for j in 1:X0+Y0+1
                if i-1 <= XS
                    ctfrom = ctdown
                else
                    ctfrom = 1.0
                end
                if i-1 <= XS-1
                    ctto = ctdown
                else
                    ctto = 1.0
                end
                if (i+j > X0+Y0+2) || (i,j)==(X0+1,Y0+1)
                    continue
                end
                pVals[i,j,t] = pVals[i,j,t-1]
                if i<X0+1 && j>1
                    pVals[i,j,t] += dt*ctto*rSI*(i)*(j-2)*pVals[i+1,j-1,t-1]
                end
                pVals[i,j,t] -= dt*(j-1)*(ctfrom*rSI*(i-1)+rIR)*pVals[i,j,t-1]
                if j<X0+Y0+1
                    pVals[i,j,t] += dt*rIR*(j)*pVals[i,j+1,t-1]
                end
            end
        end
    end
    #display(pVals[:,:,end])
    println("SUM:",sum([sum(pVals[:,j,end]) for j in 1:X0+Y0+1]))
    #=
    EIT = 0.0
    for i in 1:X0+1
        for j in 1:X0+Y0+1
            EIT += (j-1)*pVals[i,j,end]
        end
    end
    return EIT
    =#
    EPn = Float64(X0)
    for i in 1:X0+1
        for j in 1:X0+Y0+1
            EPn -= (i-1)*pVals[i,j,end]
        end
    end
    return EPn
end


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

function makeQSIR(X0::Int,Y0::Int,rSI::Float64,rIR::Float64,ct::Float64)
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
            Q[from,Setats[ifrom-1,jfrom+1]] = ct*rSI*ifrom*jfrom
            Q[from,from] -= ct*rSI*ifrom*jfrom
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


function Q6dTime(X0::Int,Y0::Int,rSI::Float64,rIR::Float64,XS ::Int,ctdown::Float64;T=500.0::Float64,M=1e2)
    if XS >= X0
        println("XS should be smaller than X0")
        return
    end
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
    Q,InitState = makeQSIR(X0,Y0,rSI,rIR,1.0)
    n = size(Q)[1]
    Pjump  = (Q-diagm(0 => diag(Q)))./-diag(Q)
    lamVec = abs.(-diag(Q))
    SwitchStates = []
    Extinct = []
    for _ in 1:M
        state  = sample(1:n,weights(InitState))
        while true
            state = sample(1:n,weights(Pjump[state,:]))
            InterpretState = States[state]
            if InterpretState[1] == XS
                #println("broken")
                break
            end
            if InterpretState[2] == 0
                push!(Extinct,X0 - InterpretState[1])
                break
            end
        end
        push!(SwitchStates,state)
        #push!(SwitchStates,States[state])
    end

    FinalInfections =zeros(length(SwitchStates))
    for m in 1:Int(length(SwitchStates))
        xs,ys = States[SwitchStates[m]]
        Q,InitState = makeQSIR(xs,ys,rSI,rIR,ctdown)
        FinalInfections[m] = X0 - States[doobGillespie(T,Q,InitState)][1]
    end
    push!(FinalInfections,Extinct...)
    return Statistics.mean(FinalInfections)
end

function Q6d(X0::Int, Y0::Int, rSI::Float64, rIR::Float64, ctdownlist::Array{Float64,1}; T=1500.0::Float64, M=1e2,dX=3::Int)
    p = plot(legend=:bottomleft)
    upper = Q6dTime(X0,Y0,rSI,rIR,X0-1,1.0,M=1e2)
    hline!([upper],label="Without Distancing",color=:black)
    hline!([upper/2],label="Half of Without Distancing",color=:black,linestyle=:dash)
    col = distinguishable_colors(length(ctdownlist),RGB(255/255,127/255,80/255))
    #for ctdown in ctdownlist
    for c in 1:length(ctdownlist)
        ctdown = ctdownlist[c]
        scatter!(1:dX:X0-1,[Q6dTime(X0,Y0,rSI,rIR,XS,ctdown,T=T,M=M) for XS in 0:dX:X0-1],label="c_t=$ctdown",color=col[c])
        plot!(1:dX:X0-1,[Q6b(X0,Y0,rSI,rIR,XS,ctdown) for XS in 1:dX:X0-1],label=:none,color=col[c])
    end
    xlabel!("XS")
    ylabel!("Number of new infectives")
end

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

function Q6HowLong(X0::Int,Y0::Int,rSI::Float64,rIR::Float64,XS::Int,ctdown::Float64;T=1000.0)
    tVals = range(0,stop=T,length=50000)
    dt = tVals[2] - tVals[1]
    pVals = zeros(X0+1,X0+Y0+1,length(tVals))
    pVals[X0+1,Y0+1,1] = 1.0
    for t in 2:length(tVals)
        pVals[X0+1,Y0+1,t] = pVals[X0+1,Y0+1,t-1] - dt*Y0*(1.0*rSI*X0+rIR)*pVals[X0+1,Y0+1,t-1]
        for i in 1:X0+1
            for j in 1:X0+Y0+1
                if i-1 <= XS
                    ctfrom = ctdown
                else
                    ctfrom = 1.0
                end
                if i-1 <= XS-1
                    ctto = ctdown
                else
                    ctto = 1.0
                end
                if (i+j > X0+Y0+2) || (i,j)==(X0+1,Y0+1)
                    continue
                end
                pVals[i,j,t] = pVals[i,j,t-1]
                if i<X0+1 && j>1
                    pVals[i,j,t] += dt*ctto*rSI*(i)*(j-2)*pVals[i+1,j-1,t-1]
                end
                pVals[i,j,t] -= dt*(j-1)*(ctfrom*rSI*(i-1)+rIR)*pVals[i,j,t-1]
                if j<X0+Y0+1
                    pVals[i,j,t] += dt*rIR*(j)*pVals[i,j+1,t-1]
                end
            end
        end
    end
    println("SUM:",sum([sum(pVals[:,j,end]) for j in 1:X0+Y0+1]))
    PExtinctList = zeros(length(tVals))
    for t in 1:length(tVals)
        PExtinct = 0.0
        for i in 1:X0+1
            PExtinct += pVals[i,1,t]
        end
        PExtinctList[t] = PExtinct
    end
    return PExtinctList
end

function Q6LookExtinctions(X0::Int,Y0::Int,rSI::Float64,rIR::Float64,XSList::Array,ctdownList::Array)
    if !(length(XSList) == length(ctdownList))
        println("They're not the same size")
        return
    end
    tVals = range(0,stop=1000.0,length=50000)
    Vals=[]
    p = plot(legend=:outertopleft)
    for i in 1:length(XSList)
        XS,ctdown = XSList[i],ctdownList[i]
        val = Q6HowLong(X0,Y0,rSI,rIR,XS,ctdown)
        plot!(tVals,val,label="XS=$XS, ct=$ctdown")
        push!(Vals,val)
        valdash = []
        for t in 2:length(tVals)
            push!(valdash,(val[t]-val[t-1]))#/(tVals[t]-tVals[t-1]))
        end
        texpect = sum([valdash[i]*tVals[i] for i in 1:length(tVals)-1])
        #texpect = tVals[2:end]''*valdash
        println("With XS=$XS and ct=$ctdown, the expected time until the infection goes extinct is $texpect")
    end
    for t in 1:length(tVals)
        probs = [Vals[i][t] for i in 1:length(Vals)]
        if minimum(probs) > 1-1e-3
            xaxis!((0,tVals[t]))
            break
        end
    end
    xlabel!("Time")
    ylabel!("Extinction Probability")
    return p
end
