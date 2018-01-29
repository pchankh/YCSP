
function rowLookIPHeuristic(H, T, n, N, SI, SP, SR, SO, Stot, contInStackToRetrieve, indicesToRetrieve, heightsInitial, posCraneInitial, costMove, costPreMove, costToGo, alpha, posteriorStacks, anteriorStacks, IOPoint, printSolver, gapOfMIP, limitOfTimeHeur, rowLook)

    neighborStacks = Dict{Int64,Array{Int64}}();
    AllneighborStacks = Array{Int64}(0);
    for s in SR
        neighborStacks[s] = union(collect(S*(max(realStack[s][1]-rowLook,1)-1)+1:S*min(realStack[s][1]+rowLook,R)),IOPoint[s]);
        append!(AllneighborStacks,neighborStacks[s]);
    end
    AllneighborStacks = unique(AllneighborStacks);
    SRLoc = SR;
    SBLoc = intersect(SB,AllneighborStacks);
    SOLoc = intersect(SO,AllneighborStacks);
    SILoc = intersect(SI,AllneighborStacks);
    SPLoc = intersect(SP,AllneighborStacks);
    StotLoc = intersect(Stot,AllneighborStacks);

    anteriorStacksLoc = Dict{Int64,Array{Int64}}();
    posteriorStacksLoc = Dict{Int64,Array{Int64}}();
    for s in SPLoc
        anteriorStacksLoc[s] = intersect(anteriorStacks[s],AllneighborStacks);
        posteriorStacksLoc[s] = intersect(posteriorStacks[s],AllneighborStacks);
    end
    for s in SOLoc
        anteriorStacksLoc[s] = intersect(anteriorStacks[s],AllneighborStacks);
    end
    SODiff = setdiff(SO,SOLoc);
    fixedHeight = Dict{Array{Int64},Int64}();
    for s in SODiff
        for h = 0:H
            if heightsInitial[realStack[s][1],realStack[s][2]] == h
                fixedHeight[[s,h]] = 1;
            else
                fixedHeight[[s,h]] = 0;
            end
        end
    end

    IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTimeHeur));
    x = @variable(IPModel, [s in SPLoc,r in posteriorStacksLoc[s],1:T], Bin);
    d = @variable(IPModel, [StotLoc,SPLoc,2:T], Bin);
    dInit = @variable(IPModel, [SPLoc], Bin);
    heights = @variable(IPModel, [SRLoc,0:H,1:T], Bin);
    finalh = @variable(IPModel, [SOLoc,0:H], Bin);
    w = @variable(IPModel, [1:n,1:T], Bin);


    @objective(IPModel, Min, sum(costMove[[s,r]]*x[s,r,t] for s in SPLoc for r in posteriorStacksLoc[s] for t=1:T) + sum(costPreMove[[s,r]]*d[s,r,t] for s in StotLoc for r in SPLoc for t = 2:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SPLoc) + costToGo * sum( alpha[h] * (sum(heights[s,h,T] for s in SRLoc) + sum(finalh[s,h] for s in SOLoc)) for h = 2:H));

    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################


    for t = 2:T
        ## Uniqueness of pre-move
        @constraint(IPModel, sum(d[s,r,t] for s in StotLoc for r in SPLoc) == 1);
    end
    for t = 1:T
        ## Uniqueness of move
        @constraint(IPModel, sum(x[s,r,t] for s in SPLoc for r in posteriorStacksLoc[s]) == 1);
    end

    @constraint(IPModel, sum(dInit[r] for r in SPLoc) == 1);

    for s in SPLoc
        @constraint(IPModel, sum(x[s,r,1] for r in posteriorStacksLoc[s]) - dInit[s] == 0);
        for t = 2:T
            ## Conservation of flow (1)
            @constraint(IPModel, sum(x[s,r,t] for r in posteriorStacksLoc[s]) - sum(d[r,s,t] for r in StotLoc) == 0);
        end
    end

    for s in StotLoc
        for t = 2:T
            ## Conservation of flow (2)
            @constraint(IPModel, sum(d[s,r,t] for r in SPLoc) - sum(x[r,s,t-1] for r in anteriorStacksLoc[s]) == 0);
        end
    end

    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################

    ## Number of stacks
    @constraint(IPModel, sum(x[s,r,t] for s in SILoc for r in posteriorStacksLoc[s] for t = 1:T) == N - n);

    for m = 1:n
        ## Number of deliveries
        @constraint(IPModel, sum(w[m,t] for t = 1:T) == 1);
    end

    for s in SRLoc
        for t = 1:T
            ## Relation between variables x and w
            @constraint(IPModel, x[s,IOPoint[s],t] - sum(w[m,t] for m in contInStackToRetrieve[s]) == 0);
        end
    end

    for m = 1:n
        if heightsInitial[realStack[indicesToRetrieve[m][1]][1],realStack[indicesToRetrieve[m][1]][2]] > indicesToRetrieve[m][2]
            ## Relation between variables w and θ
            @constraint(IPModel, w[m,1] <= 0);
        end
        for t = 2:T
            ## Relation between variables w and θ
            @constraint(IPModel, w[m,t] - heights[indicesToRetrieve[m][1],indicesToRetrieve[m][2],t-1] <= 0);
        end
        for t = 1:T
            ## Relation between variables w and θ (2)
            @constraint(IPModel, sum(heights[indicesToRetrieve[m][1],h,t] for h = 0:indicesToRetrieve[m][2]-1) - sum(w[m,u] for u = 1:t) <= 0);
        end
    end

    #####################################################################
    ########################### Heights Update ##########################
    #####################################################################

    for s in SRLoc
        ## Initial update
        @constraint(IPModel, sum(h * heights[s,h,1] for h = 1:H) - sum(x[r,s,1] for r in anteriorStacksLoc[s]) + sum(x[s,r,1] for r in posteriorStacksLoc[s]) == heightsInitial[realStack[s][1],realStack[s][2]]);
        for t = 2:T
            ## Further updates
            @constraint(IPModel, sum(h * heights[s,h,t] for h = 1:H) - sum(h * heights[s,h,t-1] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacksLoc[s]) + sum(x[s,r,t] for r in posteriorStacksLoc[s]) == 0);
        end
        for t = 1:T
            ## Uniqueness of height
            @constraint(IPModel, sum(heights[s,h,t] for h = 0:H) == 1);
        end
    end

    for s in SOLoc
        ## Final heights
        @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacksLoc[s] for t = 1:T) == heightsInitial[realStack[s][1],realStack[s][2]]);
        ## Uniqueness of final height
        @constraint(IPModel, sum(finalh[s,h] for h = 0:H) == 1);
    end

    #####################################################################
    ################################ Cuts ###############################
    #####################################################################

    for s in SILoc
        if length(anteriorStacksLoc[s]) == 0
            for r in SPLoc
                for t = 2:T
                    @constraint(IPModel, d[s,r,t] == 0);
                end
            end
        end
    end

    for s in SRLoc
        for m in contInStackToRetrieve[s]
            for t = 1:T
                @constraint(IPModel, sum(x[r,s,t] for r in anteriorStacksLoc[s]) - sum(w[m,u] for u = 1:t) <= 0);
            end
        end
    end

    println(MathProgBase.numvar(IPModel) , " variables");
    println(MathProgBase.numconstr(IPModel), " constraints");
    println("Solving ....");

    tic();
    status = solve(IPModel);
    println("Solved !");
    timeToSolveHeur = toc();

    objLoc = getobjectivevalue(IPModel);
    for s in SODiff
        for h = 0:H
            if fixedHeight[[s,h]] == 1
                objLoc += costToGo * alpha[h];
            end
        end
    end
    println("The objective value of the heuristic is: ", objLoc);

    XLoc = getvalue(x);
    DLoc = getvalue(d);
    DinitLoc = getvalue(dInit);
    newHeightsLoc = getvalue(heights);
    finalHeightsLoc = getvalue(finalh);
    WLoc = getvalue(w);

    XHeur = Dict{Array{Int64},Int64}();
    for s in SP
        for r in posteriorStacks[s]
            for t = 1:T
                if s in SPLoc && r in posteriorStacksLoc[s] && XLoc[s,r,t] > 0.9
                    XHeur[[s,r,t]] = 1;
                else
                    XHeur[[s,r,t]] = 0;
                end
            end
        end
    end
    DHeur = Dict{Array{Int64},Int64}();
    for s in Stot
        for r in SP
            for t = 1:T
                if t == 1
                    if s == posCraneInitial && r in SPLoc && DinitLoc[r] > 0.9
                        DHeur[[s,r,t]] = 1;
                    else
                        DHeur[[s,r,t]] = 0;
                    end
                else
                    if s in StotLoc && r in SPLoc && DLoc[s,r,t] > 0.9
                        DHeur[[s,r,t]] = 1;
                    else
                        DHeur[[s,r,t]] = 0;
                    end
                end
            end
        end
    end
    newHeightsHeur = Dict{Array{Int64},Int64}();
    for s in SR
        for h = 0:H
            for t = 1:T
                if newHeightsLoc[s,h,t] > 0.9
                    newHeightsHeur[[s,h,t]] = 1;
                else
                    newHeightsHeur[[s,h,t]] = 0;
                end
            end
        end
    end
    finalHeightsHeur = Dict{Array{Int64},Int64}();
    for s in SO
        for h = 0:H
            if s in SOLoc && finalHeightsLoc[s,h] > 0.9
                finalHeightsHeur[[s,h]] = 1;
            elseif s in SOLoc
                finalHeightsHeur[[s,h]] = 0;
            else
                finalHeightsHeur[[s,h]] = fixedHeight[[s,h]];
            end
        end
    end
    WHeur = Dict{Array{Int64},Int64}();
    for m = 1:n
        for t = 1:T
            if WLoc[m,t] > 0.9
                WHeur[[m,t]] = 1;
            else
                WHeur[[m,t]] = 0;
            end
        end
    end

    return (XHeur,DHeur,newHeightsHeur,finalHeightsHeur,WHeur,timeToSolveHeur);
end

function MIP(H, T, n, N, SI, SP, SR, SO, Stot, contInStackToRetrieve, indicesToRetrieve, heightsInitial, posCraneInitial, costMove, costPreMove, costToGo, alpha, posteriorStacks, anteriorStacks, IOPoint, printSolver, gapOfMIP, limitOfTime, XHeur, DHeur, newHeightsHeur, finalHeightsHeur, WHeur)

    IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));
    x = @variable(IPModel, [s in SP,r in posteriorStacks[s],1:T], Bin);
    d = @variable(IPModel, [Stot,SP,1:T], Bin);
    heights = @variable(IPModel, [SR,0:H,1:T], Bin);
    finalh = @variable(IPModel, [SO,0:H], Bin);
    w = @variable(IPModel, [1:n,1:T], Bin);

    #####################################################################
    ######################## Set Initial Solution #######################
    #####################################################################

    for s in SP
        for r in posteriorStacks[s]
            for t = 1:T
                if XHeur[[s,r,t]] == 1
                    setvalue(x[s,r,t], 1);
                else
                    setvalue(x[s,r,t], 0);
                end
            end
        end
    end
    for s in Stot
        for r in SP
            for t = 1:T
                if DHeur[[s,r,t]] == 1
                    setvalue(d[s,r,t], 1);
                else
                    setvalue(d[s,r,t], 0);
                end
            end
        end
    end
    for s in SR
        for h = 0:H
            for t = 1:T
                if newHeightsHeur[[s,h,t]] == 1
                    setvalue(heights[s,h,t], 1);
                else
                    setvalue(heights[s,h,t], 0);
                end
            end
        end
    end
    for s in SO
        for h = 0:H
            if finalHeightsHeur[[s,h]] == 1
                setvalue(finalh[s,h], 1);
            else
                setvalue(finalh[s,h], 0);
            end
        end
    end
    for m = 1:n
        for t = 1:T
            if WHeur[[m,t]] == 1
                setvalue(w[m,t], 1);
            else
                setvalue(w[m,t], 0);
            end
        end
    end

    @objective(IPModel, Min, sum( sum(costMove[[s,r]]*x[s,r,t] for s in SP for r in posteriorStacks[s]) + sum(costPreMove[[s,r]]*d[s,r,t] for s in Stot for r in SP) for t = 1:T) + costToGo * sum( alpha[h] * (sum(heights[s,h,T] for s in SR) + sum(finalh[s,h] for s in SO))  for h = 2:H));

    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################


    for t = 1:T
        ## Uniqueness of pre-move
        @constraint(IPModel, sum(d[s,r,t] for s in Stot for r in SP) == 1);
        ## Uniqueness of move
        @constraint(IPModel, sum(x[s,r,t] for s in SP for r in posteriorStacks[s]) == 1);
    end

    for s in setdiff(Stot,posCraneInitial)
        for r in SP
            ## Starting point condition
            @constraint(IPModel, d[s,r,1] == 0);
        end
    end

    for s in SP
        for t = 1:T
            ## Conservation of flow (1)
            @constraint(IPModel, sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in Stot) == 0);
        end
    end

    for s in Stot
        for t = 2:T
            ## Conservation of flow (2)
            @constraint(IPModel, sum(d[s,r,t] for r in SP) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
        end
    end

    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################

    ## Number of stacks
    @constraint(IPModel, sum(x[s,r,t] for s in SI for r in posteriorStacks[s] for t = 1:T) == N - n);

    for m = 1:n
        ## Number of deliveries
        @constraint(IPModel, sum(w[m,t] for t = 1:T) == 1);
    end

    for s in SR
        for t = 1:T
            ## Relation between variables x and w
            @constraint(IPModel, x[s,IOPoint[s],t] - sum(w[m,t] for m in contInStackToRetrieve[s]) == 0);
        end
    end

    for m = 1:n
        if heightsInitial[realStack[indicesToRetrieve[m][1]][1],realStack[indicesToRetrieve[m][1]][2]] > indicesToRetrieve[m][2]
            ## Relation between variables w and θ
            @constraint(IPModel, w[m,1] <= 0);
        end
        for t = 2:T
            ## Relation between variables w and θ
            @constraint(IPModel, w[m,t] - heights[indicesToRetrieve[m][1],indicesToRetrieve[m][2],t-1] <= 0);
        end
        for t = 1:T
            ## Relation between variables w and θ (2)
            @constraint(IPModel, sum(heights[indicesToRetrieve[m][1],h,t] for h = 0:indicesToRetrieve[m][2]-1) - sum(w[m,u] for u = 1:t) <= 0);
        end
    end

    #####################################################################
    ########################### Heights Update ##########################
    #####################################################################

    for s in SR
        ## Initial update
        @constraint(IPModel, sum(h * heights[s,h,1] for h = 1:H) - sum(x[r,s,1] for r in anteriorStacks[s]) + sum(x[s,r,1] for r in posteriorStacks[s]) == heightsInitial[realStack[s][1],realStack[s][2]]);
        for t = 2:T
            ## Further updates
            @constraint(IPModel, sum(h * heights[s,h,t] for h = 1:H) - sum(h * heights[s,h,t-1] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s]) + sum(x[s,r,t] for r in posteriorStacks[s]) == 0);
        end
        for t = 1:T
            ## Uniqueness of height
            @constraint(IPModel, sum(heights[s,h,t] for h = 0:H) == 1);
        end
    end

    for s in SO
        ## Final heights
        @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightsInitial[realStack[s][1],realStack[s][2]]);
        ## Uniqueness of final height
        @constraint(IPModel, sum(finalh[s,h] for h = 0:H) == 1);
    end

    #####################################################################
    ################################ Cuts ###############################
    #####################################################################

    @constraint(IPModel, sum(d[posCraneInitial,r,1] for r in SP) == 1);

    for s in SI
        if length(anteriorStacks[s]) == 0
            for r in SP
                for t = 2:T
                    @constraint(IPModel, d[s,r,t] == 0);
                end
            end
        end
    end

    for s in SR
        for m in contInStackToRetrieve[s]
            for t = 1:T
                @constraint(IPModel, sum(x[r,s,t] for r in anteriorStacks[s]) - sum(w[m,u] for u = 1:t) <= 0);
            end
        end
    end

    println(MathProgBase.numvar(IPModel) , " variables");
    println(MathProgBase.numconstr(IPModel), " constraints");
    println("Solving ....");

    tic();
    status = solve(IPModel);
    println("Solved !");
    timeToSolve = toc();

    X = getvalue(x);
    D = getvalue(d);
    newHeights = getvalue(heights);
    finalHeights = getvalue(finalh);
    W = getvalue(w);

    return (X,D,newHeights,finalHeights,W,timeToSolve);
end
