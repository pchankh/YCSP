

function MIP(H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime)

    IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));
    # x = @variable(IPModel, [s in SX,r in posteriorStacks[s],1:T], Bin);
    # dInit = @variable(IPModel, [SX], Bin);
    # d = @variable(IPModel, [SY,SX,2:T], Bin);
    # heights = @variable(IPModel, [SR,0:H,1:T], Bin);
    # finalh = @variable(IPModel, [SO,0:H], Bin);
    w = @variable(IPModel, [1:N,1:T], Bin);

    @variable(IPModel, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(IPModel, 0 <= dInit[SX] <= 1);
    @variable(IPModel, 0 <= d[SY,SX,2:T] <= 1);
    @variable(IPModel, 0 <= heights[SR,0:H,1:T] <= 1);
    @variable(IPModel, 0 <= finalh[SO,0:H] <= 1);
    # @variable(IPModel, 0 <= w[1:N,1:T] <= 1);

    t = 0;
    remainToStack = N - n;
    didRetrieve = true;
    for s in SR
        if didRetrieve && remainToStack > 0
            t = t + 1;
            println(t)
            setvalue(w[N - remainToStack + 1,t],1);
            remainToStack = remainToStack - 1;
            didRetrieve = false;
        end
        minHeight = heightsInitial[realStack[s][1],realStack[s][2]];
        for h = heightsInitial[realStack[s][1],realStack[s][2]]:-1:1
            for m = 1:n
                if stackOf[m] == s && heightOf[m] == h
                    t = t + (minHeight - h) + 1;
                    println(t)
                    #setvalue(w[m,t],1);
                    minHeight = heightOf[m];
                    didRetrieve = true;
                end
            end
        end
    end
    while remainToStack > 0
        t = t + 1;
        setvalue(w[N - remainToStack + 1,t],1);
        remainToStack = remainToStack - 1;
    end

    @objective(IPModel, Min, sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum( alpha[h] * (sum(heights[s,h,T] for s in SR) + sum(finalh[s,h] for s in SO))  for h = 2:H));

    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################

    ## Uniqueness of pre-move
    @constraint(IPModel, sum(dInit[r] for r in SX) == 1);
    for t = 2:T
        @constraint(IPModel, sum(d[s,r,t] for s in SY for r in SX) == 1);
    end

    for t = 1:T
        ## Uniqueness of move
        @constraint(IPModel, sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    end

    for s in SX
        @constraint(IPModel, sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
        for t = 2:T
            ## Conservation of flow (1)
            @constraint(IPModel, sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
        end
    end

    for s in SY
        for t = 2:T
            ## Conservation of flow (2)
            @constraint(IPModel, sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
        end
    end

    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################

    for m = 1:N
        ## Number of productive moves
        @constraint(IPModel, sum(w[m,t] for t = 1:T) == 1);
    end

    for t = 1:T
        ## Number of productive moves
        @constraint(IPModel, sum(w[m,t] for m = 1:N) <= 1);
        # @constraint(IPModel, sum(w[m,t] for m = 1:N) + sum(x[s,r,t] for s in SR for r in intersect(posteriorStacks[s],SB)) == 1);
    end

    # for t = 1:T
    #     @constraint(IPModel, sum(w[m,t] for m = n+1:N if toBeUnloaded[m-n]=="up") - sum(x[s,r,t] for s=intersect(R*S+1:R*S+S,SU) for r in posteriorStacks[s]) == 0);
    #     @constraint(IPModel, sum(w[m,t] for m = n+1:N if toBeUnloaded[m-n]=="down") - sum(x[s,r,t] for s=intersect((R+1)*S+1:R*S+2*S,SU) for r in posteriorStacks[s]) == 0);
    # end

    for m = n+1:N
        for t = 1:T
            # Relation between variables x and w for stackings
            @constraint(IPModel, w[m,t] - sum(x[s,r,t] for s in unloadFrom[m] for r in posteriorStacks[s]) <= 0);
        end
    end

÷

    for m = n+2:N
        for t = 2:T
            # Relation between variables x and w for retrievals
            @constraint(IPModel, w[m,t] - sum(w[m-1,u] for u = 1:t-1) <= 0);
        end
    end

    for m = 1:n
        for t = 1:T
            # Relation between variables x and w for retrievals
            @constraint(IPModel, w[m,t] - sum(x[stackOf[m],r,t] for r in loadOf[m]) <= 0);
        end
    end

    # for s in SR
    #     for r in intersect(posteriorStacks[s],IOPoints[s])
    #         for t = 1:T
    #             @constraint(IPModel, x[s,r,t] - sum(w[m,t] for m in contStackIOPoint[[s,r]]) == 0);
    #         end
    #     end
    # end

    for m = 1:n
        if heightsInitial[realStack[stackOf[m]][1],realStack[stackOf[m]][2]] > heightOf[m]
            ## Relation between variables w and θ
            @constraint(IPModel, w[m,1] == 0);
        end
        for t = 2:T
            ## Relation between variables w and θ
            @constraint(IPModel, w[m,t] - heights[stackOf[m],heightOf[m],t-1] <= 0);
        end
        for t = 1:T
            ## Relation between variables w and θ (2)
            @constraint(IPModel, sum(heights[stackOf[m],h,t] for h = 0:heightOf[m]-1) - sum(w[m,u] for u = 1:t) <= 0);
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

    for s in union(SL)
        if length(anteriorStacks[s]) == 0
            for r in SP
                for t = 2:T
                    @constraint(IPModel, d[s,r,t] == 0);
                end
            end
        end
    end

    for s in SR
        for m = 1:n
            if stackOf[m] == s
                for t = 1:T
                    @constraint(IPModel, sum(x[r,s,t] for r in anteriorStacks[s]) - sum(w[m,u] for u = 1:t) <= 0);
                end
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
    DInit = getvalue(dInit)
    D = getvalue(d);
    newHeights = getvalue(heights);
    finalHeights = getvalue(finalh);
    W = getvalue(w);

    return (X,DInit,D,newHeights,finalHeights,W,timeToSolve);
end
