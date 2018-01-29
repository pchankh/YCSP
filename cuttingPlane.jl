function LP(H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime)

    IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));

    #####################################################################
    ############################# Variables #############################
    #####################################################################


    @variable(IPModel, 0 <= w[1:T,1:T] <= 1);
    @variable(IPModel, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(IPModel, 0 <= dInit[SX] <= 1);
    @variable(IPModel, 0 <= d[SY,SX,2:T] <= 1);
    @variable(IPModel, 0 <= finalh[SB,0:H] <= 1);

    #####################################################################
    ############################# Objective #############################
    #####################################################################

    @objective(IPModel, Min, sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum( alpha[h] * sum(finalh[s,h] for s in SB)  for h = 2:H));


    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################

    for m = 1:T
        ## Uniqueness of move
        @constraint(IPModel, sum(w[m,t] for t = 1:T) == 1);
    end

    for t = 1:T
        ## Uniqueness of move (2)
        @constraint(IPModel, sum(w[m,t] for m = 1:T) == 1);
    end

    for m = 1:T
        if previousContToMove[m] != 0
            for t = 1:T
                # Relation between variables x and w for retrievals
                @constraint(IPModel, sum(w[previousContToMove[m],u] for u = 1:t-1) - w[m,t] >= 0);
            end
        end
    end

    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################

    @constraintref relXandWConst[1:T,1:T];

    for m = 1:n
        for t = 1:T
            # Relation between variables x and w for retrievals
            relXandWConst[m,t] = @constraint(IPModel, sum(x[stackOf[m],r,t] for r in loadOf[m]) - w[m,t] >= 0);
        end
    end

    for m = n+1:N
        for t = 1:T
            # Relation between variables x and w for stackings
            relXandWConst[m,t] = @constraint(IPModel, sum(x[s,r,t] for s in unloadFrom[m] for r in posteriorStacks[s])  - w[m,t] >= 0);
        end
    end

    for m = N+1:T
        for t = 1:T
            # Relation between variables x and w for relocations
            relXandWConst[m,t] = @constraint(IPModel, sum(x[stackOf[m],r,t] for r in intersect(posteriorStacks[stackOf[m]],SB)) - w[m,t] >= 0);
        end
    end

    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################

    ## Uniqueness of move without a container
    @constraint(IPModel, sum(dInit[r] for r in SX) == 1);
    for t = 2:T
        @constraint(IPModel, sum(d[s,r,t] for s in SY for r in SX) == 1);
    end

    for t = 1:T
        ## Uniqueness of move with a container
        @constraint(IPModel, sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    end

    for s in SX
        ## Conservation of flow (1)
        @constraint(IPModel, sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
        for t = 2:T
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
    ########################### Final Heights ###########################
    #####################################################################

    for s in SR
        ## Final heights
        @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightOf[contMinHeightStack[s]] - 1);
    end

    for s in SO
        ## Final heights
        @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightsInitial[realStack[s][1],realStack[s][2]]);
    end

    for s in SB
        ## Uniqueness of final height
        @constraint(IPModel, sum(finalh[s,h] for h = 0:H) == 1);
    end

    #####################################################################
    ################################ Cuts ###############################
    #####################################################################

    # for s in SR
    #     for t = 1:T
    #         @constraint(IPModel, sum(x[r,s,t] for r in anteriorStacks[s]) - sum(w[contMinHeightStack[s],u] for u = 1:t) <= 0);
    #     end
    # end
    #
    # for s in SL
    #     if length(anteriorStacks[s]) == 0
    #         for r in SX
    #             for t = 2:T
    #                 @constraint(IPModel, d[s,r,t] == 0);
    #             end
    #         end
    #     end
    # end

    # println(MathProgBase.numvar(IPModel) , " variables");
    # println(MathProgBase.numconstr(IPModel), " constraints");
    # println("Solving ....");

    # tic();
    status = solve(IPModel);
    # println("Solved !");
    # timeToSolve = toc();

    PLP = getdual(relXandWConst);

    WLP = getvalue(w);

    ZLP = 0;
    for m = 1:T
        for t = 1:T
            ZLP = ZLP + WLP[m,t]*PLP[m,t];
        end
    end

    objLP = getobjectivevalue(IPModel);

    return (ZLP,PLP,objLP);
end



function subProblem(Wgiven, H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime)

    subLP = Model(solver = GurobiSolver(OutputFlag = 0,MIPGap = gapOfMIP, TimeLimit = limitOfTime));

    #####################################################################
    ############################# Variables #############################
    #####################################################################

    @variable(subLP, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(subLP, 0 <= dInit[SX] <= 1);
    @variable(subLP, 0 <= d[SY,SX,2:T] <= 1);
    @variable(subLP, 0 <= finalh[SB,0:H] <= 1);

    #####################################################################
    ############################# Objective #############################
    #####################################################################

    @objective(subLP, Min, sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum( alpha[h] * sum(finalh[s,h] for s in SB)  for h = 2:H));

    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################

    @constraintref relXandWConst[1:T,1:T];

    for m = 1:n
        for t = 1:T
            # Relation between variables x and w for retrievals
            relXandWConst[m,t] = @constraint(subLP, sum(x[stackOf[m],r,t] for r in loadOf[m]) >= Wgiven[m,t]);
        end
    end

    for m = n+1:N
        for t = 1:T
            # Relation between variables x and w for stackings
            relXandWConst[m,t] = @constraint(subLP, sum(x[s,r,t] for s in unloadFrom[m] for r in posteriorStacks[s]) >= Wgiven[m,t]);
        end
    end

    for m = N+1:T
        for t = 1:T
            # Relation between variables x and w for relocations
            relXandWConst[m,t] = @constraint(subLP, sum(x[stackOf[m],r,t] for r in intersect(posteriorStacks[stackOf[m]],SB)) >= Wgiven[m,t]);
        end
    end

    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################

    ## Uniqueness of move without a container
    @constraint(subLP, sum(dInit[r] for r in SX) == 1);
    for t = 2:T
        @constraint(subLP, sum(d[s,r,t] for s in SY for r in SX) == 1);
    end

    for t = 1:T
        ## Uniqueness of move with a container
        @constraint(subLP, sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    end

    for s in SX
        ## Conservation of flow (1)
        @constraint(subLP, sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
        for t = 2:T
            @constraint(subLP, sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
        end
    end

    for s in SY
        for t = 2:T
            ## Conservation of flow (2)
            @constraint(subLP, sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
        end
    end

    #####################################################################
    ########################### Final Heights ###########################
    #####################################################################

    for s in SB
        ## Final heights
        if s in SR
            @constraint(subLP, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightOf[contMinHeightStack[s]] - 1);
        else
            @constraint(subLP, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightsInitial[realStack[s][1],realStack[s][2]]);
        end
        ## Uniqueness of final height
        @constraint(subLP, sum(finalh[s,h] for h = 0:H) == 1);
    end

    #####################################################################
    ################################ Cuts ###############################
    #####################################################################


    # tic();
    status = solve(subLP);
    # timeToSolveLP = toc();

    P = getdual(relXandWConst);

    Z = 0;
    for m = 1:T
        for t = 1:T
            Z = Z + Wgiven[m,t]*P[m,t];
        end
    end

    Obj = getobjectivevalue(subLP);

    return (Z,P,Obj);
end

function MIP_2(H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime)

    (ZLP,PLP,objLP) = LP(H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime)

    CutModel = Model(solver = GurobiSolver(PreCrush=1, OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));

    #####################################################################
    ############################# Variables #############################
    #####################################################################

    @variable(CutModel, w[1:T,1:T], Bin);

    @variable(CutModel, z >= 0);

    #####################################################################
    ############################# Objective #############################
    #####################################################################

    @objective(CutModel, Min, z);


    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################

    for m = 1:T
        ## Uniqueness of move
        @constraint(CutModel, sum(w[m,t] for t = 1:T) == 1);
    end

    for t = 1:T
        ## Uniqueness of move (2)
        @constraint(CutModel, sum(w[m,t] for m = 1:T) == 1);
    end

    for m = 1:T
        if previousContToMove[m] != 0
            for t = 1:T
                # Relation between variables x and w for retrievals
                @constraint(CutModel, sum(w[previousContToMove[m],u] for u = 1:t-1) - w[m,t] >= 0);
            end
        end
    end

    function mycutgenerator(cb)
        w_val = getvalue(w);
        z_val = getvalue(z);

        # println("z_val = ",z_val);

        # Allow for some impreciseness in the solution
        TOL = 1e-3

        (Z,P,Obj) = subProblem(w_val, H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime);

        if z_val < Obj - TOL
            @lazyconstraint(cb, z - sum(P[m,t]*w[m,t] for m = 1:T for t = 1:T) >= Obj - Z);
        end

    end  # End of callback function

    # Tell JuMP/Gurobi to use our callback function
    addlazycallback(CutModel, mycutgenerator, fractional=true);

    # tic();
    status = solve(CutModel);
    # println("Solved !");
    # timeToSolve = toc();

    W = getvalue(w);

    return (W);
end

function MIP_3(H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime)

    saveP = Dict{Int64,Array{Float64}}();
    saveD = Dict{Int64,Float64}();
    nCutsToBeginWith = min(10,factorial(n));
    for i = 1:nCutsToBeginWith
        Wgiven = zeros(Int64,T,T);
        t = 0;
        remainToStack = N - n;
        orderDrawn = randperm(n);
        for m in orderDrawn
            if sum(Wgiven[m,:]) == 0
                blockingCont = [m];
                while previousContToMove[blockingCont[length(blockingCont)]] != 0
                    append!(blockingCont,previousContToMove[blockingCont[length(blockingCont)]]);
                end
                for c = length(blockingCont):-1:1
                    t = t + 1;
                    Wgiven[blockingCont[c],t] = 1;
                end
                if remainToStack > 0
                    t = t + 1;
                    Wgiven[N - remainToStack + 1,t] = 1;
                    remainToStack = remainToStack - 1;
                end
            end
        end
        while remainToStack > 0
            t = t + 1;
            Wgiven[N - remainToStack + 1,t] = 1;
            remainToStack = remainToStack - 1;
        end
        (Z,P,Obj) = subProblem(Wgiven, H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime);
        saveP[i] = P;
        saveD[i] = Obj - Z;
    end
    nCut = 0;
    z_val = 0;
    TOL = 1e-3;
    Z = 1e4;
    cutMax = 30;
    bestObj = Z;
    Obj =
    bestW = zeros(Int64,(T,T));
    while (nCut < cutMax && z_val < Obj - TOL)
        CutModel = Model(solver = GurobiSolver(OutputFlag = 0,MIPGap = gapOfMIP, TimeLimit = limitOfTime));
        @variable(CutModel, w[1:T,1:T], Bin);
        @variable(CutModel, z >= 0);
        @objective(CutModel, Min, z);
        for m = 1:T
            @constraint(CutModel, sum(w[m,t] for t = 1:T) == 1);
        end
        for t = 1:T
            @constraint(CutModel, sum(w[m,t] for m = 1:T) == 1);
        end
        for m = 1:T
            if previousContToMove[m] != 0
                for t = 1:T
                    @constraint(CutModel, sum(w[previousContToMove[m],u] for u = 1:t-1) - w[m,t] >= 0);
                end
            end
        end
        for i = 1:nCut
            @constraint(CutModel, z - sum(saveP[i][m,t]*w[m,t] for m = 1:T for t = 1:T) >= saveD[i]);
        end
        solve(CutModel);
        w_val = getvalue(w);
        z_val = getvalue(z);
        (Z,P,Obj) = subProblem(w_val, H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime);
        println(z_val, " < ", Z);
        if z_val < Obj - TOL
            nCut = nCut + 1;
            saveP[nCut] = P;
            saveD[nCut] = Obj - Z;
        end
        if bestObj > Obj
            bestObj = Obj;
            bestW = w_val;
            println(round((bestObj - z_val)/bestObj*100,2),"%");
        end
    end

    return (W);
end
