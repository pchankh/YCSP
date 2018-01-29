function subProblem(WtoBeCut, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime)

    LP = Model(solver = GurobiSolver(OutputFlag = printSolver, TimeLimit = limitOfTime));

    #####################################################################
    ############################# Variables #############################
    #####################################################################

    @variable(LP, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(LP, 0 <= dInit[SX] <= 1);
    @variable(LP, 0 <= d[SY,SX,2:T] <= 1);
    @variable(LP, 0 <= finalh[SB,0:H] <= 1);

    #####################################################################
    ############################# Objective #############################
    #####################################################################

    @objective(LP, Min, sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum( alpha[h] * sum(finalh[s,h] for s in SB)  for h = 2:H));

    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################

    @constraint(LP, constW_1[m = 1:T,s in moveFrom[m],t = 1:T], sum(x[s,r,t] for r in posteriorContStacks[[m,s]]) >= WtoBeCut[m,s,t]);

    @constraint(LP, constW_2[s in SR,t = 1:T-1], sum(x[r,s,t] for r in anteriorStacks[s]) <= sum(WtoBeCut[contMinHeightStack[s],s,u] for u = 1:t-1));

    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################

    ## Uniqueness of move without a container
    @constraint(LP, sum(dInit[r] for r in SX) == 1);
    for t = 2:T
        @constraint(LP, sum(d[s,r,t] for s in SY for r in SX) == 1);
    end

    for t = 1:T
        ## Uniqueness of move with a container
        @constraint(LP, sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    end

    for s in SX
        ## Conservation of flow (1)
        @constraint(LP, sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
        for t = 2:T
            @constraint(LP, sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
        end
    end

    for s in SY
        for t = 2:T
            ## Conservation of flow (2)
            @constraint(LP, sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
        end
    end

    #####################################################################
    ########################### Final Heights ###########################
    #####################################################################

    for s in SB
        ## Final heights
        @constraint(LP, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == artificialHeights[s]);
    end

    for s in SB
        ## Uniqueness of final height
        @constraint(LP, sum(finalh[s,h] for h = 0:H) == 1);
    end

    println("Finding new cut ....");
    tic();
    status = solve(LP);
    # println("Solved !");
    timeToSolve = toc();

    Obj = getobjectivevalue(LP);

    return (Z,P,Obj)

end
