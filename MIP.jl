

function MIP(H, N, n, artificialHeights, moveFrom, toBeUnloaded, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, deadlines, printSolver, gapOfMIP, limitOfTime)

    IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));

    #####################################################################
    ############################# Variables #############################
    #####################################################################

    @variable(IPModel, w[m = 1:T, s in moveFrom[m], t = 1:T], Bin);

    @variable(IPModel, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(IPModel, 0 <= dInit[SX] <= 1);
    @variable(IPModel, 0 <= d[SY,SX,2:T] <= 1);
    @variable(IPModel, 0 <= finalh[SB,0:H] <= 1);
    @variable(IPModel, servingTimes[1:T,1:T] >= 0);


    #####################################################################
    ############################# Incumbent #############################
    #####################################################################

    # Wgiven = zeros(Int64,T,T);
    # t = 0;
    # remainToStack = N - n;
    # for m = 1:n
    #     if sum(Wgiven[m,:]) == 0
    #         blockingCont = [m];
    #         while previousContToMove[blockingCont[length(blockingCont)]] != 0
    #             append!(blockingCont,previousContToMove[blockingCont[length(blockingCont)]]);
    #         end
    #         for c = length(blockingCont):-1:1
    #             t = t + 1;
    #             Wgiven[blockingCont[c],t] = 1;
    #         end
    #         if remainToStack > 0
    #             t = t + 1;
    #             Wgiven[N - remainToStack + 1,t] = 1;
    #             remainToStack = remainToStack - 1;
    #         end
    #     end
    # end
    # while remainToStack > 0
    #     t = t + 1;
    #     Wgiven[N - remainToStack + 1,t] = 1;
    #     remainToStack = remainToStack - 1;
    # end
    # for t = 1:T
    #     for m = 1:n
    #         setvalue(w[m,moveFrom[m][1],t],Wgiven[m,t]);
    #     end
    #     for m = n+1:N
    #         setvalue(w[m,moveFrom[m][1],t],Wgiven[m,t]);
    #     end
    #     for m = N+1:T
    #         setvalue(w[m,moveFrom[m][1],t],Wgiven[m,t]);
    #     end
    # end

    M = T*(maximum(values(costMove)) + maximum(values(costPreMove)));

    #####################################################################
    ############################# Objective #############################
    #####################################################################

    @objective(IPModel, Min, sum(servingTimes[m,T] for m = 1:T) + costToGo * sum(alpha[h] * sum(finalh[s,h] for s in SB)  for h = 2:H));

    for m = 1:T
        if m <= N
            @constraint(IPModel, sum(servingTimes[m,t] for t = 1:T) <= deadlines[m]);
        end
        for t = 1:T
            @constraint(IPModel, servingTimes[m,t] - M*sum(w[m,s,t] for s in moveFrom[m]) <= 0);
        end
    end
    for t = 1:T
        @constraint(IPModel, sum(costMove[[s,r]]*x[s,r,u] for s in SX for r in posteriorStacks[s] for u = 1:t) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,u] for s in SY for r in SX for u = 2:t) -  sum(servingTimes[m,t] for m = 1:T) <= 0);
    end

    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################

    for m = 1:T
        ## Uniqueness of move
        @constraint(IPModel, sum(w[m,s,t] for t = 1:T for s in moveFrom[m]) == 1);
    end

    for t = 1:T
        ## Uniqueness of move (2)
        @constraint(IPModel, sum(w[m,s,t] for m = 1:T for s in moveFrom[m]) == 1);
    end

    for m = 1:T
        if previousContToMove[m] != 0
            for t = 1:T
                # Relation between variables x and w for retrievals
                @constraint(IPModel, sum(w[previousContToMove[m],s,u] for u = 1:t-1 for s in moveFrom[previousContToMove[m]]) - sum(w[m,s,t] for s in moveFrom[m]) >= 0);
            end
        end
    end

    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################

    for m = 1:T
        for s in moveFrom[m]
            for t =1:T
                # Relation between variables x and w for retrievals
                @constraint(IPModel, sum(x[s,r,t] for r in posteriorContStacks[[m,s]]) - w[m,s,t] >= 0);
            end
        end
    end

    for s in SR
        for t = 1:T-1
            @constraint(IPModel, sum(w[contMinHeightStack[s],s,u] for u = 1:t-1) - sum(x[r,s,t] for r in anteriorStacks[s]) >= 0);
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

    for s in SB
        ## Final heights
        @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == artificialHeights[s]);
    end

    # for s in SR
    #     ## Final heights
    #     @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightOf[contMinHeightStack[s]] - 1);
    # end
    #
    # for s in SO
    #     ## Final heights
    #     @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightsInitial[realStack[s][1],realStack[s][2]]);
    # end

    for s in SB
        ## Uniqueness of final height
        @constraint(IPModel, sum(finalh[s,h] for h = 0:H) == 1);
    end

    # println(MathProgBase.numvar(IPModel) , " variables");
    # println(MathProgBase.numconstr(IPModel), " constraints");
    # println("Solving ....");

    tic();
    status = solve(IPModel);
    # println("Solved !");
    timeToSolve = toc();

    X = getvalue(x);
    DInit = getvalue(dInit)
    D = getvalue(d);
    finalHeights = getvalue(finalh);
    W = getvalue(w);
    ServTimes = getvalue(servingTimes);

    obj = getobjectivevalue(IPModel);

    return (X,DInit,D,finalHeights,ServTimes,W,obj,timeToSolve);
end
