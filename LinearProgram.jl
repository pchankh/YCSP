function LinearProgram(Wgiven, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, limitOfTime)

    LP = Model(solver = GurobiSolver(OutputFlag = 0, TimeLimit = limitOfTime));

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

    @constraint(LP, constW_1[m = 1:T,s in moveFrom[m],t = 1:T], sum(x[s,r,t] for r in posteriorContStacks[[m,s]]) >= Wgiven[[m,s,t]]);

    @constraint(LP, constW_2[s in SR,t = 1:T], sum(x[r,s,t] for r in anteriorStacks[s]) <= sum(Wgiven[[contMinHeightStack[s],s,u]] for u = 1:t-1));

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

    status = solve(LP);
    # println("Solved !");

    # P1 = getdual(constW_1);
    # P2 = getdual(constW_2);
    # P = Dict{Array{Int64},Float64}();
    # for m = 1:T
    #     for s in moveFrom[m]
    #         for t = 1:T
    #             if round(P1[m,s,t],5) != 0
    #                 P[[m,s,t]] = round(P1[m,s,t],5);
    #             end
    #         end
    #     end
    # end
    # for s in SR
    #     for t = 1:T
    #         if round(P2[s,t],5) != 0
    #             for u = 1:t-1
    #                 if [contMinHeightStack[s],s,u] in keys(P)
    #                     P[[contMinHeightStack[s],s,u]] = round(P[[contMinHeightStack[s],s,u]] + P2[s,t],5);
    #                 else
    #                     P[[contMinHeightStack[s],s,u]] = round(P2[s,t],5);
    #                 end
    #             end
    #         end
    #     end
    # end

    Obj = getobjectivevalue(LP);

    # return (Obj,P);
    return Obj;
end

function initialSolution(T, moveFrom, n, N, previousContToMove, SI, H, posteriorStacks, heightsInitial)
    Wgiven = Dict{Array{Int64},Int64}();
    orderCont = zeros(Int64,T);
    orderTime = zeros(Int64,T);
    for m = 1:T
        for s in moveFrom[m]
            for t = 1:T
                Wgiven[[m,s,t]] = 0;
            end
        end
    end
    alreadyRetrieved = zeros(n);
    potentialStack = Dict{Int64,Int64}();
    for s in SI
        potentialStack[s] = H*length(posteriorStacks[s]);
        for r in posteriorStacks[s]
            potentialStack[s] = potentialStack[s] - heightsInitial[r];
        end
    end
    t = 0;
    remainToStack = N - n;
    for m = 1:n
        if alreadyRetrieved[m] == 0
            blockingCont = [m];
            while previousContToMove[blockingCont[length(blockingCont)]] != 0
                append!(blockingCont,previousContToMove[blockingCont[length(blockingCont)]]);
            end
            for c = length(blockingCont):-1:1
                t = t + 1;
                Wgiven[[blockingCont[c],moveFrom[m][1],t]] = 1;
                orderCont[blockingCont[c]] = t;
                orderTime[t] = blockingCont[c];
                if blockingCont[c] <= n
                    alreadyRetrieved[blockingCont[c]] = 1;
                end
            end
            # if remainToStack > 0
            #     t = t + 1;
            #     sta = rand(1:length(moveFrom[N - remainToStack + 1]));
            #     while potentialStack[moveFrom[N - remainToStack + 1][sta]] == 0
            #         sta = rand(1:length(moveFrom[N - remainToStack + 1]));
            #     end
            #     potentialStack[sta] = potentialStack[moveFrom[N - remainToStack + 1][sta]] - 1;
            #     Wgiven[[N - remainToStack + 1,moveFrom[N - remainToStack + 1][sta],t]] = 1;
            #     orderCont[N - remainToStack + 1] = t;
            #     orderTime[t] = N - remainToStack + 1;
            #     remainToStack = remainToStack - 1;
            # end
        end
    end
    while remainToStack > 0
        t = t + 1;
        sta = rand(1:length(moveFrom[N - remainToStack + 1]));
        while potentialStack[moveFrom[N - remainToStack + 1][sta]] == 0
            sta = rand(1:length(moveFrom[N - remainToStack + 1]));
        end
        potentialStack[sta] = potentialStack[moveFrom[N - remainToStack + 1][sta]] - 1;
        Wgiven[[N - remainToStack + 1,moveFrom[N - remainToStack + 1][sta],t]] = 1;
        orderCont[N - remainToStack + 1] = t;
        orderTime[t] = N - remainToStack + 1;
        remainToStack = remainToStack - 1;
    end
    positiveW = Dict{Int64,Array{Int64}}();
    i = 0;
    for m = 1:T
        for s in moveFrom[m]
            if Wgiven[[m,s,orderCont[m]]] == 1
                i = i + 1;
                positiveW[i] = [m,s,orderCont[m]];
            end
        end
    end
    return (Wgiven,orderCont,orderTime,positiveW);
end

function isExchange(indRem,indAdd,orderCont,orderTime,previousContToMove)
    isExchangable = true;
    tA = indRem[3];
    tR = orderCont[indAdd[1]];
    for t = min(tA,tR):T
        if isExchangable && t > tR && previousContToMove[indRem[1]] == orderTime[t]
            isExchangable = false;
        end
        if isExchangable && t > tA && previousContToMove[indAdd[1]] == orderTime[t]
            isExchangable = false;
        end
    end
    # if indRem[3] != indAdd[3]
    #     isExchangable = false;
    # elseif indRem[1] != indAdd[1]
    #     if indRem[3] < orderCont[indAdd[1]]
    #         t = indRem[3];
    #         while isExchangable && t <= orderCont[indAdd[1]]
    #             isExchangable = ((previousContToMove[orderTime[t]] != indRem[1]) && (previousContToMove[indAdd[1]] != orderTime[t]));
    #             t = t + 1;
    #         end
    #     else
    #         t = orderCont[indAdd[1]];
    #         while isExchangable && t <= indRem[3]
    #             isExchangable = ((previousContToMove[indRem[1]] != orderTime[t]) && (previousContToMove[orderTime[t]] != indAdd[1]));
    #             t = t + 1;
    #         end
    #     end
    # end
    return isExchangable
end

function exchangeIndices(Wgiven,T,indRem,indAdd,moveFrom,orderCont,orderTime,positiveW)
    Wnew = Dict{Array{Int64},Int64}();
    orderContnew = zeros(Int64,T);
    orderTimenew = zeros(Int64,T);
    positiveWnew = Dict{Int64,Array{Int64}}();
    i = 0;
    for m = 1:T
        orderContnew[m] = orderCont[m];
        orderTimenew[m] = orderTime[m];
        for s in moveFrom[m]
            for t = 1:T
                Wnew[[m,s,t]] = Wgiven[[m,s,t]];
                if Wgiven[[m,s,t]] == 1 && [m,s,t] != indRem
                    i = i + 1;
                    positiveWnew[i] = [m,s,t];
                end
            end
        end
    end
    Wnew[indRem] = 0;
    Wnew[indAdd] = 1;
    if indRem[1] != indAdd[1]
        for s in moveFrom[indAdd[1]]
            Wnew[[indAdd[1],s,orderCont[indAdd[1]]]] = 0;
        end
        Wnew[[indRem[1],indRem[2],orderCont[indAdd[1]]]] = 1;
        orderContnew[indRem[1]] = orderCont[indAdd[1]];
        orderContnew[indAdd[1]] = indAdd[3];
        orderTimenew[indRem[3]] = indAdd[1];
        orderTimenew[orderCont[indAdd[1]]] = indRem[1];
    end
    i = i + 1;
    positiveWnew[i] = indAdd;
    return (Wnew,orderContnew,orderTimenew,positiveWnew);
end
