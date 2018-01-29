function costFunction(W,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights)
    LP = Model(solver = GurobiSolver(OutputFlag = 0));
    #####################################################################
    ############################# Variables #############################
    #####################################################################
    @variable(LP, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(LP, 0 <= dInit[SX] <= 1);
    @variable(LP, 0 <= d[SY,SX,2:T] <= 1);
    @variable(LP, 0 <= finalh[SB,0:Z] <= 1);
    #####################################################################
    ############################# Objective #############################
    #####################################################################
    @objective(LP, Min,sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum(alpha[h]*finalh[s,h] for s in SB for h = 2:Z));
    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################
    @constraint(LP, constW_1[m = 1:T,s in moveFrom[m],t = 1:T], sum(x[s,r,t] for r in posteriorContStacks[[m,s]]) >= W[[m,s,t]]);
    @constraint(LP, constW_2[s in SR,t = 1:T], sum(x[r,s,t] for r in anteriorStacks[s]) <= sum(W[[contMinHeightStack[s],s,u]] for u = 1:t-1));
    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################
    ## Uniqueness of move without a container
    @constraint(LP, constDInit, sum(dInit[r] for r in SX) == 1);
    @constraint(LP, constDUnique[t = 2:T], sum(d[s,r,t] for s in SY for r in SX) == 1);
    ## Uniqueness of move with a container
    @constraint(LP, constXUnique[t = 1:T],sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    ## Conservation of flow (1)
    @constraint(LP, conservFlowInit_1[s in SX], sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
    @constraint(LP, conservFlow_1[s in SX,t = 2:T], sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
    ## Conservation of flow (2)
    @constraint(LP, conservFlow_2[s in SY,t = 2:T], sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
    #####################################################################
    ########################### Final Heights ###########################
    #####################################################################
    ## Final heights
    @constraint(LP, finalHeightConst[s in SB], sum(h * finalh[s,h] for h = 1:Z) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == artificialHeights[s]);
    ## Uniqueness of final height
    @constraint(LP, finalHeightUniq[s in SB], sum(finalh[s,h] for h = 0:Z) == 1);
    TT = STDOUT; # save original STDOUT stream
    redirect_stdout();
    status = solve(LP);
    redirect_stdout(TT);
    Obj = getobjectivevalue(LP);
    moveWithCont = getvalue(x);
    moveInit = getvalue(dInit);
    moveWithoutCont = getvalue(d);
    finalHeights = getvalue(finalh);
    return (Obj,moveWithCont,moveInit,moveWithoutCont,finalHeights);
end

################################################################################
################################### LowerBound #################################
################################################################################

function LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights)
    LP = Model(solver = GurobiSolver(OutputFlag = 0));
    #####################################################################
    ############################# Variables #############################
    #####################################################################
    @variable(LP, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(LP, 0 <= dInit[SX] <= 1);
    @variable(LP, 0 <= d[SY,SX,2:T] <= 1);
    @variable(LP, 0 <= finalh[SB,0:Z] <= 1);
    @variable(LP, 0 <= wStack[m = 1:T, s in moveFrom[m]] <= 1);
    #####################################################################
    ############################# Objective #############################
    #####################################################################
    @objective(LP, Min,sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum(alpha[h]*finalh[s,h] for s in SB for h = 2:Z));
    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################
    @constraint(LP, constW_1[m = 1:T,s in moveFrom[m]], sum(x[s,r,orderCont[m]] for r in posteriorContStacks[[m,s]]) >= wStack[m,s]);
    @constraint(LP, constW_2[s in SR,t = 1:orderCont[contMinHeightStack[s]]], sum(x[r,s,t] for r in anteriorStacks[s]) == 0);
    @constraint(LP, constW_3[m = 1:T], sum(wStack[m,s] for s in moveFrom[m]) == 1);
    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################
    ## Uniqueness of move without a container
    @constraint(LP, constDInit, sum(dInit[r] for r in SX) == 1);
    @constraint(LP, constDUnique[t = 2:T], sum(d[s,r,t] for s in SY for r in SX) == 1);
    ## Uniqueness of move with a container
    @constraint(LP, constXUnique[t = 1:T],sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    ## Conservation of flow (1)
    @constraint(LP, conservFlowInit_1[s in SX], sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
    @constraint(LP, conservFlow_1[s in SX,t = 2:T], sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
    ## Conservation of flow (2)
    @constraint(LP, conservFlow_2[s in SY,t = 2:T], sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
    #####################################################################
    ########################### Final Heights ###########################
    #####################################################################
    ## Final heights
    @constraint(LP, finalHeightConst[s in SB], sum(h * finalh[s,h] for h = 1:Z) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == artificialHeights[s]);
    ## Uniqueness of final height
    @constraint(LP, finalHeightUniq[s in SB], sum(finalh[s,h] for h = 0:Z) == 1);
    TT = STDOUT; # save original STDOUT stream
    redirect_stdout();
    status = solve(LP);
    redirect_stdout(TT);
    LBObj = getobjectivevalue(LP);
    Wsta = getvalue(wStack);
    moveWithCont = getvalue(x);
    moveInit = getvalue(dInit);
    moveWithoutCont = getvalue(d);
    finalHeights = getvalue(finalh);
    nonIntegralSolution = Dict{Int64,Array{Int64}}();
    integralSolution = Dict{Int64,Int64}();
    capacityStack_1 = Dict{Int64,Float64}();
    for m = 1:T
        for s in moveFrom[m]
            if round(Wsta[m,s],6) < 1 - 0.000001 &&  round(Wsta[m,s],6) > 0.000001
                for r in posteriorStacks[s]
                    if round(moveWithCont[s,r,orderCont[m]],6) > 0.000001
                        if m in keys(nonIntegralSolution)
                            append!(nonIntegralSolution[m],[r]);
                        else
                            nonIntegralSolution[m] = [r];
                        end
                        if r in keys(capacityStack_1)
                            capacityStack_1[r] = capacityStack_1[r] + round(moveWithCont[s,r,orderCont[m]],6);
                        else
                            capacityStack_1[r] = round(moveWithCont[s,r,orderCont[m]],6);
                        end
                    end
                end
            else
                for r in posteriorStacks[s]
                    if round(moveWithCont[s,r,orderCont[m]],6) > 1 - 0.000001
                        integralSolution[m] = r;
                    end
                end
            end
        end
    end
    capacityStack = Dict{Int64,Int64}();
    for r in keys(capacityStack_1)
        capacityStack[r] = ceil(capacityStack_1[r]);
    end
    return (LBObj,nonIntegralSolution,integralSolution,capacityStack,moveWithCont,moveInit,moveWithoutCont,finalHeights);
end

function computeAssignmentCost(invOrder,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove)
    previousSta = posCraneInitial;
    cost = 0;
    for o = 1:T
        Cont = invOrder[o];
        if Cont <= N
            cost = cost + costPreMove[[previousSta,intersect(anteriorStacks[StackCont[Cont]],moveFrom[Cont])[1]]] + costMove[[intersect(anteriorStacks[StackCont[Cont]],moveFrom[Cont])[1],StackCont[Cont]]];
            previousSta = StackCont[Cont];
        end
    end
    return cost;
end

function firstAssignment(T,integralSolution,nonIntegralSolution,capacityStack)
    StackCont = zeros(Int64,T);
    for m in keys(integralSolution)
        StackCont[m] = integralSolution[m];
    end
    SlackRemaining = Dict{Int64,Int64}();
    capacityStackLoc = Dict{Int64,Int64}();
    nContperStack = Dict{Int64,Int64}();
    for m in keys(nonIntegralSolution)
        SlackRemaining[m] = 0;
        for r in nonIntegralSolution[m]
            if r in keys(nContperStack)
                nContperStack[r] = nContperStack[r] + 1;
            else
                nContperStack[r] = 1;
                capacityStackLoc[r] = capacityStack[r];
            end
            SlackRemaining[m] = SlackRemaining[m] + capacityStackLoc[r];
        end
    end
    while length(collect(keys(SlackRemaining))) > 0
        minSlacks = Array{Int64}(0);
        for m in keys(SlackRemaining)
            if SlackRemaining[m] == minimum(collect(values(SlackRemaining)))
                append!(minSlacks,[m]);
            end
        end
        selectedCont = 0;
        selectedStack = 0;
        selectedRatio = 0;
        for m in minSlacks
            for r in nonIntegralSolution[m]
                if capacityStackLoc[r]/nContperStack[r] > selectedRatio
                    selectedRatio = capacityStackLoc[r]/nContperStack[r];
                    selectedStack = r;
                    selectedCont = m;
                end
            end
        end
        capacityStackLoc[selectedStack] = capacityStackLoc[selectedStack] - 1;
        for r in nonIntegralSolution[selectedCont]
            nContperStack[r] = nContperStack[r] - 1;
        end
        StackCont[selectedCont] = selectedStack;
        delete!(SlackRemaining, selectedCont);
    end
    return StackCont;
end

function addDrop(m,StackCont,r)
    StackContLoc = zeros(Int64,length(StackCont));
    for i = 1:length(StackCont)
        StackContLoc[i] = StackCont[i];
    end
    StackContLoc[m] = r;
    return StackContLoc;
end

function pairwiseExchange(m,StackCont,l)
    StackContLoc = zeros(Int64,length(StackCont));
    for i = 1:length(StackCont)
        StackContLoc[i] = StackCont[i];
    end
    StackContLoc[m], StackContLoc[l] = StackContLoc[l], StackContLoc[m];
    return StackContLoc;
end

function LocalImprovement(orderCont,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove,capacityStack,nonIntegralSolution)
    invOrder = zeros(Int64,T);
    for m = 1:T
        invOrder[orderCont[m]] = m;
    end
    currentValue = computeAssignmentCost(invOrder,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove);
    localOptimum = false;
    while !localOptimum
        localOptimum = true;
        bestNewStackCont = zeros(Int64,T);
        bestNewValue = currentValue;
        unusedCapacity = Dict{Int64,Int64}();
        for r in keys(capacityStack)
            unusedCapacity[r] = capacityStack[r];
        end
        for m in keys(nonIntegralSolution)
            unusedCapacity[StackCont[m]] = unusedCapacity[StackCont[m]] - 1;
        end
        for m in keys(nonIntegralSolution)
            for r in nonIntegralSolution[m]
                if r != StackCont[m] && unusedCapacity[r] > 0
                    StackContLoc = addDrop(m,StackCont,r);
                    locValue = computeAssignmentCost(invOrder,StackContLoc,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove);
                    if locValue < bestNewValue
                        bestNewValue = locValue;
                        bestNewStackCont = StackContLoc;
                        localOptimum = false;
                    end
                end
            end
        end
        for m in keys(nonIntegralSolution)
            for l in keys(nonIntegralSolution)
                if m != l && StackCont[m] in nonIntegralSolution[l] && StackCont[l] in nonIntegralSolution[m]
                    StackContLoc = pairwiseExchange(m,StackCont,l);
                    locValue = computeAssignmentCost(invOrder,StackContLoc,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove);
                    if locValue < bestNewValue
                        bestNewValue = locValue;
                        bestNewStackCont = StackContLoc;
                        localOptimum = false;
                    end
                end
            end
        end
        if !localOptimum
            currentValue = bestNewValue;
            for m = 1:N
                StackCont[m] = bestNewStackCont[m];
            end
        end
    end
    return StackCont;
end

function generalizedAssignmentModel(orderCont,T,integralSolution,nonIntegralSolution,capacityStack,posCraneInitial,N,costPreMove,anteriorStacks,moveFrom,costMove)
    StackCont = firstAssignment(T,integralSolution,nonIntegralSolution,capacityStack);
    StackCont = LocalImprovement(orderCont,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove,capacityStack,nonIntegralSolution);
    return StackCont;
end

function UpperBound(T,N,orderCont,nonIntegralSolution,integralSolution,capacityStack,posCraneInitial,costPreMove,anteriorStacks,moveFrom,costMove,SX,posteriorStacks,SY,Z,costToGo,alpha,SR,contMinHeightStack,SB,artificialHeights,moveWithCont,moveInit,moveWithoutCont,finalHeights,LBObj)
    if length(keys(integralSolution)) == T
        StackCont = zeros(Int64,T);
        for m in keys(integralSolution)
            StackCont[m] = integralSolution[m];
        end
        UBObj = LBObj;
    else
        StackCont = generalizedAssignmentModel(orderCont,T,integralSolution,nonIntegralSolution,capacityStack,posCraneInitial,N,costPreMove,anteriorStacks,moveFrom,costMove);
        W = Dict{Array{Int64},Int64}();
        for m = 1:T
            for s in moveFrom[m]
                for t = 1:T
                    if t == orderCont[m] && s == intersect(anteriorStacks[StackCont[m]],moveFrom[m])[1]
                        W[[m,s,t]] = 1;
                    else
                        W[[m,s,t]] = 0;
                    end
                end
            end
        end
        (UBObj,moveWithCont,moveInit,moveWithoutCont,finalHeights) = costFunction(W,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
    end
    return (UBObj,StackCont,moveWithCont,moveInit,moveWithoutCont,finalHeights);
end

function fullOrderContainers(T, N, permutationProductive, blockingCont)
    orderCont = zeros(Int64,T);
    t = 0;
    sta = 0;
    for m = 1:N
        if permutationProductive[m] in keys(blockingCont)
            for i = length(blockingCont[permutationProductive[m]]):-1:1
                t = t + 1;
                orderCont[blockingCont[permutationProductive[m]][i]] = t;
            end
        else
            t = t + 1;
            orderCont[permutationProductive[m]] = t;
        end
    end
    return orderCont;
end

function feasibleSwap(k,l,currentPermuOrder,changeOrderReal,typeOfTruck,clusterToRealOrder)
    isSwapPossible = (abs(currentPermuOrder[k] - l) <= changeOrderReal[currentPermuOrder[k]] && abs(currentPermuOrder[l] - k) <= changeOrderReal[currentPermuOrder[l]]);
    if isSwapPossible
        isSwapPossible = !(typeOfTruck[clusterToRealOrder[currentPermuOrder[k]]] == "internal" && typeOfTruck[clusterToRealOrder[currentPermuOrder[l]]] == "internal");
    end
    return isSwapPossible;
end
