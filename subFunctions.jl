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
    return Obj;
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
    WLB = Dict{Array{Int64},Float64}();
    
    nonIntegralSolution = Dict{Int64,Array{Int64}}();
    for m = 1:T
        for s in moveFrom[m]
            WLB[[m,s]] = round(Wsta[m,s],8);
            if WLB[[m,s]] < 1 - 0.0001 &&  WLB[[m,s]] > 0.0001
                for r in posteriorStacks[s]
                    if
                if m in keys(nonIntegralSolution)
                    append!(nonIntegralSolution[m],[s]);
                else
                    nonIntegralSolution[m] = [s];
                end
            end
        end
    end
    return (LBObj,WLB,nonIntegralSolution);
end

function UpperBound(LBObj,WLB,nonIntegralSolution,orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights)
    if length(nonIntegralSolution) == 0
        StackCont = zeros(Int64,T);
        for m = 1:T
            for s in moveFrom[m]
                if WLB[[m,s]] >= 1 - 0.0001
                    StackCont[m] = s;
                end
            end
        end
        UBObj = LBObj;
    else
        nPerIOPoints = Dict{Int64,Float64}();
        for m in keys(nonIntegralSolution)
            for s in nonIntegralSolution[m]
                if s in keys(nPerIOPoints)
                    nPerIOPoints[s] = nPerIOPoints[s] + WLB[[m,s]];
                else
                    nPerIOPoints[s] = WLB[[m,s]];
                end
            end
        end
        for s in keys(nPerIOPoints)
            nPerIOPoints[s] = ceil(round(nPerIOPoints[s],7));
        end
        StackCont = zeros(Int64,T);
        for m = 1:T
            if m in keys(nonIntegralSolution)
                for s in nonIntegralSolution[m]
                    if nPerIOPoints[s] != 0 && StackCont[m] == 0
                        StackCont[m] = s;
                        nPerIOPoints[s] = nPerIOPoints[s] - 1;
                    end
                end
            else
                for s in moveFrom[m]
                    if WLB[[m,s]] >= 1 - 0.0001
                        StackCont[m] = s;
                    end
                end
            end
        end
        W = Dict{Array{Int64},Int64}();
        for m = 1:T
            for s in moveFrom[m]
                for t = 1:T
                    if t == orderCont[m] && s == StackCont[m]
                        W[[m,s,t]] = 1;
                    else
                        W[[m,s,t]] = 0;
                    end
                end
            end
        end
        UBObj = costFunction(W,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
    end
    return (UBObj,StackCont);
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
