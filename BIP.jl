################################################################################
########################## BIP (Binary Integer Program) ########################
################################################################################
## This function solves the binary integer program and return all the solution
## variables and the objective value
## For testing
function BIP(limitOfTime,T,moveFrom,SX,posteriorStacks,SY,SB,Z,costMove,costPreMove,posCraneInitial,costToGo,alpha,N,realToClusterOrder,changeOfOrder,typeOfTruck,previousContToMove,posteriorContStacks,SR,contMinHeightStack,anteriorStacks,artificialHeights)

    IPModel = Model(solver = GurobiSolver(TimeLimit = limitOfTime));
    #####################################################################
    ############################# Variables #############################
    #####################################################################
    @variable(IPModel, w[m = 1:T, s in moveFrom[m], t = 1:T], Bin);
    @variable(IPModel, x[s in SX,r in posteriorStacks[s],1:T], Bin);
    @variable(IPModel, dInit[SX], Bin);
    @variable(IPModel, d[SY,SX,2:T], Bin);
    @variable(IPModel, finalh[SB,0:Z], Bin);
    # @variable(IPModel, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    # @variable(IPModel, 0 <= dInit[SX] <= 1);
    # @variable(IPModel, 0 <= d[SY,SX,2:T] <= 1);
    # @variable(IPModel, 0 <= finalh[SB,0:Z] <= 1);
    #####################################################################
    ############################# Objective #############################
    #####################################################################
    @objective(IPModel, Min, sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum( alpha[h] * sum(finalh[s,h] for s in SB)  for h = 2:Z));
    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################
    @constraint(IPModel, uniquenessMove[m = 1:T], sum(w[m,s,t] for t = 1:T for s in moveFrom[m]) == 1);
    @constraint(IPModel, uniquenessMove_2[t = 1:T], sum(w[m,s,t] for m = 1:T for s in moveFrom[m]) == 1);
    @constraint(IPModel, orderConstraint_1[o = 1:N-1,t = changeOfOrder[typeOfTruck[realToClusterOrder[o]]]+1:T], sum(w[realToClusterOrder[p],s,u] for p = o+1:N for s in moveFrom[realToClusterOrder[p]] for u = 1:t) - (t - min(changeOfOrder[typeOfTruck[realToClusterOrder[o]]],T)) * sum(w[realToClusterOrder[o],s,u] for s in moveFrom[realToClusterOrder[o]] for u = 1:t) <= min(changeOfOrder[typeOfTruck[realToClusterOrder[o]]],T));
    @constraint(IPModel, orderConstraint_2[o = 2:N,t = 1:T-changeOfOrder[typeOfTruck[realToClusterOrder[o]]]], sum(w[realToClusterOrder[p],s,u] for p = 1:o-1 for s in moveFrom[realToClusterOrder[p]] for u = t:T) - (T - t + 1 - min(changeOfOrder[typeOfTruck[realToClusterOrder[o]]],T)) * sum(w[realToClusterOrder[o],s,u] for s in moveFrom[realToClusterOrder[o]] for u = t:T) <= min(changeOfOrder[typeOfTruck[realToClusterOrder[o]]],T));
    for m = 1:T
        if (m <= n || m >= N+1) && previousContToMove[m] != 0
            for t = 2:T
                # Relation between variables x and w for retrievals
                @constraint(IPModel, sum(w[previousContToMove[m],s,t-1] for s in moveFrom[previousContToMove[m]]) - sum(w[m,s,t] for s in moveFrom[m]) == 0);
            end
        end
    end
    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################
    @constraint(IPModel, constW_1[m = 1:T,s in moveFrom[m],t = 1:T], sum(x[s,r,t] for r in posteriorContStacks[[m,s]]) - w[m,s,t] >= 0);
    @constraint(IPModel, constW_2[s in SR,t = 1:T], sum(w[contMinHeightStack[s],s,u] for u = 1:t-1) - sum(x[r,s,t] for r in anteriorStacks[s]) >= 0);
    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################
    ## Uniqueness of move without a container
    @constraint(IPModel, constDInit, sum(dInit[r] for r in SX) == 1);
    @constraint(IPModel, constDUnique[t = 2:T], sum(d[s,r,t] for s in SY for r in SX) == 1);
    ## Uniqueness of move with a container
    @constraint(IPModel, constXUnique[t = 1:T],sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    ## Conservation of flow (1)
    @constraint(IPModel, conservFlowInit_1[s in SX], sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
    @constraint(IPModel, conservFlow_1[s in SX,t = 2:T], sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
    ## Conservation of flow (2)
    @constraint(IPModel, conservFlow_2[s in SY,t = 2:T], sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
    #####################################################################
    ########################### Final Heights ###########################
    #####################################################################
    ## Final heights
    @constraint(IPModel, finalHeightConst[s in SB], sum(h * finalh[s,h] for h = 1:Z) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == artificialHeights[s]);
    ## Uniqueness of final height
    @constraint(IPModel, finalHeightUniq[s in SB], sum(finalh[s,h] for h = 0:Z) == 1);
    status = solve(IPModel);
    moveWithCont = getvalue(x);
    moveInit = getvalue(dInit);
    moveWithoutCont = getvalue(d);
    finalHeights = getvalue(finalh);
    orderContStack = getvalue(w);
    obj = getobjectivevalue(IPModel);
    return (moveWithCont,moveInit,moveWithoutCont,finalHeights,orderContStack,obj);
end
