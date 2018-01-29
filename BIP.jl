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
    @objective(IPModel, Min, sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum(sum((alpha[h] - alpha[max(artificialHeights[s],1)])*finalh[s,h] for s in SB)  for h = 2:Z));
    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################
    @constraint(IPModel, uniquenessMove[m = 1:T], sum(w[m,s,t] for t = 1:T for s in moveFrom[m]) == 1);
    @constraint(IPModel, uniquenessMove_2[t = 1:T], sum(w[m,s,t] for m = 1:T for s in moveFrom[m]) == 1);
    @constraint(IPModel, orderConstraint_1[o = 1:N,t = 1:T], sum(w[m,s,u] for m = 1:N for s in moveFrom[m] for u = 1:t-1) + (t - (o+changeOfOrder[typeOfTruck[o]])) * sum(w[realToClusterOrder[o],s,t] for s in moveFrom[realToClusterOrder[o]]) <= t-1);
    @constraint(IPModel, orderConstraint_2[o = 1:N,t = 1:T], sum(w[m,s,u] for m = 1:N for s in moveFrom[m] for u = t+1:T) + (T - t - (N-o+changeOfOrder[typeOfTruck[o]])) * sum(w[realToClusterOrder[o],s,t] for s in moveFrom[realToClusterOrder[o]]) <= T - t);
    # @constraint(IPModel, orderConstraint_1[o = 1:N-1,t = changeOfOrder[typeOfTruck[o]]+1:T], sum(w[realToClusterOrder[p],s,u] for p = o+1:N for s in moveFrom[realToClusterOrder[p]] for u = 1:t) - (t - min(changeOfOrder[typeOfTruck[o]],T)) * sum(w[realToClusterOrder[o],s,u] for s in moveFrom[realToClusterOrder[o]] for u = 1:t) <= min(changeOfOrder[typeOfTruck[o]],T));
    # @constraint(IPModel, orderConstraint_2[o = 2:N,t = 1:T-changeOfOrder[typeOfTruck[o]]], sum(w[realToClusterOrder[p],s,u] for p = 1:o-1 for s in moveFrom[realToClusterOrder[p]] for u = t:T) - (T - t + 1 - min(changeOfOrder[typeOfTruck[o]],T)) * sum(w[realToClusterOrder[o],s,u] for s in moveFrom[realToClusterOrder[o]] for u = t:T) <= min(changeOfOrder[typeOfTruck[o]],T));
    for m = 1:T
        if (m <= n || m >= N+1) && previousContToMove[m] != 0 && previousContToMove[m] >= N+1
            for t = 2:T
                # Relation between variables x and w for retrievals
                @constraint(IPModel, sum(w[previousContToMove[m],s,t-1] for s in moveFrom[previousContToMove[m]]) - sum(w[m,s,t] for s in moveFrom[m]) == 0);
            end
        end
    end
    for m = 1:T
        if (m <= n || m >= N+1) && previousContToMove[m] != 0 && previousContToMove[m] <= n
            for t = 2:T
                # Relation between variables x and w for retrievals
                @constraint(IPModel, sum(w[previousContToMove[m],s,u] for s in moveFrom[previousContToMove[m]] for u = 1:t-1) - sum(w[m,s,t] for s in moveFrom[m]) >= 0);
            end
        end
    end
    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################
    @constraint(IPModel, constW_1[m = 1:T,s in moveFrom[m],t = 1:T], sum(x[s,r,t] for r in posteriorContStacks[[m,s]]) - w[m,s,t] >= 0);
    @constraint(IPModel, constW_2[s in SR,t = 1:T], sum(w[contMinHeightStack[s],s,u] for u = 1:t) - sum(x[r,s,t] for r in anteriorStacks[s]) >= 0);
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


function printResultTest(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights)
    f = open("SolutionTest.txt","w")
    write(f,string("-------------------------------\n"));
    write(f,string("BLOCK STATE AT THE BEGINNING\n"));
    str = "";
    if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:X
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tRight");
    end
    write(f,string(str,"\n"));
    if IOPointsPosition == "up-and-down"
        write(f,string("Seaside\n"));
    end
    sta = 0;
    for r = 1:Y
        str = string("Row_",r);
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        for s = 1:X
            sta += 1;
            str = string(str,"\t",Int64(round(heightsInitial[realStack[sta][1],realStack[sta][2]])));
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        write(f,string(str,"\n"));
    end
    if IOPointsPosition == "up-and-down"
        write(f,string("Landside\n"));
    end
    lastCranePosition = Array{Any}(2,2);
    lastCranePosition[1,1] = "Row";
    lastCranePosition[2,1] = "Stack";
    lastCranePosition[1,2] = "TBD";
    lastCranePosition[2,2] = "TBD";
    for t=1:T
        write(f,string("-----------------------------\n"));
        write(f,string("Time ", t,"\n"));
        if t == 1
            for r in SX
                if moveInit[r] > 0.9
                    if r == posCraneInitial
                        if r in SB
                            write(f,string("Crane stays at ", realStack[r],"\n"));
                        else
                            write(f,string("Crane stays at ", nameIOPoint[r],"\n"));
                        end
                    else
                        if posCraneInitial in SB && r in SB
                            write(f,string("Crane moves from ", realStack[posCraneInitial], " to ", realStack[r],"\n"));
                        elseif posCraneInitial in SB && !(r in SB)
                            write(f,string("Crane moves from ", realStack[posCraneInitial], " to ", nameIOPoint[r],"\n"));
                        elseif !(posCraneInitial in SB) && r in SB
                            write(f,string("Crane moves from ", nameIOPoint[posCraneInitial], " to ", realStack[r],"\n"));
                        else
                            write(f,string("Crane moves from ", nameIOPoint[posCraneInitial], " to ", nameIOPoint[r],"\n"));
                        end
                    end
                end
            end
        else
            for s in SY
                for r in SX
                    if moveWithoutCont[s,r,t] > 0.9
                        if r == s
                            if s in SB
                                write(f,string("Crane stays at ", realStack[s],"\n"));
                            else
                                write(f,string("Crane stays at ", nameIOPoint[s],"\n"));
                            end
                        else
                            if s in SB && r in SB
                                write(f,string("Crane moves from ", realStack[s], " to ", realStack[r],"\n"));
                            elseif s in SB && !(r in SB)
                                write(f,string("Crane moves from ", realStack[s], " to ", nameIOPoint[r],"\n"));
                            elseif !(s in SB) && r in SB
                                write(f,string("Crane moves from ", nameIOPoint[s], " to ", realStack[r],"\n"));
                            else
                                write(f,string("Crane moves from ", nameIOPoint[s], " to ", nameIOPoint[r],"\n"));
                            end
                        end
                    end
                end
            end
        end
        for s in SX
            for r in posteriorStacks[s]
                if moveWithCont[s,r,t] > 0.9
                    if s in SB && r in SB
                        write(f,string("Relocate container from ", realStack[s], " to ", realStack[r],"\n"));
                    elseif s in SU && r in SB
                        for m = n+1:N
                            if s in unloadFrom[m] && orderContStack[m,s,t] > 0.9
                                write(f,string("Stack container ", clusterToRealOrder[m]," from ",nameIOPoint[s]," on ", realStack[r],"\n"));
                                if t == T
                                    lastCranePosition[:,2] = realStack[r];
                                end
                            end
                        end
                    elseif s in SB && r in SL
                        for m = 1:n
                            if stackOf[m] == s && orderContStack[m,s,t] > 0.9
                                write(f,string("Retrieve container ", clusterToRealOrder[m], " from ", realStack[s], " to ", nameIOPoint[r],"\n"));
                                if t == T
                                    lastCranePosition[:,2] = split(nameIOPoint[r]);
                                    if all(isnumber, lastCranePosition[1,2])
                                        lastCranePosition[1,2] = parse(Int64,lastCranePosition[1,2]);
                                    end
                                    if all(isnumber, lastCranePosition[2,2])
                                        lastCranePosition[2,2] = parse(Int64,lastCranePosition[2,2]);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    write(f,string("-----------------------------\n"));
    write(f,string("BLOCK STATE AT THE END\n"));
    str = "";
    if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:X
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tRight");
    end
    write(f,string(str,"\n"));
    if IOPointsPosition == "up-and-down"
        write(f,string("Seaside\n"));
    end
    lastHeight = zeros(Int64,Y,X);
    for s in SB
        h = 0;
        for z = 1:Z
            h = h + z * finalHeights[s,z];
        end
        lastHeight[realStack[s][1],realStack[s][2]] = Int64(round(h));
    end
    for r = 1:Y
        str = string("Row_",r);
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        for s = 1:X
            str = string(str,"\t",lastHeight[r,s]);
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        write(f,string(str,"\n"));
    end
    if IOPointsPosition == "up-and-down"
        write(f,string("Landside\n"));
    end
    close(f);
    writecsv("newBlockTest.csv",lastHeight);
    writecsv("newCranePositionTest.csv",lastCranePosition);
end
