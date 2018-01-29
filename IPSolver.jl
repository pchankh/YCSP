################################################################################
################################### IPSolver ###################################
################################################################################
## This function solves the binary integer program

function IPSolver(limitOfTime,NBar,L,SL,SE,SB,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,heightsTilde,restricted,blockingCont,N,minStack,delta,SR,cstVerticalCost,vZ,NR,NS,NU,E,requestsID,heightsBlock,positionCont,blockID,realStack,nameIOPoint,printSolutions,outputFolder,period)

    IPModel = Model(solver = GurobiSolver(TimeLimit = limitOfTime));
    @variable(IPModel, w[n = 1:NBar, s in L[n],1:NBar], Bin);
    @variable(IPModel, dInit[SL], Bin);
    @variable(IPModel, dEmpty[SE,SL,2:NBar], Bin);
    @variable(IPModel, dLoaded[SL,SE,1:NBar], Bin);
    @variable(IPModel, finalZ[r in intersect(SE,SB),z = heightsTilde[r]:Z], Bin);
    @objective(IPModel, Min, sum(costEmptyDrive[[posCraneInitial,s]]*dInit[s] for s in SL) + sum(costEmptyDrive[[r,s]]*dEmpty[r,s,k] for r in SE for s in SL for k = 2:NBar) + sum(costLoadedDrive[[s,r]]*dLoaded[s,r,k] for s in SL for r in SE for k = 1:NBar) + sum(beta[z]*finalZ[r,z] for r in intersect(SE,SB) for z = heightsTilde[r]+1:Z));
    @constraint(IPModel, uniquenessMove[n = 1:NBar], sum(w[n,s,k] for s in L[n] for k = 1:NBar) == 1);
    @constraint(IPModel, uniquenessMove_2[k = 1:NBar], sum(w[n,s,k] for n = 1:NBar for s in L[n]) == 1);
    for n in keys(blockingCont)
        @constraint(IPModel, sum(w[n,s,1] for s in L[n]) == 0);
        for k = 2:NBar
            if restricted && blockingCont[n] in NU
                @constraint(IPModel, sum(w[n,s,k] for s in L[n]) - sum(w[blockingCont[n],s,k-1] for s in L[blockingCont[n]]) == 0);
            else
                @constraint(IPModel, sum(w[n,s,k] for s in L[n]) - sum(w[blockingCont[n],s,l] for s in L[blockingCont[n]] for l=1:k-1) <= 0);
            end
        end
    end
    @constraint(IPModel, orderConstraint_1[n = 1:N,k = (n+delta[n,2]+1):NBar], sum(w[m,s,l] for m = 1:N for s in L[m] for l = 1:k-1) + (k-(n+delta[n,2]))*sum(w[n,s,k] for s in L[n]) <= k-1);
    @constraint(IPModel, orderConstraint_2[n = 1:N,k = 1:(NBar-(N-n+delta[n,1])-1)], sum(w[m,s,l] for m = 1:N for s in L[m] for l = k+1:NBar) + (NBar - k - (N - n + delta[n,1]))*sum(w[n,s,k] for s in L[n]) <= NBar - k);
    @constraint(IPModel, constDInit, sum(dInit[s] for s in SL) == 1);
    @constraint(IPModel, constEmptyUnique[k = 2:NBar], sum(dEmpty[r,s,k] for r in SE for s in SL) == 1);
    @constraint(IPModel, constLoadedUnique[k = 1:NBar], sum(dLoaded[s,r,k] for s in SL for r in SE) == 1);
    @constraint(IPModel, conservFlowInit_1[s in SL], sum(dLoaded[s,r,1] for r in SE) - dInit[s] == 0);
    @constraint(IPModel, conservFlow_1[s in SL,k = 2:NBar], sum(dLoaded[s,r,k] for r in SE) - sum(dEmpty[r,s,k] for r in SE) == 0);
    @constraint(IPModel, conservFlow_2[r in SE,k = 2:NBar], sum(dEmpty[r,s,k] for s in SL) - sum(dLoaded[s,r,k-1] for s in SL) == 0);
    @constraint(IPModel, finalHeightUniq[r in intersect(SE,SB)], sum(finalZ[r,z] for z = heightsTilde[r]:Z) == 1);
    @constraint(IPModel, finalHeightConst[r in intersect(SE,SB)], sum(z * finalZ[r,z] for z = heightsTilde[r]:Z) - sum(dLoaded[s,r,k] for s in SL for k = 1:NBar) == heightsTilde[r]);
    @constraint(IPModel, constW_1[n = 1:NBar,s in L[n],k = 1:NBar], sum(dLoaded[s,r,k] for r in E[n]) - w[n,s,k] >= 0);
    @constraint(IPModel, constW_2[r in intersect(SR,SE),k = 1:NBar], sum(dLoaded[s,r,k] for s in SL) - sum(w[minStack[r],r,l] for l = 1:k-1) <= 0);
    TT = STDOUT; # save original STDOUT stream
    redirect_stdout();
    status = solve(IPModel);
    redirect_stdout(TT);
    WSol = getvalue(w);
    DSolInit = getvalue(dInit);
    DSolEmpty = getvalue(dEmpty);
    DSolLoaded = getvalue(dLoaded);
    ZSolFinal = getvalue(finalZ);
    horizontalCost = sum(costEmptyDrive[[posCraneInitial,s]]*DSolInit[s] for s in SL) + sum(costEmptyDrive[[r,s]]*DSolEmpty[r,s,k] for r in SE for s in SL for k = 2:NBar) + sum(costLoadedDrive[[s,r]]*DSolLoaded[s,r,k] for s in SL for r in SE for k = 1:NBar);
    verticalCost = cstVerticalCost;
    if length(intersect(SE,SB))>0
        verticalCost = verticalCost - 1/vZ*sum(z*(z+1)*ZSolFinal[r,z] for r in intersect(SE,SB) for z = heightsTilde[r]:Z);
    end
    if printSolutions
        outputFile = joinpath(outputFolder,string(period,"_Period"));
        f = open(outputFile,"a");
        posCrane = posCraneInitial;
        for k = 1:NBar
            for n = 1:NBar
                for s in L[n]
                    if round(WSol[n,s,k]) == 1
                        write(f,string("Stage  ", k,"\n"));
                        if n in NR
                            for r in E[n]
                                if round(DSolLoaded[s,r,k]) == 1
                                    write(f,string("Empty Drive from ", realStack[posCrane], " to ", realStack[s],"\n"));
                                    write(f,string("Retrieve container ", requestsID[n], "(", n, ") from ", realStack[s], " to ", nameIOPoint[r],"\n"));
                                    posCrane = r;
                                end
                            end
                        elseif n in NS
                            for r in E[n]
                                if round(DSolLoaded[s,r,k]) == 1
                                    write(f,string("Empty Drive from ", realStack[posCrane], " to ", realStack[s],"\n"));
                                    write(f,string("Store container ", requestsID[n], "(", n, ") from ", nameIOPoint[s], " in ", realStack[r],"\n"));
                                    posCrane = r;
                                end
                            end
                        elseif n in NU
                            for r in E[n]
                                if round(DSolLoaded[s,r,k]) == 1
                                    write(f,string("Empty Drive from ", realStack[posCrane], " to ", realStack[s],"\n"));
                                    write(f,string("Relocate container ", requestsID[n], " from ", realStack[s], " in ", realStack[r],"\n"));
                                    posCrane = r;
                                end
                            end
                        end
                    end
                end
            end
            write(f,string("-----------------------------\n"));
        end
        write(f,string("Initial State of Block\n"));
        str = "";
        if IOPointsPosition == "Asian-left" || IOPointsPosition == "two-sided"
            str = string(str,"\tLeft");
        end
        for s = 1:X
            str = string(str,"\tStack_",s);
        end
        if IOPointsPosition == "Asian-right" || IOPointsPosition == "two-sided"
            str = string(str,"\tRight");
        end
        write(f,string(str,"\n"));
        if IOPointsPosition == "Euro"
            write(f,string("Landside\n"));
        end
        sta = 0;
        for r = 1:Y
            str = string("Row_",r);
            if IOPointsPosition == "Asian-left" || IOPointsPosition == "two-sided"
                str = string(str,"\t");
            end
            for s = 1:X
                sta += 1;
                str = string(str,"\t",Int64(round(heightsBlock[sta])));
            end
            if IOPointsPosition == "Asian-right" || IOPointsPosition == "two-sided"
                str = string(str,"\t");
            end
            write(f,string(str,"\n"));
        end
        if IOPointsPosition == "Euro"
            write(f,string("Seaside\n"));
        end
        close(f);
    end
    for k = 1:NBar
        for n = 1:NBar
            for s in L[n]
                if round(WSol[n,s,k]) == 1
                    if n in NR
                        for r in E[n]
                            if round(DSolLoaded[s,r,k]) == 1
                                delete!(blockID,[s,heightsBlock[s]]);
                                heightsBlock[s] = heightsBlock[s] - 1;
                                delete!(positionCont,requestsID[n]);
                                posCraneInitial = r;
                            end
                        end
                    elseif n in NS
                        for r in E[n]
                            if round(DSolLoaded[s,r,k]) == 1
                                heightsBlock[r] = heightsBlock[r] + 1;
                                positionCont[requestsID[n]] = [r,heightsBlock[r]];
                                blockID[[r,heightsBlock[r]]] = requestsID[n];
                                posCraneInitial = r;
                            end
                        end
                    elseif n in NU
                        for r in E[n]
                            if round(DSolLoaded[s,r,k]) == 1
                                delete!(blockID,[s,heightsBlock[s]]);
                                heightsBlock[s] = heightsBlock[s] - 1;
                                heightsBlock[r] = heightsBlock[r] + 1;
                                positionCont[requestsID[n]] = [r,heightsBlock[r]];
                                blockID[[r,heightsBlock[r]]] = requestsID[n];
                                posCraneInitial = r;
                            end
                        end
                    end
                end
            end
        end
    end
    return (horizontalCost,verticalCost,heightsBlock,positionCont,blockID,posCraneInitial);
end
