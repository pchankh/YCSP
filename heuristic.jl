

function lowerBound(order,NBar,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack)
    reverseOrder = Array{Int64}(NBar);
    for k = 1:NBar
        reverseOrder[order[k]] = k;
    end
    validCycles = Dict{Int64,Array{Int64}}();
    for r in union(intersect(SE,SB),SR)
        validCycles[r] = [];
        for k = 1:NBar
            if r in E[order[k]]
                push!(validCycles[r],k);
            end
        end
    end
    LBModel = Model(solver = GurobiSolver(OutputFlag=0));
    @variable(LBModel, 0 <= w[k = 1:NBar,s in L[order[k]]] <= 1);
    @variable(LBModel, 0 <= dInit[L[order[1]]] <= 1);
    @variable(LBModel, 0 <= dEmpty[k = 2:NBar,r in E[order[k-1]],s in L[order[k]]] <= 1);
    @variable(LBModel, 0 <= dLoaded[k = 1:NBar,s in L[order[k]],r in E[order[k]]] <= 1);
    @variable(LBModel, 0 <= finalZ[r in intersect(SE,SB),z = heightsTilde[r]:Z] <= 1);
    @objective(LBModel, Min, sum(costEmptyDrive[[posCraneInitial,s]]*dInit[s] for s in L[order[1]]) + sum(costEmptyDrive[[r,s]]*dEmpty[k,r,s] for k = 2:NBar for r in E[order[k-1]] for s in L[order[k]]) + sum(costLoadedDrive[[s,r]]*dLoaded[k,s,r] for k = 1:NBar for s in L[order[k]] for r in E[order[k]]) + sum(beta[z]*finalZ[r,z] for r in intersect(SE,SB) for z = heightsTilde[r]+1:Z));
    @constraint(LBModel, uniquenessMove[k = 1:NBar], sum(w[k,s] for s in L[order[k]]) == 1);
    @constraint(LBModel, constDInit, sum(dInit[s] for s in L[order[1]]) == 1);
    @constraint(LBModel, constEmptyUnique[k = 2:NBar], sum(dEmpty[k,r,s] for r in E[order[k-1]] for s in L[order[k]]) == 1);
    @constraint(LBModel, constLoadedUnique[k = 1:NBar], sum(dLoaded[k,s,r] for s in L[order[k]] for r in E[order[k]]) == 1);
    @constraint(LBModel, conservFlowInit_1[s in L[order[1]]], sum(dLoaded[1,s,r] for r in E[order[1]]) - dInit[s] == 0);
    @constraint(LBModel, conservFlow_1[k = 2:NBar, s in L[order[k]]], sum(dLoaded[k,s,r] for r in E[order[k]]) - sum(dEmpty[k,r,s] for r in E[order[k-1]]) == 0);
    @constraint(LBModel, conservFlow_2[k = 2:NBar, r in E[order[k-1]]], sum(dEmpty[k,r,s] for s in L[order[k]]) - sum(dLoaded[k-1,s,r] for s in L[order[k-1]]) == 0);
    @constraint(LBModel, finalHeightUniq[r in intersect(SE,SB)], sum(finalZ[r,z] for z = heightsTilde[r]:Z) == 1);
    @constraint(LBModel, finalHeightConst[r in intersect(SE,SB)], sum(z * finalZ[r,z] for z = heightsTilde[r]:Z) - sum(dLoaded[k,s,r] for k in validCycles[r] for s in L[order[k]]) == heightsTilde[r]);
    @constraint(LBModel, constW_1[k = 1:NBar, s in L[order[k]]], sum(dLoaded[k,s,r] for r in E[order[k]]) - w[k,s] >= 0);
    @constraint(LBModel, constW_2_changed[r in SR, k = intersect(validCycles[r],1:reverseOrder[minStack[r]]),s in L[order[k]]], dLoaded[k,s,r] == 0);
    # TT = STDOUT;
    # redirect_stdout();
    status = solve(LBModel);
    # redirect_stdout(TT);
    ObjLB = getobjectivevalue(LBModel);
    WL = getvalue(w);
    SLB = Dict{Int64,Array{Int64}}();
    WLB = Dict{Int64,Array{Float64}}();
    for k = 1:NBar
        for s in L[order[k]]
            if WL[k,s] > 0
                if k in keys(SLB)
                    push!(SLB[k],s);
                    push!(WLB[k],WL[k,s]);
                else
                    SLB[k] = [s];
                    WLB[k] = [WL[k,s]];
                end
            end
        end
    end
    return (ObjLB,SLB,WLB)
end

function upperBound(order,sigma,NBar,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack)
    reverseOrder = Array{Int64}(NBar);
    for k = 1:NBar
        reverseOrder[order[k]] = k;
    end
    EBar = Dict{Int64,Array{Int64}}();
    costBar = Dict{Array{Int64},Float64}();
    unionEBar = [];
    for k = 1:NBar
        EBar[order[k]] = [];
        for r in E[order[k]]
            if r in SR
                if reverseOrder[minStack[r]] < k
                    push!(EBar[order[k]],r);
                    if k < NBar
                        costBar[[k,r]] = costLoadedDrive[[sigma[k],r]] + costEmptyDrive[[r,sigma[k+1]]];
                    else
                        costBar[[k,r]] = costLoadedDrive[[sigma[k],r]];
                    end
                end
            else
                push!(EBar[order[k]],r);
                if k < NBar
                    costBar[[k,r]] = costLoadedDrive[[sigma[k],r]] + costEmptyDrive[[r,sigma[k+1]]];
                else
                    costBar[[k,r]] = costLoadedDrive[[sigma[k],r]];
                end
            end
        end
        unionEBar = union(unionEBar,EBar[order[k]]);
    end
    SBBar = intersect(SE,SB,unionEBar);
    KBar = Dict{Int64,Array{Int64}}();
    for r in SBBar
        KBar[r] = [];
        for k = 1:NBar
            if r in EBar[order[k]]
                push!(KBar[r],k);
            end
        end
    end
    UBModel = Model(solver = GurobiSolver(OutputFlag=0));
    @variable(UBModel, 0 <= dBar[k = 1:NBar, r in EBar[order[k]]] <= 1);
    @variable(UBModel, 0 <= finalZBar[r in SBBar,z = heightsTilde[r]:Z] <= 1);
    @objective(UBModel, Min, sum(costBar[[k,r]]*dBar[k,r] for k = 1:NBar for r in EBar[order[k]]) + sum(beta[z]*finalZBar[r,z] for r in SBBar for z = heightsTilde[r]+1:Z));
    @constraint(UBModel, constDUniqBar[k = 1:NBar], sum(dBar[k,r] for r in EBar[order[k]]) == 1);
    @constraint(UBModel, finalZUniqBar[r in SBBar], sum(finalZBar[r,z] for z = heightsTilde[r]:Z) == 1);
    @constraint(UBModel, finalZConstBar[r in SBBar], sum(z * finalZBar[r,z] for z = heightsTilde[r]:Z) - sum(dBar[k,r] for k in KBar[r]) == heightsTilde[r]);
    # TT = STDOUT;
    # redirect_stdout();
    status = solve(UBModel);
    # redirect_stdout(TT);
    ObjUB = costEmptyDrive[[posCraneInitial,sigma[1]]] + getobjectivevalue(UBModel);
    D = getvalue(dBar);
    emptyStack = Array{Int64,1}(NBar);
    for k = 1:NBar
        for r in EBar[order[k]]
            if round(D[k,r]) == 1
                emptyStack[k] = r;
            end
        end
    end
    return (ObjUB,emptyStack);
end

function orderEval(order,NSamples,bestUB,NBar,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack)
    (ObjLB,SLB,WLB) = lowerBound(order,NBar,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack);
    visitedSigma = Dict{Int64,Array{Int64}}();
    sigma = Array{Int64,1}(NBar);
    emptyStack = Array{Int64,1}(NBar);
    ObjUB = Inf;
    if ObjLB < bestUB
        NSigma = 1;
        for k = 1:NBar
            NSigma = NSigma * length(SLB[k]);
        end
        visitedSigma = Dict{Int64,Array{Int64}}();
        emptyStack = Array{Int64,1}(NBar);
        ObjUB = Inf;
        for i=1:min(NSamples,NSigma)
            sigmaLoc = Array{Int64,1}(NBar);
            for k = 1:NBar
                sigmaLoc[k] = sample(SLB[k], Weights(WLB[k]));
            end
            while sigmaLoc in values(visitedSigma)
                sigmaLoc = Array{Int64,1}(NBar);
                for k = 1:NBar
                    sigmaLoc[k] = sample(SLB[k], Weights(WLB[k]));
                end
            end
            visitedSigma[i] = sigmaLoc;
            (ObjUBLoc,emptyStackLoc) = upperBound(order,sigmaLoc,NBar,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack);
            if ObjUBLoc < ObjUB
                ObjUB = ObjUBLoc;
                sigma = sigmaLoc;
                emptyStack = emptyStackLoc;
            end
        end
    end
    return (ObjLB,ObjUB,sigma,emptyStack)
end

function fullOrder(NBar, N, orderReq, forcedReloc)
    order = zeros(Int64,NBar);
    k = 0;
    for n = 1:N
        for c in forcedReloc[orderReq[n]]
            k = k + 1;
            order[k] = c;
        end
        k = k + 1;
        order[k] = orderReq[n];
    end
    return order;
end

function feasibleOrder(N,orderReq,UnSwapReq,delta)
    isFeasible = true;
    for k = 1:N
        n = orderReq[k];
        isFeasible = (isFeasible & (n - k <= delta[n,1] && k - n <= delta[n,2]));
        for l = k+1:N
            m = orderReq[l];
            if m < n && [m,n] in UnSwapReq
                isFeasible = false;
            end
        end
    end
    return isFeasible;
end

function localSearch(currentOrderReq,startTime,limitOfTime,N,PairsSwap,UnSwapReq,delta,NBar,forcedReloc,NSamples,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack)
    currentOrder = fullOrder(NBar, N, currentOrderReq, forcedReloc);
    (currentLB,currentUB,currentSigma,currentEmptyStack) = orderEval(currentOrder,NSamples,Inf,NBar,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack);
    visitedOrder = Set{Array{Int64}}();
    push!(visitedOrder,currentOrderReq);
    minLocalFound = false;
    while !minLocalFound && time() - startTime <= limitOfTime
        minLocalFound = true;
        randInd = randperm(Int64(N*(N-1)/2));
        indSwap = 1;
        while indSwap <= N*(N-1)/2 && minLocalFound && time() - startTime <= limitOfTime
            orderReq = copy(currentOrderReq);
            orderReq[PairsSwap[randInd[indSwap]][1]], orderReq[PairsSwap[randInd[indSwap]][2]] = orderReq[PairsSwap[randInd[indSwap]][2]], orderReq[PairsSwap[randInd[indSwap]][1]];
            if feasibleOrder(N,orderReq,UnSwapReq,delta) && !(orderReq in visitedOrder)
                order = fullOrder(NBar, N, orderReq, forcedReloc);
                (ObjLB,ObjUB,sigma,emptyStack) = orderEval(order,NSamples,currentUB,NBar,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack);
                if ObjUB < currentUB
                    minLocalFound = false;
                    currentUB = ObjUB;
                    currentOrder = order;
                    currentSigma = sigma;
                    currentEmptyStack = emptyStack;
                end
                push!(visitedOrder,orderReq);
            end
            indSwap = indSwap + 1;
        end
    end
    return (currentUB,currentOrder,currentSigma,currentEmptyStack)
end

function sampleOrder(N,delta,UnSwapReq)
    randOrderReq = Array{Int64,1}(N);
    while true
        fixedReq = zeros(Int64,N);
        for k = 1:N
            reqs = [];
            for n = 1:N
                if fixedReq[n] == 0 && n - k <= delta[n,1] && k - n <= delta[n,2]
                    push!(reqs,n);
                end
            end
            if length(reqs) > 0
                randOrderReq[k] = rand(reqs);
                fixedReq[randOrderReq[k]] = 1;
            else
                break
            end
        end
        if sum(fixedReq) == N && feasibleOrder(N,randOrderReq,UnSwapReq,delta)
            break
        end
    end
    return randOrderReq;
end

function heuristic(nTotalLocal,NSamples,N,blockingCont,NU,NR,NS,NBar,limitOfTime,delta,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack,cstVerticalCost,vZ,printSolutions,outputFolder,realStack,requestsID,nameIOPoint,IOPointsPosition,X,Y,heightsBlock,positionCont,blockID,period)
    startTime = time();
    forcedReloc = Dict{Int64,Array{Int64,1}}();
    for n = 1:N
        forcedReloc[n] = [];
        if n in keys(blockingCont) && blockingCont[n] in NU
            c = blockingCont[n];
            unshift!(forcedReloc[n],c);
            while c in keys(blockingCont) && blockingCont[c] in NU
                c = blockingCont[c];
                unshift!(forcedReloc[n],c);
            end
        end
    end
    PairsSwap = Dict{Int64,Array{Int64}}();
    UnSwapReq = Set{Array{Int64}}();
    ind = 0;
    for g = 1:N-1
        for h = g+1:N
            ind = ind + 1;
            PairsSwap[ind] = [g,h];
            if g in NR && h in NR && L[g][1] == L[h][1]
                push!(UnSwapReq,[g,h]);
            end
        end
    end
    bestUB = Inf;
    bestOrder = Array{Int64}(NBar);
    bestSigma = Array{Int64}(NBar);
    bestEmptyStack = Array{Int64}(NBar);
    nLocal = 0;
    while nLocal <= nTotalLocal && time() - startTime <= limitOfTime
        if nLocal == 0
            currentOrderReq = collect(1:N);
        else
            currentOrderReq = sampleOrder(N,delta,UnSwapReq);
        end
        TT = STDOUT;
        redirect_stdout();
        (currentUB,currentOrder,currentSigma,currentEmptyStack) = localSearch(currentOrderReq,startTime,limitOfTime,N,PairsSwap,UnSwapReq,delta,NBar,forcedReloc,NSamples,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack);
        redirect_stdout(TT);
        if currentUB < bestUB
            bestUB = currentUB;
            bestOrder = currentOrder;
            bestSigma = currentSigma;
            bestEmptyStack = currentEmptyStack;
        end
        nLocal = nLocal + 1;
    end
    if printSolutions
        outputFile = joinpath(outputFolder,string(period,"_Period"));
        f = open(outputFile,"a");
        posCrane = posCraneInitial;
        for k = 1:NBar
            n = bestOrder[k];
            s = bestSigma[k];
            r = bestEmptyStack[k];
            write(f,string("Stage  ", k,"\n"));
            if n in NR
                write(f,string("Empty Drive from ", realStack[posCrane], " to ", realStack[s],"\n"));
                write(f,string("Retrieve container ", requestsID[n], "(", n, ") from ", realStack[s], " to ", nameIOPoint[r],"\n"));
                posCrane = r;
            elseif n in NS
                write(f,string("Empty Drive from ", realStack[posCrane], " to ", realStack[s],"\n"));
                write(f,string("Store container ", requestsID[n], "(", n, ") from ", nameIOPoint[s], " in ", realStack[r],"\n"));
                posCrane = r;
            elseif n in NU
                write(f,string("Empty Drive from ", realStack[posCrane], " to ", realStack[s],"\n"));
                write(f,string("Relocate container ", requestsID[n], " from ", realStack[s], " in ", realStack[r],"\n"));
                posCrane = r;
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
    horizontalCost = 0;
    for k = 1:NBar
        n = bestOrder[k];
        s = bestSigma[k];
        r = bestEmptyStack[k];
        horizontalCost = horizontalCost + costEmptyDrive[[posCraneInitial,s]] + costLoadedDrive[[s,r]];
        if n in NR
            delete!(blockID,[s,heightsBlock[s]]);
            heightsBlock[s] = heightsBlock[s] - 1;
            delete!(positionCont,requestsID[n]);
            posCraneInitial = r;
        elseif n in NS
            heightsBlock[r] = heightsBlock[r] + 1;
            positionCont[requestsID[n]] = [r,heightsBlock[r]];
            blockID[[r,heightsBlock[r]]] = requestsID[n];
            posCraneInitial = r;
        elseif n in NU
            delete!(blockID,[s,heightsBlock[s]]);
            heightsBlock[s] = heightsBlock[s] - 1;
            heightsBlock[r] = heightsBlock[r] + 1;
            positionCont[requestsID[n]] = [r,heightsBlock[r]];
            blockID[[r,heightsBlock[r]]] = requestsID[n];
            posCraneInitial = r;
        end
    end
    verticalCost = cstVerticalCost;
    for r in intersect(SE,SB)
        verticalCost = verticalCost - 1/vZ*(heightsBlock[r])*(heightsBlock[r]+1);
    end
    return (horizontalCost,verticalCost,heightsBlock,positionCont,blockID,posCraneInitial);
end
