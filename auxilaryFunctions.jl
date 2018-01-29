

function printProblem(relocCost, rowCost, stackCost, costToGo, H, R, S, heightsInitial,  posCraneInitial, N, n, toRetrieve, IOPointsPosition)

    println("-------------------------------");
    println("COSTS");
    println("Relocation (Z) : ", relocCost);
    println("Row Move (Y) : ", rowCost);
    println("Stack Move (X) : ", stackCost);
    println("Cost to Go : ", costToGo);
    println("-------------------------------");
    println("MAIN VALUES");
    println("Number of rows : ", R);
    println("Number of stack : ", S);
    println("Number of tiers : ", H);
    println("-------------------------------");
    println("INITIAL STATE OF BLOCK");
    str = "";
    if IOPointsPosition != "right-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:S
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition != "left-sided"
        str = string(str,"\tRight");
    end
    println(str);
    for r = 1:R
        str = string("Row_",r);
        if IOPointsPosition != "right-sided"
            str = string(str,"\t");
        end
        for s = 1:S
            str = string(str,"\t",Int64(round(heightsInitial[r,s])));
        end
        if IOPointsPosition != "left-sided"
            str = string(str,"\t");
        end
        println(str);
    end
    println("-------------------------------");
    println("INITIAL STATE OF CRANE");
    println("Row\tStack",);
    if posCraneInitial <= R*S
        craneInitial = Array{Float64}(1,2);
        craneInitial[1,1] = ceil(posCraneInitial/S);
        craneInitial[1,2] = posCraneInitial - (craneInitial[1,1]-1)*S;
        println(craneInitial[1,1],"\t",craneInitial[1,2]);
    else
        posCraneInitialLoc = posCraneInitial - R*S;
        if IOPointsPosition == "left-sided"
            println(posCraneInitialLoc,"\tLeft");
        elseif IOPointsPosition == "right-sided"
            println(posCraneInitialLoc,"\tRight");
        elseif IOPointsPosition == "two-sided"
            if posCraneInitialLoc <= R
                println(posCraneInitialLoc,"\tLeft");
            else
                posCraneInitialLoc = posCraneInitialLoc - R;
                println(posCraneInitialLoc,"\tRight");
            end
        end
    end
    println("-------------------------------");
    println("ACTIONS");
    println("Number of actions : ", N);
    println("Number of stacks : ", N-n);
    println("Number of retrievals : ", n);
    if n > 0
        println("\tRow\tStack\tTier",);
        for m=1:n
            println("Cont_",m,"\t",toRetrieve[m,1],"\t",toRetrieve[m,2],"\t",toRetrieve[m,3]);
        end
    end
    println("-------------------------------");
end

function randomInitialHeightsFun(R,S,H)
    heightsInitial = rand(H-2:H,R,S);
    C = sum(heightsInitial);
    while C > R * (S*H - (H-1))
        row = rand(1:R);
        stack = rand(1:S);
        while heightsInitial[row,stack] == 0
            row = rand(1:R);
            stack = rand(1:S);
        end
        heightsInitial[row,stack] -= 1;
        C -= 1;
    end
    return heightsInitial;
end

function randomInitialCraneFun(Stot)
    posCraneInitial = rand(Stot);
    return posCraneInitial;
end

function randomNumberOfActionsFun(H,heightsInitial)
    C = sum(heightsInitial);
    N = rand(0:(2*C)/H);
    n = rand(0:N);
    return (N,n);
end

function randomRetrievalFun(n,heightsInitial)
    C = sum(heightsInitial);
    toRetrieveList = sort(shuffle(collect(1:C))[1:n]);
    toRetrieve = Array{Int64}(n,3);
    c = 0;
    i = 1;
    R = size(heightsInitial)[1];
    S = size(heightsInitial)[2];
    for r = 1:R
        for s = 1:S
            for h = 1:heightsInitial[r,s]
                c += 1;
                if i <= n && c == toRetrieveList[i]
                    toRetrieve[i,1] = r;
                    toRetrieve[i,2] = s;
                    toRetrieve[i,3] = h;
                    i += 1;
                end
            end
        end
    end
    return toRetrieve;
end

function initializeMainValues(heightsInitial,IOPointsPosition)
    R = size(heightsInitial)[1];
    S = size(heightsInitial)[2];
    SB = 1:R*S;
    if IOPointsPosition == "two-sided"
        SI = R*S+1:R*S+2*R;
    elseif IOPointsPosition == "left-sided" || IOPointsPosition == "right-sided"
        SI = R*S+1:R*S+R;
    end
    Stot = union(SB,SI);
    return (SB,SI,Stot);
end

function seqStackToRealStack(heightsInitial,R,S,SB,SI,IOPointsPosition)
    realStack = Dict{Int64,Array{Int64}}();
    for s in SB
        realStack[s] = Array{Int64}(1,2);
        realStack[s][1,1] = ceil(s/S);
        realStack[s][1,2] = s - (realStack[s][1,1]-1)*S;
    end
    for s in SI
        realStack[s] = Array{Int64}(1,2);
        sLoc = s - R*S;
        if sLoc <= R
            if IOPointsPosition == "right-sided"
                realStack[s][1,1] = sLoc;
                realStack[s][1,2] = S+1;
            else
                realStack[s][1,1] = sLoc;
                realStack[s][1,2] = 0;
            end
        else
            sLoc = sLoc - R;
            realStack[s][1,1] = sLoc;
            realStack[s][1,2] = S+1;
        end
    end
    return realStack;
end

function dataToRetrieve(toRetrieve,S,n,SB,SI)
    indicesToRetrieve = Dict{Int64,Array{Int64}}();
    SR = Set{Int64}();
    s = 0;
    for m = 1:n
        indicesToRetrieve[m] = [S*(toRetrieve[m,1]-1)+toRetrieve[m,2] toRetrieve[m,3]];
        push!(SR,indicesToRetrieve[m][1]);
    end
    contInStackToRetrieve = Dict{Int64,Array{Int64}}();
    for s in SR
        contInStackToRetrieve[s] = [];
        for m = 1:n
            if indicesToRetrieve[m][1] == s
                push!(contInStackToRetrieve[s], m);
            end
        end
    end
    SO = setdiff(SB,SR);
    SP = union(SR,SI);
    return (SR,SO,SP,indicesToRetrieve,contInStackToRetrieve);
end

function computeAlphas(H)
    alpha = Array{Float64}(H,1);
    sumInv = 0;
    for h = 1:H
        sumInv += 1/h;
        alpha[h] = h - sumInv;
    end
    return alpha;
end

function defineCost(R,S,Stot,SB,SP,SR,realStack,relocCost,rowCost,stackCost)
    costMove = Dict{Array{Int64},Float64}();
    for s in SP
        for r in Stot
            costMove[[s,r]] = 2 * relocCost + rowCost * abs(realStack[s][1] - realStack[r][1]) + stackCost * abs(realStack[s][2] - realStack[r][2]);
        end
    end
    costPreMove = Dict{Array{Int64},Float64}();
    for s in Stot
        for r in SP
            costPreMove[[s,r]] = rowCost * abs(realStack[s][1] - realStack[r][1]) + stackCost * abs(realStack[s][2] - realStack[r][2]);
        end
    end
    return (costMove, costPreMove);
end

function dataIOPoints(S,SI,SB,IOPointsPosition,realStack)
    IOPoint = Dict{Int64,Int64}();
    correspondingIOPoint = Dict{Int64,Array{Int64}}();
    for s in SI
        correspondingIOPoint[s] = [];
    end
    rowIOPoint = Dict{Int64,AbstractString}();
    for s in SB
        if IOPointsPosition == "left-sided" || IOPointsPosition == "right-sided"
            IOPoint[s] = SI[realStack[s][1]];
        elseif IOPointsPosition == "two-sided"
            middle = Int64(ceil(S/2));
            if realStack[s][2] <= middle
                IOPoint[s] = SI[realStack[s][1]];
            else
                IOPoint[s] = SI[realStack[s][1]+R];
            end
        end
        append!(correspondingIOPoint[IOPoint[s]],s);
    end
    for s in SI
        if realStack[s][2] == 0
            rowIOPoint[s] = string(realStack[s][1]," left");
        else
            rowIOPoint[s] = string(realStack[s][1]," right");
        end
    end
    return (IOPoint,correspondingIOPoint,rowIOPoint);
end

function antePostStacks(SR,SB,SI,SO,IOPoint,correspondingIOPoint)
    anteriorStacks = Dict{Int64,Array{Int64}}();
    posteriorStacks = Dict{Int64,Array{Int64}}();
    for s in SI
        anteriorStacks[s] = intersect(SR,correspondingIOPoint[s]);
        posteriorStacks[s] = correspondingIOPoint[s];
    end
    for s in SR
        anteriorStacks[s] = union(setdiff(SR,s),IOPoint[s]);
        posteriorStacks[s] = union(setdiff(SB,s),IOPoint[s]);
    end
    for s in SO
        anteriorStacks[s] = union(SR,IOPoint[s]);
    end
    return (anteriorStacks,posteriorStacks);
end


function timeUpperBound(SR,H,N,n,heightsInitial,indicesToRetrieve,realStack,contInStackToRetrieve)
    T = N - n;
    for s in SR
        minHeight = heightsInitial[realStack[s][1],realStack[s][2]];
        for m in contInStackToRetrieve[s]
            minHeight = min(minHeight,indicesToRetrieve[m][2]);
        end
        T += heightsInitial[realStack[s][1],realStack[s][2]] - minHeight + 1;
    end
    return T;
end

function timeUpperBound_2(N,n,H,indicesToRetrieve)
    T = N - n;
    for m = 1:n
        T += H - indicesToRetrieve[m][2] + 1;
    end
    return T;
end

function printResult(R,S,T,n,SP,SR,SB,Stot,X,Y,Z,A,D,ACT,newHeights,realStack,rowIOPoint,indicesToRetrieve,IOPointsPosition)

    println("-------------------------------");
    println("BLOCK STATE AT TIME 0");
    str = "";
    if IOPointsPosition != "right-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:S
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition != "left-sided"
        str = string(str,"\tRight");
    end
    println(str);
    sta = 0;
    for r = 1:R
        str = string("Row_",r);
        if IOPointsPosition != "right-sided"
            str = string(str,"\t");
        end
        for s = 1:S
            sta += 1;
            str = string(str,"\t",Int64(round(newHeights[sta,0])));
        end
        if IOPointsPosition != "left-sided"
            str = string(str,"\t");
        end
        println(str);
    end
    for t=1:T
        if ACT[t] > 0.9
            println("\n");
            println("-----------------------------");
            println("Time ", t);
            for s in Stot
                for r in SP
                    if D[s,r,t] > 0.9
                        if r == s
                            if s in SB
                                println("Crane stays at ", realStack[s]);
                            else
                                println("Crane stays at ", rowIOPoint[s]);
                            end
                        else
                            if s in SB && r in SB
                                println("Crane moves from ", realStack[s], " to ", realStack[r]);
                            elseif s in SB && !(r in SB)
                                println("Crane moves from ", realStack[s], " to ", rowIOPoint[r]);
                            elseif !(s in SB) && r in SB
                                println("Crane moves from ", rowIOPoint[s], " to ", realStack[r]);
                            else
                                println("Crane moves from ", rowIOPoint[s], " to ", rowIOPoint[r]);
                            end
                        end
                    end
                end
            end
            for s in SR
                for r in SB
                    if X[s,r,t] > 0.9
                        println("Relocate container from ", realStack[s], " to ", realStack[r]);
                    end
                end
            end
            for m=1:n
                if A[m,t-1] - A[m,t] > 0.9
                    println("Retrieve container ", m, " from ", realStack[indicesToRetrieve[m][1]]);
                end
            end
            for s in SB
                if Z[s,t] > 0.9
                    println("Stack a new container on ", realStack[s]);
                end
            end
            str = "";
            if IOPointsPosition != "right-sided"
                str = string(str,"\tLeft");
            end
            for s = 1:S
                str = string(str,"\tStack_",s);
            end
            if IOPointsPosition != "left-sided"
                str = string(str,"\tRight");
            end
            println(str);
            sta = 0;
            for r = 1:R
                str = string("Row_",r);
                if IOPointsPosition != "right-sided"
                    str = string(str,"\t");
                end
                for s = 1:S
                    sta += 1;
                    str = string(str,"\t",Int64(round(newHeights[sta,t])));
                end
                if IOPointsPosition != "left-sided"
                    str = string(str,"\t");
                end
                println(str);
            end
        end
    end
end

function printResult_2(S, R, SR, SB, SI, SP, Stot, T, posteriorStacks, X, D, W, heightsInitial, IOPointsPosition, realStack, rowIOPoint)

    updatedHeights = zeros(Int64, length(SB),T+1);
    for s in SB
        updatedHeights[s,1] = heightsInitial[realStack[s][1],realStack[s][2]];
    end
    for t = 1:T
        startStack = 0;
        endStack = 0;
        for s in SP
            for r in posteriorStacks[s]
                if X[s,r,t] > 0.9
                    if s in SB
                        startStack = s;
                    end
                    if r in SB
                        endStack = r;
                    end
                end
            end
        end
        for s in SB
            if s == startStack
                updatedHeights[s,t+1] = updatedHeights[s,t] - 1;
            elseif s == endStack
                updatedHeights[s,t+1] = updatedHeights[s,t] + 1;
            else
                updatedHeights[s,t+1] = updatedHeights[s,t];
            end
        end
    end
    println("-------------------------------");
    println("BLOCK STATE AT TIME 0");
    str = "";
    if IOPointsPosition != "right-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:S
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition != "left-sided"
        str = string(str,"\tRight");
    end
    println(str);
    sta = 0;
    for r = 1:R
        str = string("Row_",r);
        if IOPointsPosition != "right-sided"
            str = string(str,"\t");
        end
        for s = 1:S
            sta += 1;
            str = string(str,"\t",Int64(round(updatedHeights[sta,1])));
        end
        if IOPointsPosition != "left-sided"
            str = string(str,"\t");
        end
        println(str);
    end
    for t=1:T
        println("\n");
        println("-----------------------------");
        println("Time ", t);
        for s in Stot
            for r in SP
                if D[s,r,t] > 0.9
                    if r == s
                        if s in SB
                            println("Crane stays at ", realStack[s]);
                        else
                            println("Crane stays at ", rowIOPoint[s]);
                        end
                    else
                        if s in SB && r in SB
                            println("Crane moves from ", realStack[s], " to ", realStack[r]);
                        elseif s in SB && !(r in SB)
                            println("Crane moves from ", realStack[s], " to ", rowIOPoint[r]);
                        elseif !(s in SB) && r in SB
                            println("Crane moves from ", rowIOPoint[s], " to ", realStack[r]);
                        else
                            println("Crane moves from ", rowIOPoint[s], " to ", rowIOPoint[r]);
                        end
                    end
                end
            end
        end
        for s in SP
            for r in posteriorStacks[s]
                if X[s,r,t] > 0.9
                    if s in SB && r in SB
                        println("Relocate container from ", realStack[s], " to ", realStack[r]);
                    elseif s in SI && r in SB
                        println("Stack a new container on ", realStack[r]);
                    elseif s in SB && r in SI
                        for m in contInStackToRetrieve[s]
                            if W[m,t] > 0.9
                                println("Retrieve container ", m, " from ", realStack[s]);
                            end
                        end
                    end
                end
            end
        end
        str = "";
        if IOPointsPosition != "right-sided"
            str = string(str,"\tLeft");
        end
        for s = 1:S
            str = string(str,"\tStack_",s);
        end
        if IOPointsPosition != "left-sided"
            str = string(str,"\tRight");
        end
        println(str);
        sta = 0;
        for r = 1:R
            str = string("Row_",r);
            if IOPointsPosition != "right-sided"
                str = string(str,"\t");
            end
            for s = 1:S
                sta += 1;
                str = string(str,"\t",Int64(round(updatedHeights[sta,t+1])));
            end
            if IOPointsPosition != "left-sided"
                str = string(str,"\t");
            end
            println(str);
        end
    end
end
