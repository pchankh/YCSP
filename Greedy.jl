include("subProblem.jl");



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
        #     Wgiven[[N - remainToStack + 1,moveFrom[N - remainToStack + 1][1],t]] = 1;
        #     orderCont[N - remainToStack + 1] = t;
        #     orderTime[t] = N - remainToStack + 1;
        #     remainToStack = remainToStack - 1;
        # end
    end
end
while remainToStack > 0
    t = t + 1;
    Wgiven[[N - remainToStack + 1,moveFrom[N - remainToStack + 1][1],t]] = 1;
    orderCont[N - remainToStack + 1] = t;
    orderTime[t] = N - remainToStack + 1;
    remainToStack = remainToStack - 1;
end

(Obj,P) = subProblem(Wgiven, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime);

revP = Dict{Float64,Dict{Int64,Vector{Array{Int64}}}}();
for k in keys(P)
    if P[k] in keys(revP)
        i = length(collect(keys(revP[P[k]])));
        revP[P[k]][i+1] = [collect(k),zeros(Int64,3)];
    else
        revP[P[k]] = Dict{Int64,Vector{Array{Int64}}}();
        revP[P[k]][1] = [collect(k),zeros(Int64,3)];
    end
    for l in keys(P)
        if P[k]-P[l] > 0
            if P[k]-P[l] in keys(revP)
                i = length(collect(keys(revP[P[k]-P[l]])));
                revP[P[k]-P[l]][i+1] = [collect(k),collect(l)];
            else
                revP[P[k]-P[l]] = Dict{Int64,Vector{Array{Int64}}}();
                revP[P[k]-P[l]][1] = [collect(k),collect(l)];
            end
        end
    end
end
sortedValues = sort(collect(keys(revP)),rev=true);

improvedSolutionFound = false;
i = 1;
indRem = zeros(3);
indAdd = zeros(3);
numNeighbor = 0;
while !improvedSolutionFound && i <= length(sortedValues)
    val = sortedValues[i];
    j = 1;
    J = length(keys(revP[val]));
    while !improvedSolutionFound && j <= J
        indRem = revP[val][j][1];
        indAdd = revP[val][j][2];
        if indAdd == zeros(3)
            numVar = 0;
            varTriplets = Dict{Int64,Array{Int64}}();
            for m = 1:T
                if m != indRem[1]
                    for s in moveFrom[m]
                        numVar = numVar + 1;
                        varTriplets[numVar] = [m,s,indRem[3]];
                    end
                else
                    for s in moveFrom[m]
                        if s != indRem[2]
                            numVar = numVar + 1;
                            varTriplets[numVar] = [m,s,indRem[3]];
                        end
                    end
                end
            end
            toExplore = randperm(numVar);
            k = 1;
            while !improvedSolutionFound && k <= numVar
                indAdd = varTriplets[toExplore[k]];
                if !(indAdd in keys(P)) && isExchange(indRem,indAdd,orderCont,orderTime,previousContToMove)
                    numNeighbor = numNeighbor + 1;
                    # improvedSolutionFound = true;
                end
                k = k + 1;
            end
        else
            if isExchange(indRem,indAdd,orderCont,orderTime,previousContToMove)
                numNeighbor = numNeighbor + 1;
                # improvedSolutionFound = true;
            end
        end
        j = j + 1;
    end
    i = i + 1;
end

if improvedSolutionFound
    addT = orderCont[indAdd[1]];
    Wgiven[indRem] = 0;
    for s in moveFrom[indAdd[1]]
        Wgiven[[indAdd[1],s,addT]] = 0;
    end
    Wgiven[indAdd] = 1;
    Wgiven[[indRem[1],indRem[2],addT]] = 1;
    orderCont[indRem[1]] = addT;
    orderCont[indAdd[1]] = indRem[3];
    orderTime[addT] = indRem[1];
    orderTime[indRem[3]] = indAdd[1];
end
