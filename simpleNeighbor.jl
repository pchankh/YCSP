include("LinearProgram.jl");

numVar = 0;
varTriplets = Dict{Int64,Array{Int64}}();
for m = 1:T
    for s in moveFrom[m]
        numVar = numVar + 1;
        varTriplets[numVar] = [m,s];
    end
end

(Wgiven,orderCont,orderTime,positiveW) = initialSolution(T, moveFrom, n, N, previousContToMove, SI, H, posteriorStacks, heightsInitial);

currentObj = LinearProgram(Wgiven, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime);

println("Starting Objective is: ", currentObj);
improvedSolutionFound = true;
while improvedSolutionFound
    improvedSolutionFound = false;
    explorePositive = randperm(T);
    i = 1;
    while !improvedSolutionFound && i <= T
        indRem = positiveW[explorePositive[i]];
        k = 1;
        exploreNull = randperm(numVar);
        while !improvedSolutionFound && k <= numVar
            indAdd = [varTriplets[exploreNull[k]][1],varTriplets[exploreNull[k]][2],indRem[3]];
            if indAdd != indRem && isExchange(indRem,indAdd,orderCont,orderTime,previousContToMove)
                (Wnew,orderContnew,orderTimenew,positiveWnew) = exchangeIndices(Wgiven,T,indRem,indAdd,moveFrom,orderCont,orderTime,positiveW);
                objNew = LinearProgram(Wnew, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime);
                if objNew < currentObj
                    improvedSolutionFound = true;
                    currentObj = objNew;
                    Wgiven = Wnew;
                    positiveW = positiveWnew;
                    orderCont = orderContnew;
                    orderTime = orderTimenew;
                end
            end
            k = k + 1;
        end
        i = i + 1;
    end
    if improvedSolutionFound
        println("Improved Objective is: ", currentObj);
    else
        print("The local optimal solution has been found");
    end
end
