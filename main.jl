using JuMP, Gurobi, AmplNLWriter, MathProgBase, DataStructures
cd(dirname(Base.source_path()));
include("auxilaryFunctions.jl");
include("subFunctions.jl");

testing = true;

(limitOfTime,changeOfOrder,IOPointsPosition,Z,stackCost,rowCost,relocCost) = loadParametersFun();

(X,Y,Z,heightsInitial) = initialHeightsFun(Z,testing);

(SB,SI,realStack) = basicSetsFun(X,Y,IOPointsPosition);

(IOPoints,innerPoints,nameIOPoint,groupIOPoint) = IOPointsFun(SI,SB,realStack);

(posCraneInitial) = cranePositionFun(X,Y,IOPointsPosition,testing,SB,SI);

(N,n,realToClusterOrder,clusterToRealOrder,typeOfTruck,toBeLoaded,toRetrieve,toBeUnloaded) = productionMovesFun(testing,X,Y,heightsInitial);

(stackOf,heightOf,loadOf,SR,SO,SL,SY) = retrievalsBasicsFun(X,n,SB,toBeLoaded,toRetrieve,groupIOPoint,IOPoints);

(T,contMinHeightStack,artificialHeights,previousContToMove,blockingCont) = reshufflesBasicsFun(N,n,SR,SO,heightsInitial,realStack,stackOf,heightOf);

(unloadFrom,SU,SX) = storageBasicsFun(N,n,toBeUnloaded,groupIOPoint);

(anteriorStacks,posteriorStacks,moveFrom,posteriorContStacks) = antePostStacks(N,T,SR,SU,SO,SL,IOPoints,innerPoints,unloadFrom,loadOf,stackOf);

(costMove,costPreMove,costToGo,alpha) = defineCosts(n,X,Y,Z,SX,SY,posCraneInitial,posteriorStacks,rowCost,stackCost,relocCost,realStack);

printProblem(X,Y,Z,IOPointsPosition,N,n,clusterToRealOrder,realToClusterOrder,typeOfTruck,toBeLoaded,toRetrieve,toBeUnloaded,stackCost,rowCost,relocCost,costToGo,heightsInitial,posCraneInitial);


startTime = now();
println("iteration 1");
permutationProductive = realToClusterOrder[1:N];
orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);
(IdLBObj,nonIntegralSolution,integralSolution,capacityStack,moveWithCont,moveInit,moveWithoutCont,finalHeights) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
(IdUBObj,StackCont,moveWithContBest,moveInitBest,moveWithoutContBest,finalHeightsBest) = UpperBound(T,N,orderCont,nonIntegralSolution,integralSolution,capacityStack,posCraneInitial,costPreMove,anteriorStacks,moveFrom,costMove,SX,posteriorStacks,SY,Z,costToGo,alpha,SR,contMinHeightStack,SB,artificialHeights,moveWithCont,moveInit,moveWithoutCont,finalHeights,IdLBObj);
bestUBObj = IdUBObj;
println("Best Solution: ", bestUBObj);
PairsSwap = Dict{Int64,Array{Int64}}();
ind = 0;
for g = 1:N-1
    for h = g+1:N
        ind = ind + 1;
        PairsSwap[ind] = [g,h];
    end
end
i = 1;
nLocal = 0;
currentPermuOrder = collect(1:N);
changeOrderReal = zeros(Int64,N);
for o = 1:N
    changeOrderReal[o] = min(changeOfOrder[typeOfTruck[o]],N);
end
currentUBObj = IdUBObj;
bestOrder = collect(1:N);
bestOrderCont = orderCont;
bestStackCont = StackCont;
while nLocal <= 10 && now() - startTime <= Dates.Second(limitOfTime)
    visitExchange = randperm(Int64(N*(N-1)/2));
    foundImprovedNeighbor = false;
    j = 1;
    while j <= N*(N-1)/2 && !foundImprovedNeighbor
        k = PairsSwap[visitExchange[j]][1];
        l = PairsSwap[visitExchange[j]][2];
        if feasibleSwap(k,l,currentPermuOrder,changeOrderReal,typeOfTruck,clusterToRealOrder)
            i = i + 1;
            println("iteration ", i);
            currentPermuOrder[k], currentPermuOrder[l] = currentPermuOrder[l], currentPermuOrder[k];
            permutationProductive = realToClusterOrder[currentPermuOrder];
            orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);
            (LBObj,nonIntegralSolution,integralSolution,capacityStack,moveWithCont,moveInit,moveWithoutCont,finalHeights) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
            if LBObj < currentUBObj
                println("Compute Upper Bound");
                (UBObj,StackCont,moveWithCont,moveInit,moveWithoutCont,finalHeights) = UpperBound(T,N,orderCont,nonIntegralSolution,integralSolution,capacityStack,posCraneInitial,costPreMove,anteriorStacks,moveFrom,costMove,SX,posteriorStacks,SY,Z,costToGo,alpha,SR,contMinHeightStack,SB,artificialHeights,moveWithCont,moveInit,moveWithoutCont,finalHeights,LBObj);
                if currentUBObj > UBObj
                    println("Improved Neighbor found");
                    foundImprovedNeighbor = true;
                    currentUBObj = UBObj;
                    if bestUBObj > currentUBObj
                        bestOrder = currentPermuOrder;
                        bestUBObj = currentUBObj;
                        bestOrderCont = orderCont;
                        bestStackCont = StackCont;
                        moveWithContBest = moveWithCont;
                        moveInitBest = moveInit;
                        moveWithoutContBest = moveWithoutCont;
                        finalHeightsBest = finalHeights;
                        println("Best Solution Updated: ", bestUBObj);
                    end
                else
                    currentPermuOrder[k], currentPermuOrder[l] = currentPermuOrder[l], currentPermuOrder[k];
                end
            else
                currentPermuOrder[k], currentPermuOrder[l] = currentPermuOrder[l], currentPermuOrder[k];
            end
        end
        j = j + 1;
    end
    if !foundImprovedNeighbor
        nLocal = nLocal + 1;
        if nLocal <= 2
            println("---------------------------------")
            println("Start the search to the beginning")
            currentPermuOrder = collect(1:N);
            currentUBObj = IdUBObj;
        end
    end
end
println("---------------------------------")
println("Best Solution ", bestUBObj);
moveWithCont = moveWithContBest;
moveInit = moveInitBest;
moveWithoutCont = moveWithoutContBest;
finalHeights = finalHeightsBest;
orderContStack = Dict{Array{Int64},Int64}();
for m = 1:T
    for s in moveFrom[m]
        for t = 1:T
            if t == bestOrderCont[m] && s in anteriorStacks[bestStackCont[m]]
                orderContStack[[m,s,t]] = 1;
            else
                orderContStack[[m,s,t]] = 0;
            end
        end
    end
end
printResult(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights);
