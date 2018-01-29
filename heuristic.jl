include("subFunctions.jl");

startTime = now();

println("iteration 1");
permutationProductive = realToClusterOrder[1:N];
orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);
(IdLBObj,WLB,nonIntegralSolution) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
(IdUBObj,stackCont) = UpperBound(IdLBObj,WLB,nonIntegralSolution,orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
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
bestStackCont = stackCont;
while nLocal <= 1 && now() - startTime <= Dates.Second(limitOfTime)
    visitExchange = randperm(Int64(N*(N-1)/2));
    foundImprovedNeighbor = false;
    j = 1;
    while j <= N*(N-1)/2 && !foundImprovedNeighbor
        k = PairsSwap[visitExchange[j]][1];
        l = PairsSwap[visitExchange[j]][2];
        if abs(currentPermuOrder[k] - l) <= changeOrderReal[currentPermuOrder[k]] && abs(currentPermuOrder[l] - k) <= changeOrderReal[currentPermuOrder[l]]
            i = i + 1;
            println("iteration ", i);
            currentPermuOrder[k], currentPermuOrder[l] = currentPermuOrder[l], currentPermuOrder[k];
            permutationProductive = realToClusterOrder[currentPermuOrder];
            orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);
            (LBObj,WLB,nonIntegralSolution) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
            if LBObj < currentUBObj
                println("Compute Upper Bound");
                (UBObj,stackCont) = UpperBound(LBObj,WLB,nonIntegralSolution,orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
                if currentUBObj > UBObj
                    println("Improved Neighbor found");
                    foundImprovedNeighbor = true;
                    currentUBObj = UBObj;
                    if bestUBObj > currentUBObj
                        bestOrder = currentPermuOrder;
                        bestUBObj = currentUBObj;
                        bestStackCont = stackCont;
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
