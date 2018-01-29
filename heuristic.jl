include("subFunctions.jl");

startTime = now();

println("iteration 1");
bestPermutationProductive = realToClusterOrder[1:N];
orderCont = fullOrderContainers(T, N, bestPermutationProductive, blockingCont);
(bestLBObj,WLB,nonIntegralSolution) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
(bestUBObj,bestStackCont) = UpperBound(bestLBObj,WLB,nonIntegralSolution,orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
println("Best Solution: ", bestLBObj);

vectorChangeOfOrder = zeros(Int64,N);
for o = 1:N
    vectorChangeOfOrder[o] = min(changeOfOrder[typeOfTruck[o]],N);
end
# nPermut = computeNumberOfOrders(N,vectorChangeOfOrder,0,sqrt(2*factorial(N)));
# visitedPerms = collect(2:nPermut);

i = 2;
while i <= 1000 && now() - startTime <= Dates.Second(limitOfTime)
    println("iteration ", i);
    permuOrder = zeros(N);
    if false#nPermut <= sqrt(2*factorial(N))
        index = rand(1:length(visitedPerms));
        encode = visitedPerms[index];
        deleteat!(visitedPerms,index);
        permuOrder = decodeRec(encode,N,vectorChangeOfOrder,1:N);
    else
        permuOrder = randperm(N);
        while isInfeasible(permuOrder,N,vectorChangeOfOrder)
            permuOrder = randperm(N);
        end
    end
    permutationProductive = realToClusterOrder[permuOrder];
    orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);
    (LBObj,WLB,nonIntegralSolution) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
    if LBObj < bestUBObj
        (UBObj,StackCont) = UpperBound(LBObj,WLB,nonIntegralSolution,orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
        if bestUBObj > UBObj
            bestUBObj = UBObj;
            bestLBObj = LBObj;
            bestPermutationProductive = permutationProductive;
            bestStackCont = StackCont;
            println("Best Solution: ", bestLBObj);
        end
    end
    i = i + 1;
end

println("Best Solution ", bestUBObj);
