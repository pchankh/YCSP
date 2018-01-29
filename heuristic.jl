include("subFunctions.jl");

blockingCont = Dict{Int64,Array{Int64}}();
for m = 1:n
    blockingCont[m] = [m];
    while previousContToMove[blockingCont[m][length(blockingCont[m])]] != 0
        append!(blockingCont[m],previousContToMove[blockingCont[m][length(blockingCont[m])]]);
    end
end

solvingTime = 10;
bestUBObj = Inf;
bestLBObj = 0;
bestPermutationProductive = zeros(N);
bestStackCont = zeros(T);
i = 1;
startTime = now();
while i <= 1000 && Int64(now() - startTime) <= limitOfTime * 1000
    println("iteration ", i);
    permutationProductive = randperm(N);
    # permutationProductive = [3 4 1 5 6 2 7];

    orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);

    (LBObj,WLB,nonIntegralSolution) = LowerBound(orderCont, H, T, artificialHeights, moveFrom, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, contMinHeightStack, costMove, costPreMove, costToGo, alpha);

    if LBObj < bestUBObj
        (UBObj,StackCont) = UpperBound(T,H,LBObj,artificialHeights,orderCont,moveFrom,WLB,nonIntegralSolution,SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, contMinHeightStack, costMove, costPreMove, costToGo, alpha);
        if bestUBObj > UBObj
            bestUBObj = UBObj;
            bestLBObj = LBObj;
            bestPermutationProductive = permutationProductive;
            bestStackCont = StackCont;
            println("Best Solution: ", UBObj);
        end
    end
    i = i + 1;
end

println("Best Solution ", bestUBObj);
