using JuMP, Clp
################################################################################
########################## SET PROBLEM PARAMETERS ##############################
################################################################################
newPackage = false;
if Pkg.installed("JuMP") == nothing
  Pkg.add("JuMP");
  newPackage = true;
end
if Pkg.installed("Clp") == nothing
  Pkg.add("Clp");
  newPackage = true;
end
## Set the directory to current directory of main.jl
cd(dirname(Base.source_path()));
## Load the side script
include("auxilaryFunctions.jl");
## Loading Packages and setting the working directory
using JuMP, Clp

testing = true;

startTime = time();

################################################################################
############################### LOAD INPUT DATA ################################
################################################################################

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

################################################################################
################################### SOLVER #####################################
################################################################################

println("Progress : 0.00 % |                                                   |");
changeOrderReal = zeros(Int64,N);
for o = 1:N
    changeOrderReal[o] = min(changeOfOrder[typeOfTruck[o]],N);
end
currentPermuOrder = collect(1:N);
currentPermuOrder = checkFeasibilityFCFS(currentPermuOrder,N,realToClusterOrder,stackOf,heightOf,changeOrderReal);
IdCurrentPermuOrder = collect(1:N);
for m = 1:N
    IdCurrentPermuOrder[m] = currentPermuOrder[m]
end
permutationProductive = realToClusterOrder[currentPermuOrder];
orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);
(IdLBObj,nonIntegralSolution,integralSolution,capacityStack,moveWithCont,moveInit,moveWithoutCont,finalHeights) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
(IdUBObj,StackCont,moveWithContBest,moveInitBest,moveWithoutContBest,finalHeightsBest) = UpperBound(T,N,orderCont,nonIntegralSolution,integralSolution,capacityStack,posCraneInitial,costPreMove,anteriorStacks,moveFrom,costMove,SX,posteriorStacks,SY,Z,costToGo,alpha,SR,contMinHeightStack,SB,artificialHeights,moveWithCont,moveInit,moveWithoutCont,finalHeights,IdLBObj);
bestUBObj = IdUBObj;
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
currentUBObj = IdUBObj;
bestOrder = collect(1:N);
bestOrderCont = orderCont;
bestStackCont = StackCont;
nTotalLocal = T;
alreadyVisited = Dict{Array{Int64},Array{Int64}}();
currentPermuOrderLoc = zeros(N);
for ic = 1:N
    currentPermuOrderLoc[ic] = currentPermuOrder[ic];
end
alreadyVisited[currentPermuOrderLoc] = [];
while nLocal <= nTotalLocal && time() - startTime <= limitOfTime
    k_1 = round(max(nLocal/nTotalLocal,(time() - startTime)/limitOfTime)*100,2);
    if k_1 <= 99.99 && k_1 >= 10.00
        k_1 = round(k_1,1);
    end
    k_2 = Int64(ceil(k_1/2));
    if k_1 < 100
        println(string("Progress : ", string(k_1), " % |", "="^k_2 ,">", " "^(50-k_2), "|"));
    else
        println(string("Progress : ", "100  % |", "="^k_2 ,">", " "^(50-k_2), "|"));
    end
    visitExchange = randperm(Int64(N*(N-1)/2));
    foundImprovedNeighbor = false;
    j = 1;
    while j <= N*(N-1)/2 && !foundImprovedNeighbor && time() - startTime <= limitOfTime
        k = PairsSwap[visitExchange[j]][1];
        l = PairsSwap[visitExchange[j]][2];
        if !foundImprovedNeighbor && feasibleSwap(k,l,currentPermuOrder,changeOrderReal,typeOfTruck,realToClusterOrder,stackOf,heightOf) && !(visitExchange[j] in alreadyVisited[currentPermuOrder])
            i = i + 1;
            currentPermuOrder[k], currentPermuOrder[l] = currentPermuOrder[l], currentPermuOrder[k];
            permutationProductive = realToClusterOrder[currentPermuOrder];
            orderCont = fullOrderContainers(T, N, permutationProductive, blockingCont);
            (LBObj,nonIntegralSolution,integralSolution,capacityStack,moveWithCont,moveInit,moveWithoutCont,finalHeights) = LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
            if LBObj < currentUBObj
                (UBObj,StackCont,moveWithCont,moveInit,moveWithoutCont,finalHeights) = UpperBound(T,N,orderCont,nonIntegralSolution,integralSolution,capacityStack,posCraneInitial,costPreMove,anteriorStacks,moveFrom,costMove,SX,posteriorStacks,SY,Z,costToGo,alpha,SR,contMinHeightStack,SB,artificialHeights,moveWithCont,moveInit,moveWithoutCont,finalHeights,LBObj);
                if currentUBObj > UBObj
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
                    end
                    if !(currentPermuOrder in keys(alreadyVisited))
                        currentPermuOrderLoc = zeros(N);
                        for ic = 1:N
                            currentPermuOrderLoc[ic] = currentPermuOrder[ic];
                        end
                        alreadyVisited[currentPermuOrderLoc] = [];
                    end
                else
                    currentPermuOrder[k], currentPermuOrder[l] = currentPermuOrder[l], currentPermuOrder[k];
                    currentPermuOrderLoc = zeros(N);
                    for ic = 1:N
                        currentPermuOrderLoc[ic] = currentPermuOrder[ic];
                    end
                    append!(alreadyVisited[currentPermuOrderLoc],visitExchange[j]);
                end
            else
                currentPermuOrder[k], currentPermuOrder[l] = currentPermuOrder[l], currentPermuOrder[k];
                currentPermuOrderLoc = zeros(N);
                for ic = 1:N
                    currentPermuOrderLoc[ic] = currentPermuOrder[ic];
                end
                append!(alreadyVisited[currentPermuOrderLoc],visitExchange[j]);
            end
        end
        j = j + 1;
    end
    if !foundImprovedNeighbor
        nLocal = nLocal + 1;
        for m = 1:N
            currentPermuOrder[m] = IdCurrentPermuOrder[m];
        end
        currentUBObj = IdUBObj;
    end
end
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

################################################################################
################################### OUTPUT #####################################
################################################################################

printResult(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights);
