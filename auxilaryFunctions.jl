
################################################################################
############################## loadParametersFun ###############################
################################################################################
## This function loads the parameters entered as inputs in Parameters.csv
function loadParametersFun()
    parametersData = readcsv("Parameters.csv", header=false);
    limitOfTime = parametersData[2,2];
    changeOfOrder = Dict{AbstractString,Int64}();
    changeOfOrder["external"] = parametersData[4,2];
    changeOfOrder["internal"] = parametersData[5,2];
    IOPointsPosition = parametersData[7,2];
    Z = parametersData[9,2];
    stackCost = parametersData[10,2];
    rowCost = parametersData[11,2];
    relocCost = parametersData[12,2];
    return (limitOfTime,changeOfOrder,IOPointsPosition,Z,stackCost,rowCost,relocCost);
end

################################################################################
############################## initialHeightsFun ###############################
################################################################################
## This function defines heightsInitial, the heights of stacks in the block.
## It loads the current configuration of the block from block.csv and compute X
## and Y accordingly.
## If Z is set lower than the heighest stack in block.csv, then Z is increased
## (with a WARNING message).
function initialHeightsFun(Z,testing)
## For non-testing
    if !testing
        heightsInitial = readcsv("block.csv", header=false, Int64);
        X = size(heightsInitial,2);
        Y = size(heightsInitial,1);
        if Z < maximum(heightsInitial)
            println("WARNING: The maximum number of tiers was set too low given the input block configuration");
            println("The number of tiers was automatically changed from ", Z, " to ", maximum(heightsInitial));
            Z = maximum(heightsInitial);
        end
    else
## For testing
        testingData = readcsv("Testing.csv", header=false);
        X = testingData[2,2];
        Y = testingData[3,2];
        heightsInitial = rand(Int64(ceil(Z/2)):Z,Y,X);
        C = sum(heightsInitial);
        while C > Y * (X*Z - (Z-1))
            row = rand(1:Y);
            stack = rand(1:X);
            while heightsInitial[row,stack] == 0
                row = rand(1:Y);
                stack = rand(1:X);
            end
            heightsInitial[row,stack] -= 1;
            C -= 1;
        end
    end
    return (X,Y,Z,heightsInitial);
end

################################################################################
################################ basicSetsFun ##################################
################################################################################
## This function initializes the basic stack sets
## SB: stacks in the block
## SI: stacks corresponding to IO-points
## This function also creates the map realStack such that
## realStack[s] is a 2 dimensional vector
## realStack[s][1] provides the actual row of stack s (Y position)
## realStack[s][2] provides the actual stack of stack s (X position)
function basicSetsFun(X,Y,IOPointsPosition)
    SB = 1:X*Y;
    SI = UnitRange{Int64};
    if IOPointsPosition == "left-sided" || IOPointsPosition == "right-sided"
        SI = X*Y+1:X*Y+Y;
    elseif IOPointsPosition == "two-sided"
        SI = X*Y+1:X*Y+2*Y;
    elseif IOPointsPosition == "up-and-down"
        SI = X*Y+1:X*Y+2*X;
    else
        error("WARNING: IOPointsPosition should either be left-sided, right-sided, two-sided or up-and-down!")
    end
    realStack = Dict{Int64,Array{Int64}}();
    for s in SB
        realStack[s] = Array{Int64}(1,2);
        realStack[s][1] = ceil(s/X);
        realStack[s][2] = s - (realStack[s][1]-1) * X;
    end
    for s in SI
        realStack[s] = Array{Int64}(1,2);
        sLoc = s - X * Y;
        if IOPointsPosition == "left-sided"
            realStack[s][1] = sLoc;
            realStack[s][2] = 0;
        elseif IOPointsPosition == "right-sided"
            realStack[s][1] = sLoc;
            realStack[s][2] = X + 1;
        elseif IOPointsPosition == "two-sided"
            if sLoc <= Y
                realStack[s][1] = sLoc;
                realStack[s][2] = 0;
            else
                realStack[s][1] = sLoc - Y;
                realStack[s][2] = X + 1;
            end
        elseif IOPointsPosition == "up-and-down"
            if sLoc <= X
                realStack[s][1] = 0;
                realStack[s][2] = sLoc;
            else
                realStack[s][1] = Y + 1;
                realStack[s][2] = sLoc - X;
            end
        end
    end
    return (SB,SI,realStack);
end

################################################################################
################################# IOPointsFun ##################################
################################################################################
## This function defines IOPoint the map such that
## IOPoint[s] is the set of IO-points which are linked to stack s
## innerPoints[r] is the set of points in the block served by IOpoints r
## nameIOPoint[s] is for the purpose of printing the solution
function IOPointsFun(SI,SB,realStack)
    IOPoints = Dict{Int64,Array{Int64}}();
    nameIOPoint = Dict{Int64,AbstractString}();
    innerPoints = Dict{Int64,Array{Int64}}();
    for s in SB
        if IOPointsPosition == "left-sided" || IOPointsPosition == "right-sided"
            IOPoints[s] = [SI[realStack[s][1]]];
        elseif IOPointsPosition == "two-sided"
            IOPoints[s] = [SI[realStack[s][1]],SI[realStack[s][1]+Y]];
        elseif IOPointsPosition == "up-and-down"
            IOPoints[s] = [SI[realStack[s][2]],SI[realStack[s][2]+X]];
        end
        for r in IOPoints[s]
            if r in keys(innerPoints)
                append!(innerPoints[r],s);
            else
                innerPoints[r] = [s];
            end
        end
    end
    for s in SI
        if IOPointsPosition == "left-sided"
            nameIOPoint[s] = string(realStack[s][1]," left");
        elseif IOPointsPosition == "right-sided"
            nameIOPoint[s] = string(realStack[s][1]," right");
        elseif IOPointsPosition == "two-sided"
            if realStack[s][2] == 0
                nameIOPoint[s] = string(realStack[s][1]," left");
            else
                nameIOPoint[s] = string(realStack[s][1]," right");
            end
        elseif IOPointsPosition == "up-and-down"
            if realStack[s][1] == 0
                nameIOPoint[s] = string("up ",realStack[s][2]);
            else
                nameIOPoint[s] = string("down ",realStack[s][2]);
            end
        end
    end
    groupIOPoint = Dict{AbstractString,Array{Int64}}();
    if IOPointsPosition == "left-sided"
        groupIOPoint["left"] = SI;
    elseif IOPointsPosition == "right-sided"
        groupIOPoint["right"] = SI;
    elseif IOPointsPosition == "two-sided"
        groupIOPoint["left"] = SI[1:Y];
        groupIOPoint["right"] = SI[Y+1:2*Y];
        groupIOPoint["left-and-right"] = SI;
    elseif IOPointsPosition == "up-and-down"
        groupIOPoint["up"] = SI[1:X];
        groupIOPoint["down"] = SI[X+1:2*X];
        groupIOPoint["up-and-down"] = SI;
    end
    return (IOPoints,innerPoints,nameIOPoint,groupIOPoint);
end

################################################################################
############################## cranePositionFun ################################
################################################################################
## This function defines the initial position of the crane in Stot
## It takes the input file cranePosition.csv and depends on the IOPointsPosition.
## If the position is not working with the configuration of the block, it throws
## an error.
function cranePositionFun(X,Y,IOPointsPosition,testing,SB,SI)
## For non-testing
    if !testing
        craneData = readcsv("cranePosition.csv", header=false);
        craneInitial = [craneData[1,2],craneData[2,2]];
        if typeof(craneInitial[1])==Int64 && typeof(craneInitial[2])==Int64
            if 1 <= craneInitial[1] && craneInitial[1] <= Y && 1 <= craneInitial[2] && craneInitial[2] <= X
                posCraneInitial = X * (craneInitial[1]-1) + craneInitial[2];
            else
                error("WARNING: The position provided in cranePosition.csv is not valid!");
            end
        else
            if IOPointsPosition == "left-sided"
                if typeof(craneInitial[1])==Int64 && typeof(craneInitial[2])==SubString{String}
                    if 1 <= craneInitial[1] && craneInitial[1] <= Y && craneInitial[2] == "left"
                        posCraneInitial = SI[craneInitial[1]];
                    else
                        error("WARNING: The position provided in cranePosition.csv is not valid!");
                    end
                else
                    error("WARNING: The position provided in cranePosition.csv is not valid!");
                end
            elseif IOPointsPosition == "right-sided"
                if typeof(craneInitial[1])==Int64 && typeof(craneInitial[2])==SubString{String}
                    if 1 <= craneInitial[1] && craneInitial[1] <= Y && craneInitial[2] == "right"
                        posCraneInitial = SI[craneInitial[1]];
                    else
                        error("WARNING: The position provided in cranePosition.csv is not valid!");
                    end
                else
                    error("WARNING: The position provided in cranePosition.csv is not valid!");
                end
            elseif IOPointsPosition == "two-sided"
                if typeof(craneInitial[1])==Int64 && typeof(craneInitial[2])==SubString{String}
                    if 1 <= craneInitial[1] && craneInitial[1] <= Y && craneInitial[2] == "left"
                        posCraneInitial = SI[craneInitial[1]];
                    elseif 1 <= craneInitial[1] && craneInitial[1] <= Y && craneInitial[2] == "right"
                        posCraneInitial = SI[craneInitial[1]+Y];
                    else
                        error("WARNING: The position provided in cranePosition.csv is not valid!");
                    end
                else
                    error("WARNING: The position provided in cranePosition.csv is not valid!");
                end
            elseif IOPointsPosition == "up-and-down"
                if typeof(craneInitial[1])==SubString{String} && typeof(craneInitial[2])==Int64
                    if craneInitial[1] == "up" && 1 <= craneInitial[2] && craneInitial[2] <= X
                        posCraneInitial = SI[craneInitial[2]];
                    elseif craneInitial[1] == "down" && 1 <= craneInitial[2] && craneInitial[2] <= X
                        posCraneInitial = SI[craneInitial[2]+X];
                    else
                        error("WARNING: The position provided in cranePosition.csv is not valid!");
                    end
                else
                    error("WARNING: The position provided in cranePosition.csv is not valid!");
                end
            end
        end
## For testing
    else
        posCraneInitial = Int64(rand(1:length(SB)+length(SI)));
    end
    return (posCraneInitial);
end

################################################################################
############################# productionMovesFun ###############################
################################################################################
## This function initialize all the features of production moves.
## The outputs are
## 1) toRetrieve which is a 3 dimensional vector
## toRetrieve[m,1] is the row of m
## toRetrieve[m,2] is the stack of m
## toRetrieve[m,3] is the height of m
## 2) toBeLoaded which is vector of string
## toBeLoaded[m] is the group of IO-points the containers can be delivered
## it can be left, right, up or down.
## For testing and in the case of left-and-right and  up-and-down, the group is
## taken based on the type of trucks for the container.
## 3)
function productionMovesFun(testing,X,Y,heightsInitial)
    productiveMovesData = Array{Any,2}(0,0);
## For non-testing
    if !testing
        productiveMovesData = readcsv("productiveMoves.csv", header=false);
        N = size(productiveMovesData,1);
        realToClusterOrder = zeros(Int64,N);
        clusterToRealOrder = zeros(Int64,N);
        n = 0;
        for m = 1:N
            if productiveMovesData[m,3] == "delivery"
                n = n + 1;
                realToClusterOrder[m] = n;
                clusterToRealOrder[n] = m;
            end
        end
        r = 0;
        for m = 1:N
            if productiveMovesData[m,3] == "storage"
                r = r + 1;
                realToClusterOrder[m] = n + r;
                clusterToRealOrder[n + r] = m;
            end
        end
## For testing
    else
        testingData = readcsv("Testing.csv", header=false);
        N = testingData[4,2];
        n = testingData[5,2];
        productiveMovesData = Array{Any,2}(N,7);
        clusterToRealOrder = randperm(N);
        realToClusterOrder = zeros(Int64,N);
        C = sum(heightsInitial);
        toRetrieveList = shuffle(collect(1:C))[1:n];
        for m = 1:N
            c = clusterToRealOrder[m];
            realToClusterOrder[c] = m;
            productiveMovesData[c,1] = c;
            if rand() < 0.5
                productiveMovesData[c,2] = "internal";
            else
                productiveMovesData[c,2] = "external";
            end
            if m <= n
                productiveMovesData[c,3] = "delivery";
            else
                productiveMovesData[c,3] = "storage";
            end
            if IOPointsPosition == "left-sided"
                productiveMovesData[c,4] = "left"
            elseif IOPointsPosition == "right-sided"
                productiveMovesData[c,4] = "right"
            elseif IOPointsPosition == "two-sided"
                if productiveMovesData[c,2] == "internal"
                    productiveMovesData[c,4] = "right"
                else
                    productiveMovesData[c,4] = "left"
                end
            elseif IOPointsPosition == "up-and-down"
                if productiveMovesData[c,2] == "internal"
                    productiveMovesData[c,4] = "up"
                else
                    productiveMovesData[c,4] = "down"
                end
            end
            if m <= n
                d = 0;
                for r = 1:Y
                    for s = 1:X
                        for h = 1:heightsInitial[r,s]
                            d += 1;
                            if toRetrieveList[m] == d
                                productiveMovesData[c,5] = r;
                                productiveMovesData[c,6] = s;
                                productiveMovesData[c,7] = h;
                            end
                        end
                    end
                end
            else
                productiveMovesData[c,5] = "";
                productiveMovesData[c,6] = "";
                productiveMovesData[c,7] = "";
            end
        end
    end
## For non-testing
    typeOfTruck = Array{AbstractString}(N,1);
    toBeLoaded = Array{AbstractString}(n,1);
    toRetrieve = Array{Int64}(n,3);
    toBeUnloaded = Array{AbstractString}(N-n,1);
    for m = 1:n
        typeOfTruck[clusterToRealOrder[m]] = productiveMovesData[clusterToRealOrder[m],2];
        toBeLoaded[m] = productiveMovesData[clusterToRealOrder[m],4];
        toRetrieve[m,1] = productiveMovesData[clusterToRealOrder[m],5];
        toRetrieve[m,2] = productiveMovesData[clusterToRealOrder[m],6];
        toRetrieve[m,3] = productiveMovesData[clusterToRealOrder[m],7];
    end
    for m = n+1:N
        typeOfTruck[clusterToRealOrder[m]] = productiveMovesData[clusterToRealOrder[m],2];
        toBeUnloaded[m-n] = productiveMovesData[clusterToRealOrder[m],4];
    end
    return (N,n,realToClusterOrder,clusterToRealOrder,typeOfTruck,toBeLoaded,toRetrieve,toBeUnloaded);
end

################################################################################
############################## retrievalsBasicsFun #############################
################################################################################
## This function outputs
## (stackOf,heightOf,loadOf) are the three inputs for retrievals for the Integer Program
## SR: the set of stacks where at least one retrieval occurs
## SO: the set of stacks in the block where no retrieval occurs
## SL: the set of I-O points in which a delivery can occur
## SY: the set of stacks corresponding to “end” points of moves with containers
## or “start” points of moves wihtout a container
function retrievalsBasicsFun(X,n,SB,toBeLoaded,toRetrieve,groupIOPoint,IOPoints)
    stackOf = Dict{Int64,Int64}();
    heightOf = Dict{Int64,Int64}();
    loadOf = Dict{Int64,Array{Int64}}();
    SR = Set{Int64}();
    SL = Set{Int64}();
    for m = 1:n
        stackOf[m] = X * (toRetrieve[m,1]-1) + toRetrieve[m,2];
        push!(SR,stackOf[m]);
        heightOf[m] = toRetrieve[m,3];
        loadOf[m] = intersect(groupIOPoint[toBeLoaded[m]],IOPoints[stackOf[m]]);
        for r in loadOf[m]
            push!(SL,r);
        end
    end
    SO = setdiff(SB,SR);
    SY = union(SB,SL);
    return (stackOf,heightOf,loadOf,SR,SO,SL,SY);
end

################################################################################
############################## reshufflesBasicsFun #############################
################################################################################
## This function computes
## contMinHeightStack[s]: container m (=1:n) needed to be retrieved and the
## bottommost of such containers.
## artificialHeights[s]: height after each retrieval has been done and no stacking
## blockingCont[m]: containers blocking m including itself
function reshufflesBasicsFun(N,n,SR,SO,heightsInitial,realStack,stackOf,heightOf)
    contMinHeightStack = Dict{Int64,Int64}();
    for m = 1:n
        if !(stackOf[m] in keys(contMinHeightStack)) || (stackOf[m] in keys(contMinHeightStack) && heightOf[contMinHeightStack[stackOf[m]]] > heightOf[m])
            contMinHeightStack[stackOf[m]] = m;
        end
    end
    artificialHeights = Dict{Int64,Int64}();
    for s in SR
        artificialHeights[s] = heightOf[contMinHeightStack[s]] - 1;
    end
    for s in SO
        artificialHeights[s] = heightsInitial[realStack[s][1],realStack[s][2]];
    end
    previousContToMove = Dict{Int64,Int64}();
    T = N;
    for s in SR
        nextCont = contMinHeightStack[s];
        for h = heightOf[contMinHeightStack[s]]+1:heightsInitial[realStack[s][1],realStack[s][2]]
            toAssign = false;
            for m = 1:n
                if stackOf[m] == s && heightOf[m] == h
                    previousContToMove[nextCont] = m;
                    nextCont = m;
                    toAssign = true;
                end
            end
            if !(toAssign)
                T = T + 1;
                stackOf[T] = s;
                previousContToMove[nextCont] = T;
                nextCont = T;
            end
        end
        previousContToMove[nextCont] = 0;
    end
    blockingCont = Dict{Int64,Array{Int64}}();
    for m = 1:n
        blockingCont[m] = [m];
        while previousContToMove[blockingCont[m][length(blockingCont[m])]] != 0 && previousContToMove[blockingCont[m][length(blockingCont[m])]] > N
            append!(blockingCont[m],previousContToMove[blockingCont[m][length(blockingCont[m])]]);
        end
    end
    return (T,contMinHeightStack,artificialHeights,previousContToMove,blockingCont);
end

################################################################################
################################ storageBasicsFun ##############################
################################################################################
## This function computes
## unloadFrom[m] is the set of IO-points from which container m can be stacked
## SU : union of unloadFrom[m]
## SX : the set of stacks corresponding to “start” points of moves with containers
## or “end” points of moves wihtout a container
function storageBasicsFun(N,n,toBeUnloaded,groupIOPoint)
    unloadFrom = Dict{Int64,Array{Int64}}();
    SU = Set{Int64}();
    for m = n+1:N
        unloadFrom[m] = groupIOPoint[toBeUnloaded[m-n]];
        for r in unloadFrom[m]
            push!(SU,r);
        end
    end
    SX = union(SR,SU);
    return (unloadFrom,SU,SX);
end

################################################################################
################################# antePostStacks ###############################
################################################################################

## This function defines
## anteriorStacks[s] for s in SY: stacks which could potentially be visited before
## an action ending at stack s
## posteriorStacks[s] for s in SX: stacks which could potentially be visited after
## an action starting at stack s
## moveFrom[m]: stacks from which m could be moved
## posteriorContStacks[[m,s]]: stacks in posteriorStacks[s] where m could be moved
## from s
function antePostStacks(N,T,SR,SU,SO,SL,IOPoints,innerPoints,unloadFrom,loadOf,stackOf)
    contStackIOPoint = Dict{Array{Int64},Array{Int64}}();
    for m = 1:n
        for r in loadOf[m]
            if [stackOf[m],r] in keys(contStackIOPoint)
                append!(contStackIOPoint[[stackOf[m],r]],m);
            else
                contStackIOPoint[[stackOf[m],r]] = [m];
            end
        end
    end
    anteriorStacks = Dict{Int64,Array{Int64}}();
    posteriorStacks = Dict{Int64,Array{Int64}}();
    for s in SR
        anteriorStacks[s] = union(setdiff(SR,s),intersect(IOPoints[s],SU));
        potentialIOPoints = Array{Int64}(0);
        for r in IOPoints[s]
            if length(contStackIOPoint[[s,r]]) > 0
                append!(potentialIOPoints,r);
            end
        end
        posteriorStacks[s] = union(setdiff(SB,s),potentialIOPoints);
    end
    for s in SO
        anteriorStacks[s] = union(SR,intersect(IOPoints[s],SU));
    end
    for s in SL
        anteriorStacks[s] = [];
        for r in intersect(SR,innerPoints[s])
            if length(contStackIOPoint[[r,s]]) > 0
                append!(anteriorStacks[s],r);
            end
        end
    end
    for s in SU
        posteriorStacks[s] = innerPoints[s];
    end
    moveFrom = Dict{Int64,Array{Int64}}();
    posteriorContStacks = Dict{Array{Int64},Array{Int64}}();
    for m = 1:n
        moveFrom[m] = [stackOf[m]];
        for s in moveFrom[m]
            posteriorContStacks[[m,s]] = loadOf[m];
        end
    end
    for m = n+1:N
        moveFrom[m] = unloadFrom[m];
        for s in moveFrom[m]
            posteriorContStacks[[m,s]] = posteriorStacks[s];
        end
    end
    for m = N+1:T
        moveFrom[m] = [stackOf[m]];
        for s in moveFrom[m]
            posteriorContStacks[[m,s]] = intersect(posteriorStacks[s],SB);
        end
    end
    return (anteriorStacks,posteriorStacks,moveFrom,posteriorContStacks);
end

################################################################################
################################### defineCosts ################################
################################################################################
## This function defines the appropriate costs for the Integer Program
## costMove[[s,r]]: cost of moving from stack s to r with a container (counting
## up and down)
## costPreMove[[s,r]]: cost of moving from stack s to r without a container (not
## counting up and down)
## costToGo: the multiplicative factor associated with the expected number of
## blocking containers
## alpha[h]: the expected number of blocking containers in a stack with h containers
function defineCosts(n,X,Y,Z,SX,SY,posCraneInitial,posteriorStacks,rowCost,stackCost,relocCost,realStack)
    costMove = Dict{Array{Int64},Float64}();
    for s in SX
        for r in posteriorStacks[s]
            costMove[[s,r]] = rowCost * abs(realStack[s][1] - realStack[r][1]) + stackCost * abs(realStack[s][2] - realStack[r][2]) + 2 * relocCost;
        end
    end
    costPreMove = Dict{Array{Int64},Float64}();
    for s in union(SY,posCraneInitial)
        for r in SX
            costPreMove[[s,r]] = rowCost * abs(realStack[s][1] - realStack[r][1]) + stackCost * abs(realStack[s][2] - realStack[r][2]);
        end
    end
    costToGo = round((n+1)/(X*Y)*(2*(relocCost+min(rowCost,stackCost))),4);
    alpha = Array{Float64}(Z,1);
    sumInv = 0;
    for h = 1:Z
        sumInv += 1/h;
        alpha[h] = h - sumInv;
    end
    return (costMove,costPreMove,costToGo,alpha);
end

################################################################################
################################## printProblem ################################
################################################################################
## This function prints the problem inputs in a file called inputProblem.txt
function printProblem(X,Y,Z,IOPointsPosition,N,n,clusterToRealOrder,realToClusterOrder,typeOfTruck,toBeLoaded,toRetrieve,toBeUnloaded,stackCost,rowCost,relocCost,costToGo,heightsInitial,posCraneInitial)
    f = open("inputProblem.txt","w")
    write(f,string("-------------------------------\n"));
    write(f,string("MAIN VALUES\n"));
    write(f,string("Number of stack (X) : ", X, "\n"));
    write(f,string("Number of rows (Y) : ", Y, "\n"));
    write(f,string("Number of tiers (Z) : ", Z, "\n"));
    write(f,string("IO-points configuration : ", IOPointsPosition, "\n"));
    write(f,string("-------------------------------", "\n"));
    write(f,string("PRODUCTIVE MOVES", "\n"));
    write(f,string("DELIVERIES", "\n"));
    write(f,string("Number of deliveries: ", n, "\n"));
    if n > 0
        write(f,string("Order\tTruckType\tTypeRequest\tLoadOn\tRow\tStack\tTier\n"));
        for o = 1:N
            m = realToClusterOrder[o];
            if m <= n
                write(f,string(clusterToRealOrder[m],"\t",typeOfTruck[clusterToRealOrder[m]],"\tdelivery\t",toBeLoaded[m],"\t",toRetrieve[m,1],"\t",toRetrieve[m,2],"\t",toRetrieve[m,3],"\n"));
            end
        end
    end
    write(f,string("STORAGE", "\n"));
    write(f,string("Number of storage : ", N-n,"\n"));
    if N - n > 0
        write(f,string("Order\tTruckType\tTypeRequest\tUnloadFrom\n"));
        for o = 1:N
            m = realToClusterOrder[o];
            if m >= n+1 && m <= N
                write(f,string(clusterToRealOrder[m],"\t",typeOfTruck[clusterToRealOrder[m]],"\tstorage  \t",toBeUnloaded[m-n],"\n"));
            end
        end
    end
    write(f,string("-------------------------------\n"));
    write(f,string("COSTS\n"));
    write(f,string("Trolley Cost (X) : ", stackCost,"\n"));
    write(f,string("Wheel Cost (Y) : ", rowCost,"\n"));
    write(f,string("Hoisting Cost (Z) : ", relocCost,"\n"));
    write(f,string("Cost to Go : ", costToGo,"\n"));
    write(f,string("-------------------------------\n"));
    write(f,string("INITIAL STATE OF BLOCK\n"));
    str = "";
    if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:X
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tRight");
    end
    write(f,string(str,"\n"));
    if IOPointsPosition == "up-and-down"
        write(f,string("Seaside\n"));
    end
    for r = 1:Y
        str = string("Row_",r);
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        for s = 1:X
            str = string(str,"\t",Int64(round(heightsInitial[r,s])));
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        write(f,string(str,"\n"));
    end
    if IOPointsPosition == "up-and-down"
        write(f,string("Landside\n"));
    end
    write(f,string("-------------------------------\n"));
    write(f,string("INITIAL STATE OF CRANE\n"));
    write(f,string("Row\tStack\n"));
    if posCraneInitial in SB
        craneInitial = Array{Int64}(1,2);
        craneInitial[1,1] = ceil(posCraneInitial/X);
        craneInitial[1,2] = posCraneInitial - (craneInitial[1,1]-1)*X;
        write(f,string(craneInitial[1,1],"\t",craneInitial[1,2],"\n"));
    else
        posCraneInitialLoc = posCraneInitial - X*Y;
        if IOPointsPosition == "left-sided"
            write(f,string(posCraneInitialLoc,"\tleft\n"));
        elseif IOPointsPosition == "right-sided"
            write(f,string(posCraneInitialLoc,"\tright\n"));
        elseif IOPointsPosition == "two-sided"
            if posCraneInitialLoc <= Y
                write(f,string(posCraneInitialLoc,"\lleft\n"));
            else
                posCraneInitialLoc = posCraneInitialLoc - Y;
                write(f,string(posCraneInitialLoc,"\tright\n"));
            end
        elseif IOPointsPosition == "up-and-down"
            if posCraneInitialLoc <= X
                write(f,string("up\t",posCraneInitialLoc,"\n"));
            else
                posCraneInitialLoc = posCraneInitialLoc - X;
                write(f,string("down\t",posCraneInitialLoc,"\n"));
            end
        end
    end
    write(f,string("-------------------------------\n"));
    close(f);
end

################################################################################
################################## printResult #################################
################################################################################

function printResult(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights)
    f = open("Solution.txt","w")
    write(f,string("-------------------------------\n"));
    write(f,string("BLOCK STATE AT THE BEGINNING\n"));
    str = "";
    if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:X
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tRight");
    end
    write(f,string(str,"\n"));
    if IOPointsPosition == "up-and-down"
        write(f,string("Seaside\n"));
    end
    sta = 0;
    for r = 1:Y
        str = string("Row_",r);
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        for s = 1:X
            sta += 1;
            str = string(str,"\t",Int64(round(heightsInitial[realStack[sta][1],realStack[sta][2]])));
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        write(f,string(str,"\n"));
    end
    if IOPointsPosition == "up-and-down"
        write(f,string("Landside\n"));
    end
    lastCranePosition = Array{Any}(2,2);
    lastCranePosition[1,1] = "Row";
    lastCranePosition[2,1] = "Stack";
    lastCranePosition[1,2] = "TBD";
    lastCranePosition[2,2] = "TBD";
    for t=1:T
        write(f,string("-----------------------------\n"));
        write(f,string("Time ", t,"\n"));
        if t == 1
            for r in SX
                if moveInit[r] > 0.9
                    if r == posCraneInitial
                        if r in SB
                            write(f,string("Crane stays at ", realStack[r],"\n"));
                        else
                            write(f,string("Crane stays at ", nameIOPoint[r],"\n"));
                        end
                    else
                        if posCraneInitial in SB && r in SB
                            write(f,string("Crane moves from ", realStack[posCraneInitial], " to ", realStack[r],"\n"));
                        elseif posCraneInitial in SB && !(r in SB)
                            write(f,string("Crane moves from ", realStack[posCraneInitial], " to ", nameIOPoint[r],"\n"));
                        elseif !(posCraneInitial in SB) && r in SB
                            write(f,string("Crane moves from ", nameIOPoint[posCraneInitial], " to ", realStack[r],"\n"));
                        else
                            write(f,string("Crane moves from ", nameIOPoint[posCraneInitial], " to ", nameIOPoint[r],"\n"));
                        end
                    end
                end
            end
        else
            for s in SY
                for r in SX
                    if moveWithoutCont[s,r,t] > 0.9
                        if r == s
                            if s in SB
                                write(f,string("Crane stays at ", realStack[s],"\n"));
                            else
                                write(f,string("Crane stays at ", nameIOPoint[s],"\n"));
                            end
                        else
                            if s in SB && r in SB
                                write(f,string("Crane moves from ", realStack[s], " to ", realStack[r],"\n"));
                            elseif s in SB && !(r in SB)
                                write(f,string("Crane moves from ", realStack[s], " to ", nameIOPoint[r],"\n"));
                            elseif !(s in SB) && r in SB
                                write(f,string("Crane moves from ", nameIOPoint[s], " to ", realStack[r],"\n"));
                            else
                                write(f,string("Crane moves from ", nameIOPoint[s], " to ", nameIOPoint[r],"\n"));
                            end
                        end
                    end
                end
            end
        end
        for s in SX
            for r in posteriorStacks[s]
                if moveWithCont[s,r,t] > 0.9
                    if s in SB && r in SB
                        write(f,string("Relocate container from ", realStack[s], " to ", realStack[r],"\n"));
                    elseif s in SU && r in SB
                        for m = n+1:N
                            if s in unloadFrom[m] && orderContStack[[m,s,t]] > 0.9
                                write(f,string("Stack container ", clusterToRealOrder[m]," from ",nameIOPoint[s]," on ", realStack[r],"\n"));
                                if t == T
                                    lastCranePosition[:,2] = realStack[r];
                                end
                            end
                        end
                    elseif s in SB && r in SL
                        for m = 1:n
                            if stackOf[m] == s && orderContStack[[m,s,t]] > 0.9
                                write(f,string("Retrieve container ", clusterToRealOrder[m], " from ", realStack[s], " to ", nameIOPoint[r],"\n"));
                                if t == T
                                    lastCranePosition[:,2] = split(nameIOPoint[r]);
                                    if all(isnumber, lastCranePosition[1,2])
                                        lastCranePosition[1,2] = parse(Int64,lastCranePosition[1,2]);
                                    end
                                    if all(isnumber, lastCranePosition[2,2])
                                        lastCranePosition[2,2] = parse(Int64,lastCranePosition[2,2]);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    write(f,string("-----------------------------\n"));
    write(f,string("BLOCK STATE AT THE END\n"));
    str = "";
    if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:X
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tRight");
    end
    write(f,string(str,"\n"));
    if IOPointsPosition == "up-and-down"
        write(f,string("Seaside\n"));
    end
    lastHeight = zeros(Int64,Y,X);
    for s in SB
        h = 0;
        for z = 1:Z
            h = h + z * finalHeights[s,z];
        end
        lastHeight[realStack[s][1],realStack[s][2]] = Int64(round(h));
    end
    for r = 1:Y
        str = string("Row_",r);
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        for s = 1:X
            str = string(str,"\t",lastHeight[r,s]);
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        write(f,string(str,"\n"));
    end
    if IOPointsPosition == "up-and-down"
        write(f,string("Landside\n"));
    end
    close(f);
    writecsv("newBlock.csv",lastHeight);
    writecsv("newCranePosition.csv",lastCranePosition);
end

function costFunction(W,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights)
    # LP = Model(solver = GurobiSolver(OutputFlag = 0));
    LP = Model(solver = ClpSolver());
    #####################################################################
    ############################# Variables #############################
    #####################################################################
    @variable(LP, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(LP, 0 <= dInit[SX] <= 1);
    @variable(LP, 0 <= d[SY,SX,2:T] <= 1);
    @variable(LP, 0 <= finalh[SB,0:Z] <= 1);
    #####################################################################
    ############################# Objective #############################
    #####################################################################
    @objective(LP, Min,sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum(alpha[h]*finalh[s,h] for s in SB for h = 2:Z));
    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################
    @constraint(LP, constW_1[m = 1:T,s in moveFrom[m],t = 1:T], sum(x[s,r,t] for r in posteriorContStacks[[m,s]]) >= W[[m,s,t]]);
    @constraint(LP, constW_2[s in SR,t = 1:T], sum(x[r,s,t] for r in anteriorStacks[s]) <= sum(W[[contMinHeightStack[s],s,u]] for u = 1:t-1));
    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################
    ## Uniqueness of move without a container
    @constraint(LP, constDInit, sum(dInit[r] for r in SX) == 1);
    @constraint(LP, constDUnique[t = 2:T], sum(d[s,r,t] for s in SY for r in SX) == 1);
    ## Uniqueness of move with a container
    @constraint(LP, constXUnique[t = 1:T],sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    ## Conservation of flow (1)
    @constraint(LP, conservFlowInit_1[s in SX], sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
    @constraint(LP, conservFlow_1[s in SX,t = 2:T], sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
    ## Conservation of flow (2)
    @constraint(LP, conservFlow_2[s in SY,t = 2:T], sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
    #####################################################################
    ########################### Final Heights ###########################
    #####################################################################
    ## Final heights
    @constraint(LP, finalHeightConst[s in SB], sum(h * finalh[s,h] for h = 1:Z) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == artificialHeights[s]);
    ## Uniqueness of final height
    @constraint(LP, finalHeightUniq[s in SB], sum(finalh[s,h] for h = 0:Z) == 1);
    TT = STDOUT; # save original STDOUT stream
    redirect_stdout();
    status = solve(LP);
    redirect_stdout(TT);
    Obj = getobjectivevalue(LP);
    moveWithCont = getvalue(x);
    moveInit = getvalue(dInit);
    moveWithoutCont = getvalue(d);
    finalHeights = getvalue(finalh);
    return (Obj,moveWithCont,moveInit,moveWithoutCont,finalHeights);
end

################################################################################
################################### LowerBound #################################
################################################################################

function LowerBound(orderCont,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights)
    # LP = Model(solver = GurobiSolver(OutputFlag = 0));
    LP = Model(solver = ClpSolver());
    #####################################################################
    ############################# Variables #############################
    #####################################################################
    @variable(LP, 0 <= x[s in SX,r in posteriorStacks[s],1:T] <= 1);
    @variable(LP, 0 <= dInit[SX] <= 1);
    @variable(LP, 0 <= d[SY,SX,2:T] <= 1);
    @variable(LP, 0 <= finalh[SB,0:Z] <= 1);
    @variable(LP, 0 <= wStack[m = 1:T, s in moveFrom[m]] <= 1);
    #####################################################################
    ############################# Objective #############################
    #####################################################################
    @objective(LP, Min,sum(costMove[[s,r]]*x[s,r,t] for s in SX for r in posteriorStacks[s] for t = 1:T) + sum(costPreMove[[posCraneInitial,r]]*dInit[r] for r in SX) + sum(costPreMove[[s,r]]*d[s,r,t] for s in SY for r in SX for t = 2:T) + costToGo * sum(alpha[h]*finalh[s,h] for s in SB for h = 2:Z));
    #####################################################################
    ################# Relation between variables x and w ################
    #####################################################################
    @constraint(LP, constW_1[m = 1:T,s in moveFrom[m]], sum(x[s,r,orderCont[m]] for r in posteriorContStacks[[m,s]]) >= wStack[m,s]);
    @constraint(LP, constW_2[s in SR,t = 1:orderCont[contMinHeightStack[s]]], sum(x[r,s,t] for r in anteriorStacks[s]) == 0);
    @constraint(LP, constW_3[m = 1:T], sum(wStack[m,s] for s in moveFrom[m]) == 1);
    #####################################################################
    ##################### Conditions on crane moves #####################
    #####################################################################
    ## Uniqueness of move without a container
    @constraint(LP, constDInit, sum(dInit[r] for r in SX) == 1);
    @constraint(LP, constDUnique[t = 2:T], sum(d[s,r,t] for s in SY for r in SX) == 1);
    ## Uniqueness of move with a container
    @constraint(LP, constXUnique[t = 1:T],sum(x[s,r,t] for s in SX for r in posteriorStacks[s]) == 1);
    ## Conservation of flow (1)
    @constraint(LP, conservFlowInit_1[s in SX], sum(x[s,r,1] for r in posteriorStacks[s]) - dInit[s] == 0);
    @constraint(LP, conservFlow_1[s in SX,t = 2:T], sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in SY) == 0);
    ## Conservation of flow (2)
    @constraint(LP, conservFlow_2[s in SY,t = 2:T], sum(d[s,r,t] for r in SX) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
    #####################################################################
    ########################### Final Heights ###########################
    #####################################################################
    ## Final heights
    @constraint(LP, finalHeightConst[s in SB], sum(h * finalh[s,h] for h = 1:Z) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == artificialHeights[s]);
    ## Uniqueness of final height
    @constraint(LP, finalHeightUniq[s in SB], sum(finalh[s,h] for h = 0:Z) == 1);
    TT = STDOUT; # save original STDOUT stream
    redirect_stdout();
    status = solve(LP);
    redirect_stdout(TT);
    LBObj = getobjectivevalue(LP);
    Wsta = getvalue(wStack);
    moveWithCont = getvalue(x);
    moveInit = getvalue(dInit);
    moveWithoutCont = getvalue(d);
    finalHeights = getvalue(finalh);
    nonIntegralSolution = Dict{Int64,Array{Int64}}();
    integralSolution = Dict{Int64,Int64}();
    capacityStack_1 = Dict{Int64,Float64}();
    for m = 1:T
        for s in moveFrom[m]
            if round(Wsta[m,s],6) < 1 - 0.000001 &&  round(Wsta[m,s],6) > 0.000001
                for r in posteriorStacks[s]
                    if round(moveWithCont[s,r,orderCont[m]],6) > 0.000001
                        if m in keys(nonIntegralSolution)
                            append!(nonIntegralSolution[m],[r]);
                        else
                            nonIntegralSolution[m] = [r];
                        end
                        if r in keys(capacityStack_1)
                            capacityStack_1[r] = capacityStack_1[r] + round(moveWithCont[s,r,orderCont[m]],6);
                        else
                            capacityStack_1[r] = round(moveWithCont[s,r,orderCont[m]],6);
                        end
                    end
                end
            else
                for r in posteriorStacks[s]
                    if round(moveWithCont[s,r,orderCont[m]],6) > 1 - 0.000001
                        integralSolution[m] = r;
                    end
                end
            end
        end
    end
    capacityStack = Dict{Int64,Int64}();
    for r in keys(capacityStack_1)
        capacityStack[r] = ceil(capacityStack_1[r]);
    end
    return (LBObj,nonIntegralSolution,integralSolution,capacityStack,moveWithCont,moveInit,moveWithoutCont,finalHeights);
end

function computeAssignmentCost(invOrder,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove)
    previousSta = posCraneInitial;
    cost = 0;
    for o = 1:T
        Cont = invOrder[o];
        if Cont <= N
            cost = cost + costPreMove[[previousSta,intersect(anteriorStacks[StackCont[Cont]],moveFrom[Cont])[1]]] + costMove[[intersect(anteriorStacks[StackCont[Cont]],moveFrom[Cont])[1],StackCont[Cont]]];
            previousSta = StackCont[Cont];
        end
    end
    return cost;
end

function firstAssignment(T,integralSolution,nonIntegralSolution,capacityStack)
    StackCont = zeros(Int64,T);
    for m in keys(integralSolution)
        StackCont[m] = integralSolution[m];
    end
    SlackRemaining = Dict{Int64,Int64}();
    capacityStackLoc = Dict{Int64,Int64}();
    nContperStack = Dict{Int64,Int64}();
    for m in keys(nonIntegralSolution)
        SlackRemaining[m] = 0;
        for r in nonIntegralSolution[m]
            if r in keys(nContperStack)
                nContperStack[r] = nContperStack[r] + 1;
            else
                nContperStack[r] = 1;
                capacityStackLoc[r] = capacityStack[r];
            end
            SlackRemaining[m] = SlackRemaining[m] + capacityStackLoc[r];
        end
    end
    while length(collect(keys(SlackRemaining))) > 0
        minSlacks = Array{Int64}(0);
        for m in keys(SlackRemaining)
            if SlackRemaining[m] == minimum(collect(values(SlackRemaining)))
                append!(minSlacks,[m]);
            end
        end
        selectedCont = 0;
        selectedStack = 0;
        selectedRatio = 0;
        for m in minSlacks
            for r in nonIntegralSolution[m]
                if capacityStackLoc[r]/nContperStack[r] > selectedRatio
                    selectedRatio = capacityStackLoc[r]/nContperStack[r];
                    selectedStack = r;
                    selectedCont = m;
                end
            end
        end
        capacityStackLoc[selectedStack] = capacityStackLoc[selectedStack] - 1;
        for r in nonIntegralSolution[selectedCont]
            nContperStack[r] = nContperStack[r] - 1;
        end
        StackCont[selectedCont] = selectedStack;
        delete!(SlackRemaining, selectedCont);
    end
    return StackCont;
end

function addDrop(m,StackCont,r)
    StackContLoc = zeros(Int64,length(StackCont));
    for i = 1:length(StackCont)
        StackContLoc[i] = StackCont[i];
    end
    StackContLoc[m] = r;
    return StackContLoc;
end

function pairwiseExchange(m,StackCont,l)
    StackContLoc = zeros(Int64,length(StackCont));
    for i = 1:length(StackCont)
        StackContLoc[i] = StackCont[i];
    end
    StackContLoc[m], StackContLoc[l] = StackContLoc[l], StackContLoc[m];
    return StackContLoc;
end

function LocalImprovement(orderCont,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove,capacityStack,nonIntegralSolution)
    invOrder = zeros(Int64,T);
    for m = 1:T
        invOrder[orderCont[m]] = m;
    end
    currentValue = computeAssignmentCost(invOrder,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove);
    localOptimum = false;
    while !localOptimum
        localOptimum = true;
        bestNewStackCont = zeros(Int64,T);
        bestNewValue = currentValue;
        unusedCapacity = Dict{Int64,Int64}();
        for r in keys(capacityStack)
            unusedCapacity[r] = capacityStack[r];
        end
        for m in keys(nonIntegralSolution)
            unusedCapacity[StackCont[m]] = unusedCapacity[StackCont[m]] - 1;
        end
        for m in keys(nonIntegralSolution)
            for r in nonIntegralSolution[m]
                if r != StackCont[m] && unusedCapacity[r] > 0
                    StackContLoc = addDrop(m,StackCont,r);
                    locValue = computeAssignmentCost(invOrder,StackContLoc,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove);
                    if locValue < bestNewValue
                        bestNewValue = locValue;
                        bestNewStackCont = StackContLoc;
                        localOptimum = false;
                    end
                end
            end
        end
        for m in keys(nonIntegralSolution)
            for l in keys(nonIntegralSolution)
                if m != l && StackCont[m] in nonIntegralSolution[l] && StackCont[l] in nonIntegralSolution[m]
                    StackContLoc = pairwiseExchange(m,StackCont,l);
                    locValue = computeAssignmentCost(invOrder,StackContLoc,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove);
                    if locValue < bestNewValue
                        bestNewValue = locValue;
                        bestNewStackCont = StackContLoc;
                        localOptimum = false;
                    end
                end
            end
        end
        if !localOptimum
            currentValue = bestNewValue;
            for m = 1:N
                StackCont[m] = bestNewStackCont[m];
            end
        end
    end
    return StackCont;
end

function generalizedAssignmentModel(orderCont,T,integralSolution,nonIntegralSolution,capacityStack,posCraneInitial,N,costPreMove,anteriorStacks,moveFrom,costMove)
    StackCont = firstAssignment(T,integralSolution,nonIntegralSolution,capacityStack);
    StackCont = LocalImprovement(orderCont,StackCont,posCraneInitial,T,N,costPreMove,anteriorStacks,moveFrom,costMove,capacityStack,nonIntegralSolution);
    return StackCont;
end

function UpperBound(T,N,orderCont,nonIntegralSolution,integralSolution,capacityStack,posCraneInitial,costPreMove,anteriorStacks,moveFrom,costMove,SX,posteriorStacks,SY,Z,costToGo,alpha,SR,contMinHeightStack,SB,artificialHeights,moveWithCont,moveInit,moveWithoutCont,finalHeights,LBObj)
    if length(keys(integralSolution)) == T
        StackCont = zeros(Int64,T);
        for m in keys(integralSolution)
            StackCont[m] = integralSolution[m];
        end
        UBObj = LBObj;
    else
        StackCont = generalizedAssignmentModel(orderCont,T,integralSolution,nonIntegralSolution,capacityStack,posCraneInitial,N,costPreMove,anteriorStacks,moveFrom,costMove);
        W = Dict{Array{Int64},Int64}();
        for m = 1:T
            for s in moveFrom[m]
                for t = 1:T
                    if t == orderCont[m] && s == intersect(anteriorStacks[StackCont[m]],moveFrom[m])[1]
                        W[[m,s,t]] = 1;
                    else
                        W[[m,s,t]] = 0;
                    end
                end
            end
        end
        (UBObj,moveWithCont,moveInit,moveWithoutCont,finalHeights) = costFunction(W,SX,posteriorStacks,T,SY,Z,moveFrom,costMove,costPreMove,posCraneInitial,costToGo,alpha,SR,contMinHeightStack,anteriorStacks,SB,artificialHeights);
    end
    return (UBObj,StackCont,moveWithCont,moveInit,moveWithoutCont,finalHeights);
end

function fullOrderContainers(T, N, permutationProductive, blockingCont)
    orderCont = zeros(Int64,T);
    t = 0;
    sta = 0;
    for m = 1:N
        if permutationProductive[m] in keys(blockingCont)
            for i = length(blockingCont[permutationProductive[m]]):-1:1
                t = t + 1;
                orderCont[blockingCont[permutationProductive[m]][i]] = t;
            end
        else
            t = t + 1;
            orderCont[permutationProductive[m]] = t;
        end
    end
    return orderCont;
end

function feasibleSwap(k,l,currentPermuOrder,changeOrderReal,typeOfTruck,clusterToRealOrder)
    isSwapPossible = (abs(currentPermuOrder[k] - l) <= changeOrderReal[currentPermuOrder[k]] && abs(currentPermuOrder[l] - k) <= changeOrderReal[currentPermuOrder[l]]);
    if isSwapPossible
        isSwapPossible = !(typeOfTruck[clusterToRealOrder[currentPermuOrder[k]]] == "internal" && typeOfTruck[clusterToRealOrder[currentPermuOrder[l]]] == "internal");
    end
    return isSwapPossible;
end
