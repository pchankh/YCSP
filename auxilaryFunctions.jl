
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
        typeOfTruck[m] = productiveMovesData[clusterToRealOrder[m],2];
        toBeLoaded[m] = productiveMovesData[clusterToRealOrder[m],4];
        toRetrieve[m,1] = productiveMovesData[clusterToRealOrder[m],5];
        toRetrieve[m,2] = productiveMovesData[clusterToRealOrder[m],6];
        toRetrieve[m,3] = productiveMovesData[clusterToRealOrder[m],7];
    end
    for m = n+1:N
        typeOfTruck[m] = productiveMovesData[clusterToRealOrder[m],2];
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
        while previousContToMove[blockingCont[m][length(blockingCont[m])]] != 0
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
function printProblem(X,Y,Z,IOPointsPosition,N,n,realToClusterOrder,typeOfTruck,toBeLoaded,toRetrieve,toBeUnloaded,stackCost,rowCost,relocCost,costToGo,heightsInitial,posCraneInitial)
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
                write(f,string(clusterToRealOrder[m],"\t",typeOfTruck[m],"\tdelivery\t",toBeLoaded[m],"\t",toRetrieve[m,1],"\t",toRetrieve[m,2],"\t",toRetrieve[m,3],"\n"));
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
                write(f,string(clusterToRealOrder[m],"\t",typeOfTruck[m],"\tstorage  \t",toBeUnloaded[m-n],"\n"));
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
                            if s in unloadFrom[m] && orderContStack[m,s,t] > 0.9
                                write(f,string("Stack container ", clusterToRealOrder[m]," from ",nameIOPoint[s]," on ", realStack[r],"\n"));
                                if t == T
                                    lastCranePosition[:,2] = realStack[r];
                                end
                            end
                        end
                    elseif s in SB && r in SL
                        for m = 1:n
                            if stackOf[m] == s && orderContStack[m,s,t] > 0.9
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
