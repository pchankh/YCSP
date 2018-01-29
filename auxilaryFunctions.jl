
################################################################################
########################### randomInitialHeightsFun ############################
################################################################################

## This function is used by initializeInitialHeights defined below
## This functions generates random heights for the whole block given R, S and H
## Note that each height is drawn uniformely at random from H-2 to H (to keep the
## utilization rate of the block high)
## Finally, if the number of containers if too high i.e. C > R * (S*H - (H-1))
## we retrieve a container from a random column until C <= R * (S*H - (H-1))
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

################################################################################
########################### initializeInitialHeights ###########################
################################################################################

## This function defines heightsInitial, the heights of stacks in the block.
## If the input is not random, R, S, H are redifined using the input file.
function initializeInitialHeights(randomInitialHeights, R, S, H)
    if randomInitialHeights
        heightsInitial = randomInitialHeightsFun(R,S,H);
    else
        heightsInitial = readcsv("block.csv", header=false, Int64);
        if R != size(heightsInitial,1)
            println("The number of rows was automatically changed from ", R, " to ", size(heightsInitial,1));
            R = size(heightsInitial,1);
        end
        if S != size(heightsInitial,2)
            println("The number of stacks was automatically changed from ", S, " to ", size(heightsInitial,2));
            S = size(heightsInitial,2);
        end
        if H != maximum(heightsInitial)
            println("The number of tiers was automatically changed from ", H, " to ", maximum(heightsInitial));
            H = maximum(heightsInitial);
        end
    end
    return (heightsInitial, R, S, H);
end

################################################################################
############################# initializeMainValues #############################
################################################################################

## This function initializes the basic stack sets
## SB: stacks in the block
## SI: stacks corresponding to IO-points
## Stot: All potential stacks
## This function also creates the map realStack such that
## realStack[s] is a 2 dimensional vector
## realStack[s][1] provides the actual row of stack s
## realStack[s][2] provides the actual stack of stack s
## Finally, it defines IOPoint the map such that
## IOPoint[s] is the set of IO-points which are linked to stack s
## nameIOPoint[s] is for the purpose of printing the solution
function initializeMainValues(R, S, IOPointsPosition)
    SB = 1:R*S;
    if IOPointsPosition == "left-sided" || IOPointsPosition == "right-sided"
        SI = R*S+1:R*S+R;
    elseif IOPointsPosition == "two-sided"
        SI = R*S+1:R*S+2*R;
    elseif IOPointsPosition == "up-and-down"
        SI = R*S+1:R*S+2*S;
    end
    Stot = union(SB,SI);
    realStack = Dict{Int64,Array{Int64}}();
    for s in SB
        realStack[s] = Array{Int64}(1,2);
        realStack[s][1] = ceil(s/S);
        realStack[s][2] = s - (realStack[s][1]-1)*S;
    end
    for s in SI
        realStack[s] = Array{Int64}(1,2);
        sLoc = s - R*S;
        if IOPointsPosition == "left-sided"
            realStack[s][1] = sLoc;
            realStack[s][2] = 0;
        elseif IOPointsPosition == "right-sided"
            realStack[s][1] = sLoc;
            realStack[s][2] = S+1;
        elseif IOPointsPosition == "two-sided"
            if sLoc <= R
                realStack[s][1] = sLoc;
                realStack[s][2] = 0;
            else
                realStack[s][1] = sLoc - R;
                realStack[s][2] = S+1;
            end
        elseif IOPointsPosition == "up-and-down"
            if sLoc <= S
                realStack[s][1] = 0;
                realStack[s][2] = sLoc;
            else
                realStack[s][1] = R + 1;
                realStack[s][2] = sLoc - S;
            end
        end
    end
    IOPoints = Dict{Int64,Array{Int64}}();
    nameIOPoint = Dict{Int64,AbstractString}();
    innerPoints = Dict{Int64,Array{Int64}}();
    for s in SI
        innerPoints[s] = [];
    end
    for s in SB
        if IOPointsPosition == "left-sided" || IOPointsPosition == "right-sided"
            IOPoints[s] = [SI[realStack[s][1]]];
        elseif IOPointsPosition == "two-sided"
            IOPoints[s] = [SI[realStack[s][1]],SI[realStack[s][1]+R]];
        elseif IOPointsPosition == "up-and-down"
            IOPoints[s] = [SI[realStack[s][2]],SI[realStack[s][2]+S]];
        end
        for r in IOPoints[s]
            append!(innerPoints[r],s);
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
                nameIOPoint[s] = string(realStack[s][2]," up");
            else
                nameIOPoint[s] = string(realStack[s][2]," down");
            end
        end
    end
    groupIOPoint = Dict{AbstractString,Array{Int64}}();
    if IOPointsPosition == "left-sided"
        groupIOPoint["left"] = SI;
    elseif IOPointsPosition == "right-sided"
        groupIOPoint["right"] = SI;
    elseif IOPointsPosition == "two-sided"
        groupIOPoint["left"] = SI[1:R];
        groupIOPoint["right"] = SI[R+1:2*R];
        groupIOPoint["left-and-right"] = SI;
    elseif IOPointsPosition == "up-and-down"
        groupIOPoint["up"] = SI[1:S];
        groupIOPoint["down"] = SI[S+1:2*S];
        groupIOPoint["up-and-down"] = SI;
    end
    return (SB, SI, Stot, realStack, IOPoints, innerPoints, nameIOPoint, groupIOPoint);
end

################################################################################
############################ initializeInitialCrane ############################
################################################################################

## This function defines the initial position of the crane in Stot
## If it is randomly taken then it uses randomInitialCraneFun
## Otherwise, it takes the input file and depends on the IOPointsPosition.
## If the position is not working with the configuration of the block, it throws
## a message and it is taken randomly
function initializeInitialCrane(randomInitialCrane, R, S, IOPointsPosition, SI, Stot)
    if randomInitialCrane
        posCraneInitial = Int64(rand(Stot));
    else
        craneInitial = readcsv("cranePosition.csv", header=false, Int64);
        if 1 <= craneInitial[1] && craneInitial[1] <= R && 1 <= craneInitial[2] && craneInitial[2] <= S
            posCraneInitial = S*(craneInitial[1]-1)+craneInitial[2];
        else
            if IOPointsPosition == "left-sided"
                if craneInitial[1] != 0 || craneInitial[2] > R
                    println("The input crane position is not possible with the structure of IO-points");
                    println("We set the position randomly")
                    if craneInitial[1] != 0
                        println("The first component should be equal to 0");
                    end
                    if craneInitial[2] > R
                        println("The second component should be less or equal than", R);
                    end
                    posCraneInitial = Int64(rand(Stot));
                else
                    posCraneInitial = SI[craneInitial[2]];
                end
            elseif IOPointsPosition == "right-sided"
                if craneInitial[1] != S+1 || craneInitial[2] > R
                    println("The input crane position is not possible with the structure of IO-points");
                    println("We set the position randomly")
                    if craneInitial[1] != S+1
                        println("The first component should be equal to ", S+1);
                    end
                    if craneInitial[2] > R
                        println("The second component should be less or equal than", R);
                    end
                    posCraneInitial = Int64(rand(Stot));
                else
                    posCraneInitial = SI[craneInitial[2]];
                end
            elseif IOPointsPosition == "two-sided"
                if (craneInitial[1] != 0 && craneInitial[1] != S+1) || craneInitial[2] > R
                    println("The input crane position is not possible with the structure of IO-points");
                    println("We set the position randomly")
                    if (craneInitial[1] != 0 && craneInitial[1] != S+1)
                        println("The first component should be equal to ", 0, " or ", S+1);
                    end
                    if craneInitial[2] > R
                        println("The second component should be less or equal than", R);
                    end
                    posCraneInitial = Int64(rand(Stot));
                else
                    if craneInitial[1] == 0
                        posCraneInitial = SI[craneInitial[2]];
                    else
                        posCraneInitial = SI[R+craneInitial[2]];
                    end
                end
            elseif IOPointsPosition == "up-and-down"
                if (craneInitial[2] != 0 && craneInitial[2] != R+1) || craneInitial[1] > S
                    println("The input crane position is not possible with the structure of IO-points");
                    println("We set the position randomly")
                    if craneInitial[1] > S
                        println("The first component should be less or equal than", S);
                    end
                    if (craneInitial[2] != 0 && craneInitial[2] != R+1)
                        println("The second component should be equal to ", 0, " or ", R+1);
                    end
                    posCraneInitial = Int64(rand(Stot));
                else
                    if craneInitial[2] == 0
                        posCraneInitial = SI[craneInitial[1]];
                    else
                        posCraneInitial = SI[S+craneInitial[1]];
                    end
                end
            end
        end
    end
    return posCraneInitial;
end

################################################################################
############################## randomRetrievalFun ##############################
################################################################################

## This function generates n retrievals taken at random given the status of the
## block. It also generates on which IO-point the containers can be delivered
## The outputs are toRetrieve which is a 3 dimensional vector
## toRetrieve[m,1] is the row of m
## toRetrieve[m,2] is the stack of m
## toRetrieve[m,3] is the height of m
## and toBeLoaded which is vector of string
## toBeLoaded[m] is the group of IO-points the containers can be delivered
## it can be left, right, left-and-right, up or down.
## In the case of up-and-down, the group is taken at random.
function randomRetrievalFun(n,heightsInitial,IOPointsPosition)
    C = sum(heightsInitial);
    toRetrieveList = sort(shuffle(collect(1:C))[1:n]);
    toRetrieve = Array{Int64}(n,3);
    toBeLoaded = Array{AbstractString}(n,1);
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
    for m = 1:n
        if IOPointsPosition == "left-sided"
            toBeLoaded[m] = "left"
        elseif IOPointsPosition == "right-sided"
            toBeLoaded[m] = "right"
        elseif IOPointsPosition == "two-sided"
            toBeLoaded[m] = "left-and-right"
        elseif IOPointsPosition == "up-and-down"
            if rand(1)[1] <= 0.5
                toBeLoaded[m] = "up"
            else
                toBeLoaded[m] = "down"
            end
        end
    end
    return (toRetrieve, toBeLoaded);
end

################################################################################
############################# initializeRetrievals #############################
################################################################################

## This function outputs
## SR: the set of stacks where at least one retrieval occurs
## SO: the set of stacks in the block where no retrieval occurs
## SL: the set of I-O points in which a delivery can occur
## SY: the set of stacks corresponding to “end” points of moves with containers
## or “start” points of moves wihtout a container
## (stackOf,heightOf,loadOf) are the three inputs for retrievals for the Integer Program
## contStackIOPoint[[r,s]] is the list of containers which are deliverable from s to r
function initializeRetrievals(n, S, SB, heightsInitial, IOPoints, groupIOPoint, IOPointsPosition)
    if randomRetrieval || n == 0
        (toRetrieve, toBeLoaded) = randomRetrievalFun(n,heightsInitial,IOPointsPosition);
    else
        toRetrieveData = readcsv("retrievals.csv", header=false);
        toRetrieve = toRetrieveData[:,1:3];
        toRetrieve = Array{Int64}(toRetrieve);
        toBeLoaded = toRetrieveData[:,4];
        toBeLoaded = Array{AbstractString}(toBeLoaded);
    end
    stackOf = Dict{Int64,Int64}();
    heightOf = Dict{Int64,Int64}();
    loadOf = Dict{Int64,Array{Int64}}();
    SR = Set{Int64}();
    SL = Set{Int64}();
    for m = 1:n
        stackOf[m] = S*(toRetrieve[m,1]-1)+toRetrieve[m,2];
        push!(SR,stackOf[m]);
        heightOf[m] = toRetrieve[m,3];
        loadOf[m] = intersect(groupIOPoint[toBeLoaded[m]],IOPoints[stackOf[m]]);
        for r in loadOf[m]
            push!(SL,r);
        end
    end
    SO = setdiff(SB,SR);
    SY = union(SB,SL);
    contStackIOPoint = Dict{Array{Int64},Array{Int64}}();
    contMinHeightStack = Dict{Int64,Int64}();
    for s in SR
        contMinHeightStack[s] = 0;
        for r in IOPoints[s]
            contStackIOPoint[[s,r]] = [];
        end
    end
    for m = 1:n
        if contMinHeightStack[stackOf[m]] == 0 || heightOf[contMinHeightStack[stackOf[m]]] > heightOf[m]
            contMinHeightStack[stackOf[m]] = m;
        end
        for r in loadOf[m]
            append!(contStackIOPoint[[stackOf[m],r]],m);
        end
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
    previousCont = 0;
    for m = n+1:N
        previousContToMove[m] = previousCont;
        previousCont = m;
    end

    return (T, SR, SO, SL, SY, stackOf, heightOf, loadOf, contStackIOPoint, contMinHeightStack, previousContToMove, toRetrieve, toBeLoaded);
end

################################################################################
############################## initializeStackings #############################
################################################################################

## This function defines the inputs from the stacking moves
## If randomly selected, the selection is made as for the retrievals
## unloadFrom[m]: is the set of IOPoints from which m can be stacked
## SU: I-O points from which a stack can occur
## SX: stacks corresponding to “start” points of moves with containers or “end”
## points of moves wihtout a container
## contStackableIOPoint[s]: the set of indices of containers which can be stacked
## from IO-point s
function initializeStackings(N,n,randomStacking,IOPointsPosition,SR,groupIOPoint)
    toBeUnloaded = Array{AbstractString}(N-n,1);
    if randomStacking
        for m = 1:N-n
            if IOPointsPosition == "left-sided"
                toBeUnloaded[m] = "left"
            elseif IOPointsPosition == "right-sided"
                toBeUnloaded[m] = "right"
            elseif IOPointsPosition == "two-sided"
                toBeUnloaded[m] = "left-and-right"
            elseif IOPointsPosition == "up-and-down"
                if rand(1)[1] <= 0.5
                    toBeUnloaded[m] = "up"
                else
                    toBeUnloaded[m] = "down"
                end
            end
        end
    else
        toBeUnloaded = readcsv("stackings.csv", header=false);
    end
    unloadFrom = Dict{Int64,Array{Int64}}();
    SU = Set{Int64}();
    for m = n+1:N
        unloadFrom[m] = groupIOPoint[toBeUnloaded[m-n]];
        for r in unloadFrom[m]
            push!(SU,r);
        end
    end
    SX = union(SR,SU);
    contStackableIOPoint = Dict{Int64,Array{Int64}}();
    for s in SU
        contStackableIOPoint[s] = [];
    end
    for m = n+1:N
        for r in unloadFrom[m]
            append!(contStackableIOPoint[r],m);
        end
    end
    return (unloadFrom, SU, SX, contStackableIOPoint, toBeUnloaded);
end

################################################################################
################################# antePostStacks ###############################
################################################################################

## This function defines
## anteriorStacks[s] for s in SY: stacks which could potentially be visited before
## an action ending at stack s
## posteriorStacks[s] for s in SX: stacks which could potentially be visited after
## an action starting at stack s
function antePostStacks(SR,SU,SO,SL,IOPoints,innerPoints,contStackIOPoint)
    #SR,SB,SI,SO,IOPoint,correspondingIOPoint)
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
    return (anteriorStacks,posteriorStacks);
end

################################################################################
############################### defineMoveFromStack ############################
################################################################################

## This function defines the moveFrom structure which is defined by
## moveFrom[m] is the stacks from which m could be moved

function defineMoveFromStack(n,N,T,unloadFrom,stackOf)
    moveFrom = Dict{Int64,Array{Int64}}();
    for m = 1:n
        moveFrom[m] = [stackOf[m]];
    end
    for m = n+1:N
        moveFrom[m] = unloadFrom[m];
    end
    for m = N+1:T
        moveFrom[m] = [stackOf[m]];
    end
    return moveFrom;
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
function defineCosts(N, R, S, H, SX, SY, posCraneInitial, posteriorStacks, rowCost, stackCost, relocCost, realStack)
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
    # costToGo = 1/(R*S);
    costToGo = sqrt(N)/(R*S);
    # costToGo = sqrt(T/(R*S));
    # costToGo = T/sqrt(R*S);
    alpha = Array{Float64}(H,1);
    sumInv = 0;
    for h = 1:H
        sumInv += 1/h;
        alpha[h] = h - sumInv;
    end
    return (costMove, costPreMove, costToGo, alpha);
end

################################################################################
################################## printProblem ################################
################################################################################

## This function prints the problem inputs
function printProblem(R,S,H,IOPointsPosition,N,n,toRetrieve,toBeLoaded,toBeUnloaded,stackCost,rowCost,relocCost,costToGo,heightsInitial,posCraneInitial)

    println("-------------------------------");
    println("MAIN VALUES");
    println("Number of rows : ", R);
    println("Number of stack : ", S);
    println("Number of tiers : ", H);
    println("IO-points configuration : ", IOPointsPosition);
    println("-------------------------------");
    println("ACTIONS");
    println("Number of actions : ", N);
    println("Number of retrievals : ", n);
    if n > 0
        println("\tRow\tStack\tTier\tLoadOn",);
        for m=1:n
            println("Cont_",m,"\t",toRetrieve[m,1],"\t",toRetrieve[m,2],"\t",toRetrieve[m,3],"\t",toBeLoaded[m]);
        end
    end
    if N - n > 0
        println("\tUnloadFrom",);
        for m=n+1:N
            println("Cont_",m,"\t"toBeUnloaded[m-n]);
        end
    end
    println("-------------------------------");
    println("COSTS");
    println("Trolley Cost (X) : ", stackCost);
    println("Wheel Cost (Y) : ", rowCost);
    println("Hoisting Cost (Z) : ", relocCost);
    println("Cost to Go : ", costToGo);
    println("-------------------------------");
    println("INITIAL STATE OF BLOCK");
    str = "";
    if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:S
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tRight");
    end
    println(str);
    if IOPointsPosition == "up-and-down"
        println("Seaside");
    end
    for r = 1:R
        str = string("Row_",r);
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        for s = 1:S
            str = string(str,"\t",Int64(round(heightsInitial[r,s])));
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        println(str);
    end
    if IOPointsPosition == "up-and-down"
        println("Landside");
    end
    println("-------------------------------");
    println("INITIAL STATE OF CRANE");
    println("Row\tStack",);
    if posCraneInitial in SB
        craneInitial = Array{Int64}(1,2);
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
        elseif IOPointsPosition == "up-and-down"
            if posCraneInitialLoc <= S
                println("Seaside\t",posCraneInitialLoc);
            else
                posCraneInitialLoc = posCraneInitialLoc - S;
                println("Landside\t",posCraneInitialLoc);
            end
        end
    end
    println("-------------------------------");
end

################################################################################
################################## printResult ################################
################################################################################

function printResult(S,R,N,n,posCraneInitial,IOPointsPosition,heightsInitial,realStack,stackOf,unloadFrom,SB,SU,SL,SX,SY,posteriorStacks,nameIOPoint,T,X,DInit,D,W)

    updatedHeights = zeros(Int64, length(SB),T+1);
    for s in SB
        updatedHeights[s,1] = heightsInitial[realStack[s][1],realStack[s][2]];
    end
    for t = 1:T
        startStack = 0;
        endStack = 0;
        for s in SX
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
    if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tLeft");
    end
    for s = 1:S
        str = string(str,"\tStack_",s);
    end
    if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
        str = string(str,"\tRight");
    end
    println(str);
    if IOPointsPosition == "up-and-down"
        println("Seaside");
    end
    sta = 0;
    for r = 1:R
        str = string("Row_",r);
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        for s = 1:S
            sta += 1;
            str = string(str,"\t",Int64(round(updatedHeights[sta,1])));
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\t");
        end
        println(str);
    end
    if IOPointsPosition == "up-and-down"
        println("Landside");
    end
    for t=1:T
        println("\n");
        println("-----------------------------");
        println("Time ", t);
        if t == 1
            for r in SX
                if DInit[r] > 0.9
                    if r == posCraneInitial
                        if r in SB
                            println("Crane stays at ", realStack[r]);
                        else
                            println("Crane stays at ", nameIOPoint[r]);
                        end
                    else
                        if posCraneInitial in SB && r in SB
                            println("Crane moves from ", realStack[posCraneInitial], " to ", realStack[r]);
                        elseif posCraneInitial in SB && !(r in SB)
                            println("Crane moves from ", realStack[posCraneInitial], " to ", nameIOPoint[r]);
                        elseif !(posCraneInitial in SB) && r in SB
                            println("Crane moves from ", nameIOPoint[posCraneInitial], " to ", realStack[r]);
                        else
                            println("Crane moves from ", nameIOPoint[posCraneInitial], " to ", nameIOPoint[r]);
                        end
                    end
                end
            end
        else
            for s in SY
                for r in SX
                    if D[s,r,t] > 0.9
                        if r == s
                            if s in SB
                                println("Crane stays at ", realStack[s]);
                            else
                                println("Crane stays at ", nameIOPoint[s]);
                            end
                        else
                            if s in SB && r in SB
                                println("Crane moves from ", realStack[s], " to ", realStack[r]);
                            elseif s in SB && !(r in SB)
                                println("Crane moves from ", realStack[s], " to ", nameIOPoint[r]);
                            elseif !(s in SB) && r in SB
                                println("Crane moves from ", nameIOPoint[s], " to ", realStack[r]);
                            else
                                println("Crane moves from ", nameIOPoint[s], " to ", nameIOPoint[r]);
                            end
                        end
                    end
                end
            end
        end
        for s in SX
            for r in posteriorStacks[s]
                if X[s,r,t] > 0.9
                    if s in SB && r in SB
                        println("Relocate container from ", realStack[s], " to ", realStack[r]);
                    elseif s in SU && r in SB
                        for m = n+1:N
                            if s in unloadFrom[m] && W[m,s,t] > 0.9
                                println("Stack container ", m," from ",nameIOPoint[s]," on ", realStack[r]);
                            end
                        end
                    elseif s in SB && r in SL
                        for m = 1:n
                            if stackOf[m] == s && W[m,s,t] > 0.9
                                println("Retrieve container ", m, " from ", realStack[s], " to ", nameIOPoint[r]);
                            end
                        end
                    end
                end
            end
        end
        str = "";
        if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\tLeft");
        end
        for s = 1:S
            str = string(str,"\tStack_",s);
        end
        if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
            str = string(str,"\tRight");
        end
        println(str);
        if IOPointsPosition == "up-and-down"
            println("Seaside");
        end
        sta = 0;
        for r = 1:R
            str = string("Row_",r);
            if IOPointsPosition == "left-sided" || IOPointsPosition == "two-sided"
                str = string(str,"\t");
            end
            for s = 1:S
                sta += 1;
                str = string(str,"\t",Int64(round(updatedHeights[sta,t+1])));
            end
            if IOPointsPosition == "right-sided" || IOPointsPosition == "two-sided"
                str = string(str,"\t");
            end
            println(str);
        end
        if IOPointsPosition == "up-and-down"
            println("Landside");
        end
    end
end
