################################################################################
################################ readInstances ################################
################################################################################
## This function reads the basic geometry of the problem and defines:
## heightsInitial: the heights of stacks in the block of dimension X,Y,Z. The
## fill rate determines the number of containers C in the block with the the
## following relation C = floor(fillRate*(Y * (X*Z - (Z-1))))
## heightsBlock is the vectorized version of heightsInitial
## SB: stacks in the block
## SI: stacks corresponding to IO-points
## This function also creates the map realStack such that
## realStack[s] is a 2 dimensional vector
## realStack[s][1] provides the actual row of stack s (Y position)
## realStack[s][2] provides the actual stack of stack s (X position)
## nameIOPoint[s] is for the purpose of printing the solution
## groupIOPoint[string] is the sets of IOpoint corresponding to the side string
## posCraneInitial: the initial position of the crane in Stot
## costEmptyDrive[[s,r]] is the cost of having a crane empty drive from s to r
## costLoadedDrive[[s,r]] is the cost of having a crane loaded drive from s to r
## beta[z] is the term define in the paper

function readInstances(nameFolder,instanceNumber,X,Y,Z,fillRate,IOPointsPosition,vXEmpty,vXLoaded,vYEmpty,vYLoaded,vZEmpty,vZLoaded,timeHandling,gamma,N,nPeriods)
    heightsInitial = readcsv(joinpath(nameFolder,string(instanceNumber,"_Block.csv")));
    C = sum(heightsInitial);
    positionCont = Dict{Int64,Array{Int64}}();
    blockID = Dict{Array{Int64},Int64}();
    id = 0;
    for row = 1:Y
        for stack = 1:X
            for z = 1:heightsInitial[row,stack]
                id = id + 1;
                positionCont[id] = [(row-1)*X+stack,z];
                blockID[[(row-1)*X+stack,z]] = id;
            end
        end
    end
    SB = 1:X*Y;
    SI = UnitRange{Int64};
    if IOPointsPosition == "Asian-left" || IOPointsPosition == "Asian-right"
        SI = X*Y+1:X*Y+Y;
    elseif IOPointsPosition == "two-sided"
        SI = X*Y+1:X*Y+2*Y;
    elseif IOPointsPosition == "Euro"
        SI = X*Y+1:X*Y+2*X;
    else
        error("WARNING: IOPointsPosition should either be Asian-left, Asian-right, two-sided or Euro!")
    end
    heightsBlock = Array{Int64}(X*Y);
    realStack = Dict{Int64,Array{Int64}}();
    for s in SB
        realStack[s] = Array{Int64}(1,2);
        realStack[s][1] = ceil(s/X);
        realStack[s][2] = s - (realStack[s][1]-1) * X;
        heightsBlock[s] = heightsInitial[realStack[s][1],realStack[s][2]];
    end
    for s in SI
        realStack[s] = Array{Int64}(1,2);
        sLoc = s - X * Y;
        if IOPointsPosition == "Asian-left"
            realStack[s][1] = sLoc;
            realStack[s][2] = 0;
        elseif IOPointsPosition == "Asian-right"
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
        elseif IOPointsPosition == "Euro"
            if sLoc <= X
                realStack[s][1] = 0;
                realStack[s][2] = sLoc;
            else
                realStack[s][1] = Y + 1;
                realStack[s][2] = sLoc - X;
            end
        end
    end
    nameIOPoint = Dict{Int64,AbstractString}();
    for s in SI
        if IOPointsPosition == "Asian-left"
            nameIOPoint[s] = string(realStack[s][1]," left");
        elseif IOPointsPosition == "Asian-right"
            nameIOPoint[s] = string(realStack[s][1]," right");
        elseif IOPointsPosition == "two-sided"
            if realStack[s][2] == 0
                nameIOPoint[s] = string(realStack[s][1]," left");
            else
                nameIOPoint[s] = string(realStack[s][1]," right");
            end
        elseif IOPointsPosition == "Euro"
            if realStack[s][1] == 0
                nameIOPoint[s] = string("Landside ",realStack[s][2]);
            else
                nameIOPoint[s] = string("Seaside ",realStack[s][2]);
            end
        end
    end
    groupIOPoint = Dict{AbstractString,Array{Int64}}();
    if IOPointsPosition == "Asian-left"
        groupIOPoint["left"] = SI;
    elseif IOPointsPosition == "Asian-right"
        groupIOPoint["right"] = SI;
    elseif IOPointsPosition == "two-sided"
        groupIOPoint["left"] = SI[1:Y];
        groupIOPoint["right"] = SI[Y+1:2*Y];
        groupIOPoint["left-and-right"] = SI;
    elseif IOPointsPosition == "Euro"
        groupIOPoint["lanside"] = SI[1:X];
        groupIOPoint["seaside"] = SI[X+1:2*X];
    end
    posCraneInitial = Int64(floor(Y/2)*X+floor(X/2));
    costEmptyDrive = Dict{Array{Int64},Float64}();
    costLoadedDrive = Dict{Array{Int64},Float64}();
    costVerticalDrive = Dict{Int64,Float64}();
    for s in union(SB,SI)
        for r in union(SB,SI)
            costEmptyDrive[[s,r]] = max(abs(realStack[s][1] - realStack[r][1])/vYEmpty, abs(realStack[s][2] - realStack[r][2])/vXEmpty);
            costLoadedDrive[[s,r]] = max(abs(realStack[s][1] - realStack[r][1])/vYLoaded, abs(realStack[s][2] - realStack[r][2])/vXLoaded);
        end
    end
    vZ = 2/(1/vZEmpty+1/vZLoaded);
    for z = 0:Z
        costVerticalDrive[z] = 2*(Z+1-z)/vZ + timeHandling;
    end
    alpha = Dict{Int64,Float64}();
    alpha[0] = 0;
    sumInv = 0;
    for z = 1:Z
        sumInv += 1/z;
        alpha[z] = z - sumInv;
    end
    beta = Dict{Int64,Float64}();
    for z = 0:Z
        beta[z] = gamma * alpha[z] - z*(z+1)/vZ;
    end
    scenarioMatrix = readcsv(joinpath(nameFolder,string(instanceNumber,"_Requests.csv")));
    scenario = Dict{Int64,Array{Any,2}}();
    for p = 1:nPeriods
        scenario[p] = scenarioMatrix[(p-1)*N+1:p*N,:];
    end
    return (heightsBlock,positionCont,blockID,SB,SI,realStack,nameIOPoint,groupIOPoint,posCraneInitial,costEmptyDrive,costLoadedDrive,vZ,costVerticalDrive,beta,scenario);
end
