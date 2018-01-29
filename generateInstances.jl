
################################################################################
######################## GENERATE AND SAVE SCENARII ############################
################################################################################

## Set how many instances you want to generate
nInstances = 100;

## Set parameters for each instance
X = 7;
Y = 30;
Z = 4;
fillRate = 0.67;
IOPointsPosition = "Asian-right";
nRequests = 1000;
gap = X*Y;








################################################################################
########################### DO NOT CHANGE FROM HERE ############################
################################################################################
## Set the directory to current directory
cd(dirname(Base.source_path()));
## using DataFrames package
using DataFrames

## Saving parameters
nameFolder = joinpath("Inputs",string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap"));
if !ispath(nameFolder)
    mkpath(nameFolder);
end
for instance = 1:nInstances
## Generate Block Heights
    heightsInitial = fill(Z,(Y,X));
    C = X*Y*Z;
    while C > floor(fillRate*(Y * (X*Z - (Z-1))))
        row = rand(1:Y);
        stack = rand(1:X);
        while heightsInitial[row,stack] == 0
            row = rand(1:Y);
            stack = rand(1:X);
        end
        heightsInitial[row,stack] -= 1;
        C -= 1;
    end
## Save Block Heights
    writecsv(joinpath(nameFolder,string(instance,"_Block.csv")), heightsInitial);
## Generate Requests
    scenarioMatrix = Array{Any}(nRequests,4);
    containersInBlock = collect(1:C);
    arrivalReq = Dict{Int64,Int64}();
    for c = 1:C
        arrivalReq[c] = -gap;
    end
    ind = Int64(C);
    for r = 1:nRequests
        if rand() < 1/2
            scenarioMatrix[r,1] = "storage";
            ind = ind + 1;
            arrivalReq[ind] = r;
            push!(containersInBlock,Int64(ind));
            scenarioMatrix[r,4] = Int64(ind);
        else
            scenarioMatrix[r,1] = "retrieval";
            c = rand(1:length(containersInBlock));
            while r - arrivalReq[containersInBlock[c]] < gap
                c = rand(1:length(containersInBlock));
            end
            delete!(arrivalReq,containersInBlock[c]);
            scenarioMatrix[r,4] = containersInBlock[c];
            deleteat!(containersInBlock,c);
        end
        if scenarioMatrix[r,1] == "storage"
            scenarioMatrix[r,2] = "internal";
        else
            scenarioMatrix[r,2] = "external";
        end
        if IOPointsPosition == "Asian-left"
            scenarioMatrix[r,3] = "left"
        elseif IOPointsPosition == "Asian-right"
            scenarioMatrix[r,3] = "right"
        elseif IOPointsPosition == "two-sided"
            if scenarioMatrix[r,2] == "internal"
                scenarioMatrix[r,3] = "right"
            else
                scenarioMatrix[r,3] = "left"
            end
        elseif IOPointsPosition == "Euro"
            if scenarioMatrix[r,2] == "external"
                scenarioMatrix[r,3] = "landside"
            else
                scenarioMatrix[r,3] = "seaside"
            end
        end
    end
## Save Requests
    writecsv(joinpath(nameFolder,string(instance,"_Requests.csv")), scenarioMatrix);
end
