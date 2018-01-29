
################################################################################
########################## SET PROBLEM PARAMETERS ##############################
################################################################################
## Set the directory to current directory of main.jl
cd(dirname(Base.source_path()));
## Loading Packages and setting the working directory
using JuMP, Gurobi, DataFrames, StatsBase
include("readInstances.jl");
include("simulateScenario.jl");
include("IPSolver.jl");
include("heuristic.jl");

## Set the crane parameters Speeds must be in containers/s and timeHandling in sec
vXEmpty = 0.5;
vXLoaded = 0.5;
vYEmpty = 0.37;
vYLoaded = 0.2;
vZEmpty = 0.39;
vZLoaded = 0.20;
timeHandling = 20;

## Set parameters for each instance
X = 4;
Y = 4;
Z = 4;
fillRate = 0.67;
IOPointsPosition = "Asian-right";
nRequests = 60;
gap = X*Y;

## Instance number you want to consider
instanceNumber = 10;

## Say which method you want to use the problem to solver
## Choose limitOfTime,  Select between
## 1) FCFS_Myopic_Restricted (Set restricted = true, gamma = 0 and flexibilityRequests = 0)
## 2) Flexible_Myopic_Restricted (Set restricted = true, gamma = 0 and choose flexibilityRequests)
## 3) Flexible_Restricted (Set restricted = true and choose gamma and flexibilityRequests)
## 4) Flexible (Set restricted = false and choose gamma and flexibilityRequests)
## 5) Heuristic (Choose gamma and flexibilityRequests, NSamples and nTotalLocal)
method = "Heuristic";
restricted = true;
gamma = ceil(Z*(Z-1)*(1/vZEmpty+1/vZLoaded));
flexibilityRequests = Dict{AbstractString,Array{Int64}}();
flexibilityRequests["internal"] = [0 2];
flexibilityRequests["external"] = [2 0];
NSamples = 30;
nTotalLocal = 40;
if method == "FCFS_Myopic_Restricted"
    gamma = 0;
    flexibilityRequests["internal"] = [0 0];
    flexibilityRequests["external"] = [0 0];
elseif method == "Flexible_Myopic_Restricted"
    gamma = 0;
elseif method == "Flexible"
    retricted = false;
end

## Solver time limit
limitOfTime = 60;

## Set the simulation parameters
N = 5;
nPeriods = Int64(floor(nRequests/N));

## Saving parameters
printSolutions = true;


################################################################################
############################ SETTING UP SIMULATION #############################
################################################################################
subFolder = string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap");
nameFolder = joinpath("Inputs",subFolder);
outputFolder = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma));
if !ispath(outputFolder)
    mkpath(outputFolder);
end
(heightsBlock,positionCont,blockID,SB,SI,realStack,nameIOPoint,groupIOPoint,posCraneInitial,costEmptyDrive,costLoadedDrive,vZ,costVerticalDrive,beta,scenario) = readInstances(nameFolder,instanceNumber,X,Y,Z,fillRate,IOPointsPosition,vXEmpty,vXLoaded,vYEmpty,vYLoaded,vZEmpty,vZLoaded,timeHandling,gamma,N);
performanceMetrics = DataFrame(Total_Cost = Array{Float64}(nPeriods+1), Ratio_Reloc_Prod = Array{Float64}(nPeriods+1), Horiz_Cost = Array{Float64}(nPeriods+1),Vert_Cost = Array{Float64}(nPeriods+1),nReloc = Array{Int64}(nPeriods+1));
for period = 1:nPeriods
    println(period,"/",nPeriods);
    (NS,NR,delta,L,E,SR,minStack,heightsTilde,NBar,NU,blockingCont,requestsID,SL,SE,cstVerticalCost) = simulateScenario(period,scenario,N,SB,heightsBlock,realStack,groupIOPoint,positionCont,Z,flexibilityRequests,blockID,costVerticalDrive,vZ,printSolutions,outputFolder);
    if method == "Heuristic"
        (horizontalCost,verticalCost,heightsBlock,positionCont,blockID,posCraneInitial) = heuristic(nTotalLocal,NSamples,N,blockingCont,NU,NR,NS,NBar,limitOfTime,delta,L,SL,SE,SB,heightsTilde,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,E,SR,minStack,cstVerticalCost,vZ,printSolutions,outputFolder,realStack,requestsID,nameIOPoint,IOPointsPosition,X,Y,heightsBlock,positionCont,blockID,period);
    else
        (horizontalCost,verticalCost,heightsBlock,positionCont,blockID,posCraneInitial) = IPSolver(limitOfTime,NBar,L,SL,SE,SB,Z,costEmptyDrive,posCraneInitial,costLoadedDrive,beta,heightsTilde,restricted,blockingCont,N,minStack,delta,SR,cstVerticalCost,vZ,NR,NS,NU,E,requestsID,heightsBlock,positionCont,blockID,realStack,nameIOPoint,printSolutions,outputFolder,period);
    end
    performanceMetrics[:Total_Cost][period+1] = round(horizontalCost + verticalCost,2);
    performanceMetrics[:Ratio_Reloc_Prod][period+1] = round((NBar-N)/N,4);
    performanceMetrics[:Horiz_Cost][period+1] = round(horizontalCost,2);
    performanceMetrics[:Vert_Cost][period+1] = round(verticalCost,2);
    performanceMetrics[:nReloc][period+1] = NBar-N;
end
performanceMetrics[:Total_Cost][1] = round(sum(performanceMetrics[:Total_Cost][2:nPeriods+1]),2);
performanceMetrics[:Ratio_Reloc_Prod][1] = round(mean(performanceMetrics[:Ratio_Reloc_Prod][2:nPeriods+1]),4);
performanceMetrics[:Horiz_Cost][1] = round(sum(performanceMetrics[:Horiz_Cost][2:nPeriods+1]),2);
performanceMetrics[:Vert_Cost][1] = round(sum(performanceMetrics[:Vert_Cost][2:nPeriods+1]),2);
performanceMetrics[:nReloc][1] = sum(performanceMetrics[:nReloc][2:nPeriods+1]);
writetable(joinpath(outputFolder,"0_TotalResults"), performanceMetrics, separator = ',', header = true);
