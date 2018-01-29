using JuMP, Gurobi, MathProgBase
cd(dirname(Base.source_path()));
include("auxilaryFunctions.jl");
include("MIP.jl");

rowLook = 0;
limitOfTimeHeur = 30;

parametersData = readcsv("Parameters.csv", header=false);

N = parametersData[2,2];
n = N - parametersData[3,2];

limitOfTime = parametersData[5,2];
gapOfMIP = parametersData[6,2];
printSolver = parametersData[7,2];

randomInitialCrane = parametersData[9,2];
randomInitialHeights = parametersData[10,2];
randomNumberOfActions = parametersData[11,2];
randomRetrieval = parametersData[12,2];

relocCost = parametersData[14,2];
rowCost = parametersData[15,2];
stackCost = parametersData[16,2];


H = parametersData[18,2];

if randomInitialHeights
    R = parametersData[19,2];
    S = parametersData[20,2];
    heightsInitial = randomInitialHeightsFun(R,S,H);
else
    heightsInitial = readcsv("block.csv", header=false, Int64);
    R = size(heightsInitial,1);
    S = size(heightsInitial,2);
end

IOPointsPosition = parametersData[21,2];
(SB,SI,Stot) = initializeMainValues(heightsInitial,IOPointsPosition);

if randomInitialCrane
    posCraneInitial = randomInitialCraneFun(Stot);
else
    craneInitial = readcsv("cranePosition.csv", header=false, Int64);
    if (craneInitial[1] == 0 || (craneInitial[1] == S+1 && IOPointsPosition=="right-sided"))
        posCraneInitial = SI[craneInitial[2]];
    elseif craneInitial[1] == S+1
        posCraneInitial = SI[R+craneInitial[2]];
    else
        posCraneInitial = S*(craneInitial[1]-1)+craneInitial[2];
    end
end
if randomNumberOfActions
    (N,n) = randomNumberOfActionsFun(H,heightsInitial);
end
if randomRetrieval
    toRetrieve = randomRetrievalFun(n,heightsInitial);
else
    toRetrieve = readcsv("retrievals.csv", header=false, Int64);
end

realStack = seqStackToRealStack(heightsInitial,R,S,SB,SI,IOPointsPosition);

(SR,SO,SP,indicesToRetrieve,contInStackToRetrieve) = dataToRetrieve(toRetrieve,S,n,SB,SI);

alpha = computeAlphas(H);

(costMove, costPreMove) = defineCost(R,S,Stot,SB,SP,SR,realStack,relocCost,rowCost,stackCost);

(IOPoint,correspondingIOPoint,rowIOPoint) = dataIOPoints(S,SI,SB,IOPointsPosition,realStack);

(anteriorStacks,posteriorStacks) = antePostStacks(SR,SB,SI,SO,IOPoint,correspondingIOPoint);

T = timeUpperBound(SR,H,N,n,heightsInitial,indicesToRetrieve,realStack,contInStackToRetrieve);
#T = timeUpperBound_2(N,n,H,indicesToRetrieve);

# costToGo = 1/(R*S);
costToGo = sqrt(T)/(R*S);
# costToGo = sqrt(T/(R*S));
costToGo = round(costToGo,4);

printProblem(relocCost, rowCost, stackCost, costToGo, H, R, S, heightsInitial,  posCraneInitial, N, n, toRetrieve, IOPointsPosition);

(XHeur,DHeur,newHeightsHeur,finalHeightsHeur,WHeur,timeToSolveHeur) = rowLookIPHeuristic(H, T, n, N, SI, SP, SR, SO, Stot, contInStackToRetrieve, indicesToRetrieve, heightsInitial, posCraneInitial, costMove, costPreMove, costToGo, alpha, posteriorStacks, anteriorStacks, IOPoint, printSolver, gapOfMIP, limitOfTimeHeur, rowLook);

(X,D,newHeights,finalHeights,W, timeToSolve) = MIP(H, T, n, N, SI, SP, SR, SO, Stot, contInStackToRetrieve, indicesToRetrieve, heightsInitial, posCraneInitial, costMove, costPreMove, costToGo, alpha, posteriorStacks, anteriorStacks, IOPoint, printSolver, gapOfMIP, limitOfTime, XHeur,DHeur,newHeightsHeur,finalHeightsHeur,WHeur);

printResult_2(S, R, SR, SB, SI, SP, Stot, T, posteriorStacks, X, D, W, heightsInitial, IOPointsPosition, realStack, rowIOPoint);
