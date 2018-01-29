using JuMP, Gurobi, AmplNLWriter, MathProgBase, DataStructures, Permutations
cd(dirname(Base.source_path()));
include("auxilaryFunctions.jl");
include("MIP.jl");
include("BIP.jl");

parametersData = readcsv("Parameters.csv", header=false);

R = parametersData[2,2];
S = parametersData[3,2];
H = parametersData[4,2];
IOPointsPosition = parametersData[5,2];

N = parametersData[7,2];
n = parametersData[8,2];

stackCost = parametersData[10,2];
rowCost = parametersData[11,2];
relocCost = parametersData[12,2];

limitOfTime = parametersData[14,2];
gapOfMIP = parametersData[15,2];
printSolver = parametersData[16,2];

randomInitialCrane = parametersData[18,2];
randomInitialHeights = parametersData[19,2];
randomRetrieval = parametersData[20,2];
randomStacking = parametersData[21,2];

(heightsInitial, R, S, H) = initializeInitialHeights(randomInitialHeights, R, S, H);

(SB,SI,Stot,realStack,IOPoints,innerPoints,nameIOPoint,groupIOPoint) = initializeMainValues(R,S,IOPointsPosition);

posCraneInitial = initializeInitialCrane(randomInitialCrane, R, S, IOPointsPosition, SI, Stot);

(T,SR,SO,SL,SY,stackOf,heightOf,loadOf,contStackIOPoint,contMinHeightStack,previousContToMove,toRetrieve,toBeLoaded) = initializeRetrievals(n,S,SB,heightsInitial,IOPoints, groupIOPoint,IOPointsPosition);

(unloadFrom, SU, SX, contStackableIOPoint, toBeUnloaded) = initializeStackings(N,n,randomStacking,IOPointsPosition,SR,groupIOPoint);

(anteriorStacks,posteriorStacks) = antePostStacks(SR,SU,SO,SL,IOPoints,innerPoints,contStackIOPoint);

(costMove, costPreMove, costToGo, alpha) = defineCosts(N, R, S, H, SX, SY, posCraneInitial, posteriorStacks, rowCost, stackCost, relocCost, realStack);

printProblem(R,S,H,IOPointsPosition,N,n,toRetrieve,toBeLoaded,toBeUnloaded,stackCost,rowCost,relocCost,costToGo,heightsInitial,posCraneInitial);

(moveFrom,posteriorContStacks,artificialHeights) = defineMoveFromStack(n,N,T,unloadFrom,loadOf,stackOf,heightOf,SB,SR,SO,contMinHeightStack,heightsInitial,realStack);

(X,DInit,D,finalHeights,W,obj,timeToSolve) = MIP(H, N, n, artificialHeights, moveFrom, toBeUnloaded, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime);

(X_b,DInit_b,D_b,finalHeights_b,W_b,obj_b,timeToSolve_b) = BIP(H, N, n, heightsInitial, moveFrom, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime);

# printResult(S,R,N,n,posCraneInitial,IOPointsPosition,heightsInitial,realStack,stackOf,unloadFrom,SB,SU,SL,SX,SY,posteriorStacks,nameIOPoint,T,X,DInit,D,W);
