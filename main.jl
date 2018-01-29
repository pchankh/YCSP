using JuMP, Gurobi, AmplNLWriter, MathProgBase, DataStructures
cd(dirname(Base.source_path()));
include("auxilaryFunctions.jl");
include("BIP.jl");

testing = true;

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

printProblem(X,Y,Z,IOPointsPosition,N,n,realToClusterOrder,typeOfTruck,toBeLoaded,toRetrieve,toBeUnloaded,stackCost,rowCost,relocCost,costToGo,heightsInitial,posCraneInitial);

(moveWithCont,moveInit,moveWithoutCont,finalHeights,orderContStack,obj) = BIP(limitOfTime,T,moveFrom,SX,posteriorStacks,SY,SB,Z,costMove,costPreMove,posCraneInitial,costToGo,alpha,N,realToClusterOrder,changeOfOrder,typeOfTruck,previousContToMove,posteriorContStacks,SR,contMinHeightStack,anteriorStacks,artificialHeights);

printResult(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights);
