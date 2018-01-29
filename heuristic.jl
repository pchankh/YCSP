include("BIP.jl");

(moveWithCont,moveInit,moveWithoutCont,finalHeights,orderContStack,obj) = BIP(limitOfTime,T,moveFrom,SX,posteriorStacks,SY,SB,Z,costMove,costPreMove,posCraneInitial,costToGo,alpha,N,realToClusterOrder,changeOfOrder,typeOfTruck,previousContToMove,posteriorContStacks,SR,contMinHeightStack,anteriorStacks,artificialHeights);

printResultTest(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights);
