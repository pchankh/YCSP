include("BIP.jl");

using Gurobi

(moveWithCont,moveInit,moveWithoutCont,finalHeights,orderContStack,obj) = BIP(limitOfTime,T,moveFrom,SX,posteriorStacks,SY,SB,Z,costMove,costPreMove,posCraneInitial,costToGo,alpha,N,realToClusterOrder,changeOfOrder,typeOfTruck,previousContToMove,posteriorContStacks,SR,contMinHeightStack,anteriorStacks,artificialHeights);

printResultTest(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights);

function printResultTest(IOPointsPosition,X,Y,heightsInitial,realStack,T,moveInit,posCraneInitial,nameIOPoint,SB,SX,SY,moveWithoutCont,posteriorStacks,moveWithCont,N,n,SL,SU,orderContStack,stackOf,unloadFrom,clusterToRealOrder,finalHeights)
    f = open("SolutionTest.txt","w")
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
    writecsv("newBlockTest.csv",lastHeight);
    writecsv("newCranePositionTest.csv",lastCranePosition);
end
