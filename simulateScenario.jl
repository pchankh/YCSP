################################################################################
############################## simulateScenario ################################
################################################################################
## This function read a scenario and an initial configuration for a given period
## It generates all the necessary elements to solve the problem for this period
## NS: indices of storage requests
## NR: indices of retrieval requests

function simulateScenario(period,scenario,N,SB,heightsBlock,realStack,groupIOPoint,positionCont,Z,flexibilityRequests,blockID,costVerticalDrive,vZ,printSolutions,outputFolder)
    dataPeriod = scenario[period];
    requestsID = zeros(Int64,0);
    NS = Set{Int64}();
    NR = Set{Int64}();
    delta = Array{Int64,2}(N,2);
    L = Dict{Int64,Array{Int64}}();
    E = Dict{Int64,Array{Int64}}();
    SR = Set{Int64}();
    minStack = Dict{Int64,Int64}();
    contInStack = Dict{Int64,Array{Int64}}();
    heightsTilde = Array{Int64}(length(SB),1);
    for s in SB
        heightsTilde[s] = heightsBlock[s];
    end
    for n = 1:N
        push!(requestsID, dataPeriod[n,4]);
        if dataPeriod[n,1] == "storage"
            push!(NS,n);
            L[n] = groupIOPoint[dataPeriod[n,3]];
            E[n] = SB;
        else
            push!(NR,n);
            s = positionCont[dataPeriod[n,4]][1];
            z = positionCont[dataPeriod[n,4]][2];
            L[n] = [s];
            if s in SR
                if heightsTilde[s] >= z
                    minStack[s] = n;
                    heightsTilde[s] = z-1;
                end
                contInStack[s][z] = n;
            else
                push!(SR,s);
                minStack[s] = n;
                contInStack[s] = fill(-1,heightsTilde[s]);
                contInStack[s][z] = n;
                heightsTilde[s] = z-1;
            end
            E[n] = groupIOPoint[dataPeriod[n,3]];
        end
        delta[n,:] = flexibilityRequests[dataPeriod[n,2]];
    end
    for s in SR
        if length(contInStack[s]) > 1
            for i = length(contInStack[s]):-1:2
                if contInStack[s][i] != -1
                    minCont = contInStack[s][i]
                    for j = 1:i-1
                        if contInStack[s][j] != -1 && contInStack[s][j] < minCont
                            minCont = contInStack[s][j];
                        end
                    end
                    if minCont < contInStack[s][i]
                        for k = contInStack[s][i]:-1:minCont+1
                            if k-1 in NS
                                delete!(NS,k-1);
                                push!(NS,k);
                                delete!(NR,k);
                                push!(NR,k-1);
                            else
                                s = L[k-1][1];
                                ind = findfirst(contInStack[s], k-1);
                                contInStack[s][ind] = contInStack[s][ind] + 1;
                                if k-1 == minStack[s]
                                    minStack[s] = k;
                                end
                            end
                            dataPeriod[k-1,:], dataPeriod[k,:] = dataPeriod[k,:], dataPeriod[k-1,:];
                            requestsID[k-1], requestsID[k] = requestsID[k], requestsID[k-1];
                            delta[k-1,:], delta[k,:] = delta[k,:], delta[k-1,:];
                            L[k-1], L[k] = L[k], L[k-1];
                            E[k-1], E[k] = E[k], E[k-1];
                        end
                        contInStack[s][i] = minCont;
                    end
                end
            end
        end
    end
    NBar = N;
    NU = Set{Int64}();
    blockingCont = Dict{Int64,Int64}();
    for s in SR
        h = 1;
        while contInStack[s][h] != minStack[s]
            h = h + 1;
        end
        h = h + 1;
        while h <= length(contInStack[s])
            if contInStack[s][h] == -1
                NBar = NBar + 1;
                push!(requestsID, blockID[[s,h]]);
                contInStack[s][h] = NBar;
                push!(NU,NBar);
                L[NBar] = [s];
                E[NBar] = setdiff(SB,s);
                dataPeriod = [dataPeriod; ["relocation" "" ""  blockID[[s,h]]]];
            end
            blockingCont[contInStack[s][h-1]] = contInStack[s][h];
            h = h + 1;
        end
    end
    SL = Set{Int64}();
    SE = Set{Int64}();
    for n = 1:NBar
        union!(SL,Set(L[n]));
        union!(SE,Set(E[n]));
    end
    cstVerticalCost = length(union(NR,NS))*costVerticalDrive[1] + length(union(NS,NU))*costVerticalDrive[0];
    for n = union(NR,NU)
        cstVerticalCost = cstVerticalCost + costVerticalDrive[positionCont[dataPeriod[n,4]][2]];
    end
    for r in intersect(SE,SB)
        cstVerticalCost = cstVerticalCost + heightsTilde[r]*(heightsTilde[r]+1)/vZ;
    end

    if printSolutions
        outputFile = joinpath(outputFolder,string(period,"_Period"));
        f = open(outputFile,"w");
        write(f,string("-------------------------------", "\n"));
        write(f,string("PRODUCTIVE MOVES", "\n"));
        write(f,string("Number of storage requests: ", length(NS), "\n"));
        write(f,string("Number of retrieval requests: ", length(NR),"\n"));
        write(f,string("Number of relocation requests: ", length(NU),"\n"));
        write(f,string("-------------------------------", "\n"));
        write(f,string("Order\tCont_ID\tRequestType\tTruckType\tIOPoint\tRow\tStack\tTier\n"));
        order = 0;
        for n = 1:N
            if n in NS
                order = order + 1;
                write(f,string(order,"\t",requestsID[n],"\t",dataPeriod[n,1],"    \t",dataPeriod[n,2],"\t",dataPeriod[n,3],"\n"));
            else
                toPrint = [n];
                while toPrint[1] in keys(blockingCont) && !(blockingCont[toPrint[1]] in NR)
                    toPrint = [blockingCont[toPrint[1]];toPrint];
                end
                for i = 1:length(toPrint)
                    c = toPrint[i];
                    order = order + 1;
                    if c in NU
                        write(f,string(order,"\t",requestsID[c],"\t",dataPeriod[c,1],"\t\t\t\t",realStack[L[c][1]][1],"\t",realStack[L[c][1]][2],"\t",positionCont[dataPeriod[c,4]][2],"\n"));
                    else
                        write(f,string(order,"\t",requestsID[c],"\t",dataPeriod[c,1],"\t",dataPeriod[c,2],"\t",dataPeriod[c,3],"\t",realStack[L[c][1]][1],"\t",realStack[L[c][1]][2],"\t",positionCont[dataPeriod[c,4]][2],"\n"));
                    end
                end
            end
        end
        write(f,string("-------------------------------\n"));
        close(f);
    end
    return (NS,NR,delta,L,E,SR,minStack,heightsTilde,NBar,NU,blockingCont,requestsID,SL,SE,cstVerticalCost)
end
