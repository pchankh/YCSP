include("LinearProgram.jl");

(Wgiven,orderCont,orderTime) = initialSolution(T, moveFrom, n, N, previousContToMove);

(currentObj,P) = LinearProgram(Wgiven, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime);

println("Starting Objective is: ", currentObj);

revP = Dict{Float64,Dict{Int64,Vector{Int64}}}();
for k in keys(P)
    if P[k] in keys(revP)
        i = length(collect(keys(revP[abs(P[k])])));
        revP[abs(P[k])][i+1] = collect(k);
    else
        revP[abs(P[k])] = Dict{Int64,Vector{Array{Int64}}}();
        revP[abs(P[k])][1] = collect(k);
    end
end
sortedValues = sort(collect(keys(revP)),rev=true);

indRem = zeros(3);
indAdd = zeros(3);
improvedSolutionFound = false;
i = 1;
while !improvedSolutionFound && i <= length(sortedValues)
    j = 1;
    while !improvedSolutionFound && j <= length(keys(revP[sortedValues[i]]))
        if P[revP[sortedValues[i]][j]] > 0
            indRem = revP[sortedValues[i]][j];
        else
            indAdd = revP[sortedValues[i]][j];
        end
        numVar = 0;
        varTriplets = Dict{Int64,Array{Int64}}();
        for m = 1:T
            if m != revP[sortedValues[i]][j][1]
                for s in moveFrom[m]
                    numVar = numVar + 1;
                    varTriplets[numVar] = [m,s,revP[sortedValues[i]][j][3]];
                end
            else
                for s in moveFrom[m]
                    if s != revP[sortedValues[i]][j][2]
                        numVar = numVar + 1;
                        varTriplets[numVar] = [m,s,revP[sortedValues[i]][j][3]];
                    end
                end
            end
        end
        k = 1;
        while !improvedSolutionFound && k <= numVar
            if P[revP[sortedValues[i]][j]] > 0
                indAdd = varTriplets[k];
            else
                indRem = varTriplets[k];
            end
            if isExchange(indRem,indAdd,orderCont,orderTime,previousContToMove)
                (Wnew,orderContnew,orderTimenew) = exchangeIndices(Wgiven,T,indRem,indAdd,moveFrom,orderCont,orderTime);
                (objNew,newP) = LinearProgram(Wnew, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime);
                if objNew < currentObj
                    improvedSolutionFound = true;
                    currentObj = objNew;
                    Wgiven = Wnew;
                    orderCont = orderContnew;
                    orderTime = orderTimenew;
                    P = newP;
                end
            end
            k = k + 1;
        end
        j = j + 1;
    end
    println(i);
    i = i + 1;
end

if improvedSolutionFound
    println("Improved Objective is: ", currentObj);
else
    print("The optimal solution has been found");
end

m = 1;
while !improvedSolutionFound && m <= T
    st = 1;
    while !improvedSolutionFound && st <= length(moveFrom[m])
        s = moveFrom[m][st];
        t = 1;
        while !improvedSolutionFound && t <= T
            if Wgiven[[m,s,t]] == 1
                indRem = [m,s,t];
            else
                indAdd = [m,s,t];
            end
            numVar = 0;
            varTriplets = Dict{Int64,Array{Int64}}();
            for mloc = 1:T
                if mloc != m
                    for sloc in moveFrom[mloc]
                        numVar = numVar + 1;
                        varTriplets[numVar] = [mloc,sloc,t];
                    end
                else
                    for sloc in moveFrom[m]
                        if sloc != s
                            numVar = numVar + 1;
                            varTriplets[numVar] = [mloc,sloc,t];
                        end
                    end
                end
            end
            k = 1;
            while !improvedSolutionFound && k <= numVar
                if Wgiven[[m,s,t]] == 1
                    indAdd = varTriplets[k];
                else
                    indRem = varTriplets[k];
                end
                if isExchange(indRem,indAdd,orderCont,orderTime,previousContToMove)
                    (Wnew,orderContnew,orderTimenew) = exchangeIndices(Wgiven,T,indRem,indAdd,moveFrom,orderCont,orderTime);
                    (objNew,newP) = LinearProgram(Wnew, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime);
                    if objNew < currentObj
                        improvedSolutionFound = true;
                        currentObj = objNew;
                        Wgiven = Wnew;
                        orderCont = orderContnew;
                        orderTime = orderTimenew;
                        P = newP;
                    end
                end
                k = k + 1;
            end
            t = t + 1;
        end
        st = st + 1;
    end
    m = m + 1;
end





# for k in keys(P)
#     println("----")
#     println(P[k]);
#     println(Wgiven[k]);
# end

# for m = 1:T
#     if previousContToMove[m] != 0
#         for t = 1:T
#             # Relation between variables x and w for retrievals
#             firstTerm = 0;
#             for u = 1:t-1
#                 for s in moveFrom[previousContToMove[m]]
#                     firstTerm = firstTerm + Wnew[[previousContToMove[m],s,u]];
#                 end
#             end
#             secondTerm = 0;
#             for s in moveFrom[m]
#                 secondTerm = secondTerm + Wnew[[m,s,t]];
#             end
#             println(m, " and ", t)
#             println(firstTerm >= secondTerm)
#         end
#     end
# end
