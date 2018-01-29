include("LinearProgram.jl");

(W_old,orderCont_old,orderTime_old,positiveW_old) = initialSolution(T, moveFrom, n, N, previousContToMove, SI, H, posteriorStacks, heightsInitial);

cost_old = LinearProgram(W_old, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, limitOfTime);

println("Starting Objective is: ", cost_old);

temperature = 1.0;
temperatureMin = 0.0001;
temperatureDiscount = 0.8;
# while T > T_min
#     T = T*alpha
# end

numVar = 0;
varTriplets = Dict{Int64,Array{Int64}}();
for m = 1:T
    for s in moveFrom[m]
        numVar = numVar + 1;
        varTriplets[numVar] = [m,s];
    end
end

while temperature >= temperatureMin
    numInStep = 0;
    numPerStep = 20;
    while numInStep <= numPerStep
        indRem = positiveW_old[rand(1:T)];
        indpartial = varTriplets[rand(1:numVar)];
        while indpartial == [indRem[1],indRem[2]]
            indpartial = varTriplets[rand(1:numVar)];
        end
        indAdd = [indpartial[1],indpartial[2],indRem[3]];
        if isExchange(indRem,indAdd,orderCont_old,orderTime_old,previousContToMove)
            (W_new,orderCont_new,orderTime_new,positiveW_new) = exchangeIndices(W_old,T,indRem,indAdd,moveFrom,orderCont_old,orderTime_old,positiveW_old);
            cost_new = LinearProgram(W_new, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, printSolver, limitOfTime);
            if !isnan(cost_new)
                numInStep = numInStep + 1;
                acceptanceProbability = exp((cost_old - cost_new)/temperature);
                if acceptanceProbability > rand()
                    cost_old = cost_new;
                    W_old = W_new;
                    orderCont_old = orderCont_new;
                    orderTime_old = orderTime_new;
                    positiveW_old = positiveW_new;
                    println(cost_old);
                end
            end
        end
    end
    temperature = temperatureDiscount * temperature;
end

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
