count = 0;
for s in SX
    for r in posteriorStacks[s]
        count = count + T;
    end
end

count = 0;
for t = 1:T
    for s in SX
        for r in posteriorStacks[s]
            if X[s,r,t] > 0.00001
                count = count + 1;
                println(s, " , ", r, " , ", t);
                println(X[s,r,t]);
            end
        end
    end
end
println(count, " should be ", T);

count = 0;
for m =1:T
    for t = 1:T
        if W[m,t] > 0.9
            println(W[m,t])
            println(WLP[m,t])
            println(P[m,t])
            println("--------")
        end
    end
end

for m = 1:T
    count = 0;
    for t = 1:T
        if WLP[m,t] > 0.0000000001
            # println(m, " , ",t);
            # println(WLP[m,t]);
            count = count + 1
        end
    end
    # println(m," : ", count);
    countB = 0;
    for t = 1:T
        if WLP[m,t] > 2/count
            # println(m, " , ",t);
            # println(WLP[m,t]);
            countB = countB + 1
        end
    end
    println(m," : ", countB);
end


count = 0;
for s in SB
    for h = 0:H
        if finalHeights[s,h] > 0.0000000001
            println(s, " , ",h);
            println(finalHeights[s,h]);
            count = count + 1
        end
    end
end
println(count, " should be ", R*S);

Wgiven = zeros(Int64,T,T);
t = 0;
remainToStack = N - n;
for m = 1:n
    if sum(Wgiven[m,:]) == 0
        blockingCont = [m];
        while previousContToMove[blockingCont[length(blockingCont)]] != 0
            append!(blockingCont,previousContToMove[blockingCont[length(blockingCont)]]);
        end
        for c = length(blockingCont):-1:1
            t = t + 1;
            Wgiven[blockingCont[c],t] = 1;
        end
        if remainToStack > 0
            t = t + 1;
            Wgiven[N - remainToStack + 1,t] = 1;
            remainToStack = remainToStack - 1;
        end
    end
end
while remainToStack > 03
    t = t + 1;
    Wgiven[N - remainToStack + 1,t] = 1;
    remainToStack = remainToStack - 1;
end
