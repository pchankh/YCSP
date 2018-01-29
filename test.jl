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
for m = 1:N
    for t = 1:T
        if W[m,t] > 0.0000000001
            println(m, " , ",t);
            println(W[m,t]);
            count = count + 1
        end
    end
end
println(count, " should be ", N);

for h = 0:H
    for t = 1:T
        if newHeights[32,h,t] != 0
            println(t, " and ", h);
            println(newHeights[32,h,t]);
        end
    end
end
