
function MIP(costToGo, heightsInitial, H, N, n, R, S, Stot, posCraneInitial, realStack, indicesToRetrieve, stacksToRetrieve, contInStackToRetrieve, alpha, costRetrieveOrStack, costRelocate, costMove, correspondingIOPoint, T, printSolver, gapOfMIP, limitOfTime)

    startPoint = union(stacksToRetrieve,collect(Stot+1:Stot+2*R));

    IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));
    @variable(IPModel, x[stacksToRetrieve,1:Stot,1:T],Bin);
    @variable(IPModel, y[stacksToRetrieve,1:T]);
    @variable(IPModel, z[1:Stot,1:T],Bin);
    @variable(IPModel, a[1:n,0:T],Bin);
    @variable(IPModel, d[1:Stot+2*R,startPoint,1:T], Bin);
    @variable(IPModel, 0 <= heights[1:Stot,0:T] <= H);
    @variable(IPModel, active[0:T] <= 1);
    @variable(IPModel, finalh[1:Stot,0:H], Bin);

    @objective(IPModel, Min, sum{sum{costRetrieveOrStack[realStack[s][2]]*z[s,t],s=1:Stot} +sum{costRelocate[[s,r]]*x[s,r,t],s in stacksToRetrieve, r=1:Stot},t=1:T} + sum{costMove[[s,r]]*d[s,r,t],s=1:Stot+2*R,r in startPoint,t=1:T}+ costToGo * sum{sum{alpha[h]*finalh[s,h],h=1:H},s=1:Stot});

    for t=1:T
        @constraint(IPModel, sum{y[s,t], s in stacksToRetrieve} + sum{z[r,t], r = 1:Stot} + sum{x[s,r,t], s in stacksToRetrieve, r=1:Stot} - active[t] == 0);
        @constraint(IPModel, active[t] - active[t-1] <= 0);
    end

    for t = 1:T
        for s = 1:Stot
            if s in stacksToRetrieve
                @constraint(IPModel, x[s,s,t] == 0);
                @constraint(IPModel, heights[s,t] - heights[s,t-1] - z[s,t] + y[s,t] - sum{x[r,s,t], r in stacksToRetrieve} + sum{x[s,r,t],r=1:Stot} == 0);
            else
                @constraint(IPModel, heights[s,t] - heights[s,t-1] - z[s,t] - sum{x[r,s,t], r in stacksToRetrieve} == 0);
            end
        end
    end

    for t = 1:T
        for m = 1:n
            @constraint(IPModel, heights[indicesToRetrieve[m][1],t-1] - indicesToRetrieve[m][2] - H * (1 + a[m,t] - a[m,t-1]) <= 0);
            @constraint(IPModel, a[m,t] - a[m,t-1] <= 0);
            @constraint(IPModel, indicesToRetrieve[m][2] - heights[indicesToRetrieve[m][1],t] + H * a[m,t] <= H);
        end
        for s in stacksToRetrieve
            @constraint(IPModel, y[s,t] - sum{a[m,t-1] - a[m,t], m = contInStackToRetrieve[s]} == 0);
        end
    end

    for s = 1:Stot+2*R
        for r in startPoint
            if s != posCraneInitial
                @constraint(IPModel, d[s,r,1] == 0);
            else
                if r in stacksToRetrieve
                    @constraint(IPModel, sum{x[r,u,1],u=1:Stot} + y[r,1] - d[s,r,1] <= 0);
                else
                    @constraint(IPModel, sum{z[u,1],u=correspondingIOPoint[r]} - d[s,r,1] <= 0);
                end
            end
        end
    end
    @constraint(IPModel, sum{d[s,r,1],s = 1:Stot+2*R, r in startPoint} == 1);
    for t = 2:T
        @constraint(IPModel, sum{d[s,r,t],s = 1:Stot + 2*R, r in startPoint} - active[t] == 0);
        for s = 1:Stot
            for r in stacksToRetrieve
                if s in stacksToRetrieve
                    @constraint(IPModel, sum{x[s,u,t-1] + x[r,u,t],u=1:Stot} + z[s,t-1] + y[r,t] - d[s,r,t] <= 1);
                else
                    @constraint(IPModel, z[s,t-1] + sum{x[r,u,t],u=1:Stot} + y[r,t] - d[s,r,t] <= 1);
                end
            end
        end
        for s = 1:Stot
            for r = Stot+1:Stot+2*R
                if s in stacksToRetrieve
                    @constraint(IPModel, sum{x[s,u,t-1],u=1:Stot} + z[s,t-1] + sum{z[u,t],u=correspondingIOPoint[r]} - d[s,r,t] <= 1);
                else
                    @constraint(IPModel, z[s,t-1] + sum{z[u,t],u=correspondingIOPoint[r]} - d[s,r,t] <= 1);
                end
            end
        end
        for s = Stot+1:Stot+2*R
            for r in stacksToRetrieve
                @constraint(IPModel, sum{y[u,t-1],u=intersect(stacksToRetrieve,correspondingIOPoint[s])} + sum{x[r,u,t],u=1:Stot} + y[r,t] - d[s,r,t] <= 1);
            end
        end
        for s = Stot+1:Stot+2*R
            for r = Stot+1:Stot+2*R
                @constraint(IPModel, sum{y[u,t-1],u=intersect(stacksToRetrieve,correspondingIOPoint[s])} + sum{z[u,t],u=correspondingIOPoint[r]} - d[s,r,t] <= 1);
            end
        end
    end

    @constraint(IPModel, active[0] == 1);
    @constraint(IPModel, sum{z[s,t],s = 1:Stot, t=1:T} == N - n);
    @constraint(IPModel, sum{y[s,t],s in stacksToRetrieve, t=1:T} == n);
    for m = 1:n
        @constraint(IPModel, a[m,T] == 0);
        @constraint(IPModel, a[m,0] == 1);
    end
    for s = 1:Stot
        @constraint(IPModel, heights[s,0] == heightsInitial[realStack[s][1],realStack[s][2]]);
    end

    for s=1:Stot
        @constraint(IPModel, sum{h*finalh[s,h],h=1:H} - heights[s,T] == 0);
        @constraint(IPModel, sum{finalh[s,h],h=0:H} == 1);
    end

## Precedence constraint
    # for m=1:n-1
    #     for p = m+1:n
    #         if indicesToRetrieve[m][1] != indicesToRetrieve[p][1] || (indicesToRetrieve[m][1] == indicesToRetrieve[p][1] &&
    #         indicesToRetrieve[m][2] > indicesToRetrieve[p][2])
    #             @constraint(IPModel, sum{t*((a[m,t-1] - a[m,t]) - (a[p,t-1] - a[p,t])), t = 1:T} <= 0);
    #         end
    #     end
    # end

    println(MathProgBase.numvar(IPModel) , " variables");
    println(MathProgBase.numconstr(IPModel), " constraints");
    println("Solving ....");

    tic();
    status = solve(IPModel);
    println("Solved !");
    timeToSolve = toc();

    X = getvalue(x);
    Y = getvalue(y);
    Z = getvalue(z);
    A = getvalue(a);
    D = getvalue(d);
    ACT = getvalue(active);
    newHeights = getvalue(heights);

    return (X,Y,Z,A,D,ACT,newHeights,timeToSolve);
end
