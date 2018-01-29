function MIP_2(H, T, n, N, SI, SP, SR, SO, Stot, contInStackToRetrieve, indicesToRetrieve, heightsInitial, posCraneInitial, costMove, costPreMove, costToGo, alpha, posteriorStacks, anteriorStacks, IOPoint, printSolver, gapOfMIP, limitOfTime)

   IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));
   x = @variable(IPModel, [s in SP,r in posteriorStacks[s],1:T], Bin);
   heights = @variable(IPModel, [SR,0:H,1:T], Bin);
   d = @variable(IPModel, [Stot,SP,1:T], Bin);
   finalh = @variable(IPModel, [SO,0:H], Bin);
   w = @variable(IPModel, [1:n,1:T], Bin);
   #active = @variable(IPModel, [1:T], Bin);


   @objective(IPModel, Min, sum( sum(costMove[[s,r]]*x[s,r,t] for s in SP for r in posteriorStacks[s]) + sum(costPreMove[[s,r]]*d[s,r,t] for s in Stot for r in SP) for t = 1:T) + costToGo * sum( alpha[h] * (sum(heights[s,h,T] for s in SR) + sum(finalh[s,h] for s in SO))  for h = 2:H));

   for t = 1:T
       # @constraint(IPModel, sum(d[s,r,t] for s in Stot for r in SP) - active[t] == 0);
       # @constraint(IPModel, sum(x[s,r,t] for s in SP for r in posteriorStacks[s]) - active[t] == 0);
       @constraint(IPModel, sum(d[s,r,t] for s in Stot for r in SP) == 1);
       @constraint(IPModel, sum(x[s,r,t] for s in SP for r in posteriorStacks[s]) == 1);
   end

   # for t = 2:T
   #     @constraint(IPModel, active[t] - active[t-1] <= 0);
   # end

   @constraint(IPModel, sum(x[s,r,t] for s in SI for r in posteriorStacks[s] for t = 1:T) == N - n);

   for m = 1:n
       @constraint(IPModel, sum(w[m,t] for t = 1:T) == 1);
   end

   for s in SR
       for t = 1:T
           @constraint(IPModel, x[s,IOPoint[s],t] - sum(w[m,t] for m in contInStackToRetrieve[s]) == 0);
       end
   end

   for m = 1:n
       if heightsInitial[realStack[indicesToRetrieve[m][1]][1],realStack[indicesToRetrieve[m][1]][2]] > indicesToRetrieve[m][2]
           @constraint(IPModel, w[m,1] == 0);
       end
       for t = 2:T
           @constraint(IPModel, w[m,t] - heights[indicesToRetrieve[m][1],indicesToRetrieve[m][2],t-1] <= 0);
       end
       for t = 1:T
           @constraint(IPModel, sum(heights[indicesToRetrieve[m][1],h,t] for h = 0:indicesToRetrieve[m][2]-1) - sum(w[m,u] for u = 1:t) <= 0);
       end
   end

   for s in SR
       @constraint(IPModel, sum(h * heights[s,h,1] for h = 0:H) - sum(x[r,s,1] for r in anteriorStacks[s]) + sum(x[s,r,1] for r in posteriorStacks[s]) == heightsInitial[realStack[s][1],realStack[s][2]]);
       for t = 2:T
           @constraint(IPModel, sum(h * heights[s,h,t] for h = 0:H) - sum(h * heights[s,h,t-1] for h = 0:H) - sum(x[r,s,t] for r in anteriorStacks[s]) + sum(x[s,r,t] for r in posteriorStacks[s]) == 0);
       end
       for t = 1:T
           @constraint(IPModel, sum(heights[s,h,t] for h = 0:H) == 1);
       end
   end

   for s in SO
       @constraint(IPModel, sum(h * finalh[s,h] for h = 0:H) - sum(x[r,s,t] for r in anteriorStacks[s] for t = 1:T) == heightsInitial[realStack[s][1],realStack[s][2]]);
       @constraint(IPModel, sum(finalh[s,h] for h = 0:H) == 1);
   end

   for s in setdiff(Stot,posCraneInitial)
       for r in SP
           @constraint(IPModel, d[s,r,1] == 0);
       end
   end

   for t = 1:T
       for s in SP
           @constraint(IPModel, sum(x[s,r,t] for r in posteriorStacks[s]) - sum(d[r,s,t] for r in Stot) == 0);
       end
   end

   for t = 2:T
       for s in Stot
           @constraint(IPModel, sum(d[s,r,t] for r in SP) - sum(x[r,s,t-1] for r in anteriorStacks[s]) == 0);
       end
   end


   for s in SI
       if length(anteriorStacks[s]) == 0
           for t = 2:T
               for r in SP
                   @constraint(IPModel, d[s,r,t] == 0);
               end
           end
       end
   end

   # for s in SR
   #     for t = 1:T
   #         for m in contInStackToRetrieve[s]
   #             @constraint(IPModel, sum(x[r,s,t] for r in anteriorStacks[s]) - w[m,t] <= 0);
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
   newHeights = getvalue(heights);
   D = getvalue(d);
   finalHeights = getvalue(finalh);
   #ACT = getvalue(active);
   W = getvalue(w);

   return (X,newHeights,D,finalHeights,W, timeToSolve);
end




function MIP_1(costToGo, heightsInitial, H, N, n, R, S, SR, SB, SP, SO, Stot, posCraneInitial, realStack, indicesToRetrieve, contInStackToRetrieve, alpha, costMove, costPreMove, IOPoint, correspondingIOPoint, T, printSolver, gapOfMIP, limitOfTime);

   IPModel = Model(solver = GurobiSolver(OutputFlag = printSolver,MIPGap = gapOfMIP, TimeLimit = limitOfTime));
   x = @variable(IPModel, [SR,SB,1:T], Bin);
   y = @variable(IPModel, [SR,1:T]);
   z = @variable(IPModel, [SB,1:T], Bin);
   a = @variable(IPModel, [1:n,0:T], Bin);
   d = @variable(IPModel, [Stot,SP,1:T], Bin);
   heights = @variable(IPModel, [SB,0:T], lowerbound=0, upperbound=H);
   finalh = @variable(IPModel, [SB,0:H], Bin);
   active = @variable(IPModel, [1:T], upperbound=1);

   @objective(IPModel, Min, sum(costPreMove[[s,r]]*d[s,r,t] for s in Stot for r in SP for t = 1:T) + sum(sum(costMove[[s,r]]*x[s,r,t] for s in SR for r in SB) + sum(costMove[[IOPoint[s],s]] * z[s,t] for s in SB) + sum(costMove[[s,IOPoint[s]]]*y[s,t] for s in SR) for t = 1:T) + costToGo * sum(alpha[h]*finalh[s,h] for h = 2:H for s in SB));

   # Actions
   # N-n stackings
   @constraint(IPModel, sum(z[s,t] for s in SB for t = 1:T) == N - n);
   # n retrievals
   @constraint(IPModel, sum(y[s,t] for s in SR for t = 1:T) == n);

   # If active then 1 action is possible
   for t = 1:T
       @constraint(IPModel, sum(y[s,t] for s in SR) + sum(z[r,t] for r in SB) + sum(x[s,r,t] for s in SR for r in SB) - active[t] == 0);
   end
   # if unactive at time t-1, then remains inactive at time t
   for t = 2:T
       @constraint(IPModel, active[t] - active[t-1] <= 0);
   end
   # No relocation from s to s is possible at any time
   for t = 1:T
       for s in SR
           @constraint(IPModel, x[s,s,t] == 0);
       end
   end

   # Enforce retrievals
   for t = 1:T
       for m = 1:n
           # If container m was already retrieved at time t-1, it is sill true at time t
           @constraint(IPModel, a[m,t] - a[m,t-1] <= 0);
           # If heights[indicesToRetrieve[m][1],t-1] > indicesToRetrieve[m][2] then a[m,t] = a[m,t-1]
           @constraint(IPModel, heights[indicesToRetrieve[m][1],t-1] - indicesToRetrieve[m][2] - H * (1 + a[m,t] - a[m,t-1]) <= 0);
           # If heights[indicesToRetrieve[m][1],t-1] < indicesToRetrieve[m][2] then a[m,t] = 0
           @constraint(IPModel, indicesToRetrieve[m][2] - heights[indicesToRetrieve[m][1],t] + H * (a[m,t] - 1) <= 0);
       end
       # relation between y and a
       for s in SR
           @constraint(IPModel, y[s,t] - sum(a[m,t - 1] - a[m,t] for m = contInStackToRetrieve[s]) == 0);
       end
   end

   # Update heights
   for t = 1:T
       for s in SR
           @constraint(IPModel, heights[s,t] - heights[s,t-1] - z[s,t] + y[s,t] - sum(x[r,s,t] for r in SR) + sum(x[s,r,t] for r in SB) == 0);
       end
       for s in SO
           @constraint(IPModel, heights[s,t] - heights[s,t-1] - z[s,t] - sum(x[r,s,t] for r in SR) == 0);
       end
   end

   # Final heights
   for s in SB
       @constraint(IPModel, sum(h * finalh[s,h] for h = 1:H) - heights[s,T] == 0);
       @constraint(IPModel, sum(finalh[s,h] for h = 0:H) == 1);
   end

   # Initialization
   for s in SB
       @constraint(IPModel, heights[s,0] == heightsInitial[realStack[s][1],realStack[s][2]]);
   end
   for m = 1:n
       @constraint(IPModel, a[m,0] == 1);
       @constraint(IPModel, a[m,T] == 0);
   end

   # Pre-moves and Actions
   # 1st time step
   # crane starts at posCraneInitial
   for s in setdiff(Stot,posCraneInitial)
       for r in SP
           @constraint(IPModel, d[s,r,1] == 0);
       end
   end
   # if destination is r in SR then reloc or retrieve from r
   for r in SR
       @constraint(IPModel, sum(x[r,u,1] for u in SB) + y[r,1] - d[posCraneInitial,r,1] <= 0);
   end
   # if destination is r in SI then stack from r
   for r in SI
       @constraint(IPModel, sum(z[u,1] for u in posteriorStacks[r]) - d[posCraneInitial,r,1] <= 0);
   end
   # Uniqness of move
   @constraint(IPModel, sum(d[posCraneInitial,r,1] for r in SP) - active[1] == 0);

   # Other time steps
   for t = 2:T
       # Uniqness of move
       @constraint(IPModel, sum(d[s,r,t] for s in Stot for r in SP) - active[t] == 0);
       # if destination is r in SR then reloc or retrieve from r
       for r in SR
           # if origin is s in SR then reloc or stack to s
           for s in SB
               @constraint(IPModel, sum(x[u,s,t-1] for u in SR) + z[s,t-1] + sum(x[r,u,t] for u in SB) + y[r,t] - d[s,r,t] <= 1);
           end
           # if origin is s in SI s.t. anteriorStacks[s] is not empty then retrieve to s
           for s in SI
               if length(anteriorStacks[s]) != 0
                   @constraint(IPModel, sum(y[u,t - 1] for u = anteriorStacks[s]) + sum(x[r,u,t] for u in SB) + y[r,t] - d[s,r,t] <= 1);
               end
           end
       end
       # if destination is r in SI then stack from r
       for r in SI
           # if origin is s in SB then reloc or stack to s
           for s in SB
               @constraint(IPModel, sum(x[u,s,t-1] for u in SR) + z[s,t-1] + sum(z[u,t] for u = posteriorStacks[r]) - d[s,r,t] <= 1);
           end
           # if origin is s in SI s.t. anteriorStacks[s] is not empty then retrieve to s
           for s in SI
               if length(anteriorStacks[s]) != 0
                   @constraint(IPModel, sum(y[u,t - 1] for u = intersect(SR,correspondingIOPoint[s])) + sum(z[u,t] for u = correspondingIOPoint[r]) - d[s,r,t] <= 1);
               end
           end
       end

       for s in SI
           if length(anteriorStacks[s]) == 0
               for r in SP
                   @constraint(IPModel, d[s,r,t] == 0);
               end
           end
       end
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

# (X,Y,Z,A,D,ACT,newHeights,timeToSolve) = MIP_1(costToGo, heightsInitial, H, N, n, R, S, SR, SB, SP, SO, Stot, posCraneInitial, realStack, indicesToRetrieve, contInStackToRetrieve, alpha, costMove, costPreMove, IOPoint, correspondingIOPoint, T, printSolver, gapOfMIP, limitOfTime);
#
# printResult(R,S,T,n,SP,SR,SB,Stot,X,Y,Z,A,D,ACT,newHeights,realStack,rowIOPoint,indicesToRetrieve,IOPointsPosition);

# (X,newHeights,D,finalHeights,W, timeToSolve) = MIP_2(H, T, n, N, SI, SP, SR, SO, Stot, contInStackToRetrieve, indicesToRetrieve, heightsInitial, posCraneInitial, costMove, costPreMove, costToGo, alpha, posteriorStacks, anteriorStacks, IOPoint, printSolver, gapOfMIP, limitOfTime);

# printResult_2(S, R, SR, SB, SI, SP, Stot, T, posteriorStacks, X, D, W, heightsInitial, IOPointsPosition, realStack, rowIOPoint);
