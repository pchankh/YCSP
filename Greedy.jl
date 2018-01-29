# Discrete Space Hill Climbing Algorithm

function Greedy(nIter,tol,T,N,n,previousContToMove,moveFrom)

    ## Create the start node
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
    while remainToStack > 0
        t = t + 1;
        Wgiven[N - remainToStack + 1,t] = 1;
        remainToStack = remainToStack - 1;
    end
    for t = 1:T
        for m = 1:n
            setvalue(w[m,moveFrom[m][1],t],Wgiven[m,t]);
        end
        for m = n+1:N
            setvalue(w[m,moveFrom[m][1],t],Wgiven[m,t]);
        end
        for m = N+1:T
            setvalue(w[m,moveFrom[m][1],t],Wgiven[m,t]);
        end
    end



end

   loop do
      L = NEIGHBORS(currentNode);
      nextEval = -INF;
      nextNode = NULL;
      for all x in L
         if (EVAL(x) > nextEval)
              nextNode = x;
              nextEval = EVAL(x);
      if nextEval <= EVAL(currentNode)
         //Return current node since no better neighbors exist
         return currentNode;
      currentNode = nextNode;
