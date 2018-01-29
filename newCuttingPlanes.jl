
include("subProblem.jl");

function Solver(H, N, n, heightsInitial, stackOf, heightOf, loadOf, toBeUnloaded, realStack, IOPoints, SR, SB, SO, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove,  costMove, costPreMove, costToGo, alpha, printSolver, gapOfMIP, limitOfTime)

    CutModel = Model(solver = GurobiSolver(PreCrush=1, OutputFlag = printSolver, TimeLimit = limitOfTime));

    #####################################################################
    ############################# Variables #############################
    #####################################################################

    @variable(CutModel, w[m = 1:T, s in moveFrom[m], t = 1:T], Bin);
    @variable(CutModel, z >= 0);

    #####################################################################
    ############################# Objective #############################
    #####################################################################

    @objective(CutModel, Min, z);

    #####################################################################
    ################### Conditions on productive moves ##################
    #####################################################################

    for m = 1:T
        ## Uniqueness of move
        @constraint(CutModel, sum(w[m,s,t] for t = 1:T for s in moveFrom[m]) == 1);
    end

    for t = 1:T
        ## Uniqueness of move (2)
        @constraint(CutModel, sum(w[m,s,t] for m = 1:T for s in moveFrom[m]) == 1);
    end

    for m = 1:T
        if previousContToMove[m] != 0
            for t = 1:T
                # Relation between variables x and w for retrievals
                @constraint(CutModel, sum(w[previousContToMove[m],s,u] for u = 1:t-1 for s in moveFrom[previousContToMove[m]]) - sum(w[m,s,t] for s in moveFrom[m]) >= 0);
            end
        end
    end

    function mycutgenerator(cb)
        WtoBeCut = getvalue(w);
        z_val = getvalue(z);

        # Allow for some impreciseness in the solution
        TOL = 1e-3;

        (Z,P,Obj) = subProblem(WtoBeCut, H, artificialHeights, moveFrom, IOPoints, SR, SB, SX, SY, anteriorStacks, posteriorStacks, posCraneInitial, T, contMinHeightStack, previousContToMove, costMove, costPreMove, costToGo, alpha, 0, limitOfTime);

        if z_val < Obj - TOL
            @lazyconstraint(cb, z - sum(P[m,s,t]*w[m,s,t] for m = 1:T for s in moveFrom[m] for t = 1:T) >= Obj - Z);
        end

    end  # End of callback function

    # Tell JuMP/Gurobi to use our callback function
    addlazycallback(CutModel, mycutgenerator);

    # tic();
    status = solve(CutModel);
    # println("Solved !");
    # timeToSolve = toc();

    W = getvalue(w);

    return (W);
end
