using DataFrames, Gadfly, Cairo, Fontconfig

nInstances = 30;
Gammas = [0,25,50,75,100,200,500,1000];
X = 7;
Y = 30;
Z = 4;
fillRate = 0.67;
IOPointsPosition = "Asian-right";
nRequests = 1500;
gap = X*Y;
N = 5;
vZEmpty = 0.39;
vZLoaded = 0.20;

nPeriods = nRequests/N;
minGam = Int64(ceil(Z*(Z-1)*(1/vZEmpty+1/vZLoaded)));
subFolder = string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap");
nMethods = 2 * length(Gammas);
Results = DataFrame(Baseline = Array{Float64}(nInstances));
for i in 1:nInstances
    outputFile = joinpath("Outputs",subFolder,string(i),string(N),"FCFS_Myopic_Restricted","0","0_TotalResults");
    data = readtable(outputFile, nrows = 1);
    Results[:Baseline][i] = data[:Total_Cost][1]/nRequests;
end
for gamma in Gammas
    Results[Symbol(string("IP_","$gamma"))] = Array{Float64}(nInstances);
    for i in 1:nInstances
        outputFile = joinpath("Outputs",subFolder,string(i),string(N),"Flexible_Restricted",string(gamma),"0_TotalResults");
        data = readtable(outputFile, nrows = 1);
        Results[Symbol(string("IP_","$gamma"))][i] = data[:Total_Cost][1]/nRequests;
    end
    Results[Symbol(string("Heur_","$gamma"))] = Array{Float64}(nInstances);
    for i in 1:nInstances
        outputFile = joinpath("Outputs",subFolder,string(i),string(N),"Heuristic",string(gamma),"0_TotalResults");
        data = readtable(outputFile, nrows = 1);
        Results[Symbol(string("Heur_","$gamma"))][i] = data[:Total_Cost][1]/nRequests;
    end
end
writetable("ResultsExpe1.csv",Results);
Table_1 = DataFrame(Values = Array{AbstractString}(3));
Table_1[:Values][1] = "Average_Percent_Diff_Best";
Table_1[:Values][2] = "Std_Percent_Diff_Best";
Table_1[:Values][3] = "Ttest_Percent_Diff_Best";
bestMethod = Symbol();
bestAvg = Inf;
for method in names(Results)
    avg = 0;
    for i in 1:nInstances
        avg = avg + Results[method][i]/nInstances;
    end
    if avg < bestAvg
        bestMethod = method;
        bestAvg = avg;
    end
end
for method in names(Results)
    Table_1[method] = zeros(3);
    for i in 1:nInstances
        Table_1[method][1] = Table_1[method][1] + (Results[method][i] - Results[bestMethod][i])/(Results[bestMethod][i])*100;
    end
    Table_1[method][1] = Table_1[method][1]/nInstances;
    for i in 1:nInstances
        Table_1[method][2] = Table_1[method][2] + (((Results[method][i] - Results[bestMethod][i])/Results[bestMethod][i]*100) - Table_1[method][1])^2;
    end
    Table_1[method][2] = sqrt(Table_1[method][2]/(nInstances-1));
    Table_1[method][3] = Table_1[method][1] * sqrt(nInstances)/Table_1[method][2];
    Table_1[method][1] = round(Table_1[method][1],2);
    Table_1[method][2] = round(Table_1[method][2],2);
    Table_1[method][3] = round(Table_1[method][3],3);
end
writetable("Table_Alg.csv",Table_1);
