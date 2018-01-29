using DataFrames, Gadfly, Cairo, Fontconfig

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
minGam = Int64(ceil(Z*(Z-1)*(1/vZEmpty+1/vZLoaded)));

instancesToPlot = [1,2];
Gammas = [100,200,500];
methodnonGamma = ["FCFS_Myopic_Restricted"];
methodGamma = ["Heuristic","Flexible_Restricted"];
nPoints = length(instancesToPlot) * (length(methodnonGamma) + length(methodGamma) * length(Gammas));
nPointsError = length(methodnonGamma) + length(methodGamma) * length(Gammas);
# df = DataFrame(instance = Array{Int64}(nPoints), method = Array{AbstractString}(nPoints), gamma = Array{Float64}(nPoints), Total_Cost = Array{Float64}(nPoints));
# p = 0;
# subFolder = string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap");
# for instanceNumber in instancesToPlot
#     for method in methodnonGamma
#         gamma = 0;
#         outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),"0_TotalResults");
#         data = readtable(outputFile, nrows = 1);
#         p = p + 1;
#         df[:instance][p] = instanceNumber;
#         df[:method][p] = method;
#         df[:gamma][p] = gamma;
#         df[:Total_Cost][p] = data[:Total_Cost][1];
#     end
#     for method in methodGamma
#         for gamma in Gammas
#             outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),"0_TotalResults");
#             data = readtable(outputFile, nrows = 1);
#             p = p + 1;
#             df[:instance][p] = instanceNumber;
#             df[:method][p] = method;
#             df[:gamma][p] = gamma;
#             df[:Total_Cost][p] = data[:Total_Cost][1];
#         end
#     end
# end
# plot(df, x=:gamma, y=:Total_Cost, color=:method, Geom.point, shape=:instance)

dfError = DataFrame(Algorithms = Array{AbstractString}(nPointsError), gamma = Array{Int64}(nPointsError), mean_Cost = Array{Float64}(nPointsError), conf_down = Array{Float64}(nPointsError), conf_up = Array{Float64}(nPointsError));
nPointsError = length(methodnonGamma) + length(methodGamma) * length(Gammas);
subFolder = string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap");
p = 0;
baseLine = 0;
for method in methodnonGamma
    gamma = 0;
    dataMethod = Array{Float64}(0);
    for instanceNumber in instancesToPlot
        outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),"0_TotalResults");
        data = readtable(outputFile, nrows = 1);
        push!(dataMethod,data[:Total_Cost][1]/nRequests);
    end
    p = p + 1;
    dfError[:Algorithms][p] = method;
    dfError[:gamma][p] = gamma;
    dfError[:mean_Cost][p] = mean(dataMethod);
    baseLine = mean(dataMethod);
    dfError[:conf_down][p] = mean(dataMethod)+std(dataMethod);
    dfError[:conf_up][p] = mean(dataMethod)-std(dataMethod);
end
for method in methodGamma
    for gamma in Gammas
        dataMethod = Array{Float64}(0);
        for instanceNumber in instancesToPlot
            outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),"0_TotalResults");
            data = readtable(outputFile, nrows = 1);
            push!(dataMethod,data[:Total_Cost][1]/nRequests);
        end
        p = p + 1;
        dfError[:Algorithms][p] = method;
        dfError[:gamma][p] = gamma;
        dfError[:mean_Cost][p] = mean(dataMethod);
        dfError[:conf_down][p] = mean(dataMethod)+std(dataMethod);
        dfError[:conf_up][p] = mean(dataMethod)-std(dataMethod);
    end
end

myplot = plot(dfError, x=:gamma, y=:mean_Cost, ymin=:conf_down, ymax=:conf_up, color=:Algorithms, Geom.point, Geom.errorbar, Guide.xlabel("γ"), Guide.ylabel("Crane travel time per productive requests"), Guide.xticks(ticks = Gammas), Theme(background_color="white",key_position=:right),Guide.manual_color_key("Algorithms", ["Baseline","NRLP","IP"], ["red","green","purple"]),Scale.color_discrete_manual("red","green","purple"),xintercept=[minGam],Geom.vline(style=[[1mm,1mm]]),yintercept=[baseLine],Geom.hline(style=[[1mm,1mm]],color="red"));

draw(PDF("myplot.pdf", 10cm, 10cm), myplot)
