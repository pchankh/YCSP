using DataFrames, Gadfly, Cairo, Fontconfig

X = 7;
Y = 30;
Z = 4;
fillRate = 0.67;
IOPointsPosition = "Asian-right";
nRequests = 1500;
gap = X*Y;
N = 5;
nPeriods = nRequests/N;
vZEmpty = 0.39;
vZLoaded = 0.20;
minGam = Int64(ceil(Z*(Z-1)*(1/vZEmpty+1/vZLoaded)));

instancesToPlot = collect(1);
Gammas = [0,25,50,75,100,200,500,1000];
methodnonGamma = ["FCFS_Myopic_Restricted"];
methodGamma = ["Heuristic","Flexible_Restricted"];
nPointsError = length(methodnonGamma) + length(methodGamma) * length(Gammas);

dfError = DataFrame(Algorithms = Array{AbstractString}(nPointsError), gamma = Array{Int64}(nPointsError), mean_Cost = Array{Float64}(nPointsError), conf_down = Array{Float64}(nPointsError), conf_up = Array{Float64}(nPointsError), non_Opt = Array{AbstractString}(nPointsError));
subFolder = string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap");
p = 0;
baseLine = 0;
for method in methodnonGamma
    gamma = 0;
    dataMethod = Array{Float64}(0);
    nonOpt = 0;
    for instanceNumber in instancesToPlot
        outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),"0_TotalResults");
        data = readtable(outputFile, nrows = 1);
        push!(dataMethod,data[:Total_Cost][1]/nRequests);
        nonOpt = nonOpt + data[:non_Optimal][1]/(nPeriods*length(instancesToPlot));
    end
    p = p + 1;
    if method == "FCFS_Myopic_Restricted"
        dfError[:Algorithms][p] = "Baseline";
    elseif method == "FCFS_Myopic"
        dfError[:Algorithms][p] = "Baseline_2";
    elseif method == "Heuristic"
        dfError[:Algorithms][p] = "Heuristic";
    elseif method == "Flexible_Restricted"
        dfError[:Algorithms][p] = "IP";
        dfError[:solved][p] = data[]
    elseif method == "Flexible"
        dfError[:Algorithms][p] = "IP_relax";
    end
    dfError[:gamma][p] = gamma;
    dfError[:mean_Cost][p] = mean(dataMethod);
    baseLine = mean(dataMethod);
    dfError[:conf_down][p] = mean(dataMethod)+1.645std(dataMethod);
    dfError[:conf_up][p] = mean(dataMethod)-1.645std(dataMethod);
    if dfError[:Algorithms][p] == "IP" || dfError[:Algorithms][p] == "Baseline"
        dfError[:non_Opt][p] = string(round(nonOpt,2));
    else
        dfError[:non_Opt][p] = "";
    end
end
for method in methodGamma
    for gamma in Gammas
        dataMethod = Array{Float64}(0);
        nonOpt = 0;
        for instanceNumber in instancesToPlot
            if method == "Heuristic"
                outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),"0_TotalResults");
            elseif method == "Flexible_Restricted"
                outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),"0_VeryLongResults.csv");
            end
            data = readtable(outputFile, nrows = 1);
            push!(dataMethod,data[:Total_Cost][1]/nRequests);
            nonOpt = nonOpt + data[:non_Optimal][1]/(nPeriods*length(instancesToPlot));
        end
        p = p + 1;
        if method == "FCFS_Myopic_Restricted"
            dfError[:Algorithms][p] = "Baseline";
        elseif method == "Heuristic"
            dfError[:Algorithms][p] = "Heuristic";
        elseif method == "Flexible_Restricted"
            dfError[:Algorithms][p] = "IP";
        elseif method == "Flexible"
            dfError[:Algorithms][p] = "IP_relax";
        end
        dfError[:gamma][p] = gamma;
        dfError[:mean_Cost][p] = mean(dataMethod);
        dfError[:conf_down][p] = mean(dataMethod)+1.645std(dataMethod);
        dfError[:conf_up][p] = mean(dataMethod)-1.645std(dataMethod);
        if dfError[:Algorithms][p] == "IP" || dfError[:Algorithms][p] == "Baseline"
            dfError[:non_Opt][p] = string(round(nonOpt,2));
        else
            dfError[:non_Opt][p] = "";
        end
    end
end

myplot = plot(dfError, x=:gamma, y=:mean_Cost, ymin=:conf_down, ymax=:conf_up, color=:Algorithms, Geom.point, Geom.errorbar, Guide.xlabel("γ"), Guide.title("Performance of algorithms"), Guide.ylabel("Crane travel time per productive request",orientation=:vertical), Guide.xticks(ticks = Gammas,orientation=:horizontal), Theme(point_size = 3pt,background_color="white",grid_line_width=0.5mm,key_position=:right,key_title_font_size=14pt,key_label_font_size=14pt,major_label_font_size=14pt,minor_label_font_size=11pt), Scale.color_discrete_manual("red","green","purple"), xintercept=[minGam], Geom.vline(style=[[1mm,1mm]],color=["blue"]), yintercept=[baseLine], Geom.hline(style=[[1mm,1mm]],color="red"));

myplot

draw(PDF("myplot_Alg.pdf", 27cm, 13cm), myplot)
