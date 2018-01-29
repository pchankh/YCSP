using DataFrames, Gadfly, Cairo, Fontconfig

X = 7;
Y = 30;
Z = 4;
fillRate = 0.67;
IOPointsPosition = "Asian-right";
nRequests = 1500;
gap = X*Y;

gamma = Int64(50);
method = "Heuristic"

instancesToPlot = collect(1:30);
NList = [2,4,5,6,10];
nN = length(NList);

dfError = DataFrame(N = Array{Int64}(nN), mean_Cost = Array{Float64}(nN), conf_down = Array{Float64}(nN), conf_up = Array{Float64}(nN));
subFolder = string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap");
p = 0;
for N in NList
    dataN = Array{Float64}(0);
    flex = Int64(floor(N/2));
    for instanceNumber in instancesToPlot
        outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),string("0_",flex,"Flex_TotalResults"));
        data = readtable(outputFile, nrows = 1);
        push!(dataN,data[:Total_Cost][1]/nRequests);
    end
    p = p + 1;
    dfError[:N][p] = N;
    dfError[:mean_Cost][p] = mean(dataN);
    dfError[:conf_down][p] = mean(dataN)+1.645std(dataN);
    dfError[:conf_up][p] = mean(dataN)-1.645std(dataN);
end


myplot = plot(dfError, x=:N, y=:mean_Cost, ymin=:conf_down, ymax=:conf_up, Geom.point, Geom.errorbar, Guide.title("Impact of N with variable flexibility"),Guide.xlabel("N"), Guide.ylabel("Crane travel time per productive requests",orientation=:vertical),Guide.yticks(ticks = collect(140:10:180)), Guide.xticks(ticks = NList,orientation=:horizontal), Theme(background_color="white",key_position=:right,default_color = "green"));

myplot

draw(PDF("myplot_N_2.pdf", 28cm, 12cm), myplot);
