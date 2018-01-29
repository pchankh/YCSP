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
    for instanceNumber in instancesToPlot
        outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),string("0_TotalResults"));
        data = readtable(outputFile, nrows = 1);
        push!(dataN,data[:Total_Cost][1]/nRequests);
    end
    p = p + 1;
    dfError[:N][p] = N;
    dfError[:mean_Cost][p] = mean(dataN);
    dfError[:conf_down][p] = mean(dataN)+1.645std(dataN);
    dfError[:conf_up][p] = mean(dataN)-1.645std(dataN);
end


myplot = plot(dfError, x=:N, y=:mean_Cost, ymin=:conf_down, ymax=:conf_up, Geom.point, Geom.errorbar, Guide.title("Impact of N with fixed flexibility"),Guide.xlabel("N",orientation=:horizontal), Guide.ylabel("Crane travel time per productive request",orientation=:vertical),Guide.yticks(ticks = collect(140:10:180)), Guide.xticks(ticks = NList,orientation=:horizontal), Theme(point_size = 3pt,background_color="white",grid_line_width=0.5mm,key_position=:right,key_title_font_size=14pt,key_label_font_size=14pt,major_label_font_size=14pt,minor_label_font_size=11pt,default_color = "green"));

draw(PDF("myplot_N.pdf", 20cm, 12cm), myplot);
