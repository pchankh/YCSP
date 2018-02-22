using DataFrames, Gadfly, Cairo, Fontconfig

X = 7;
Y = 30;
Z = 4;
fillRate = 0.67;
IOPointsPosition = "Asian-right";
nRequests = 1500;
gap = X*Y;
N = 5;
gamma = Int64(50);
method = "Heuristic"

# instancesToPlot = setdiff(collect(1:30),[19,20,22,23,29]);
instancesToPlot = collect(1:30);
Deltas = collect(0:5);
nDeltas = length(Deltas)

dfError = DataFrame(delta = Array{Int64}(nDeltas), mean_Cost = Array{Float64}(nDeltas), conf_down = Array{Float64}(nDeltas), conf_up = Array{Float64}(nDeltas));
subFolder = string(X,"X_",Y,"Y_",Z,"Z_",fillRate,"fill_",IOPointsPosition,"_",nRequests,"req_",gap,"gap");
p = 0;
for d in Deltas
    dataDelta = Array{Float64}(0);
    for instanceNumber in instancesToPlot
        outputFile = joinpath("Outputs",subFolder,string(instanceNumber),string(N),method,string(gamma),string("0_",d,"Flex_TotalResults"));
        data = readtable(outputFile, nrows = 1);
        push!(dataDelta,data[:Total_Cost][1]/nRequests);
    end
    dfError[:delta][d+1] = d;
    dfError[:mean_Cost][d+1] = mean(dataDelta);
    dfError[:conf_down][d+1] = mean(dataDelta)+0.310*std(dataDelta);
    dfError[:conf_up][d+1] = mean(dataDelta)-0.310*std(dataDelta);
end


myplot = plot(dfError, x=:delta, y=:mean_Cost, ymin=:conf_down, ymax=:conf_up, Geom.point, Geom.errorbar, Guide.title("Impact of flexibility level"),Guide.xlabel("δ",orientation=:horizontal), Guide.ylabel("Crane travel time per productive request",orientation=:vertical),Guide.yticks(ticks = collect(140:10:180)), Guide.xticks(ticks = Deltas,orientation=:horizontal), Theme(point_size = 3pt,background_color="white",grid_line_width=0.5mm,key_position=:right,key_title_font_size=14pt,key_label_font_size=14pt,major_label_font_size=14pt,minor_label_font_size=11pt,default_color = "green"));

draw(PDF("myplot_Delta.pdf", 20cm, 12cm), myplot);
