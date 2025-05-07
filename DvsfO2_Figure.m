% import data from Excel file
filepath = "Source_DvsfO2_May2025.xlsx";
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(filepath);

% set up figure 
figure("Name", "log(DU) vs fO2", 'Position', [1 1 1100 500])
t = tiledlayout(1,2);
nexttile; nexttile;

for output = ["logfO2", "deltaQFM"]


% pick one, comment the other
%output = "deltaQFM";
%output = "logfO2";

% set dataset and its regression bounds
if output == "logfO2"
    itile = 1;
    fO2_upperBound = -6;
    fO2_lowerBound = -15;
    fO2 = dataTable.logfO2;
elseif output == "deltaQFM"
    itile = 2;
    fO2_upperBound = 5;
    fO2_lowerBound = -Inf;
    fO2 = dataTable.deltaQFM;
else, disp("unknown output format")
end

DU_lowerBound = -Inf;
DU_upperBound = Inf;

% data sources and correpsonding rows in dataTable
sourceLabels = [
    "this study, EPMA"; 
    "Burnham and Berry (2012), SIMS"; 
    "Luo and Ayers (2009), LA-ICPMS"; 
    "Rubatto and Hermann (2007), LA-ICPMS"; 
    "Ayers and Peters (2018), LA-ICPMS"];
sourceRows = {1:20; 21:34; 35:46; 47:52; 53:57};
sourceMarkerShape = ["o", "square", "o", "o", "o"];
sourceMarkerColor = ['k', 'k', "#A9A9A9", 'b', 'r'];


%% Data handling

% extract DU from data table
DU_all = dataTable.DUMCWGlsStdErr;
DU_1sAbs = dataTable.DUErrprWGlsStdErr;

% take logs
logDU = log(DU_all);
logDU_1sAbs = DU_1sAbs ./ DU_all; % approximation

% find rows without missing data
hasAllUData = all(~isnan([fO2, DU_all, DU_1sAbs]), 2);

% locate data with reasonable fO2 for crust, zircon formation
infO2Bounds = (fO2_lowerBound < fO2) & (fO2 < fO2_upperBound);
inDUBounds = (DU_lowerBound < DU_all) & (DU_all < DU_upperBound);

inROI = hasAllUData & infO2Bounds & inDUBounds;
nUdata_ROI = sum(inROI);

% toss NaNs, enforce bounds, subset data for regression
logDU_regression = logDU(inROI);
logDU_1sAbs_regression = logDU_1sAbs(inROI);
fO2_regression = fO2(inROI);

% set up weighted least squares regression
% y = ax + b
designMatrix = [fO2_regression ones(nUdata_ROI,1)];
weightVector = logDU_1sAbs_regression.^-2;
dataCovMat = diag(logDU_1sAbs_regression.^2);

% calculate excess variance needed to bring mse = 1
odfun = @(od) lscov_od(designMatrix, logDU_regression, dataCovMat, od);
od = fsolve(odfun, 1); % overdispersion
dataCov_od = dataCovMat + od*eye(nUdata_ROI); 

% perform linear regression with overdispersion term
[x, stdx, mse, S] = lscov(designMatrix, logDU_regression, dataCov_od);
logDU_1sAbs_od = sqrt(diag(dataCov_od));


%% plot U results

ax = t.Children(itile);
ax.NextPlot = 'add';

nSources = length(sourceRows);
marks = gobjects([nSources, 1]);
for iSource = 1:length(sourceRows)

    label = sourceLabels(iSource);
    rows = sourceRows{iSource};
    markershape = sourceMarkerShape(iSource);
    markercolor = sourceMarkerColor(iSource);
    sourcefO2 = fO2(rows);
    sourceLogDU = logDU(rows);
    sourceLogDU_1sAbs = logDU_1sAbs(rows);

    marks(iSource) = plot(ax, sourcefO2, logDU(rows), ...
    "Marker", markershape, ...
    "MarkerFaceColor", markercolor, ...
    "MarkerEdgeColor", 'k', ...
    "LineStyle", "none", ...
    "MarkerSize", 9);

    line(ax, [sourcefO2'; sourcefO2'], ...
        [sourceLogDU'-2*sourceLogDU_1sAbs'; sourceLogDU'+2*sourceLogDU_1sAbs'], ...
        'Color', markercolor, 'LineWidth', 1.5)

end

% Add regression and uncertainty envelope
xbuf = 0.5; npts = 500;
xlimits = [min(fO2_regression) - xbuf max(fO2_regression) + xbuf];
xvector = linspace(xlimits(1), xlimits(2), npts)';
yvector = xvector * x(1) + x(2);
positionVector = [xvector ones(npts, 1)];

uncty = zeros(npts, 1);
for idx = 1:npts
    uncty(idx) = sqrt( positionVector(idx,:) * S * positionVector(idx,:)');
end

plot(ax, xvector, yvector, '-', "Color", rgb('Dark Pastel Green'), ...
    "LineWidth", 2)
plot(ax, xvector, yvector + 2*uncty, '-', "Color", rgb('Pine'), "LineWidth", 1.5)
plot(ax, xvector, yvector - 2*uncty, '-', "Color", rgb('Pine'), "LineWidth", 1.5)

% Axes and scale formatting
set(ax, "FontSize", 18)
ylabel(ax, "log(DU)", "FontSize", 24)

currentXLim = xlim(ax);
newXLim(1) = max([fO2_lowerBound, currentXLim(1), xlimits(1)]);
ax.XTickMode = "auto";

if output == "logfO2"
    xlim(ax, [newXLim(1) 0]) % for plotting fO2
    xlabel(ax, "log({\itf}O2)", "FontSize", 24)
    annotation_left = 0.55;
elseif output == "deltaQFM"
    xlim(ax, [newXLim(1) 7]) % for plotting fO2-deltaQFM
    xlabel(ax, "log({\itf}O2) \DeltaQFM", "FontSize", 24)
    annotation_left = 0.1;
else, disp("did not recognize output format")
end

annotation("textbox", [annotation_left 0.11 0.3, 0.1], ...
    "String", ...
    "Slope = " + round(x(1),2) + " ± " + round(2*stdx(1),2) + " (2\sigma)", ...
    "FontSize", 20, ...
    "HorizontalAlignment", "center", ...
    "VerticalAlignment","middle", ...
    "LineStyle", "none")

end % for output = [QFM, logfO2]

%% add legend
ax = t.Children(1);
lgd = legend(ax, marks, sourceLabels);
lgd.Position = lgd.Position + [0.08 0.07 0 0];

set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'DvsfO2_side-by-side.png')

%% Local functions

% use in nonlinear equation solver to force mse = 1 for fit by adding
% overdispersion (excess variance) of od to each analytical uncertainty term
function mseMinusOne = lscov_od(designMatrix, logDU, dataCovMat, od)
    
    n = length(logDU);
    dataCovOD = dataCovMat + od*eye(n);
    [~, ~, mse, ~] = lscov(designMatrix, logDU, dataCovOD);
    mseMinusOne = mse - 1;

end

function rgbTriple = rgb(name)

    [~, rgbTriple] = colornames('xkcd', name);

end
