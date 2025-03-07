%% Weighted least squares regression of DU, DTh vs. fO2

% import data from Excel file
filepath = "Source_DvsfO2_Feb2025.xlsx";
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(filepath);

fO2_upperBound = -6;
fO2_lowerBound = -Inf;

DU_lowerBound = -Inf;
DU_upperBound = Inf;

fO2 = dataTable.fO2;

%% DU vs. fO2

% extract DU from data table
DU = dataTable.DUMCWGlsStdErr;
DU_1sAbs = dataTable.DUErrprWGlsStdErr;

% find rows without missing data
hasAllUData = all(~isnan([fO2, DU, DU_1sAbs]), 2);

% locate data with reasonable fO2 for crust, zircon formation
infO2Bounds = (fO2_lowerBound < fO2) & (fO2 < fO2_upperBound);
inDUBounds = (DU_lowerBound < DU) & (DU < DU_upperBound);

nUdata_ROI = sum(hasAllUData & infO2Bounds & inDUBounds);

% toss NaNs, enforce bounds
DU = DU(hasAllUData & infO2Bounds & inDUBounds);
DU_1sAbs = DU_1sAbs(hasAllUData & infO2Bounds & inDUBounds);
fO2_U = fO2(hasAllUData & infO2Bounds & inDUBounds);

% take logs
logDU = log(DU);
logDU_1sAbs = DU_1sAbs ./ DU; % approximation

% set up weighted least squares regression
% y = ax + b
designMatrix = [fO2_U ones(nUdata_ROI,1)];
weightVector = logDU_1sAbs.^-2;
dataCovMat = diag(logDU_1sAbs.^2);

% calculate excess variance needed to bring mse = 1
odfun = @(od) lscov_od(designMatrix, logDU, dataCovMat, od);
od = fsolve(odfun, 1); % overdispersion
dataCov_od = dataCovMat + od*eye(nUdata_ROI); 

% perform linear regression with overdispersion term
[x, stdx, mse, S] = lscov(designMatrix, logDU, dataCov_od);
logDU_1sAbs_od = sqrt(diag(dataCov_od));

% re-import out-of-range DU measurements for plotting
DU_all = dataTable.DUMCWGlsStdErr;
DU_1sAbs_all = dataTable.DUErrprWGlsStdErr;
logDU_highfO2 = log(DU_all(hasAllUData & ~infO2Bounds & inDUBounds));
logDU_1sAbs_highfO2 = DU_1sAbs_all(hasAllUData & ~infO2Bounds & inDUBounds)...
                   ./ exp(logDU_highfO2);
fO2_highfO2 = fO2(hasAllUData & ~infO2Bounds & inDUBounds);


%% plot U results

figure("Name", "log(DU) vs fO2")
t = tiledlayout(1,1);
ax1 = axes(t);

plot(ax1, fO2_U, logDU, '.', 'MarkerSize', 15)
hold on
line(ax1, [fO2_U'; fO2_U'], [logDU' - 2*logDU_1sAbs_od'; logDU' + 2*logDU_1sAbs_od'], ...
    'Color', rgb('Barney'), 'LineWidth', 2)
line(ax1, [fO2_U'; fO2_U'], [logDU' - 2*logDU_1sAbs'; logDU' + 2*logDU_1sAbs'], ...
    'Color', rgb('Vibrant Blue'), 'LineWidth', 1.5)

xbuf = 0.5; npts = 500;
xlimits = [min(fO2_U) - xbuf max(fO2_U) + xbuf];
xvector = linspace(xlimits(1), xlimits(2), npts)';
yvector = xvector * x(1) + x(2);
positionVector = [xvector ones(npts, 1)];

uncty = zeros(npts, 1);
for idx = 1:npts
    uncty(idx) = sqrt( positionVector(idx,:) * S * positionVector(idx,:)');
end

plot(ax1, xvector, yvector, '-', "Color", rgb('Dark Pastel Green'), ...
    "LineWidth", 2)
plot(ax1, xvector, yvector + 2*uncty, '-', "Color", rgb('Pine'), "LineWidth", 1.5)
plot(ax1, xvector, yvector - 2*uncty, '-', "Color", rgb('Pine'), "LineWidth", 1.5)

set(ax1, "FontSize", 18)
xlabel(ax1, "log({\itf}O2)", "FontSize", 24)
ylabel(ax1, "log(DU)", "FontSize", 24)

currentXLim = xlim(gca);
newXLim(1) = max([fO2_lowerBound, currentXLim(1), xlimits(1)]);
newXLim(2) = min([fO2_upperBound, currentXLim(2), xlimits(2)]) + 0.2;
xlim(ax1, newXLim)
ax1.XTickMode = "auto";

% % uncomment if plotting high-fO2 points
% plot(fO2_highfO2, logDU_highfO2, '.b', 'MarkerSize', 15)
% line([fO2_highfO2'; fO2_highfO2'], ...
%     [logDU_highfO2' - 2*logDU_1sAbs_highfO2'; logDU_highfO2' + 2*logDU_1sAbs_highfO2'], ...
%     'Color', rgb('black'), 'LineWidth', 2)
% xlim([newXLim(1) 0])

annotation("textbox", [0.4 0.14 0.5, 0.1], ...
    "String", ...
    "Slope = " + round(x(1),2) + " ± " + round(2*stdx(1),2) + " (2\sigma)", ...
    "FontSize", 20, ...
    "HorizontalAlignment", "center", ...
    "VerticalAlignment","middle", ...
    "LineStyle", "none")

ax2 = axes(t);
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
set(ax2, "FontSize", 17)

ax1xlim = xlim(ax1);
T = 1000 + 273.15; % assume a temperature in K for calcs
FMQ = 9 - 25738/T; % per Mike
ax2xlim = ax1xlim - FMQ;
xlim(ax2, ax2xlim)
xlabel(ax2, "log({\itf}O2), \DeltaQFM", "FontSize", 23)
ylabel(ax2, "")
ax2.XTickMode = "auto";
ax2.YTickLabel = [];
ax2.YTick = [];

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