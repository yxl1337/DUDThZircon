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


%% plot U results

figure("Name", "log(DU) vs fO2")
plot(fO2_U, logDU, '.', 'MarkerSize', 15)
hold on
line([fO2_U'; fO2_U'], [logDU' - 2*logDU_1sAbs_od'; logDU' + 2*logDU_1sAbs_od'], ...
    'Color', rgb('Barney'), 'LineWidth', 2)
line([fO2_U'; fO2_U'], [logDU' - 2*logDU_1sAbs'; logDU' + 2*logDU_1sAbs'], ...
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

plot(xvector, yvector, '-', "Color", rgb('Dark Pastel Green'), ...
    "LineWidth", 2)
plot(xvector, yvector + 2*uncty, '-', "Color", rgb('Pine'), "LineWidth", 1.5)
plot(xvector, yvector - 2*uncty, '-', "Color", rgb('Pine'), "LineWidth", 1.5)

set(gca, "FontSize", 18)
xlabel("fO2", "FontSize", 24)
ylabel("log(DU)", "FontSize", 24)

currentXLim = xlim(gca);
newXLim(1) = max([fO2_lowerBound, currentXLim(1), xlimits(1)]);
newXLim(2) = min([fO2_upperBound, currentXLim(2), xlimits(2)]);
xlim(newXLim)

annotation("textbox", [0.4 0.12 0.5, 0.1], ...
    "String", "Slope = " + round(x(1),2) + " ± " + round(2*stdx(1),2), ...
    "FontSize", 18, ...
    "HorizontalAlignment", "center", ...
    "VerticalAlignment","middle", ...
    "LineStyle", "none")


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