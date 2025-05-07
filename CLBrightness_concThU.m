%% Plot DU and Dth measurements of experimental zircon CL light vs dark zones

%% read in data
opts = spreadsheetImportOptions; 
opts.VariableNames = {'mixName', 'CLBrightness', 'DU', 'DTh'};
opts.VariableTypes = {'string', 'categorical', 'double', 'double'};
opts.Sheet = 'Tidy Data';
opts.DataRange = 'E2:H269';
opts.MissingRule = 'omitrow';

data = readtable("CL Brightness DUvsDTh mol%.xlsx", opts);

%% calculate DTh/DU

data.DThDU = data.DTh ./ data.DU;
mean_by_mix_and_CL = groupsummary( ...
    data, ...
    ["mixName", "CLBrightness"], ...
    "all", ...
    ["DU", "DTh", "DThDU"]...
    );

vars_keep = ["mixName", "CLBrightness", "GroupCount", ...
    "mean_DU", "std_DU", "mean_DTh", "std_DTh", "mean_DThDU", "std_DThDU"];
stats = mean_by_mix_and_CL(:, vars_keep);

stderr = @(col) stats(:,col) ./ sqrt(stats.GroupCount);
stats(:,["ste_DU", "ste_DTh", "ste_DThDU"]) = stderr(["std_DU", "std_DTh", "std_DThDU"]);

dark = stats(1:2:end,:);
light = stats(2:2:end,:);

meanColumns = ["mean_DU", "mean_DTh", "mean_DThDU"];
unctColumns = ["ste_DU", "ste_DTh", "ste_DThDU"];
diffs = dark(:,meanColumns) - light(:,meanColumns);
uncts = sqrt( dark(:,unctColumns).^2 + light(:,unctColumns).^2 );

%% set up a plot

figure
errorbar(diffs.mean_DU, 2*uncts.ste_DU, '.', 'LineWidth', 5)
yline(0)

figure
errorbar(diffs.mean_DU, 2*uncts.ste_DTh, '.', 'LineWidth', 5)
yline(0)

figure
errorbar(diffs.mean_DU, 2*uncts.ste_DThDU, '.', 'LineWidth', 5)
yline(0)