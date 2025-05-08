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

unctlinewidth = 10;
ylinewidth = 2;

figure
errorbar(diffs.mean_DU, 2*uncts.ste_DU, '.', 'LineWidth', ...
    unctlinewidth, 'CapSize', 0, 'MarkerSize', 35, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
set(gca, "FontSize", 16)
xlim([0, 12])
xticks(1:11)
xticklabels(dark.mixName)
yline(0, 'LineWidth', ylinewidth)
ylabel("DU_{dark} - DU_{light}", 'FontSize', 20)

figure
errorbar(diffs.mean_DTh, 2*uncts.ste_DTh, '.', 'LineWidth', ...
    unctlinewidth, 'CapSize', 0, 'MarkerSize', 35, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
set(gca, "FontSize", 16)
xlim([0, 12])
xticks(1:11)
xticklabels(dark.mixName)
yline(0, 'LineWidth', ylinewidth)
ylabel("DTh_{dark} - DTh_{light}", 'FontSize', 20)

figure
errorbar(diffs.mean_DThDU, 2*uncts.ste_DThDU, '.', 'LineWidth', ...
    unctlinewidth, 'CapSize', 0, 'MarkerSize', 35, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
set(gca, "FontSize", 16)
xlim([0, 12])
xticks(1:11)
xticklabels(dark.mixName)
yline(0, 'LineWidth', ylinewidth)
ylabel("DTh/DU_{dark} - DTh/DU_{light}", 'FontSize', 20)
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'DThDU_sector_zoning.png')