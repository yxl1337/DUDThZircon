%% Plot DU and Dth measurements of experimental zircon CL light vs dark zones

%% read in data
opts = spreadsheetImportOptions; 
opts.VariableNames = {'mixName', 'CLBrightness', 'DU', 'DTh'};
opts.VariableTypes = {'string', 'categorical', 'double', 'double'};
opts.Sheet = 'Tidy Data';
opts.DataRange = 'E2:H288';
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

% calculate mean and overdispersion for DTh/DU differences
function L = loglik(data, uncts, MLE)
    mu = MLE(1);
    xi = MLE(2);
    varxi = uncts.^2 + xi;
    logliki = (data - mu).^2 ./ varxi + log(varxi);
    L = sum(logliki);
end % function loglik
% Perform maximum likelihood estimation for the log-likelihood function
initialMLE = [0, 1]; % Initial guesses for mu and xi
options = optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', 1e6, 'MaxFunEval', 1e6);
MLE = fminunc(@(params) loglik(diffs.mean_DThDU, uncts.ste_DThDU, params), initialMLE, options);


%% set up a plot

unctlinewidth = 10;
ylinewidth = 2;

figure
errorbar(diffs.mean_DU, 2*uncts.ste_DU, '.', 'LineWidth', ...
    unctlinewidth, 'CapSize', 0, 'MarkerSize', 35, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
set(gca, "FontSize", 16)
xlim([0, 13])
xticks(1:12)
xticklabels(dark.mixName)
yline(0, 'LineWidth', ylinewidth)
ylabel("DU_{dark} - DU_{light}", 'FontSize', 20)

figure
errorbar(diffs.mean_DTh, 2*uncts.ste_DTh, '.', 'LineWidth', ...
    unctlinewidth, 'CapSize', 0, 'MarkerSize', 35, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
set(gca, "FontSize", 16)
xlim([0, 13])
xticks(1:12)
xticklabels(dark.mixName)
yline(0, 'LineWidth', ylinewidth)
ylabel("DTh_{dark} - DTh_{light}", 'FontSize', 20)

%%
figure
hold on
line([0, 13], [0, 0], 'LineWidth', ylinewidth, 'Color', 'k')
line([0, 13], [MLE(1), MLE(1)], 'LineWidth', 2, 'Color', '#006400')
eb = errorbar(diffs.mean_DThDU, 2*uncts.ste_DThDU, '_', 'LineWidth', ...
    unctlinewidth, 'CapSize', 0, 'MarkerSize', 8, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
eb.Color = 'b';
set(gca, "FontSize", 16)
xlim([0, 13])
yregion(MLE(1)+2*sqrt(MLE(2)), MLE(1)-2*sqrt(MLE(2)), 'LineWidth', 1, 'FaceColor', '#77DD77')
xticks(1:12)
xticklabels(dark.mixName)
ylabel("DTh/DU_{dark} - DTh/DU_{light}", 'FontSize', 20)
exportgraphics(gca, 'DThDU_sector_zoning_updated.png', 'BackgroundColor', 'white')