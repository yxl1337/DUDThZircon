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
    "mean", ...
    ["DU", "DTh", "DThDU"]...
    );
