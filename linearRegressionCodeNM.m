% linear regression for logD vs 1/T fits, plots

%addpath('../McLeanLinearRegression')

% Run Yuanyuan's code:
importDataFromSpreadsheet;
% adjust temperatures to soak/start/saturation 
lnDVSTxal;
% add uncertainty column (1s absolute)
InputMultipleRegression;

% Results in DTh, DU, and DThDU. Column 1 value, col 2 is 1s abs

% remove outliers, per WUSTL
remSamples = 54:55;
T(remSamples,:) = [];
DU(remSamples,:) = [];
DTh(remSamples,:) = [];
DThDU(remSamples,:) = [];
% see also volcanic range on line 31


%% set plot properties here

ourExperiments = 1:20;
ourExperimentsColor = 'k';
otherExperiments = 21:35;
otherExperimentsColor = 'b';
plutonic = 36:40;
plutonicColor = [11, 163, 59]/255; % dark green
volcanic = 41:64; %41:66; % removed two outliers
volcanicColor = 'r';

fitLineWidth = 2;
fitLineColor = 'b';
envelopeLineWidth = 1;
envelopeColor = 'g';
envelopeSigmaLevel = 1;

dataMarkerSize = 25;
unctBarLineWidth = 1;
unctBarColor = 'k';
unctBarSigmaLevel = 1;

figureFontSize = 16;
axisLabelFontSize = 20;

% common T axis parameters
nT = 500;
Trange = max(T(:,1)) - min(T(:,1));
Tvec = linspace(min(T(:,1)) - 0.2*Trange, ...
                max(T(:,1)) + 0.2*Trange, nT)';
TvecC = Tvec - 273.15;

minmaxTC = [600, 1450]; % for axes limits on plots vs T (C)


%% Do regression for ln(DU) vs. 1/T

nSamples = size(T,1);
dataunct = zeros(nSamples,5);

dataunct(:,1) = 1./T(:,1);
dataunct(:,2) = T(:,2);
dataunct(:,3) = log(DU(:,1));
dataunct(:,4) = DU(:,2)./DU(:,1);

skipv = ones(nSamples,1);
a1 = 0; v1 = 1; abspct = 1;

McLeanLinearRegression;

figure("Name", "log(DU) vs. 1/T")
hold on

% plot data points
ourE = plot(1./T(ourExperiments,1), log(DU(ourExperiments,1)), '.', ...
    'MarkerEdgeColor', ourExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
othE = plot(1./T(otherExperiments,1), log(DU(otherExperiments,1)), '.', ...
    'MarkerEdgeColor', otherExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
plut = plot(1./T(plutonic,1), log(DU(plutonic,1)), '.', ...
    'MarkerEdgeColor', plutonicColor, ...
    'MarkerSize', dataMarkerSize);
volc = plot(1./T(volcanic,1), log(DU(volcanic,1)), '.', ...
    'MarkerEdgeColor', volcanicColor, ...
    'MarkerSize', dataMarkerSize);

% plot uncertainties
logDU1s = DU(:,2)./DU(:,1);
line(1./([T(:,1)'; T(:,1)']), ...
     [log(DU(:,1))' - unctBarSigmaLevel*logDU1s'; ...
      log(DU(:,1))' + unctBarSigmaLevel*logDU1s'], ...
     'Color', 'k')
line([1./T(:,1)' - unctBarSigmaLevel*T(:,2)'; ...
      1./T(:,1)' + unctBarSigmaLevel*T(:,2)'], ...
     [log(DU(:,1))'; ...
     log(DU(:,1))'], ...
     'Color', unctBarColor, 'LineWidth', unctBarLineWidth)

% plot best fit line
S_DU = MSWD*Sav;
av_DU(1) = a(2);
av_DU(2) = v(2);
MSWD_DU = MSWD;

plot(1./Tvec, av_DU(1) + av_DU(2)./Tvec, '-', ...
    'Color', fitLineColor, 'LineWidth', fitLineWidth)
residualVariance = var(log(DU(:,1)) - (av_DU(1)+av_DU(2)./T(:,1)));

for iT = 1:nT
    logDU_fits2(iT) = [1 1/Tvec(iT,1)]*S_DU*[1; 1/Tvec(iT,1)];
end

logDU_fit1s = sqrt(logDU_fits2 + residualVariance);

plot(1./Tvec, av_DU(1) + av_DU(2)./Tvec - ...
     envelopeSigmaLevel*logDU_fit1s, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)
plot(1./Tvec, av_DU(1) + av_DU(2)./Tvec + ...
     envelopeSigmaLevel*logDU_fit1s, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)
  
xlim([min(1./Tvec), max(1./Tvec)])
set(gca, 'FontSize', figureFontSize)
xlabel("1/T  (K^{-1})", 'FontSize', axisLabelFontSize)
ylabel("log(DU)", 'FontSize', axisLabelFontSize)

legend([ourE othE plut volc], ...
    {'Our Experiments', 'Other Experiments', ...
     'Plutonic Samples', 'Volcanic Samples'}, 'Location', 'NW')
annDU = annotation('textbox', [0.63 0.1 0.3 0.1], 'String', ...
    "all uncertainties \pm" + num2str(unctBarSigmaLevel) + "\sigma");
annDU.FontSize = figureFontSize;
annDU.EdgeColor = "none";


%% Do regression for ln(DTh) vs. 1/T

dataunct = zeros(nSamples,5);

dataunct(:,1) = 1./T(:,1);
dataunct(:,2) = T(:,2);
dataunct(:,3) = log(DTh(:,1));
dataunct(:,4) = DTh(:,2)./DTh(:,1);

skipv = ones(nSamples,1);
a1 = 0; v1 = 1; abspct = 1;

McLeanLinearRegression;

figure("Name", "log(DTh) vs. 1/T")
hold on

% plot data points
ourE = plot(1./T(ourExperiments,1), log(DTh(ourExperiments,1)), '.', ...
    'MarkerEdgeColor', ourExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
othE = plot(1./T(otherExperiments,1), log(DTh(otherExperiments,1)), '.', ...
    'MarkerEdgeColor', otherExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
plut = plot(1./T(plutonic,1), log(DTh(plutonic,1)), '.', ...
    'MarkerEdgeColor', plutonicColor, ...
    'MarkerSize', dataMarkerSize);
volc = plot(1./T(volcanic,1), log(DTh(volcanic,1)), '.', ...
    'MarkerEdgeColor', volcanicColor, ...
    'MarkerSize', dataMarkerSize);

% plot uncertainties
logDTh1s = DTh(:,2)./DTh(:,1);
line(1./([T(:,1)'; T(:,1)']), ...
     [log(DTh(:,1))' - unctBarSigmaLevel*logDTh1s'; ...
      log(DTh(:,1))' + unctBarSigmaLevel*logDTh1s'], ...
     'Color', 'k')
line([1./T(:,1)' - unctBarSigmaLevel*T(:,2)'; ...
      1./T(:,1)' + unctBarSigmaLevel*T(:,2)'], ...
     [log(DTh(:,1))'; ...
     log(DTh(:,1))'], ...
     'Color', unctBarColor, 'LineWidth', unctBarLineWidth)

% plot best fit line
S_DTh = MSWD*Sav;
av_DTh(1) = a(2);
av_DTh(2) = v(2);
MSWD_DTh = MSWD;

plot(1./Tvec, av_DTh(1) + av_DTh(2)./Tvec, '-', ...
    'Color', fitLineColor, 'LineWidth', fitLineWidth)
residualVariance = var(log(DTh(:,1)) - (av_DTh(1)+av_DTh(2)./T(:,1)));

for iT = 1:nT
    logDTh_fits2(iT) = [1 1/Tvec(iT,1)]*S_DTh*[1; 1/Tvec(iT,1)];
end

logDTh_fit1s = sqrt(logDTh_fits2 + residualVariance);

plot(1./Tvec, av_DTh(1) + av_DTh(2)./Tvec - ...
     envelopeSigmaLevel*logDTh_fit1s, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)
plot(1./Tvec, av_DTh(1) + av_DTh(2)./Tvec + ...
     envelopeSigmaLevel*logDTh_fit1s, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)

xlim([min(1./Tvec), max(1./Tvec)])
set(gca, 'FontSize', figureFontSize)
xlabel("1/T (K^{-1})", 'FontSize', axisLabelFontSize)
ylabel("log(DTh)", 'FontSize', axisLabelFontSize)

legend([ourE othE plut volc], ...
    {'Our Experiments', 'Other Experiments', ...
     'Plutonic Samples', 'Volcanic Samples'}, 'Location', 'NW')
annDTh = annotation('textbox', [0.63 0.1 0.3 0.1], 'String', ...
    "all uncertainties \pm" + num2str(unctBarSigmaLevel) + "\sigma");
annDTh.FontSize = figureFontSize;
annDTh.EdgeColor = "none";



%% Do regression for ln(DTh/DU) vs. 1/T

dataunct = zeros(nSamples,5);

dataunct(:,1) = 1./T(:,1);
dataunct(:,2) = T(:,2);
dataunct(:,3) = log(DThDU(:,1));
dataunct(:,4) = DThDU(:,2)./DThDU(:,1);

skipv = ones(nSamples,1);
a1 = 0; v1 = 1; abspct = 1;

McLeanLinearRegression;

figure("Name", "log(DTh/DU) vs. 1/T")
hold on

% plot data points
ourE = plot(1./T(ourExperiments,1), log(DThDU(ourExperiments,1)), '.', ...
    'MarkerEdgeColor', ourExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
othE = plot(1./T(otherExperiments,1), log(DThDU(otherExperiments,1)), '.', ...
    'MarkerEdgeColor', otherExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
plut = plot(1./T(plutonic,1), log(DThDU(plutonic,1)), '.', ...
    'MarkerEdgeColor', plutonicColor, ...
    'MarkerSize', dataMarkerSize);
volc = plot(1./T(volcanic,1), log(DThDU(volcanic,1)), '.', ...
    'MarkerEdgeColor', volcanicColor, ...
    'MarkerSize', dataMarkerSize);

% plot uncertainties
DThDU1s = DThDU(:,2)./DThDU(:,1);
line(1./([T(:,1)'; T(:,1)']), ...
     [log(DThDU(:,1))' - unctBarSigmaLevel*DThDU1s'; ...
      log(DThDU(:,1))' + unctBarSigmaLevel*DThDU1s'], ...
     'Color', 'k',  'LineWidth', unctBarLineWidth)
line([1./T(:,1)' - unctBarSigmaLevel*T(:,2)'; ...
      1./T(:,1)' + unctBarSigmaLevel*T(:,2)'], ...
     [log(DThDU(:,1))'; ...
     log(DThDU(:,1))'], ...
     'Color', unctBarColor, 'LineWidth', unctBarLineWidth)

% plot best fit line
S_DThDU = MSWD*Sav;
av_DThDU(1) = a(2);
av_DThDU(2) = v(2);
MSWD_DThDU = MSWD;

plot(1./Tvec, av_DThDU(1) + av_DThDU(2)./Tvec, '-', ...
    'Color', fitLineColor, 'LineWidth', fitLineWidth)
%residualVariance = var(log(DThDU(:,1)) - (av_DThDU(1)+av_DThDU(2)./T(:,1)));
resid = log(DThDU(:,1)) - (av_DThDU(1)+av_DThDU(2)./T(:,1));
vardiff = max(resid.^2 - DThDU(:,2).^2, 0);
residualVariance = mean(vardiff); % pretty rough estimate

for iT = 1:nT
    logDThDU_fits2(iT) = [1 1/Tvec(iT,1)]*S_DThDU*[1; 1/Tvec(iT,1)];
end

logDThDU_fit1s = sqrt(logDThDU_fits2 + residualVariance);

plot(1./Tvec, av_DThDU(1) + av_DThDU(2)./Tvec - ...
     envelopeSigmaLevel*logDThDU_fit1s, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)
plot(1./Tvec, av_DThDU(1) + av_DThDU(2)./Tvec + ...
     envelopeSigmaLevel*logDThDU_fit1s, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)

xlim([min(1./Tvec), max(1./Tvec)])
set(gca, 'FontSize', figureFontSize)
xlabel("1/T  (K^{-1})", 'FontSize', axisLabelFontSize)
ylabel("log(DTh/DU)", 'FontSize', axisLabelFontSize)

legend([ourE othE plut volc], ...
    {'Our Experiments', 'Other Experiments', ...
     'Plutonic Samples', 'Volcanic Samples'}, 'Location', 'south')
annDThDU = annotation('textbox', [0.63 0.82 0.3 0.1], 'String', ...
    "all uncertainties \pm" + num2str(unctBarSigmaLevel) + "\sigma");
annDThDU.FontSize = figureFontSize;
annDThDU.EdgeColor = "none";


%% now convert to DUDTh vs T


T1s = T(:,2)./(1./T(:,1)).*T(:,1);
DThDU_fit = exp( av_DThDU(1) + av_DThDU(2)./Tvec );

DThDU_fit1s = zeros(nT,1);
for iT = 1:nT

    dfda = exp( av_DThDU(1) + av_DThDU(2)/Tvec(iT) );
    dfdv = exp( av_DThDU(1) + av_DThDU(2)/Tvec(iT) ) / Tvec(iT);

    s2i = [dfda dfdv]*S_DThDU*[dfda; dfdv] +  dfda*residualVariance;
    DThDU_fit1s(iT) = sqrt(s2i);
end


figure("Name", "DTh/DU vs. T")
hold on

% plot data points
ourE = plot(T(ourExperiments,1)-273.15, DThDU(ourExperiments,1), '.', ...
    'MarkerEdgeColor', ourExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
othE = plot(T(otherExperiments,1)-273.15, DThDU(otherExperiments,1), '.', ...
    'MarkerEdgeColor', otherExperimentsColor, ...
    'MarkerSize', dataMarkerSize);
plut = plot(T(plutonic,1)-273.15, DThDU(plutonic,1), '.', ...
    'MarkerEdgeColor', plutonicColor, ...
    'MarkerSize', dataMarkerSize);
volc = plot(T(volcanic,1)-273.15, DThDU(volcanic,1), '.', ...
    'MarkerEdgeColor', volcanicColor, ...
    'MarkerSize', dataMarkerSize);

% plot uncertainties
DThDU1s = DThDU(:,2);
line(([T(:,1)'; T(:,1)']-273.15), ...
     max(0,[DThDU(:,1)' - unctBarSigmaLevel*DThDU1s'; ...
      DThDU(:,1)' + unctBarSigmaLevel*DThDU1s']), ...
     'Color', 'k',  'LineWidth', unctBarLineWidth)
line([T(:,1)' - unctBarSigmaLevel*T1s'; ...
      T(:,1)' + unctBarSigmaLevel*T1s']-273.15, ...
     [DThDU(:,1)'; ...
      DThDU(:,1)'], ...
     'Color', unctBarColor, 'LineWidth', unctBarLineWidth)

% plot best fit line
plot(TvecC, DThDU_fit, '-', ...
    'Color', fitLineColor, 'LineWidth', fitLineWidth)

plot(TvecC, max(0,DThDU_fit - DThDU_fit1s), ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)
plot(TvecC, DThDU_fit + DThDU_fit1s, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)

xlim(minmaxTC) % min and max T (C)
set(gca, 'FontSize', figureFontSize)
xlabel("T (C)", 'FontSize', axisLabelFontSize)
ylabel("DTh/DU", 'FontSize', axisLabelFontSize)

legend([ourE othE plut volc], ...
    {'Our Experiments', 'Other Experiments', ...
     'Plutonic Samples', 'Volcanic Samples'}, 'Location', 'NE')
annDThDU = annotation('textbox', [0.13 0.83 0.3 0.1], 'String', ...
    "all uncertainties \pm" + num2str(unctBarSigmaLevel) + "\sigma");
annDThDU.FontSize = figureFontSize;
annDThDU.EdgeColor = "none";



%% and delta-t vs. T

lambda238 = 1.55125e-10; % year^-1
lambda230 = 9.1705e-6;   % year^-1

age = 50*10^6; % years
r206238 = exp(lambda238*age)-1;

t1 = 1/lambda238 * log(r206238 + 1 - ...
                       lambda238/lambda230*(DThDU_fit - 1));
t2 = 1/lambda238 * log(r206238 + 1);
deltat_fit = t1 - t2;

% rough in uncertainties for asymmetric unct prop

DThDUupper = exp(av_DThDU(1) + av_DThDU(2)./Tvec + ...
     envelopeSigmaLevel*logDThDU_fit1s);
DThDUlower = exp(av_DThDU(1) + av_DThDU(2)./Tvec - ...
     envelopeSigmaLevel*logDThDU_fit1s);
t1Upper = 1/lambda238 * log(r206238 + 1 - ...
                       lambda238/lambda230*(DThDUupper - 1));
t1Lower = 1/lambda238 * log(r206238 + 1 - ...
                       lambda238/lambda230*(DThDUlower - 1));

deltatUpper = t1Upper - t2;
deltatLower = t1Lower - t2;



figure
hold on
plot(TvecC, deltatUpper/1000, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)
plot(TvecC, deltatLower/1000, ...
     'Color', envelopeColor, 'LineWidth', envelopeLineWidth)
plot(TvecC, deltat_fit/1000, '-', ...
    'Color', fitLineColor, 'LineWidth', fitLineWidth)

xlim(minmaxTC);
set(gca, 'FontSize', figureFontSize)
xlabel('T (C)', 'FontSize', axisLabelFontSize)
ylabel('\Deltat (kyr)', 'FontSize', axisLabelFontSize)
annDThDU = annotation('textbox', [0.17 0.13 0.3 0.1], 'String', ...
    "uncertainties \pm" + num2str(unctBarSigmaLevel) + "\sigma");
annDThDU.FontSize = figureFontSize;
annDThDU.EdgeColor = "none";



%% symbolic math for delta-t
% 
% syms t r68 dthdu lambda238 lambda230
% 
% t1 = 1/lambda238 * log( r68 + 1 - lambda238/lambda230*(dthdu - 1));
% t2 = 1/lambda238 * log(r68 + 1);
% deltat = t2 - t1;
% simplify(diff(deltat, dthdu),100)

%% old DThDU vs T plotting

% figure
% plot(TvecC, DThDU_fit, '-b', 'LineWidth', 2)
% hold on
% 
% plot(TvecC, DThDU_fit - DThDU_fit1s, '-g', 'LineWidth', 1)
% plot(TvecC, DThDU_fit + DThDU_fit1s, '-g', 'LineWidth', 1)
% 
% plot(T(:,1)-273.15, DThDU(:,1), '.k', 'MarkerSize', 25)
% line([T(:,1)'; T(:,1)']-273.15, ...
%      [DThDU(:,1)'-DThDU(:,2)'; DThDU(:,1)'+DThDU(:,2)'], ...
%      'Color', unctBarColor, 'LineWidth', unctBarLineWidth)
% line([T(:,1)'-T1s'; T(:,1)'+T1s']-273.15, ...
%      [DThDU(:,1)'; DThDU(:,1)'], ...
%      'Color', unctBarColor, 'LineWidth', unctBarLineWidth)
% 
% xlabel('T (C)', 'FontSize', 16)
% ylabel('DTh/DU', 'FontSize', 16)
