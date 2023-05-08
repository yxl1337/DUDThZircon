%% Define T, SiO2, and DU for multipleRegression_TSiO2DU

% Tzir has 2 columns, the first one is T zircon crystallization, which is the same as Txtal2 
% The second column is an arbitary number as 1-sigma absolute uncertainty
% DU, DTh are in mol%, DU and DTh are from Monte Carlo simulation with U
% and Th (mol%) in zircon and glass and their uncertainties 

T = zeros(height(DvsTSiO2source),3);
T(:,1) = Txtal2(:,1);
T(:,2) = ((1./Txtal2(:,1) - 1./(Txtal2(:,1)+Txtal2(:,2)))); %neg side for higher Temp
%Txtal2 (:,2) is the one std dev for T, but to plot error bar, have to calculate the length of 1 std 
T(:,3) = ((1./(Txtal2(:,1)-Txtal2(:,2)) - 1./Txtal2(:,1))); %pos side for lower tmep

%SiO2(:,2) = SiO2 .* DvsTSiO2source{:,27}/100;

% DU(1:20,1) = DvsTSiO2source{1:20,3};
% DU(21:26,1) = UmolZir(21:26)./UmolGls(21:26);
% DU(26:38,1) = DvsTSiO2source{26:38,3};
% DU(39:74,1) = UmolZir(39:74)./UmolGls(39:74);
% 
% 
% DTh(1:20,1) = DvsTSiO2source{1:20,5};
% DTh(21:26,1) = ThmolZir(21:26)./ThmolGls(21:26);
% DTh(26:38,1) = DvsTSiO2source{26:38,5};
% DTh(39:74,1) = ThmolZir(39:74)./ThmolGls(39:74);

DU = zeros(height(DvsTSiO2source),4);
DU(1:20,1) = DvsTSiO2source{1:20,3};
DU(21:66,1) = DUMC_mean(21:66,1);
DU(1:20,2) = DvsTSiO2source{1:20,4}; %real value for 1std for DU
DU(21:66,2) = DUMC_std(21:66,1);
DU(:,3) = log(DU(:,1)+DU(:,2)) - log(DU(:,1)); %legth for DU error bar pos side for larger DU value
DU(:,4) = log(DU(:,1)) - log(DU(:,1)-DU(:,2)); %legth for DU error bar neg side

DTh = zeros(height(DvsTSiO2source),4);
DTh(1:20,1) = DvsTSiO2source{1:20,5};
DTh(21:66,1) = DThMC_mean(21:66,1);
DTh(1:20,2) = DvsTSiO2source{1:20,6};
DTh(21:66,2) = DThMC_std(21:66,1);
DTh(:,3) = log(DTh(:,1)+DTh(:,2)) - log(DTh(:,1)); %legth for DTh error bar pos side
DTh(:,4) = log(DTh(:,1)) - log(DTh(:,1)-DTh(:,2)); %legth for DTh error bar neg side

DThDU = zeros(height(DvsTSiO2source),4);
DThDU(1:20,1) = DvsTSiO2source{1:20,5} ./ DvsTSiO2source{1:20,3};
DThDU(21:66,1) = DThDUMC_mean(21:66,1);
DThDU(1:20,2) = sqrt((1./DU(1:20,1).^2).*DTh(1:20,2).^2 + ...
                  (DTh(1:20,1).^2./DU(1:20,1).^4).*DU(1:20,2).^2);
DThDU(21:66,2) = DThDUMC_std(21:66,1);
DThDU(:,3) = log(DThDU(:,1)+DThDU(:,2)) - log(DThDU(:,1)); %legth for DTh/DU error bar pos side
DThDU(:,4) = log(DThDU(:,1)) - log(DThDU(:,1)-DThDU(:,2)); %legth for DTh/DU error bar neg side