% Main file that runs all simulations for mechanism 1. 
% The results of this file are in Info_Mech1_Final.mat
% Random seed for reproducibility
rng(0);
nTrial = 8;
protsperTrial = 6;
% In order: TEN, TNE, ETN, NTE
Kmvals = [1 10 1 100];
Kcvals = [100 100 1 50];
Information = cell(nTrial,4);
for iRegime=1:4
    for iTrial = 1:nTrial
        PullingCellFiberNetwork;
        % Save variables we need
        Information{iTrial,iRegime} = {jPtslist, nucpos, AcorError, AnucError, ...
            CorAR};
    end
end
load('Info_Mech1_Final.mat')
% Compute the statistics
CorAs = zeros(nTrial,4);
NucAs = zeros(nTrial,4);
for iRegime=1:4
    NucARs{iRegime}=[];
    CorARs{iRegime}=[];
    for iTrial=1:nTrial
        nucds = Information{iTrial,iRegime}{1,2};
        normds(iTrial)=norm(nucds(end,:));
        thisprotds=nucds-[0 0; nucds(1:length(ppoints)-1,:)];
        protds(:,iTrial)=sqrt(thisprotds(2:end,1).^2+thisprotds(2:end,2).^2);
        CorAs(iTrial,iRegime)=mean(abs(Information{iTrial,iRegime}{1,3}))/(pi*rc^2);
        NucAs(iTrial,iRegime)=mean(abs(Information{iTrial,iRegime}{1,4}))/(pi*rm^2);
        CorARs{iRegime}=[CorARs{iRegime} Information{iTrial,iRegime}{1,5}];
    end
    allCARs{iRegime}=max(CorARs{iRegime},1./CorARs{iRegime});
    allNARs{iRegime}=max(NucARs{iRegime},1./NucARs{iRegime});
    allnds(:,iRegime)=normds(:);
    allprots(:,iRegime)=protds(:);
end
% Plot the results
figure;
errorbar(mean(allprots),[4 3 1 2],std(allprots)/sqrt(nTrial*protsperTrial),...
    'horizontal','o','LineWidth',2.0,'MarkerSize',8)
xlabel('Distance/protrusion')
yticks([1:4])
ylim([0.5 4.5])
yticklabels({'$E > T > N$','$N > T > E$','$T > N > E$','$T > E > N$'});
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',16)
figure;
errorbar(mean(allnds),[4 3 1 2],std(allnds)/sqrt(nTrial*protsperTrial),...
    'horizontal','o','LineWidth',2.0,'MarkerSize',8)
xlabel('Total distance (6 cycles)')
yticks([1:4])
ylim([0.5 4.5])
yticklabels({'$E > T > N$','$N > T > E$','$T > N > E$','$T > E > N$'});
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',16)
figure;
errorbar([mean(allCARs{1}) mean(allCARs{2}) mean(allCARs{3}) mean(allCARs{4})],...
    [4 3 1 2],[std(allCARs{1}) std(allCARs{2}) std(allCARs{3}) std(allCARs{4})]/sqrt(nTrial),...
    'horizontal','o','LineWidth',2.0,'MarkerSize',8)
xlabel('Aspect ratio')
yticks([1:4])
ylim([0.5 4.5])
yticklabels({'$E > T > N$','$N > T > E$','$T > N > E$','$T > E > N$'});
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',16)
% Compute mean areas
meanCAerror = mean(CorAs(:))
meanNAerror = mean(NucAs(:))


