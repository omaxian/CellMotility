% Main file that runs all simulations for mechanism 2. 
% The results of this file are in Info_Mech2_Final.mat
maxTrial = 8;
protsperTrial = 4;
% In order: TEN, TNE, ETN, NTE
Kmvals = [0.1 10 0.1 1000];
Kcvals = [100 100 10 100];
Information = cell(maxTrial,4);
for iRegime=1:4
    % Random seed for reproducibility (mechanism by mechanism this time)
    rng(iRegime-1);
    nTrial=0;
    % In mechanism 2, there is a chance that there won't be 2 bindings over
    % 4 cycles. For this reason, we only accept the cycles that move a
    % certain distance (0.35) which is far enough to ensure there was at
    % least 1 binding. 
    while (nTrial < maxTrial)
        RammingCellFiberNetwork;
        % Save variables we need (only when cell actually moves)
        if (norm(nucpos(end,:)) > 0.35)
            nTrial=nTrial+1;
            Information{nTrial,iRegime} = {jPtslist, nucpos, AcorError, AnucError, ...
                CorAR, movieframes};
        end
    end
end
load('Info_Mech2_Final.mat');
% Compute the statistics
CorAs = zeros(maxTrial,4);
NucAs = zeros(maxTrial,4);
for iRegime=1:4
    NucARs{iRegime}=[];
    CorARs{iRegime}=[];
    for iTrial=1:maxTrial
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
xlabel('Total distance (4 cycles)')
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


