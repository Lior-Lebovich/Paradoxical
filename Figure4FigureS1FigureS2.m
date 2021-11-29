%% Figure 4, Supplementary figure 1, Supplementary figure 2:


% Note: Because we were interested in the effect of backgrounds on both the 
% signal and noise and as we used a custom code, in which the diffusion 
% variance is set to 1 (see Pythin code), the fitted parameters are scaled 
% to quantify the backgrounds effect on the drift and diffusion variance 
% (see Method and the code below). 


% Figure 4: Average drift (A) and diffusion (B) as a function of the number 
% of background odorants.
figure4 = figure;

% Supplementary figure 1 - Comparison of experimental and DDM predicted 
% decisions (A), and decision times (B).
figureS1A = figure;
figureS1B = figure;

% Supplementary figure 2 - Drift (A) and Diffusion (B) extracted from the 
% DDMs fit to individual mice.
figureS2 = figure;



% Read qunatile optimization output:

cdName = 'D:\all hddm analysis\mice3\quantiles';

nMice = 6;
targetNames = {'A','B','N'};
minTarN = [1,1,0];
targeLongNames = {'target A','target B','no target'};
targetColors = [0 0 1; .49 .18 .56;1 0 0];
fitTypes = {'chi','ML'};
modelNames = {'regular','fixedThresh','fixedDrift'};


% plot parameters for the chi-square DDM fit: 

fitType = 'chi'; % to compare with ML change to 'ML'.
modelName = 'regular'; 

tableParams = readtable( [cdName '\' modelName '\subj_params_miceq_' ...
    fitType '.csv'] );

param_tar_n_m.c2New = nan(3,6,6);
param_tar_n_m.vNew = nan(3,6,6);

for tar = 1:length(targetNames)
    targetName = targetNames{tar};
    for m = 1:nMice
        paramsSub = tableParams(m,:);
        for n = 1:6
            nDist = n - minTarN(tar);
            aOld = paramsSub.(['a_' targetName num2str(nDist) '_']);
            vOld = paramsSub.(['v_' targetName num2str(nDist) '_']);
            param_tar_n_m.c2New(tar,n,m) = (1 ./ aOld) ^ 2;
            param_tar_n_m.vNew(tar,n,m) = vOld ./ aOld;
        end
        % Supplementary figure 2 
        % plot by-mouse DDM parameters:
        figure(figureS2);
        subplot(2,6,m+6);
        plot( 1-minTarN(tar):6-minTarN(tar), ...
            param_tar_n_m.c2New(tar,:,m), 'Color', targetColors(tar,:), ...
            'Marker', '.', 'lineWidth', 1 ); hold on;
        subplot(2,6,m);
        plot( 1-minTarN(tar):6-minTarN(tar), ...
            param_tar_n_m.vNew(tar,:,m), 'Color', targetColors(tar,:), ...
            'Marker', '.', 'lineWidth', 1 ); hold on;
    end
end


% Figure 4
% Plot DDM parameters, averaged over mice:
figure(figure4);
byMiceParam.c2.on = [reshape( param_tar_n_m.c2New(1,:,:), [6,6,1] )'; ...
    reshape( param_tar_n_m.c2New(2,:,:), [6,6,1] )'];
byMiceParam.c2.off = reshape( param_tar_n_m.c2New(3,:,:), [6,6,1] )';
byMiceParam.v.on = [reshape( param_tar_n_m.vNew(1,:,:), [6,6,1] )'; ...
    reshape( param_tar_n_m.vNew(2,:,:), [6,6,1] )'];
byMiceParam.v.off = reshape( param_tar_n_m.vNew(3,:,:), [6,6,1] )';
tO_Names = {'on', 'off'};
tars = [1,3];
for tO = 1:length(tO_Names)
    tar = tars(tO);
    tO_Name = tO_Names{tO};
    subplot(1,2,2);
    c2New_avgMice = mean( byMiceParam.c2.(tO_Name) );
    c2New_semMice = std( byMiceParam.c2.(tO_Name) ) ./ ...
        sqrt( length(byMiceParam.c2.(tO_Name)) );
    errorbar( 1-minTarN(tar):6-minTarN(tar), c2New_avgMice, c2New_semMice, ...
        'Color', targetColors(tar,:), ...
        'lineWidth', 1, 'CapSize', 0.1, 'Marker', 'o', 'MarkerSize', 7,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', targetColors(tar,:) );
    hold on;
    subplot(1,2,1);
    vNew_avgMice = mean( byMiceParam.v.(tO_Name) );
    vNew_semMice = std( byMiceParam.v.(tO_Name) ) ./ ...
        sqrt( length(byMiceParam.v.(tO_Name)) );
    errorbar( 1-minTarN(tar):6-minTarN(tar), vNew_avgMice, vNew_semMice, ...
        'Color', targetColors(tar,:), ...
        'lineWidth', 1, 'CapSize', 0.1, 'Marker', 'o', 'MarkerSize', 7,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', targetColors(tar,:) );
    hold on;
end
subplot(1,2,1); 
plot( [0,6], [0,0], 'k--' );
ylabel('drift'); xlabel('#BG odorant');
xticks(0:6); box off; ylim([-2,2]); xlim([-.1,6.1]);
subplot(1,2,2); 
ylabel('diffusion (c^2)'); xlabel('#BG odorant');
xticks(0:6); box off; ylim([0,1]); xlim([-.1,6.1]);


% Supplementary figure 2
figure(figureS2);
subplot(2,6,7); ylabel('diffusion (c^2)'); 
subplot(2,6,1); ylabel('drift');
subplot(2,6,3); xlabel('#BG odorant');
subplot(2,6,9); xlabel('#BG odorant');
for k = 1:12
    subplot(2,6,k); xticks(0:6); xlim([-.1,6.1]); box off;
    if k <= 6
        plot( [0,6], [0,0], 'k--' );
    end
end



% Supplementary figure 1
% compare model predictions (100 sims per mouse-target-#BG) to observed p 
% and mRT:


% load/create csv table of behavioral data:
isBehavTable = exist('goodAllMice_id_stim_rtNorm_choiceQ_0.csv'); 
if isBehavTable == 2
    tableObs = readtable('goodAllMice_id_stim_rtNorm_choiceQ_0.csv');
else 
    prepareRawData
    tableObs = readtable('goodAllMice_id_stim_rtNorm_choiceQ_0.csv');
end

% load DDM-based simulated data:
tableExp = readtable( ['DDM\simData_miceq_' fitType '.csv'] );

p.Obs.p = nan(3,6,6);
p.Exp.p = nan(3,6,6);
p.Obs.sem = nan(3,6,6);
p.Exp.sem = nan(3,6,6);
RT.Obs.med = nan(3,6,6);
RT.Exp.med = nan(3,6,6);
RT.Obs.madmed = nan(3,6,6);
RT.Exp.madmed = nan(3,6,6);


% Read and plot by-mouse RT obs/sim: 

for tar = 1:length(targetNames)
    targetName = targetNames{tar};
    for m = 1:nMice
        paramsSub = tableParams(m,:);
        for n = 1:6
            nDist = n - minTarN(tar);
            
            % Read observed data:
            tableObs_tarMouse = tableObs( tableObs.subj_idx == m-1 & ...
                strcmp( tableObs.stim, [targetName num2str(nDist)] ), : );
            p.Obs.p(tar,n,m) = mean( tableObs_tarMouse.response );
            p.Obs.sem(tar,n,m) = sqrt( p.Obs.p(tar,n,m) * ...
                (1-p.Obs.p(tar,n,m)) / size(tableObs_tarMouse,1) );
            RT.Obs.med(tar,n,m) = median( tableObs_tarMouse.rt );
            RT.Obs.madmed(tar,n,m) = mad( tableObs_tarMouse.rt, 1 );
            
            % Read simulated data:
            tableExp_tarMouse = tableExp( ...
                strcmp( tableExp.condition, ['subj_' num2str(m-1) ...
                '__stim_' targetName num2str(nDist)] ), : );
            p.Exp.p(tar,n,m) = mean( tableExp_tarMouse.response );
            p.Exp.sem(tar,n,m) = sqrt( p.Exp.p(tar,n,m) * ...
                (1-p.Exp.p(tar,n,m)) / size(tableExp_tarMouse,1) );
            RT.Exp.med(tar,n,m) = median( tableExp_tarMouse.rt );
            RT.Exp.madmed(tar,n,m) = mad( tableExp_tarMouse.rt, 1 );
            
        end
        
        % Plot by-mouse P obs/sim:
        figure(figureS1A);
        subplot(3,6,m+6*(tar-1));
        patch( [(1-minTarN(tar)):(6-minTarN(tar)), ...
            flip((1-minTarN(tar)):(6-minTarN(tar)))], ...
            [p.Obs.p(tar,:,m)+p.Obs.sem(tar,:,m), ...
            flip(p.Obs.p(tar,:,m)-p.Obs.sem(tar,:,m))], targetColors(tar,:), ...
            'EdgeColor', 'none', 'FaceAlpha', .3 ); hold on;
        patch( [(1-minTarN(tar)):(6-minTarN(tar)), ...
            flip((1-minTarN(tar)):(6-minTarN(tar)))], ...
            [p.Exp.p(tar,:,m)+p.Exp.sem(tar,:,m), ...
            flip(p.Exp.p(tar,:,m)-p.Exp.sem(tar,:,m))], ...
            'k', 'EdgeColor', 'none', 'FaceAlpha', .3 ); 
        hold on;
        plot( (1-minTarN(tar)):(6-minTarN(tar)), p.Obs.p(tar,:,m), 'Color', targetColors(tar,:), ...
            'lineWidth', 1, 'Marker', '.' ); hold on;
        plot( (1-minTarN(tar)):(6-minTarN(tar)), p.Exp.p(tar,:,m), ...
            'color', 'k', 'lineWidth', 1, 'Marker', '.' );
        hold on;
        xticks(0:6); xlim([-.1,6.1]); ylim([0,1]); box off;
        if tar == 1
            title(['m = ' num2str(m)]);
        end
        
        % Plot by-mouse RT obs/sim:
        figure(figureS1B);
        subplot(3,6,m+6*(tar-1));
        patch( [(1-minTarN(tar)):(6-minTarN(tar)), ...
            flip((1-minTarN(tar)):(6-minTarN(tar)))], ...
            [RT.Obs.med(tar,:,m)+RT.Obs.madmed(tar,:,m), ...
            flip(RT.Obs.med(tar,:,m)-RT.Obs.madmed(tar,:,m))], targetColors(tar,:), ...
            'EdgeColor', 'none', 'FaceAlpha', .3 ); hold on;
        patch( [(1-minTarN(tar)):(6-minTarN(tar)), ...
            flip((1-minTarN(tar)):(6-minTarN(tar)))], ...
            [RT.Exp.med(tar,:,m)+RT.Exp.madmed(tar,:,m), ...
            flip(RT.Exp.med(tar,:,m)-RT.Exp.madmed(tar,:,m))], ...
            'k', 'EdgeColor', 'none', 'FaceAlpha', .3 ); 
        hold on;
        plot( (1-minTarN(tar)):(6-minTarN(tar)), RT.Obs.med(tar,:,m), 'color', targetColors(tar,:), ...
            'lineWidth', 1, 'Marker', '.' ); hold on;
        plot( (1-minTarN(tar)):(6-minTarN(tar)), RT.Exp.med(tar,:,m), ...
            'k', 'lineWidth', 1, 'Marker', '.' );
        hold on;
        xticks(0:6); xlim([-.1,6.1]); box off; %ylim([0,1.1]);
        if tar == 1
            title(['m = ' num2str(m)]);
        end
    end
end

figure(figureS1A); 
subplot(3,6,7); ylabel('p');
subplot(3,6,15); xlabel('#BG odorants');
figure(figureS1B);
subplot(3,6,7); ylabel('median RT (s)');
subplot(3,6,15); xlabel('#BG odorants');

