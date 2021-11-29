%% Figure 2, Figure 3, Supplementary figure 3:



% Figure 2 - The effect of background odorants on decision probabilities:
% A-C. The probability of reporting that the target odorant is present (p) 
% as a function of the number of background odorants, for the three 
% interference hypotheses: false signal (A), signal reduction (B), and 
% noise boost (C). D. The fraction of trials that mice reported that the 
% target is present as a function of the number of background odorants. 
figure2 = figure;

% Figure 3 - The effect of background odorants on decision times: 
% A. Mean decision time (see equation 5, Methods) under the “noise-boost” 
% mechanism as a function of the number of background odorants. B. The mean
% mouse reaction time as a function of the number of background odorants. 
figure3 = figure;

% Supplementary figure 3 - The number of target-on and target-off trials 
% performed by each mouse for each number of background odorants.
figureS3 = figure;



% DDM predictions: 


% Define p and mean(DT) function over DDM moel parameters (solutions of 
% first passage problem, see [1] and [2]):

% [1] Bogacz, R., Brown, E., Moehlis, J., Holmes, P. and Cohen, J.D., 2006. 
% The physics of optimal decision making: a formal analysis of models of 
% performance in two-alternative forced-choice tasks. Psychological review, 
% 113(4), p.700.
% [2] Srivastava, V., Holmes, P. and Simen, P., 2016. Explicit moments of 
% decision times for single-and double-threshold drift-diffusion processes. 
% Journal of Mathematical Psychology, 75, pp.96-109.

fun_p = @(drift,diffusion,threshold,SP)(...
    max([(drift~=0) .* (...
    ( exp( 2 * drift .* threshold ./ diffusion ) - ...
    exp( -2 * SP .* drift ./ diffusion ) ) ./ ...
    ( exp( 2 * drift .* threshold ./ diffusion ) - ...
    exp( -2 * drift .* threshold ./ diffusion ) ) ); ...
    (drift==0) .* (...
    ones( size(diffusion) ) * (threshold+SP) ./ (2*threshold)  )]) ...
    );

fun_EDT =  @(drift,diffusion,threshold,SP)(...
    max([(drift~=0) .* (...
    ( (threshold ./ drift) .* ...
    coth( 2 * drift .* threshold ./ diffusion ) ) - ...
    ( (threshold ./ drift) .* ...
    exp( -2 * SP .* drift ./ diffusion ) .* ...
    csch( 2 * drift .* threshold ./ diffusion ) ) - ...
    (SP ./ drift) ); ...
    (drift==0) .* ...
    ( ( (threshold.^2) - (SP.^2) ) ./ diffusion )]) ...
    );

% baseline parameters:
thresh = 1; 
SP = 0 * thresh; 
blDiffusion = 3;
blDrifts = [-2.5, 2];

% target and distractors vect:
targetCols = {'r', 'b'};


% Plot models predictions (Figure 2A-C, Figure 3A):

for tar = 1:length(blDrifts)
    blDrift = blDrifts(tar);
    if blDrift > 0
        dist_vect = 0:.001:5; 
        dist_vect6 = 0:5;
    else
        dist_vect = 1:.001:6; 
        dist_vect6 = 1:6;
    end
    
    % *** model 1: False signal ***
    % model 1: drift = target + driftbias + #BG-odors
    drift = blDrift + dist_vect;
    drift6 = blDrift + dist_vect6;
    diffusion = blDiffusion;
    p = fun_p(drift,diffusion,thresh,SP);
    p6 = fun_p(drift6,diffusion,thresh,SP);
    model1.P.(['tar' num2str(tar-1)]) = p;
    % plot p as a function of #BG-odors:
    figure(figure2); subplot(1,4,1);
    plot( dist_vect, p, ...
        targetCols{tar}, 'lineWidth', 2, ...
        'lineStyle', '-' ); hold on;
    plot( dist_vect6, p6, targetCols{tar}, 'marker', '.', ...
        'markerSize', 10,'lineStyle', 'none' ); hold on;
    xlabel('#BG odorants'); ylabel('p'); xlim([0,7]); ylim([0,1]);
    yticks(0:.25:1); xticks(0:2:6);
    title('False signal');
    
    % *** model 2: Signal reduction ***
    % model 2: drift = target + driftbias - #BG-odors
    drift = blDrift - dist_vect;
    drift6 = blDrift - dist_vect6;
    diffusion = blDiffusion;
    p = fun_p(drift,diffusion,thresh,SP);
    p6 = fun_p(drift6,diffusion,thresh,SP);
    model2.P.(['tar' num2str(tar-1)]) = p;
    % plot p as a function of #BG-odors:
    figure(figure2); subplot(1,4,2);
    plot( dist_vect, p, ...
        targetCols{tar}, 'lineWidth', 2, ...
        'lineStyle', '-' ); hold on;
    plot( dist_vect6, p6, targetCols{tar}, 'marker', '.', ...
        'markerSize', 10,'lineStyle', 'none' ); hold on;
    xlabel('#BG odorants'); ylabel('p'); xlim([0,7]); ylim([0,1]);
    yticks(0:.25:1); xticks(0:2:6);
    title('Signal reduction');
    
    % *** model 3: Noise boost ***
    % model 3: drift = target + driftbias
    % model 3: diffusion = 1 + #BG-odors
    m = 3;
    drift = blDrift;
    diffusion = blDiffusion + dist_vect; 
    diffusion10 = blDiffusion + dist_vect6;
    p = fun_p(drift,diffusion,thresh,SP);
    p6 = fun_p(drift,diffusion10,thresh,SP);
    dt = fun_EDT(drift,diffusion,thresh,SP);
    dt6 = fun_EDT(drift,diffusion10,thresh,SP);
    model3.P.(['tar' num2str(tar-1)]) = p;
    model3.EDT.(['tar' num2str(tar-1)]) = dt;
    % plot p as a function of #BG-odors:
    figure(figure2); subplot(1,4,3);
    plot( dist_vect, p, ...
        targetCols{tar}, 'lineWidth', 2, ...
        'lineStyle', '-' ); hold on;
    plot( dist_vect6, p6, targetCols{tar}, 'marker', '.', ...
        'markerSize', 10,'lineStyle', 'none' ); hold on;
    xlabel('#BG odorants'); ylabel('p'); xlim([0,7]); ylim([0,1]);
    yticks(0:.25:1); xticks(0:2:6);
    title('Noise boost');
    % plot mean RT as a function of #BG-odors:
    figure(figure3); subplot(1,2,1);
    plot( dist_vect, 1000 * dt, ...
        targetCols{tar}, 'lineWidth', 2, ...
        'lineStyle', '-' ); hold on;
    plot( dist_vect6, 1000 * dt6, targetCols{tar}, 'marker', '.', ...
        'markerSize', 10,'lineStyle', 'none' ); hold on;
    xlabel('#BG odorants'); ylabel('mean decision time (a.u.)'); xticks(0:2:6);
    xlim([0,7]); 
    title('Noise boost');

end



% Behavioral data:


clearvars -except figure2 figure3 figureS3;

% load/create csv table of behavioral data:
isBehavTable = exist('goodAllMice_id_stim_rtNorm_choiceQ_0.csv'); 
if isBehavTable == 2
    tableObs = readtable('goodAllMice_id_stim_rtNorm_choiceQ_0.csv');
else 
    prepareRawData
    tableObs = readtable('goodAllMice_id_stim_rtNorm_choiceQ_0.csv');
end


targetNames = {'A','B','N'};
minTarN = [1,1,0];
nMice = 6;
targetColors = [0 0 1; .49 .18 .56;1 0 0];

NNN = nan(3,6,6);
p.Obs.p = nan(3,6,6);
RT.Obs.mean = nan(3,6,6);


% compute p and RT for each mouse in each of the 6x3 conditions:
for tar = 1:length(targetNames)
    targetName = targetNames{tar};
    for m = 1:nMice
        for n = 1:6
            nDist = n - minTarN(tar);
            % observed data:
            tableObs_tarMouse = tableObs( tableObs.subj_idx == m-1 & ...
                strcmp( tableObs.stim, [targetName num2str(nDist)] ), : );
            NNN(tar,n,m) = length( tableObs_tarMouse.response );
            p.Obs.p(tar,n,m) = mean( tableObs_tarMouse.response );
            RT.Obs.mean(tar,n,m) = mean( tableObs_tarMouse.rt );
        end
    end
end


% plot p and RT for each mouse in each of the 6x2 conditions (Target On vs 
% Tarfet off):
for tar = [1,3]
    pPlt = nan(6,6);
    mRtPlt = nan(6,6);
    targetName = targetNames{tar};
    tarVect = (1-minTarN(tar)):(6-minTarN(tar));
    for m = 1:nMice
        if tar == 1
            pPlt(m,:) = mean( p.Obs.p(1:2,:,m), 1 );
            mRtPlt(m,:) = mean( RT.Obs.mean(1:2,:,m), 1 );
            nTrialsPlt = sum( NNN(1:2,:,m), 1 );
        elseif tar == 3
            pPlt(m,:) = p.Obs.p(3,:,m);
            mRtPlt(m,:) = RT.Obs.mean(3,:,m);
            nTrialsPlt = NNN(3,:,m);
        end
        % plot by-mouse P obs:
        figure(figure2); subplot(1,4,4);
        plot( tarVect, pPlt(m,:), 'Color', targetColors(tar,:) ); hold on;
        % plot by-mouse mean RT obs:
        figure(figure3); subplot(1,2,2);
        plot( tarVect, mRtPlt(m,:), 'color', targetColors(tar,:) ); hold on;
        % plot by-mouse #trials:
        figure(figureS3); subplot(2,3,m);
        plot( tarVect, nTrialsPlt, 'color', targetColors(tar,:) ); hold on;
        xticks(0:2:6); xlim([0,7]); box off;
        xlabel('#BG odorants'); ylabel('# trials'); 
        title(['mouse ' num2str(m)])
    end
   % plot P obs averaged over mice:
    figure(figure2); subplot(1,4,4);
    plot( tarVect, mean( pPlt ), 'Color', targetColors(tar,:), ...
        'lineWidth', 2, 'Marker', '.', 'MarkerSize', 20 ); hold on;
    xticks(0:2:6); xlim([0,7]); ylim([0,1]); box off;
    xlabel('#BG odorants'); ylabel('p'); 
    title('Experiment');
    % plot mean RT obs averaged over mice:
    figure(figure3); subplot(1,2,2);
    plot( tarVect, mean( mRtPlt ), 'color', targetColors(tar,:), ...
        'lineWidth', 2, 'Marker', '.', 'MarkerSize', 20 ); hold on;
    xticks(0:2:6); xlim([0,7]); box off;
    xlabel('#BG odorants'); ylabel('mean decision time (ms)'); 
    title('Experiment');
end

