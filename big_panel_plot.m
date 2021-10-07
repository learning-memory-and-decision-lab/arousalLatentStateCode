
%% just hit "Run and Advance" for this cell 

clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                set path (directory) and params                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%$$$$ pay attention to the path changes

% change whichComp accordingly!!
whichComp=1;


showPlots=true;

if whichComp==1 % iMac in lab manager office
    baseDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    sharedFuncDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities';
elseif whichComp==2 % lab Macbook 
    baseDir='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    sharedFuncDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities';
else % all other computers 
    % change the local path to own username 
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities';
end

addpath(genpath(baseDir))
addpath(genpath(sharedFuncDir))  


%% Big panel plot learning behavior
%clear;
defaultPlotParameters;
vars = shared_variables();

vars.replace = false;
vars.nTrials = 30;
vars.H = .15;
dataCP = learningAndPerception(vars);
dataOB = oddBall_Inference(vars);
close all;

cp_kl = zscore(dataCP.kl');
cp_surprise = zscore(dataCP.surprise');
ob_kl = zscore(dataOB.kl');
ob_surprise = zscore(dataOB.surprise');
both_surprise = zscore([dataCP.surprise, dataOB.surprise]');
both_kl = zscore([ dataCP.kl dataOB.kl]');


b_cp = regress( cp_surprise, [ones(size(dataCP.kl')), cp_kl]);
b_ob = regress( ob_surprise, [ones(size(dataOB.kl')), ob_kl]);
b_both = regress(both_surprise, [ones(size(both_kl)), both_kl]);
max_val = 5;
min_val = -2;

vars.average = false; 

%both = both_models_function(vars);


[lr, up, pe] = difference(deg2rad(dataCP.X), deg2rad(dataCP.maxLikePostMu(1:end-1)), 'polar');


new_trials = zeros(size(dataCP.X));
new_trials(1) = 1;
% [lr, up, pe] = computeLR(dataCP.X, dataCP.maxLikePostMu(1:end-1),new_trials , 'polar');

num        = 1;
wid        = 17; % total width
hts        = [8, 6.5, 6.5 ]; % height of each row
cols       = {[.50 .50], [.50 .50], [.50 .50]}; % width of columns
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1.5, 2.5, [], ''); % etc
set(axs,'Units','normalized');
set(gcf,'color','white')
fig_.Position = [10 10 30 40];
left_color = [0 0 255] / 256;
right_color = [255 0 0] / 256;
set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% draw in each panel, one at a time


lw=1.5
lw2=3
exSub=12

learningCP=dataCP.surprise.*ones(1,length(dataCP.surprise));
learningOB=dataOB.surprise.*(ones(1,length(dataOB.surprise))*-1);


for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
     
        set(gca, 'visible', 'off')
        
    elseif xx==2

        set(gca, 'visible', 'off')
        
        
        
    elseif xx==3
        plot(1:dataCP.nTrials, dataCP.X, 'o', 'markerSize', 8, 'color', 'black');
        plot(1:dataCP.nTrials, dataCP.maxLikePostX, 'x', 'markerSize', 10, 'color', 'blue');
        ylim([0, 360]);
         xlim([5 30])
        ylabel('Color Value');
        xlabel('Trial');
        set(gca, 'box', 'off')
    plot(dataCP.maxLikePostMu, '-', 'color', 'r');   
     
     
        
         
        
        
    elseif xx==4
        plot(1:dataOB.nTrials, dataOB.B, 'o', 'markerSize', 8, 'color', 'black');
        plot(1:dataCP.nTrials, dataOB.maxLikePostB, 'x', 'markerSize', 10, 'color', 'blue');
        ylim([0, 360]);
         xlim([5 30])
        ylabel('Color Value');
        xlabel('Trial');
        set(gca, 'box', 'off')
      plot(dataOB.maxLikePostC, '-', 'color', 'r'); 
     
        
        
    elseif xx==5
%       kl = dataCP.kl / max(dataCP.kl);
      kl = dataCP.kl;
      
      
      
      maxSurprise = max(dataCP.surprise);
      yyaxis left;
      ylim([0 1])
      %ylabel('KL Divergence');
      plot(dataCP.muUpdate./dataCP.predError, '-', 'markerSize', 8, 'color', 'blue');
      plot(dataCP.surprise, '-', 'markerSize', 8, 'color', 'red');
%       legend('KL divergence', 'Surprise', 'location', 'northeast');
      plot(learningCP, '-', 'markerSize', 8, 'color', 'green');
      
      yyaxis right
      ylim([4 6.5])
      plot(dataCP.entropy,'-', 'color', 'm')
      
set(gca, 'box', 'off')
        hold on
    set(gca, 'box', 'off')
    xlabel('Trial');
    %ylabel('Surprise');
    xlim([5 30])
        
             
    elseif xx==6
%         kl = dataOB.kl / max(dataOB.kl);

        kl = dataOB.kl;
        yyaxis left;
        ylim([0 1])
        %ylabel('KL Divergence');
        plot(dataOB.cUpdate./dataOB.predError, '-', 'markerSize', 8, 'color', 'blue');
        
      plot(dataOB.surprise, '-', 'markerSize', 8, 'color', 'red');
%       legend('KL divergence', 'Surprise', 'location', 'northeast');
      plot(learningOB, '-', 'markerSize', 8, 'color', 'm  ');
    yyaxis right 
    ylim([4 6.5])
    plot(dataOB.entropy,'-', 'color', 'green')
set(gca, 'box', 'off')
        hold on
 set(gca, 'box', 'off')
 xlabel('Trial');
 xlim([5 30])
    %ylabel('Surprise');

    
        
               
    end
    
    
    setPLOT_panelLabel(gca, xx);
end

 kk=annotation('textbox')
 set(kk, 'string', 'Nassar et al 2009 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
 


saveas(gcf,  "~/Dropbox (Brown)/arousalLearningPerception/registered report/figures/model", 'eps')
% saveas(gcf,  'figure1new.eps', 'epsc2')
%close(gcf)



%% Big panel plot
clear;
defaultPlotParameters;
vars = shared_variables();
vars.replace = false;
vars.nTrials = 30;
vars.H = .15;
dataCP = learningAndPerception(vars);
dataOB = oddBall_Inference(vars);
close all;

[lr, up, pe] = difference(deg2rad(dataCP.X), deg2rad(dataCP.maxLikePostMu(1:end-1)), 'polar');


new_trials = zeros(size(dataCP.X));
new_trials(1) = 1;
% [lr, up, pe] = computeLR(dataCP.X, dataCP.maxLikePostMu(1:end-1),new_trials , 'polar');

num        = 1;
wid        = 17; % total width
hts        = [8, 6.5, 6.5 ]; % height of each row
cols       = {[.50 .50], [.50 .50], [.50 .50]}; % width of columns
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1.5, 2.5, [], ''); % etc
set(axs,'Units','normalized');
set(gcf,'color','white')
fig_.Position = [10 10 30 40];
left_color = [0 0 255] / 256;
right_color = [255 0 0] / 256;
set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% draw in each panel, one at a time


lw=1.5
lw2=3
exSub=12





for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
     
        set(gca, 'visible', 'off')
        
    elseif xx==2

        set(gca, 'visible', 'off')
        
        
        
    elseif xx==3
        plot(1:dataCP.nTrials, dataCP.X, 'o', 'markerSize', 8, 'color', 'black');
        plot(1:dataCP.nTrials, dataCP.maxLikePostX, 'x', 'markerSize', 10, 'color', 'blue');
        ylim([0, 360]);
        xlim([0, dataCP.nTrials]);
        ylabel('Color Value');
        xlabel('Trial');
        set(gca, 'box', 'off')
    plot(dataCP.maxLikePostMu, '-', 'color', 'r');   
     
     
        
         
        
        
    elseif xx==4
        plot(1:dataOB.nTrials, dataOB.B, 'o', 'markerSize', 8, 'color', 'black');
        plot(1:dataCP.nTrials, dataOB.maxLikePostB, 'x', 'markerSize', 10, 'color', 'blue');
        ylim([0, 360]);
        xlim([0, dataCP.nTrials]);
        ylabel('Color Value');
        xlabel('Trial');
        set(gca, 'box', 'off')
      plot(dataOB.maxLikePostC, '-', 'color', 'r'); 
     
        
        
    elseif xx==5
%       kl = dataCP.kl / max(dataCP.kl);
      kl = dataCP.kl;
      maxSurprise = max(dataCP.surprise);
      %yyaxis left;
      %ylabel('KL Divergence');
      plot(kl, '-', 'markerSize', 8, 'color', 'blue');
      
      %yyaxis right
      plot(dataCP.surprise, '-', 'markerSize', 8, 'color', 'red');
%       legend('KL divergence', 'Surprise', 'location', 'northeast');
      set(gca, 'box', 'off')
        hold on
    set(gca, 'box', 'off')
    xlabel('Trial');
    
    %ylabel('Surprise');
        
             
    elseif xx==6
%         kl = dataOB.kl / max(dataOB.kl);

        kl = dataOB.kl;
        %yyaxis left;
        %ylabel('KL Divergence');
        plot(kl, '-', 'markerSize', 8, 'color', 'blue');
        %yyaxis right
      plot(dataOB.surprise, '-', 'markerSize', 8, 'color', 'red');
%       legend('KL divergence', 'Surprise', 'location', 'northeast');
     set(gca, 'box', 'off')
        hold on
 set(gca, 'box', 'off')
 xlabel('Trial');
    %ylabel('Surprise');
        
               
    end
    
    
    setPLOT_panelLabel(gca, xx);
end

 kk=annotation('textbox')
 set(kk, 'string', 'Nassar et al 2009 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
 


saveas(gcf,  "~/Dropbox (Brown)/arousalLearningPerception/registered report/figures/model", 'eps')
% saveas(gcf,  'figure1new.eps', 'epsc2')
%close(gcf)


%% KL and Surprise Correlation
%clear;
vars = shared_variables();
vars.nTrials = 90;
vars.H = .15;
dataCP = learningAndPerception(vars);
dataOB = oddBall_Inference(vars);


cp_kl = zscore(dataCP.kl');
cp_surprise = zscore(dataCP.surprise');
ob_kl = zscore(dataOB.kl');
ob_surprise = zscore(dataOB.surprise');
both_surprise = zscore([dataCP.surprise, dataOB.surprise]');
both_kl = zscore([ dataCP.kl dataOB.kl]');


b_cp = regress( cp_surprise, [ones(size(dataCP.kl')), cp_kl]);
b_ob = regress( ob_surprise, [ones(size(dataOB.kl')), ob_kl]);
b_both = regress(both_surprise, [ones(size(both_kl)), both_kl]);
max_val = 5;
min_val = -2;

vars.average = false; 
vars.replace = false;
vars.nTrials = 300;
both = both_models_function(vars);
close all;
% figure;
% hold on;
% plot(cp_kl,cp_surprise, 'o', 'color', 'r');
% plot(ob_kl, ob_surprise, 'o', 'color', 'blue');
% plot([min_val, max_val], b_cp(1) + b_cp(2)*[min_val, max_val], '-', 'color', 'r');
% plot([min_val, max_val], b_ob(1) + b_ob(2)*[min_val, max_val], '-', 'color', 'blue');
% plot([min_val, max_val], b_both(1) + b_both(2)*[min_val, max_val], '-', 'color', 'black');
% legend('Changepoint value', 'Oddball values', 'Changepoint line of best fit', 'Oddball line of best fit'...
%         ,'Overall line of best fit', 'Location', 'northwest');
%     
% legend boxoff
% xlabel('KL divergence');
% ylabel('Surprise');
% set(gca, 'box', 'off');


num        = 1;
wid        = 17; % total width
hts        = [5.5, 5.5, 5.5, 5.5 ]; % height of each row
cols       = {[.5, .5], [.5 .5],  1, 1}; % width of columns
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 1.5, [], ''); % etc
set(axs,'Units','normalized');
set(gcf,'color','white')
fig_.Position = [10 10 30 40];
left_color = [0 0 255] / 256;
right_color = [255 0 0] / 256;
set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% draw in each panel, one at a time


lw=1.5
lw2=3
exSub=12

threshold = 0.5;

labels = {'\beta_0';'\beta_1';'\beta_2';'\beta_3';'\beta_4';'\beta_5'};
string = sprintf("H = "+ vars.H + '\nsigma_{mu} = ' + vars.sigma_c + '\nsigma_x = '...
         + vars.sigma_x + '\nsigma_y = ' + vars.sigma_y + '\nTrials = ' +  vars.nTrials);

     

     
for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
%         small_error = abs(dataCP.predError) < threshold;
%         plot([-pi, pi], [-pi, pi], '--k');
%         plot([-pi, pi], [0, 0], '--k');
%         small_b = regress(dataCP.muUpdate(small_error)',[ones(size(dataCP.predError(small_error)))' dataCP.predError(small_error)']);
%         large_b = regress(dataCP.muUpdate(~small_error)',[ones(size(dataCP.predError(~small_error)))' dataCP.predError(~small_error)']);
%         plot(dataCP.predError(small_error) , dataCP.muUpdate(small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'blue');
%         plot(dataCP.predError(~small_error) , dataCP.muUpdate(~small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'red');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         xlabel('Prediction error')
%         set(gca, 'box', 'off')
    elseif xx == 2
%         small_error = abs(dataOB.predError) < threshold;
%         plot([-pi, pi], [-pi, pi], '--k');
%         plot([-pi, pi], [0, 0], '--k');
%         small_b = regress(dataOB.cUpdate(small_error)',[ones(size(dataOB.predError(small_error)))' dataOB.predError(small_error)']);
%         large_b = regress(dataOB.cUpdate(~small_error)',[ones(size(dataOB.predError(~small_error)))' dataOB.predError(~small_error)']);
%         plot(dataOB.predError(small_error) , dataOB.cUpdate(small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'blue');
%         plot(dataOB.predError(~small_error) , dataOB.cUpdate(~small_error), 'o', 'markerSize', 8,  'markerEdgeColor', 'red');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         xlabel('Prediction error')
%         set(gca, 'box', 'off')
    elseif xx == 3
%         small_error = abs(dataCP.predictionErrorOnX) < threshold;
%         small_b = regress(dataCP.perceptualErrorOnX(small_error)',[ones(size(dataCP.predictionErrorOnX(small_error)))' dataCP.predictionErrorOnX(small_error)']);
%         large_b = regress(dataCP.perceptualErrorOnX(~small_error)',[ones(size(dataCP.predictionErrorOnX(~small_error)))' dataCP.predictionErrorOnX(~small_error)']);
%         plot(dataCP.predictionErrorOnX(small_error), dataCP.perceptualErrorOnX(small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'blue');
%         plot(dataCP.predictionErrorOnX(~small_error), dataCP.perceptualErrorOnX(~small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'red');
%         plot([-pi, pi], [-pi, pi], '--k');
%         plot([-pi, pi], [0, 0], '--k');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         set(gca, 'box', 'off')
%         ylabel('Perceptual Error')
%         xlabel('Prediction Error')
    elseif xx == 4
%         small_error = abs(dataOB.predictionErrorOnB) < threshold;
%         small_b = regress(dataOB.perceptualErrorOnB(small_error)',[ones(size(dataOB.predictionErrorOnB(small_error)))' dataOB.predictionErrorOnB(small_error)']);
%         large_b = regress(dataOB.perceptualErrorOnB(~small_error)',[ones(size(dataOB.predictionErrorOnB(~small_error)))' dataOB.predictionErrorOnB(~small_error)']);
%         plot(dataOB.predictionErrorOnB(small_error), dataOB.perceptualErrorOnB(small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'blue');
%         plot(dataOB.predictionErrorOnB(~small_error), dataOB.perceptualErrorOnB(~small_error) , 'o', 'markerSize', 8, 'markerEdgeColor', 'red');
%         plot([-pi, pi], [-pi, pi], '--k');
%         plot([-pi, pi], [0, 0], '--k');
%         plot([-pi, pi], small_b(1) + small_b(2)*[-pi, pi], '-', 'color', 'blue');
%         plot([-pi, pi], large_b(1) + large_b(2)*[-pi, pi], '-', 'color', 'red');
%         ylabel('Update')
%         set(gca, 'box', 'off')
%         ylabel('Perceptual Bias')
%         xlabel('Prediction error')
    elseif xx == 5
        
        plot([.1:.1:.6 ; .1:.1:.6], [both.bint_update(:,1)'; both.bint_update(:,2)'], '-', 'Color', 'black')
        plot(.1:.1:.6, both.b_update, 'o', 'MarkerFaceColor', 'Cyan', 'MarkerEdgeColor', 'black')

        % plot([.1:.1:.4 ; .1:.1:.4], [bint_update(:,1)'; bint_update(:,2)'], 's', 'markerSize', 8);
        plot([0, .6], [0, 0], '--')
        ylim([-0.3 0.9])
        xlim([0 0.45])
        set(gca,'XTick',[.1:.1:.6], 'xticklabel', labels);
        ylabel("Value");
        set(gca, 'box', 'off');
    elseif xx == 6 
        hold on
        plot([.1:.1:.6 ; .1:.1:.6], [both.bint_perceptual(:,1)'; both.bint_perceptual(:,2)'], '-', 'Color', 'black')
        plot(.1:.1:.6, both.b_perceptual, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'black')

        % plot([.1:.1:.4 ; .1:.1:.4], [bint_perceptual(:,1)'; bint_perceptual(:,2)'], '.', 'markerSize', 8)
        plot([0, .6], [0, 0], '--')
        ylim([-0.2 0.6])
        xlim([0 0.45])
        set(gca,'XTick',[.1:.1:.6], 'xticklabel', labels);
        ylabel("Value");
        set(gca, 'box', 'off');
    end
    setPLOT_panelLabel(gca, xx);

end



%% Eye Regression
close all;
clear;
result = eye_analysis_martin('~/Dropbox (Brown)/arousalLearningPerception/vwm_task/ET_data/Meera022620.mat',200, 4000, 150);

eye_analysis_martin('~/Dropbox (Brown)/arousalLearningPerception/vwm_task/ET_data/2001.mat',200, 4000, 150);

%eye_analysis_martin('~/Dropbox (Brown)/arousalLearningPerception/vwm_task/ET_data/Meera022620.mat',200, 4000, 150);

%% 