
whichComp=1;


showPlots=true;

if whichComp==1 % iMac in lab manager office
    baseDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    sharedFuncDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities';
elseif whichComp==2 % lab Macbook
    baseDir='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    sharedFuncDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities';
elseif whichComp==3% all other computers
    % change the local path to own username
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities';
elseif whichComp==4
    % change the local path to own username
    baseDir='~/Dropbox/arousalLearningPerception/vwm_task';
    sharedFuncDir='~/Dropbox/sharedMatlabUtilities';
    
end

addpath(genpath(baseDir))
addpath(genpath(sharedFuncDir))


fileName=sprintf('behaveData/subCombined/%s_allBlockData.mat','2003');
    
DP = fullfile(baseDir, fileName);


vars = shared_variables(DP);
vars.replace = false;
vars.nTrials = 500;
vars.H = .1;
dataCP = learningAndPerception(vars);
dataOB = oddBall_Inference(vars);
close all;
num        = 1;
wid        = 17; % total width
hts        = [7, 2.5, 2.5, 2.5 2.5]; % height of each row
cols       = {[.50 .50], [.50 .50], [.50 .50], [.50 .50],[.50 .50]}; % width of columns
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1.5, 2, [], ''); % etc
set(axs,'Units','normalized');
set(gcf,'color','white')
fig_.Position = [10 10 30 40];

% draw in each panel, one at a time


lw=1.5
lw2=3
exSub=12





for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        imageOld = (imread("~/Desktop/neuro/generative_changepoint.jpg"));
        image(flip(rot90(imageOld, 2), 2));        
        set(gca, 'visible', 'off')
        
    elseif xx==2
        imageOld = (imread("~/Desktop/neuro/generative_oddball.jpg"));
        image(flip(rot90(imageOld, 2), 2));
        set(gca, 'visible', 'off')
        
        
        
    elseif xx==3
    plot(1:dataCP.nTrials, dataCP.X, '--', 'color', 'black');
        plot(1:dataCP.nTrials, dataCP.X, 'o', 'markerSize', 8, 'color', 'black');
        ylim([0, 360]);
        xlim([0, dataCP.nTrials]);
        set(gca,'XTick',[]);
%         set(gca, 'box', 'off')
    plot(dataCP.maxLikePostMu, 'o', 'color', 'r');   
     
     
        
         
        
        
    elseif xx==4
      plot(1:dataOB.nTrials, dataOB.B, '--', 'color', 'black');
        plot(1:dataOB.nTrials, dataOB.B, 'o', 'markerSize', 8, 'color', 'black');
        ylim([0, 360]);
        xlim([0, dataCP.nTrials]);
        set(gca,'XTick',[]);
        set(gca, 'box', 'off')
      plot(dataOB.maxLikePostB, 'o', 'color', 'r'); 
     
        
          
    elseif xx==5
      
        
      imagesc(dataCP.posterior')
     ylim([0, 360]);
     xlim([0, dataCP.nTrials]);
     set(gca, 'box', 'off')
     
    elseif xx == 6
        imagesc(dataOB.posterior')
     ylim([0, 360]);
     xlim([0, dataOB.nTrials]);
     set(gca, 'box', 'off')
    
    elseif xx==7
        
     imagesc(dataCP.posteriorPercept')
     ylim([0, 360]);
     xlim([0, dataCP.nTrials]);
     set(gca, 'box', 'off')  
        
    elseif xx==8
        
     imagesc(dataOB.posteriorPercept')
     ylim([0, 360]);
     xlim([0, dataCP.nTrials]);
     set(gca, 'box', 'off') 
        
    elseif xx==9
      kl = dataCP.kl / max(dataCP.kl);
      maxSurprise = max(dataCP.surprise);
      plot(kl, '.');
      plot(dataCP.surprise, '.');
      legend('KL divergence', 'Surprise', 'location', 'northeast');
      set(gca, 'box', 'off')
        hold on
    set(gca, 'box', 'off')
        
             
    elseif xx==10
        kl = dataOB.kl / max(dataOB.kl);
        plot(kl, '.');
      plot(dataOB.surprise, '.');
      legend('KL divergence', 'Surprise', 'location', 'northeast');
     set(gca, 'box', 'off')
        hold on
 set(gca, 'box', 'off')
        
               
    end
    
    
    setPLOT_panelLabel(gca, xx);
end

 kk=annotation('textbox')
 set(kk, 'string', 'Nassar et al 2009 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
 


saveas(gcf,  "~/Dropbox (Brown)/arousalLearningPerception/registered report/figures/.fig", 'fig')
% saveas(gcf,  'figure1new.eps', 'epsc2')
close(gcf)