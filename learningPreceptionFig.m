% makeFig1
addpath("~/Dropbox (Brown)/arousalLearningPerception");
vars = shared_variables;
dataSubject = both_models_function(vars);
vars.replace = false;
vars.nTrials = 500;
dataModel = both_models_function(vars);
close all;
num        = 10;
wid        = 17; % Totalw width
hts        = [6, 6]; % height of each row
cols       = {[.5 .5], [.5 .5]}; % width of columns
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, 5, [], [], ''); % etc
set(axs,'Units','normalized');
set(gcf,'color','white')
fig_.Position = [10 10 30 20];
fig_.Name = 'Learning and Perceptual bias regression coefficients';
% draw in each panel, one at a time


% labels = {"Intercept"; "Prediction Error"; "PredError \times KL"; "PredError \times Surprise"};
labels = {"\beta_0"; "\beta_1"; "\beta_2"; "\beta_3"};

string = sprintf("H = "+ vars.H + '\nsigma_{mu} = ' + vars.sigma_c + '\nsigma_x = '...
         + vars.sigma_x + '\nsigma_y = ' + vars.sigma_y + '\nTrials = ' +  vars.nTrials);


for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        plot([.1:.1:.4 ; .1:.1:.4], [dataModel.bint_update(:,1)'; dataModel.bint_update(:,2)'], '-', 'color', 'black')
        plot([.1:.1:.4 ; .1:.1:.4], [dataModel.bint_update(:,1)'; dataModel.bint_update(:,2)'], '.', 'markerSize', 8, 'color', 'black')
        plot([0, .6], [0, 0], ':')
        set(gca,'XTick',[.1:.1:.4],'xticklabel', labels, 'FontSize', 12);
        set(gca, 'box', 'off')
        
    elseif xx==2
        plot([.1:.1:.4 ; .1:.1:.4], [dataModel.bint_perceptual(:,1)'; dataModel.bint_perceptual(:,2)'], '-', 'color', 'black')
        plot([.1:.1:.4 ; .1:.1:.4], [dataModel.bint_perceptual(:,1)'; dataModel.bint_perceptual(:,2)'], '.', 'markerSize', 8, 'color', 'black')
        plot([0, .6], [0, 0], ':')
        set(gca,'XTick',[.1:.1:.4],'xticklabel', labels, 'FontSize', 12);
        
        set(gca, 'box', 'off')
        
        
        
        
        
    elseif xx==3
 
     
    plot([.1:.1:.4 ; .1:.1:.4], [dataSubject.bint_update(:,1)'; dataSubject.bint_update(:,2)'], '-', 'color', 'black')
    plot([.1:.1:.4 ; .1:.1:.4], [dataSubject.bint_update(:,1)'; dataSubject.bint_update(:,2)'], '.', 'markerSize', 8, 'color', 'black')
    plot([0, .6], [0, 0], ':')
    set(gca,'XTick',[.1:.1:.4],'xticklabel', labels, 'FontSize', 12);
     set(gca, 'box', 'off')
        
           
        
        
    elseif xx==4
      plot([.1:.1:.4 ; .1:.1:.4], [dataSubject.bint_perceptual(:,1)'; dataSubject.bint_perceptual(:,2)'], '-', 'color', 'black')
        plot([.1:.1:.4 ; .1:.1:.4], [dataSubject.bint_perceptual(:,1)'; dataSubject.bint_perceptual(:,2)'], '.', 'markerSize', 8, 'color', 'black')
        plot([0, .6], [0, 0], ':')
        set(gca,'XTick',[.1:.1:.4],'xticklabel', labels, 'FontSize', 12);
      set(gca, 'box', 'off')
        
          
    elseif xx==5
     set(gca, 'box', 'off')
        
           
        
    elseif xx==6
        
        hold on
    set(gca, 'box', 'off')
        
             
    elseif xx==7
        hold on
 set(gca, 'box', 'off')
        
               
    end
    
    
    setPLOT_panelLabel(gca, xx);
end

 kk=annotation('textbox')
 set(kk, 'string', 'Nassar et al 2009 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
 


saveas(gcf,  "~/Dropbox (Brown)/arousalLearningPerception/registered report/figures/learningAndPerBiasRegression.eps", 'eps')
% saveas(gcf,  'figure1new.eps', 'epsc2')
close(gcf)