%use this
function params  = both_models_function(vars)
%Change-point condition data
% vars = shared_variables;
addpath(genpath("~/Dropbox (Brown)/sharedMatlabUtilities"));
suffix = '_50';
if vars.generate_data
    dataCP = learningAndPerception(vars);
    
    %Oddball condition datap
    dataOB = oddBall_Inference(vars);
else
%     "specify path to data from both conditions"
%     dataCP = importdata(vars.path_CP); % will replace this with path_CP
%     dataOB = importdata(vars.path_OB);
    dataCP = learningAndPerception(vars);
    
    %Oddball condition datap
    dataOB = oddBall_Inference(vars);


% dataPath='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined/2001_allBlockData.mat';
% %     load(dataPath);
% %     
% %     sel=alldata.block==3;
% %     dataCP = selBehav(alldata, sel);
% %     
% %     sel=alldata.block==4;
% %     dataOB = selBehav(alldata, sel);
% %     
%     
%     
%     dataCP1 = learningAndPerception(shared_variables(dataPath),1);
%     dataOB1 = oddBall_Inference(shared_variables(dataPath),1);
%     close all
%     
%     dataCP2 = learningAndPerception(shared_variables(dataPath),2);
%     dataOB2 = oddBall_Inference(shared_variables(dataPath),2);
%     close all
%     
%     dataCP1straight=straightStruct(dataCP1);
%     dataCP2straight=straightStruct(dataCP2);
%     dataOB1straight=straightStruct(dataOB1);
%     dataOB2straight=straightStruct(dataOB2);
%     
%     f=fieldnames(dataCP1straight);
%     for n=1:length(f)
%         dataCP.(f{n})=[dataCP1straight.(f{n});dataCP2straight.(f{n})];
%     end
%     f=fieldnames(dataOB1straight);
%     for n=1:length(f)
%         dataOB.(f{n})=[dataOB1straight.(f{n});dataOB2straight.(f{n})];
%     end
    
    
end

dataCP=straightStruct(dataCP);
dataOB=straightStruct(dataOB);


if vars.generate_data
    disp(8888)
    nTrials = length(dataCP.X);
    OBpredError = dataOB.predError';
    CPpredError = dataCP.predError';
    OBkl = dataOB.kl;
    CPkl = dataCP.kl;
    OBsurprise = dataOB.surprise';
    CPsurprise = dataCP.surprise';
    OBupdate = dataOB.cUpdate';
    CPupdate = dataCP.muUpdate';
    OBentropy = dataOB.entropy';
    CPentropy = dataCP.entropy';
    predictionErrorOnB = dataOB.predictionErrorOnB';
    predictionErrorOnX = dataCP.predictionErrorOnX';
    perceptualErrorOnX = dataCP.perceptualErrorOnX';
    perceptualErrorOnB = dataOB.perceptualErrorOnB';
    
else
    nTrials = length(dataCP.block)*.5;
    
%     OBpredError = nanmean(dataOB.predictErr(dataOB.block==4,:),2);
%     OBpredError(isnan(OBpredError))=0;
%     CPpredError = nanmean(dataCP.predictErr(dataCP.block==3,:),2);
%     CPpredError(isnan(CPpredError))=0;
%     
%     OBsurprise = dataOB.surprise(dataOB.block==4,:);
%     CPsurprise = dataCP.surprise(dataCP.block==3,:);
%     
%     OBupdate = nanmean(dataOB.predUpdate(dataOB.block==4,:),2);
%     OBupdate(isnan(OBupdate))=0;
%     CPupdate = nanmean(dataCP.predUpdate(dataOB.block==3,:),2);
%     CPupdate(isnan(CPupdate))=0;
%     
%     OBentropy = dataOB.entropy(dataOB.block==4,:);
%     CPentropy = dataCP.entropy(dataCP.block==3,:);
%     
%     perceptualErrorOnX = nanmean(dataCP.estErr(dataCP.block==3,:),2);
%     perceptualErrorOnX(isnan(perceptualErrorOnX))=0;
%     perceptualErrorOnB = nanmean(dataOB.estErr(dataOB.block==4,:),2);
%     perceptualErrorOnB(isnan(perceptualErrorOnB))=0;
    
    OBpredError = [dataOB.predictErr(dataOB.block==3,1);dataOB.predictErr(dataOB.block==3,2)];
    OBpredError(isnan(OBpredError))=0;
    CPpredError = [dataCP.predictErr(dataCP.block==4,1);dataCP.predictErr(dataCP.block==4,2)];
    CPpredError(isnan(CPpredError))=0;
    
    OBpredErrorObj = [dataOB.objPredErr(dataOB.block==3,1);dataOB.objPredErr(dataOB.block==3,2)];
    OBpredErrorObj(isnan(OBpredErrorObj))=0;
    CPpredErrorObj = [dataCP.objPredErr(dataCP.block==4,1);dataCP.objPredErr(dataCP.block==4,2)];
    CPpredErrorObj(isnan(CPpredErrorObj))=0;
    
    
    OBsurprise = dataOB.surprise;
    CPsurprise = dataCP.surprise;
    
    OBupdate = [dataOB.predUpdate(dataOB.block==3,1);dataOB.predUpdate(dataOB.block==3,2)];
    OBupdate(isnan(OBupdate))=0;
    CPupdate = [dataCP.predUpdate(dataOB.block==4,1);dataCP.predUpdate(dataOB.block==4,2)];
    CPupdate(isnan(CPupdate))=0;
    
    
    OBentropy = dataOB.entropy;
    CPentropy = dataCP.entropy;
    
    perceptualErrorOnX = [dataCP.estErr(dataCP.block==4,1);dataCP.estErr(dataCP.block==4,2)];
    perceptualErrorOnX(isnan(perceptualErrorOnX))=0;
    perceptualErrorOnB = [dataOB.estErr(dataOB.block==3,1);dataOB.estErr(dataOB.block==3,1)];
    perceptualErrorOnB(isnan(perceptualErrorOnB))=0;
    
    
end




if vars.average && vars.generate_data
    for i = 1:vars.average_number-1 
        dataCP = learningAndPerception(vars);
        dataOB = oddBall_Inference(vars);
        close all;
        OBpredError = OBpredError + dataOB.predError;
        CPpredError = CPpredError + dataCP.predError;
        OBkl = OBkl + dataOB.kl;
        CPkl = CPkl + dataCP.kl;
        OBsurprise = OBsurprise + dataOB.surprise;
        CPsurprise = CPsurprise + dataCP.surprise;
        OBupdate = OBupdate+  dataOB.cUpdate;
        CPupdate = CPupdate + dataCP.muUpdate;
        predictionErrorOnB = predictionErrorOnB + dataOB.predictionErrorOnB;
        predictionErrorOnX = predictionErrorOnX + dataCP.predictionErrorOnX;
        perceptualErrorOnX = perceptualErrorOnX + dataCP.perceptualErrorOnX;
        perceptualErrorOnB = perceptualErrorOnB + dataOB.perceptualErrorOnB;
        
    end
        OBpredError = OBpredError / vars.average_number;
        CPpredError = CPpredError / vars.average_number;
        OBkl = OBkl / vars.average_number;
        CPkl = CPkl / vars.average_number;
        OBsurprise = OBsurprise / vars.average_number;
        CPsurprise = CPsurprise / vars.average_number;
        OBupdate = OBupdate/ vars.average_number;
        CPupdate = CPupdate / vars.average_number;
        predictionErrorOnB = predictionErrorOnB / vars.average_number;
        predictionErrorOnX = predictionErrorOnX / vars.average_number;
        perceptualErrorOnX = perceptualErrorOnX / vars.average_number;
        perceptualErrorOnB = perceptualErrorOnB / vars.average_number;
end






close all;
%% Regression for the Update
medians = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_medians", suffix, '.mat'], ''));
medians_kl = medians.medians_kl;
params = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_params", suffix, '.mat'], ''));
params_kl = params.params_kl;

condNumCP=ones(length(CPsurprise),1);
condNumOB=ones(length(OBsurprise),1)*-1;

condNum=[condNumOB;condNumCP];

% xes = [ ones(4*nTrials,1), ...
%   ([CPpredError;OBpredError]),...
%   ([CPpredError;OBpredError] .* nanzscore([CPsurprise; OBsurprise].* condNum)),...
%   ([CPpredError;OBpredError] .* nanzscore([CPsurprise; OBsurprise])),...
%   ([CPpredError;OBpredError] .* condNum), ...
%   ([CPpredError;OBpredError].* nanzscore([CPentropy;OBentropy]))];

xes = [ ones(4*nTrials,1), ...
  ([OBpredError;CPpredError]),...
  ([OBpredError;CPpredError] .* nanzscore([OBsurprise; CPsurprise].* condNum)),...
  ([OBpredError;CPpredError] .* nanzscore([OBsurprise; CPsurprise])),...
  ([OBpredError;CPpredError] .* condNum), ...
  ([OBpredError;CPpredError].* nanzscore([OBentropy; CPentropy]))];


xes=zScoreX(xes);

data.Y=[perceptualErrorOnX; perceptualErrorOnB];
data.X=xes;
data.startPoint=[5,zeros(1,size(data.X,2)),0.05];
data.whichParams=logical([1,1,1,1,1,1,1,0]);
data.includeUniform=1;

[paramsUpdate, negLogLikeUpdate]=fitLinearModWCircErrs(data);

% xes = [ ones(2*nTrials,1), ...
%   ([CPpredError';OBpredError']),...
%   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'].* condNum)),...
%   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'])),...
%   ([CPpredError';OBpredError'] .* condNum), ([CPpredError';OBpredError'].* nanzscore([CPentropy';OBentropy']))];
%   
[b_update, bint_update] = regress([OBupdate;CPupdate], xes);
% change to other circular regression

%% Regression for Perceptual Error
medians = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_medians", suffix, '.mat'], ''));
medians_kl = medians.medians_kl;
params = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_params", suffix, '.mat'], ''));
params_kl = params.params_kl;

% xes = [ ones(4*nTrials,1), ...
%   ([CPpredErrorObj;OBpredErrorObj]),...
%   ([CPpredErrorObj;OBpredErrorObj] .* nanzscore([CPsurprise; OBsurprise].* condNum)),...
%   ([CPpredErrorObj;OBpredErrorObj] .* nanzscore([CPsurprise; OBsurprise])),...
%   ([CPpredErrorObj;OBpredErrorObj] .* condNum), ([CPpredErrorObj;OBpredErrorObj].* nanzscore([CPentropy;OBentropy]))];


xes = [ ones(4*nTrials,1), ...
  ([OBpredErrorObj;CPpredErrorObj]),...
  ([OBpredErrorObj;CPpredErrorObj] .* nanzscore([OBsurprise; CPsurprise].* condNum)),...
  ([OBpredErrorObj;CPpredErrorObj] .* nanzscore([OBsurprise; CPsurprise])),...
  ([OBpredErrorObj;CPpredErrorObj] .* condNum), ...
  ([OBpredErrorObj;CPpredErrorObj].* nanzscore([OBentropy; CPentropy]))];



xes=zScoreX(xes);

data.Y=[perceptualErrorOnX; perceptualErrorOnB];
data.X=xes;
data.startPoint=[5,zeros(1,size(data.X,2)),0.05];
data.whichParams=logical([1,1,1,1,1,1,1,0]);
data.includeUniform=1;

[paramsBias, negLogLikeBias]=fitLinearModWCircErrs(data);


% xes = [ ones(2*nTrials,1), ...
%   ([CPpredError';OBpredError']),...
%   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'].* condNum)),...
%   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'])),...
%   ([CPpredError';OBpredError'] .* condNum), ([CPpredError';OBpredError'].* nanzscore([CPentropy';OBentropy']))];

[b_perceptual, bint_perceptual] = regress([perceptualErrorOnB;perceptualErrorOnX], xes);


%% Saving coefficient values
params.bint_perceptual = bint_perceptual;
params.bint_update = bint_update;
params.b_perceptual = b_perceptual;
params.b_update = b_update;

params.circPerceptual=paramsBias;
params.circUpdate=paramsUpdate;
params.neglogLikePerceptual=negLogLikeBias;
params.neglogLikeUpdate=negLogLikeUpdate;

%% Coefficient plots
% Comparing paramters for Learning regression and Perceptual bias
% regression
% defaultPlotParameters;
% close all
% labels = {"Intercept"; "Prediction Error"; "PredError \times KL"; "PredError \times Surprise"};
% string = sprintf("H = "+ vars.H + '\nsigma_{mu} = ' + vars.sigma_c + '\nsigma_x = '...
%          + vars.sigma_x + '\nsigma_y = ' + vars.sigma_y + '\nTrials = ' +  vars.nTrials);
% fig = figure;
% 
% subplot(2, 1, 1);
% hold on
% 
% plot([.1:.1:.4 ; .1:.1:.4], [bint_update(:,1)'; bint_update(:,2)'], '-', 'Color', 'black')
% plot(.1:.1:.4, b_update, 'o', 'MarkerFaceColor', 'Cyan', 'MarkerEdgeColor', 'black')
% 
% % plot([.1:.1:.4 ; .1:.1:.4], [bint_update(:,1)'; bint_update(:,2)'], 's', 'markerSize', 8);
% plot([0, .6], [0, 0], '--')
% set(gca,'XTick',[.1:.1:.4],'xticklabel', labels);
% ylabel("Value");
% title("Regression on Learning");
% set(gca, 'box', 'off');
% 
% subplot(2, 1, 2);
% hold on
% plot([.1:.1:.4 ; .1:.1:.4], [bint_perceptual(:,1)'; bint_perceptual(:,2)'], '-', 'Color', 'black')
% plot(.1:.1:.4, b_perceptual, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'black')
% 
% % plot([.1:.1:.4 ; .1:.1:.4], [bint_perceptual(:,1)'; bint_perceptual(:,2)'], '.', 'markerSize', 8)
% plot([0, .6], [0, 0], '--')
% set(gca,'XTick',[.1:.1:.4],'xticklabel', labels);
% ylabel("Value");
% title("Regression on Perceptual Bias");
% annotation('textbox', [.7, .75, 0.1, 0.1], 'String', string, 'Interpreter','tex');
% set(gca, 'box', 'off');
% fig.Position = [10 800 1200 700];


%% Other plots

% corr([CPkl' CPsurprise'])
% corr([OBkl' OBsurprise'])
% corr([OBkl' OBsurprise'; CPkl' CPsurprise'])




% figure
% subplot(2, 1, 1)
% plot([-360,360], [-360,360], '--k');
% plot([-360,360], [0,0], '--k');
% plot([CPpredError, OBpredError], [CPupdate,OBupdate], 'or', 'markerSize', 6);
% ylabel('update');
% xlabel('prediction error');
% 
% 
% 
% figure;
% subplot(2,1,1);
% plot([CPpredError, OBpredError], [[CPpredError, OBpredError] .* [zscore(dataCP.relative_entropy_Mu), zscore(dataOB.relative_entropy_C)]], ".");
% xlabel("predError");
% ylabel("predError * KL");
% subplot(2,1,2);
% plot([CPpredError, OBpredError], [[CPpredError, OBpredError] .* [zscore(dataCP.shannon), zscore(dataOB.shannon)]], "." );
% xlabel("predError");
% ylabel("predError * shannon");
% 
% 
% figure
% hold on
% 
% plot([predictionErrorOnX, predictionErrorOnB], [perceptualErrorOnX, perceptualErrorOnB], 'o', 'markerSize', 6, 'markerFaceColor', 'b' );
% plot([predictionErrorOnX, predictionErrorOnB].* small, [perceptualErrorOnX, perceptualErrorOnB] .* small, 'o', 'markerSize', 6, 'markerFaceColor', 'b' );
% plot([predictionErrorOnX, predictionErrorOnB].* ~small, [perceptualErrorOnX, perceptualErrorOnB] .* ~small, 'o', 'markerSize', 6, 'markerFaceColor', 'r' );
% plot(0, 0,'o', 'markerSize', 6, 'markerFaceColor', 'b'  );
% plot([-2*pi, 2*pi], B(1) + B(2).*[-2*pi, 2*pi], '--k');
% plot([-2*pi, 2*pi], B(1) + B(3).*[-2*pi, 2*pi], '--k');
% 
% % plot([-2*pi, 2*pi], [0, 0], '--k');
% 
% xlabel('prediction error')
% ylabel('perceptual error')

end 






















%%


% function params  = both_models_function(vars)
% %Change-point condition data
% % vars = shared_variables;
% addpath(genpath("~/Dropbox (Brown)/sharedMatlabUtilities"));
% suffix = '_50';
% if vars.generate_data
%     dataCP = learningAndPerception(vars);
%     
%     %Oddball condition datap
%     dataOB = oddBall_Inference(vars);
% else
%     %     "specify path to data from both conditions"
%     %     dataCP = importdata(vars.path_CP); % will replace this with path_CP
%     %     dataOB = importdata(vars.path_OB);
%     dataPath='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined/2001_allBlockData.mat';
%     dataCP1 = learningAndPerception(shared_variables(dataPath),1);
%     dataOB1 = oddBall_Inference(shared_variables(dataPath),1);
%     close all
%     
%     dataCP2 = learningAndPerception(shared_variables(dataPath),2);
%     dataOB2 = oddBall_Inference(shared_variables(dataPath),2);
%     close all
%     
%     dataCP1straight=straightStruct(dataCP1);
%     dataCP2straight=straightStruct(dataCP2);
%     dataOB1straight=straightStruct(dataOB1);
%     dataOB2straight=straightStruct(dataOB2);
%     
%     f=fieldnames(dataCP1straight);
%     for n=1:length(f)
%         dataCP.(f{n})=[dataCP1straight.(f{n});dataCP2straight.(f{n})];
%     end
%     f=fieldnames(dataOB1straight);
%     for n=1:length(f)
%         dataOB.(f{n})=[dataOB1straight.(f{n});dataOB2straight.(f{n})];
%     end
%     
% end
% nTrials = length(dataCP.X);
% OBpredError = dataOB.predError;
% CPpredError = dataCP.predError;
% OBkl = dataOB.kl;
% CPkl = dataCP.kl;
% OBsurprise = dataOB.surprise;
% CPsurprise = dataCP.surprise;
% OBupdate = dataOB.cUpdate;
% CPupdate = dataCP.muUpdate;
% predictionErrorOnB = dataOB.predictionErrorOnB;
% predictionErrorOnX = dataCP.predictionErrorOnX;
% perceptualErrorOnX = dataCP.perceptualErrorOnX;
% perceptualErrorOnB = dataOB.perceptualErrorOnB;
% if vars.average && vars.generate_data
%     for i = 1:vars.average_number-1 
%         dataCP = learningAndPerception(vars);
%         dataOB = oddBall_Inference(vars);
%         close all;
%         OBpredError = OBpredError + dataOB.predError;
%         CPpredError = CPpredError + dataCP.predError;
%         OBkl = OBkl + dataOB.kl;
%         CPkl = CPkl + dataCP.kl;
%         OBsurprise = OBsurprise + dataOB.surprise;
%         CPsurprise = CPsurprise + dataCP.surprise;
%         OBupdate = OBupdate+  dataOB.cUpdate;
%         CPupdate = CPupdate + dataCP.muUpdate;
%         predictionErrorOnB = predictionErrorOnB + dataOB.predictionErrorOnB;
%         predictionErrorOnX = predictionErrorOnX + dataCP.predictionErrorOnX;
%         perceptualErrorOnX = perceptualErrorOnX + dataCP.perceptualErrorOnX;
%         perceptualErrorOnB = perceptualErrorOnB + dataOB.perceptualErrorOnB;
%     end
%         OBpredError = OBpredError / vars.average_number;
%         CPpredError = CPpredError / vars.average_number;
%         OBkl = OBkl / vars.average_number;
%         CPkl = CPkl / vars.average_number;
%         OBsurprise = OBsurprise / vars.average_number;
%         CPsurprise = CPsurprise / vars.average_number;
%         OBupdate = OBupdate/ vars.average_number;
%         CPupdate = CPupdate / vars.average_number;
%         predictionErrorOnB = predictionErrorOnB / vars.average_number;
%         predictionErrorOnX = predictionErrorOnX / vars.average_number;
%         perceptualErrorOnX = perceptualErrorOnX / vars.average_number;
%         perceptualErrorOnB = perceptualErrorOnB / vars.average_number;
% end
% close all;
% %% Regression for the Update
% 
% medians = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_medians", suffix, '.mat'], ''));
% medians_kl = medians.medians_kl;
% params = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_params", suffix, '.mat'], ''));
% params_kl = params.params_kl;
% 
% condNumCP=ones(length(CPsurprise),1);
% condNumOB=ones(length(OBsurprise),1)*-1;
% condNum=[condNumCP;condNumOB];
% 
% % xes = [ ones(2*nTrials,1), ...
% %   ([CPpredError;OBpredError]),...
% %   ([[CPpredError;OBpredError] .* nanzscore(interp1(medians_kl, params_kl, [CPkl; OBkl]))]),...
% %   ([([CPpredError;OBpredError] .* zscore(([CPsurprise; OBsurprise])))])];
% 
% xes = [ ones(2*nTrials,1), ...
%   ([CPpredError;OBpredError]),...
%   ([[CPpredError;OBpredError] .* zscore(([CPsurprise; OBsurprise])).*condNum]),...
%   ([([CPpredError;OBpredError] .* zscore(([CPsurprise; OBsurprise])))]),...
%   condNum,[zscore(([CPsurprise; OBsurprise])).*condNum],[zscore(([CPsurprise; OBsurprise]))]];
% 
% xes(isnan(xes(:,3)), 3) = max(xes(:, 3));
% 
% [b_update, bint_update] = regress([CPupdate;OBupdate], xes);
% 
% 
% %% Regression for Perceptual Error
% medians = load("/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/perceptual_regular_kl_medians_100");
% medians_kl = medians.medians_kl;
% params = load("/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/perceptual_regular_kl_params_100");
% params_kl = params.params_kl';
% 
% % xes = [ ones(2*nTrials,1), ...
% %       ([predictionErrorOnX; predictionErrorOnB]),...
% %       ([[predictionErrorOnX; predictionErrorOnB] .* nanzscore(interp1(medians_kl, params_kl, [CPkl; OBkl]))]),...
% %       ([([predictionErrorOnX; predictionErrorOnB] .* zscore([CPsurprise; OBsurprise]))])];
% 
%   xes = [ ones(2*nTrials,1), ...
%       ([predictionErrorOnX; predictionErrorOnB]),...
%       ([[predictionErrorOnX; predictionErrorOnB] .* zscore(([CPsurprise; OBsurprise])).*condNum]),...
%       ([([predictionErrorOnX; predictionErrorOnB] .* zscore([CPsurprise; OBsurprise]))]),...
%       condNum,[zscore(([CPsurprise; OBsurprise])).*condNum],[zscore(([CPsurprise; OBsurprise]))]];
% 
%   
% xes(isnan(xes(:,3)), 3) = interp1(medians_kl, params_kl, max(medians_kl));
% 
% 
% [b_perceptual, bint_perceptual] = regress([perceptualErrorOnX; perceptualErrorOnB], xes);
% 
% 
% %% Saving coefficient values
% params.bint_perceptual = bint_perceptual;
% params.bint_update = bint_update;
% params.b_perceptual = b_perceptual;
% params.b_update = b_update;
% 
% %% Coefficient plots
% % Comparing paramters for Learning regression and Perceptual bias
% % regression
% % defaultPlotParameters;
% % close all
% % labels = {"Intercept"; "Prediction Error"; "PredError \times KL"; "PredError \times Surprise";'condition'};
% % string = sprintf("H = "+ vars.H + '\nsigma_{mu} = ' + vars.sigma_c + '\nsigma_x = '...
% %          + vars.sigma_x + '\nsigma_y = ' + vars.sigma_y + '\nTrials = ' +  vars.nTrials);
% % fig = figure;
% % 
% % subplot(2, 1, 1);
% % hold on
% % 
% % plot([.1:.1:.5 ; .1:.1:.5], [bint_update(:,1)'; bint_update(:,2)'], '-', 'Color', 'black')
% % plot(.1:.1:.5, b_update, 'o', 'MarkerFaceColor', 'Cyan', 'MarkerEdgeColor', 'black')
% % 
% % % plot([.1:.1:.4 ; .1:.1:.4], [bint_update(:,1)'; bint_update(:,2)'], 's', 'markerSize', 8);
% % plot([0, .6], [0, 0], '--')
% % set(gca,'XTick',[.1:.1:.4],'xticklabel', labels);
% % ylabel("Value");
% % title("Regression on Learning");
% % set(gca, 'box', 'off');
% % 
% % subplot(2, 1, 2);
% % hold on
% % plot([.1:.1:.4 ; .1:.1:.4], [bint_perceptual(:,1)'; bint_perceptual(:,2)'], '-', 'Color', 'black')
% % plot(.1:.1:.4, b_perceptual, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'black')
% % 
% % % plot([.1:.1:.4 ; .1:.1:.4], [bint_perceptual(:,1)'; bint_perceptual(:,2)'], '.', 'markerSize', 8)
% % plot([0, .6], [0, 0], '--')
% % set(gca,'XTick',[.1:.1:.4],'xticklabel', labels);
% % ylabel("Value");
% % title("Regression on Perceptual Bias");
% % annotation('textbox', [.7, .75, 0.1, 0.1], 'String', string, 'Interpreter','tex');
% % set(gca, 'box', 'off');
% % fig.Position = [10 800 1200 700];
% 
% 
% %% Other plots
% 
% corr([CPkl' CPsurprise'])
% corr([OBkl' OBsurprise'])
% corr([OBkl' OBsurprise'; CPkl' CPsurprise'])
% 
% 
% 
% 
% % figure
% % subplot(2, 1, 1)
% % plot([-360,360], [-360,360], '--k');
% % plot([-360,360], [0,0], '--k');
% % plot([CPpredError, OBpredError], [CPupdate,OBupdate], 'or', 'markerSize', 6);
% % ylabel('update');
% % xlabel('prediction error');
% % 
% % 
% % 
% % figure;
% % subplot(2,1,1);
% % plot([CPpredError, OBpredError], [[CPpredError, OBpredError] .* [zscore(dataCP.relative_entropy_Mu), zscore(dataOB.relative_entropy_C)]], ".");
% % xlabel("predError");
% % ylabel("predError * KL");
% % subplot(2,1,2);
% % plot([CPpredError, OBpredError], [[CPpredError, OBpredError] .* [zscore(dataCP.shannon), zscore(dataOB.shannon)]], "." );
% % xlabel("predError");
% % ylabel("predError * shannon");
% % 
% % 
% % figure
% % hold on
% % 
% % plot([predictionErrorOnX, predictionErrorOnB], [perceptualErrorOnX, perceptualErrorOnB], 'o', 'markerSize', 6, 'markerFaceColor', 'b' );
% % plot([predictionErrorOnX, predictionErrorOnB].* small, [perceptualErrorOnX, perceptualErrorOnB] .* small, 'o', 'markerSize', 6, 'markerFaceColor', 'b' );
% % plot([predictionErrorOnX, predictionErrorOnB].* ~small, [perceptualErrorOnX, perceptualErrorOnB] .* ~small, 'o', 'markerSize', 6, 'markerFaceColor', 'r' );
% % plot(0, 0,'o', 'markerSize', 6, 'markerFaceColor', 'b'  );
% % plot([-2*pi, 2*pi], B(1) + B(2).*[-2*pi, 2*pi], '--k');
% % plot([-2*pi, 2*pi], B(1) + B(3).*[-2*pi, 2*pi], '--k');
% % 
% % % plot([-2*pi, 2*pi], [0, 0], '--k');
% % 
% % xlabel('prediction error')
% % ylabel('perceptual error')
% 
% end 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
