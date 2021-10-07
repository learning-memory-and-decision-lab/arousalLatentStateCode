

whichComp=2;

if whichComp==1
    basePath='/Users/ttli/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/ttli/Dropbox (Brown)/sharedMatlabUtilities/';
    modelDataDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined';
elseif whichComp==2
    basePath='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/';
elseif whichComp==4
    basePath='/Users/martin/Dropbox/arousalLearningPerception/vwm_task';
    toolPath='/Users/martin/Dropbox/sharedMatlabUtilities/';
elseif whichComp==5
    basePath='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/';
else
    basePath='/Users/CHANGE_HERE/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/CHANGE_HERE/Dropbox (Brown)/sharedMatlabUtilities/';
    
end


%%

subList=[2003,2004,2005];
subCondList=[1,2,2,2];

mastStructCP=[];
mastStructOB=[];

for i=1:length(subList)
    subNum=subList(i);
    subStr=num2str(subNum);
    
    subCond=subCondList(i);
    
    DP=['/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined/',subStr,'_allBlockData.mat'];
    %DP="/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/Meera1.mat";
    vars = shared_variables(DP);
    vars.H = .15;
    
    dataPath=DP;
    
    
    
   
    
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
        
        
        if subCond==1
            OBpredError = [dataOB.predictErr(dataOB.block==4,1);dataOB.predictErr(dataOB.block==4,2)];
            OBpredError(isnan(OBpredError))=0;
            CPpredError = [dataCP.predictErr(dataCP.block==3,1);dataCP.predictErr(dataCP.block==3,2)];
            CPpredError(isnan(CPpredError))=0;
            
            OBpredErrorObj = [dataOB.objPredErr(dataOB.block==4,1);dataOB.objPredErr(dataOB.block==4,2)];
            OBpredErrorObj(isnan(OBpredErrorObj))=0;
            CPpredErrorObj = [dataCP.objPredErr(dataCP.block==3,1);dataCP.objPredErr(dataCP.block==3,2)];
            CPpredErrorObj(isnan(CPpredErrorObj))=0;
            
            
            OBsurprise = dataOB.surprise;
            CPsurprise = dataCP.surprise;
            
            OBupdate = [dataOB.predUpdate(dataOB.block==4,1);dataOB.predUpdate(dataOB.block==4,2)];
            OBupdate(isnan(OBupdate))=0;
            CPupdate = [dataCP.predUpdate(dataOB.block==3,1);dataCP.predUpdate(dataOB.block==3,2)];
            CPupdate(isnan(CPupdate))=0;
            
            
            OBentropy = dataOB.entropy;
            CPentropy = dataCP.entropy;
            
            perceptualErrorOnX = [dataCP.estErr(dataCP.block==3,1);dataCP.estErr(dataCP.block==3,2)];
            perceptualErrorOnX(isnan(perceptualErrorOnX))=0;
            perceptualErrorOnB = [dataOB.estErr(dataOB.block==4,1);dataOB.estErr(dataOB.block==4,1)];
            perceptualErrorOnB(isnan(perceptualErrorOnB))=0;
            
        else
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
    
    dataOB.OBpredError=OBpredError;
    dataCP.CPpredError=CPpredError;      
           
    dataOB.OBpredErrorObj=OBpredErrorObj;
    dataCP.CPpredErrorObj=CPpredErrorObj;
    
    dataOB.OBsurprise=OBsurprise;
    dataCP.CPsurprise=CPsurprise;
       
    dataOB.OBupdate=OBupdate;
    dataCP.CPupdate=CPupdate;
    
    dataOB.OBentropy=OBentropy;
    dataCP.CPentropy=CPentropy;
       
    dataOB.perceptualErrorOnB=perceptualErrorOnB;        
    dataCP.perceptualErrorOnX=perceptualErrorOnX;    
   
    
    saveDir=[basePath, '/dataForModel'];
    
    fnCP=fullfile(saveDir,[subStr,'_allBlockDataModelCP.mat']);
    fnOB=fullfile(saveDir,[subStr,'_allBlockDataModelOB.mat']);

    
    save(fnCP,'dataCP')
    save(fnOB,'dataOB')
    
    close all;
    
    if i==1
        mastStructCP=dataCP;
        mastStructOB=dataOB;
    else
        
        mastStructCP=catBehav(dataCP,mastStructCP,0);
        mastStructOB=catBehav(dataOB,mastStructOB,0);
    end
end



%%


figure()
plot(mastStructOB.OBpredError,mastStructOB.OBupdate,'rx')
hold on
plot(mastStructCP.CPpredError,mastStructCP.CPupdate,'bx')
xlabel('subPredErr')
ylabel('update')
hold off

figure()
plot(mastStructOB.OBpredErrorObj,mastStructOB.perceptualErrorOnB,'ro')
hold on
plot(mastStructCP.CPpredErrorObj,mastStructCP.perceptualErrorOnX,'bo')
xlabel('objPredErr')
ylabel('perceptual error')
hold off

%%


OBpredError=mastStructOB.OBpredError;
CPpredError=mastStructCP.CPpredError;

OBpredErrorObj=mastStructOB.OBpredErrorObj;
CPpredErrorObj=mastStructCP.CPpredErrorObj;

OBsurprise=mastStructOB.OBsurprise;
CPsurprise=mastStructCP.CPsurprise;

OBentropy=mastStructOB.OBentropy;
CPentropy=mastStructCP.CPentropy;

OBupdate=mastStructOB.OBupdate;
CPupdate=mastStructCP.CPupdate;

perceptualErrorOnB=mastStructOB.perceptualErrorOnB;
perceptualErrorOnX=mastStructCP.perceptualErrorOnX;



%% Regression for the Update
medians = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_medians", suffix, '.mat'], ''));
medians_kl = medians.medians_kl;
params = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_params", suffix, '.mat'], ''));
params_kl = params.params_kl;

condNumCP=ones(length(CPsurprise),1);
condNumOB=ones(length(OBsurprise),1)*-1;

condNum=[condNumOB;condNumCP;];

% xes = [ ones(4*nTrials,1), ...
%   ([CPpredError;OBpredError]),...
%   ([CPpredError;OBpredError] .* nanzscore([CPsurprise; OBsurprise].* condNum)),...
%   ([CPpredError;OBpredError] .* nanzscore([CPsurprise; OBsurprise])),...
%   ([CPpredError;OBpredError] .* condNum), ...
%   ([CPpredError;OBpredError].* nanzscore([CPentropy;OBentropy]))];

xes = [ ones(4*length(subList)*nTrials,1), ...
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













