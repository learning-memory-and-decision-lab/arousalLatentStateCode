
function result=eeg_analysis(subNum)

% addpath(genpath("~/Dropbox (Brown)/sharedMatlabUtilities"));
% adddpath(genpath("~/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis"));
% adddpath(genpath("~/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined"));



whichComp=1;
showPlots=true;

if whichComp==1 % iMac in lab manager office
    
    baseDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    modelDataDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined';
elseif whichComp==2 % lab Macbook
    baseDir='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    
elseif whichComp==3 % Meera's Macbook
    baseDir='/Users/meerasingh/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/meerasingh/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/meerasingh/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    
elseif whichComp==4 %  Kaitlyn's Macbook
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    
else % all other computers
    % change the local path to own username
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities';
    
end

% if whichComp>2
%     datapath=fullfile(baseDir,'EEG_data'); %eeg data folder
% else
%     datapath=fullfile(baseDir,'eeg_data');
% end

% Add entire set of subdirectories to the path
addpath(genpath(baseDir))
addpath(genpath(eeglabDir))
addpath(genpath(sharedFuncDir))

eegDataDir=fullfile(baseDir, 'eegData');

behaveDataDir=fullfile(baseDir, 'behaveData');

% OPTIONS FOR ANALYSIS:

exChannel='CPz';
sampNum=4;  % this number specifies how many EEG samples we use to get 1 datapoint for analysis
doRegress=1;

% Cheat sheet for subject data:

% presented color: what color were they trying to remember
% chosenTarget: which target were they asked to remember color of
% isFixateEEG: did we register a potential horizontal saccade and abort trial?
% HEOGdif:      continuous measure that we used to ID horizontal saccades.
% tarloc :      which side of screen is cued for the working memory task?

nBlock=4;
nStim=2;

stimOn=0;
stimOff=500; % ish (might be 548?)ALSO maskOn
maskOff=1500; % mask on for 1s
cue1on=2400; % retention interval for 900ms

cd(baseDir)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     1) Load EEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegDP=fullfile(eegDataDir,['/',subNum,'_FILT_STIM','.mat']);
eegDat = load(eegDP);


% EEG data lives in: eegDat.EEG.data
channelLabels={eegDat.EEG.chanlocs(1:63).labels};
sel=strcmp(channelLabels, exChannel);


% 1) LOAD BEHAVIORAL DATA

behaveFiles=dir(fullfile(behaveDataDir, ['vwm_multiStim_color_', subNum, '*']));
behavFNs={behaveFiles.name};
%keyboard
for block=1:length(behavFNs)
    
    load(fullfile(behaveDataDir, behavFNs{block}))
    %         fileName=sprintf('vwm_multiStim_color_%s%s.mat',sn_str,num2str(block));
    %         fileDir=fullfile(dataDir,fileName);
    %         load(fileDir)
    
    nTrial=size(data.reportedValue,1);
    
    estErr=nan(nTrial,nStim);
    est=nan(nTrial,nStim);
    
    predictErr=nan(nTrial,nStim);
    pred=nan(nTrial,nStim);
    predDeg=nan(nTrial,nStim);
    
    
    trueEst=nan(nTrial,nStim);
    truePred=nan(nTrial,nStim);
    
    predEstDiff=nan(nTrial,nStim);
    predUpdate=nan(nTrial,nStim);
    
    surpriseTrial=nan(nTrial,nStim);
    
    
    %        if ~data.prefs.doFeedback
    for r=1:nTrial
        for c=1:nStim
            %if data.chosenTarg(r,c)==1
            if data.chosenTargEst(r,c)==1
                estDeg(r,c)=data.reportedValue(r,1);
                est(r,c)=deg2rad(data.reportedValue(r,1));
                estErr(r,c)=circ_dist(deg2rad(data.reportedValue(r,1)),deg2rad(data.presentedColor(r,1)));
                trueEst(r,c)=deg2rad(data.presentedColor(r,1));
            else
                estDeg(r,c)=data.reportedValue(r,2);
                est(r,c)=deg2rad(data.reportedValue(r,2));
                estErr(r,c)=circ_dist(deg2rad(data.reportedValue(r,2)),deg2rad(data.presentedColor(r,2)));
                trueEst(r,c)=deg2rad(data.presentedColor(r,2));
            end
        end
    end
    %         else
    %             for r=1:nTrial
    %                 for c=1:nStim
    %
    %                 end
    %             end
    %end
    
    if data.prefs.doPredict
        predV=data.predictedValue;
        predVUpdated=cat(1,predV(2:end,:),[0 0]);
        data.chosenTargPredict=cat(1,[nan nan],data.chosenTargPredict);
        for r=2:nTrial
            for c=1:nStim
                %if data.chosenTarg(r,c)==1
                if data.chosenTargPredict(r,c)==1
                    predDeg(r,c)=data.predictedValue(r-1,1);
                    pred(r,c)=deg2rad(data.predictedValue(r-1,1));
                    predictErr(r,c)=circ_dist(deg2rad(data.predictedValue(r-1,1)),deg2rad(data.toPredict(r-1,1)));
                    truePred(r,c)=deg2rad(data.toPredict(r-1,1));
                    predEstDiff(r,c)=circ_dist(deg2rad(data.reportedValue(r-1,1)),deg2rad(data.predictedValue(r-1,1)));
                    predUpdate(r,c)=circ_dist(deg2rad(predVUpdated(r-1,1)),deg2rad(predV(r-1,1)));
                    surpriseTrial(r,c)=data.surprise(r-1,1);
                else
                    predDeg(r,c)=data.predictedValue(r-1,2);
                    pred(r,c)=deg2rad(data.predictedValue(r-1,2));
                    predictErr(r,c)=circ_dist(deg2rad(data.predictedValue(r-1,2)),deg2rad(data.toPredict(r-1,2)));
                    truePred(r,c)=deg2rad(data.toPredict(r-1,2));
                    predEstDiff(r,c)=circ_dist(deg2rad(data.reportedValue(r-1,2)),deg2rad(data.predictedValue(r-1,2)));
                    predUpdate(r,c)=circ_dist(deg2rad(predVUpdated(r-1,2)),deg2rad(predV(r-1,2)));
                    surpriseTrial(r,c)=data.surprise(r-1,2);
                end
            end
        end
    else
        predVUpdated=nan(nTrial,2);
        data.chosenTargPredict=nan(nTrial,2);
        predDeg=nan(nTrial,2);
        pred=nan(nTrial,2);
        predictErr=nan(nTrial,2);
        truePred=nan(nTrial,2);
        predEstDiff=nan(nTrial,2);
        predUpdate=nan(nTrial,2);
        surpriseTrial=nan(nTrial,2);
    end
    
    blockData=data;
    
    blockData = rmfield(blockData, 'prefs');
    blockData = rmfield(blockData, 'totScore');
    blockData = rmfield(blockData, 'initRT');
    if data.prefs.doPredict
        blockData = rmfield(blockData, 'predictedValue');
        blockData = rmfield(blockData, 'toPredict');
    end
    
    blockData.block=ones(length(blockData.colorArray),1).*block;
    blockData.trueEst=trueEst;
    blockData.truePred=truePred;
    blockData.est=est;
    blockData.estErr=estErr;
    blockData.pred=pred;
    blockData.predictErr=predictErr; %objective
    blockData.predEstDiff=predEstDiff; %subjective
    blockData.predUpdate=predUpdate;
    blockData.surpriseTrial=surpriseTrial;
    blockData.predDeg=predDeg;
    blockData.estDeg=estDeg;
    
    blockData=straightStruct(blockData);
    
    if ~exist('alldata')
        alldata=blockData;
    else
        alldata=catBehav(blockData, alldata, true);
    end
    
    %fn=fullfile(behaveDataDir,'subCombined',[subNum,'_allBlockData.mat']);
    %save(fn,'alldata');
end

%keyboard

practiceEndTrial=find(alldata.block==3,1)-1;

dataPath=[modelDataDir,'/',subNum,'_allBlockData.mat'];
%keyboard
load(dataPath);

%keyboard
%getting KL and surprise from behavioral data
dataCP1 = learningAndPerception(shared_variables(dataPath),1);
dataOB1 = oddBall_Inference(shared_variables(dataPath),1);
close all

dataCP2 = learningAndPerception(shared_variables(dataPath),2);
dataOB2 = oddBall_Inference(shared_variables(dataPath),2);
close all

klCP=cat(1,dataCP1.kl,dataCP2.kl);
klOB=cat(1,dataOB1.kl,dataOB2.kl);

surpriseCP=cat(1,dataCP1.surprise,dataCP2.surprise);
surpriseOB=cat(1,dataOB1.surprise,dataOB2.surprise);

kl = [mean(klCP,1),mean(klOB,1)]';
surprise = [mean(surpriseCP,1),mean(surpriseOB,1)]';

condCP=length(klCP);
condOB=length(klOB);


% 2)
% Find trials where EEG data was good:
isGoodEEG=false(length(alldata.chosenTargEst), 1);
ind_OBCPstart=find(eegDat.epochNumbers>practiceEndTrial,1); %practice 5, random 20
epoch_OBCP=eegDat.epochNumbers(ind_OBCPstart:end);

epoch_OBCP=(epoch_OBCP-practiceEndTrial)';
isGoodEEG(epoch_OBCP)=true;

alldata=rmfield(alldata, 'condition');
%keyboard
goodData=selBehav(alldata, epoch_OBCP);



condNumCP=sum(isGoodEEG(1:condCP));
condNumOB=sum(isGoodEEG(condCP+1:end));

cpTrial=ones(condNumCP,1);
obTrial=ones(condNumOB,1).*-1;


condNum=[cpTrial;obTrial];

%create a new variable that indicate whether first cued stim was target
%on the left or right
goodData.firstEstLoc=goodData.chosenTargEst(:,1);
goodData.surprise=nansum(goodData.surpriseTrial,2);
goodData.update=nansum(abs(goodData.predUpdate),2);



% Run a subject level regression to see whether any electrode at any
% timepont has trial-to-trial responses that relate to
if doRegress
    t1=find(eegDat.EEG.times>-200, 1,'first');
    t2=find(eegDat.EEG.times<4000, 1,'last');
    
    %downSampData=eegDat.EEG.data(:, 1:sampNum:end,:);
    downSampData=eegDat.EEG.data(:, t1:t2,ind_OBCPstart:end);
    %downSampTimes=eegDat.EEG.times(1:sampNum:end);
    downSampTimes=eegDat.EEG.times(t1:t2);
    
    
    kl=kl(epoch_OBCP,1);
    surprise=surprise(epoch_OBCP,1);
    
    x0=ones(length(goodData.colorArray), 1);
    x1=zscore(kl);
    x2=zscore(surprise);
    x3=zscore(surprise.*condNum);
    x4=condNum;
    % control for baseline
    x5=squeeze(nanmean(nanmean(eegDat.EEG.data(:,eegDat.EEG.times<0&eegDat.EEG.times>-500,ind_OBCPstart:end),1),2));
    
    X=[x0,x1,x2,x3,x4,x5];
    
    %preallocate matrix for each subject
    mat(:,:,:)=nan(length({eegDat.EEG.chanlocs.labels}), length(downSampTimes), size( X, 2));
    for chan=1:length({eegDat.EEG.chanlocs.labels})
        for t=1:length(downSampTimes)
            Y=squeeze(downSampData(chan,t,:));
            coef=regress(Y,X);
            mat(chan,t,:) = coef;
        end
    end
end



%% save results useful for plot 

result.downSampTimes=downSampTimes;
result.mat=mat;


end
