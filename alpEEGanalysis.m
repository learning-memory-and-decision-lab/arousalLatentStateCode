





clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     0) Set path and analysis choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Get an example channel to perform time frequency analysis on:

fileStruct=dir(fullfile(eegDataDir, '*.mat'))
fns={fileStruct.name};


% Loop through EEG data files:

% preallocate mat:

for sub=1:length(fns)
    fn=fns{sub};
    sn_str=fn(1:4);
    
    eegDat=load(fullfile(eegDataDir,fn));
    
    
    % EEG data lives in: eegDat.EEG.data
    channelLabels={eegDat.EEG.chanlocs(1:63).labels};
    sel=strcmp(channelLabels, exChannel);
    
    
    % 1) LOAD BEHAVIORAL DATA
    
    behaveFiles=dir(fullfile(behaveDataDir, ['vwm_multiStim_color_', sn_str, '*']));
    behavFNs={behaveFiles.name};
    
    for block=3:length(behavFNs)
        
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
        
        fn=fullfile(behaveDataDir,'subCombined',[sn_str,'_allBlockData.mat']);
        save(fn,'alldata');
    end
    
    dataPath=[modelDataDir,'/',sn_str,'_allBlockData.mat'];
    
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
    
    % 2) 
     % Find trials where EEG data was good:
    isGoodEEG=false(length(alldata.chosenTargEst), 1);
    ind_OBCPstart=find(eegDat.epochNumbers>25,1); %practice 5, random 20
    epoch_OBCP=eegDat.epochNumbers(ind_OBCPstart:end);
    epoch_OBCP=(epoch_OBCP-25)';
    isGoodEEG(epoch_OBCP)=true;
    goodData=selBehav(alldata, epoch_OBCP);
    
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
        downSampTimeMat(sub,:,:)=downSampTimes;
       
        kl=kl(epoch_OBCP,1);
        surprise=surprise(epoch_OBCP,1);
      
        x0=ones(length(goodData.colorArray), 1); 
        x1=zscore(kl); 
        x2=zscore(surprise); 
        % control for baseline
        x3=squeeze(nanmean(nanmean(eegDat.EEG.data(:,eegDat.EEG.times<0&eegDat.EEG.times>-500,ind_OBCPstart:end),1),2));
     
        X=[x0,x1,x2,x3];
        
        %preallocate matrix for each subject 
        mat(sub,:,:,:)=nan(length({eegDat.EEG.chanlocs.labels}), length(downSampTimes), size( X, 2));
        for chan=1:length({eegDat.EEG.chanlocs.labels})
            for t=1:length(downSampTimes)
                Y=squeeze(downSampData(chan,t,:));
                coef=regress(Y,X);
                mat(sub,chan,t,:) = coef;
            end
        end
    end
    
    
    time2plot=2730;
    
    time_ind=find(eegDat.EEG.times==time2plot,1);

    %[minval,minidx]=min(abs(eegDat.EEG.times-time2plot));
    
    data2plot=squeeze(nanmean(eegDat.EEG.data(:,time_ind,:)));
    
    figure
    topoplot(data2plot,eegDat.EEG.chanlocs)
    
    
    
    t1=find(eegDat.EEG.times>-300, 1,'first');
    t2=find(eegDat.EEG.times<5000, 1,'last');
    
    downSampData=eegDat.EEG.data(:, t1:t2,:);
    downSampTimes=eegDat.EEG.times(t1:t2);
    
    % Get relevant electrodes:
    chanLabs={eegDat.EEG.chanlocs.labels};
   
    % Select channels to look at -- for P1/N1 maybe O1/O2?
    cpz_ind=find(strcmp( 'CPz', chanLabs));
    cpz=squeeze(eegDat.EEG.data(cpz_ind,:,:));
    figure
    plot(eegDat.EEG.times,nanmean(cpz, 2), 'r');
    
    o1_ind=find(strcmp( 'O1', chanLabs));
    o2_ind=find(strcmp( 'O2', chanLabs));
    o1=squeeze(eegDat.EEG.data(o1_ind,:,:));
    o2=squeeze(eegDat.EEG.data(o2_ind,:,:));
    
    bn_o1=o1-nanmean(o1(eegDat.EEG.times<0&eegDat.EEG.times>-500,:));
    bn_o2=o2-nanmean(o2(eegDat.EEG.times<0&eegDat.EEG.times>-500,:));
    
    %o1_ERP=squeeze(nanmean(eegDat.EEG.data(o1_ind,:,:),3));
    
    pz_ind=find(strcmp( 'Pz', chanLabs));
    pz=squeeze(eegDat.EEG.data(pz_ind,:,:));
    
    bn_pz=pz-nanmean(pz(eegDat.EEG.times<0&eegDat.EEG.times>-500,:));
    
    figure
    plot(eegDat.EEG.times,nanmean(bn_o1, 2), 'r');
    
    figure
    plot(eegDat.EEG.times,nanmean(o2, 2), 'b');
    
    figure
    plot(eegDat.EEG.times,nanmean(pz, 2), 'm');
    
    figure
    imagesc(o1')
    
    pz_ind=find(strcmp( 'Pz', chanLabs));
    pz=squeeze(eegDat.EEG.data(pz_ind,:,:));
    
    bn_pz=pz-nanmean(pz(eegDat.EEG.times<0&eegDat.EEG.times>-500,:));

    
    chan_ind=find(strcmp( 'Pz', chanLabs));
    erp=nanmean(eegDat.EEG.data(chan_ind,:,:),3);
    
    plot(eegDat.EEG.times,nanmean(eegDat.EEG.data,3))
    
    
    
    
    % time frequency decomposition of EEG data to look at changes in power across different frequencies. 
    
    % parameters:
    tf_params.minTime=-200;
    tf_params.maxTime=4000;
    tf_params.sites  =1:eegDat.EEG.nbchan;
    tf_params.numcycles=4;
    tf_params.num_freqs=6; %6 frequencies
    tf_params.baseTime=[-100 0];
    tf_params.low=8;
    tf_params.high=12;%baseline before stim
    
    % Decompose EEG signal into the power in 6 frequencies
    tf_struct=TF_Script(eegDat.EEG, tf_params); %freq, time, trial, electrode
    
    %chanidx=find(strcmp( 'O1', chanLabs));
    
  
    %total power(average across trials and channels? or do i pick a channel?)
    trialAlphaPower=squeeze(tf_struct.trialPowerMap(4,:,:,:)); %8-12
    alphaPower=nanmean(trialAlphaPower,2);
    alphaPower=nanmean(alphaPower,3);
     
    %ERP power?
    
    %how to plot?
    figure
    imagesc(alphaPower') 
    colorbar
    
    
    

end
%% plot individual subject regression

for sub=1:length(fns)
    subplot_n=1;
    figure()
    hold on
    for coef=1:(size(X,2)-1)
        coefV=zscore(mat(sub,63,:,coef));
        %coefVchanMean=nanmean(coefV,2);
        downSampTimes=downSampTimeMat(sub,1,:);
        
        subplot(3,1,subplot_n)
        subplot_n=subplot_n+1;
        
        sterr =  std(squeeze(coefV)) / (size(coefV,3)^(.5));
        errorAll=zeros(size(coefV,3),1)+sterr;
        
        regressErrorPlot = shadedErrorBar(squeeze(downSampTimes),squeeze(coefV),errorAll);
        set(regressErrorPlot.mainLine,'HandleVisibility','off');
        set(regressErrorPlot.patch,'HandleVisibility','off');
        set(regressErrorPlot.edge,'HandleVisibility','off');
        
        
        
        plot(squeeze(downSampTimes),squeeze(coefV))
        xlim([-200 1000])
        xline(0,'--')
        xline(530,'--')
        %xline(1500,'--')
%         xline(2520,'--')
%         xline(2600,'--')
        yline(0,'--')
        if coef==1
            title('intercept')
        elseif coef==2
            title('KL')
        elseif coef==3
            title('surprise')
        elseif coef==4
            title('control:baseline')
        end
        
%         time2plot=400;
%     
%         time_ind=find(downSampTimes==time2plot,1);
%         
%         data2plot=mat(sub,:,time_ind,coef);
%         
%         
% %         subplot(3,2,subplot_n)
% %         subplot_n=subplot_n+1;
%         figure()
%         topoplot(data2plot,eegDat.EEG.chanlocs)
        
%         colorbar
%         if coef==1
%             title('intercept')
%         elseif coef==2
%             title('KL')
%         elseif coef==3
%             title('surprise')
%         elseif coef==4
%             title('control:baseline')
%         end
%         
%         coefS=num2str(coef);
%         vidObj=VideoWriter([sn_str,'_',coefS,'topoplot_stimsAverage.avi']);
%         %currFrame=getframe;
%         %writeVideo(vidObj,currFrame);
%         
%         open(vidObj);
%         
%         for i=1:length(downSampTimes)
%             
%             data2plot=mat(sub,:,i,coef);
%             
%             currPlot=topoplot(data2plot, eegDat.EEG.chanlocs);
%             
%             time=downSampTimes(i);
%             title(num2str(time))
%             
%             currFrame=getframe(gcf);
%             
%             writeVideo(vidObj,currFrame);
%             clf
%             delete(currPlot)
%         end
%         
%         close(vidObj)
% %     
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot heat maps for individual subject regression    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subN=2;

for k=1:size(mat,4)
%select the coefficient 

figure
imagesc(downSampTimes, [],squeeze(mat(subN,:,:,k)))
colorbar
set(gca, 'clim', [-5, 5])
if k==1
    title('Intercept')
elseif k==2
    title('First Cue Location')
elseif k==3
    title('Block Type')
elseif k==4
    title('Surprise')
elseif k==5
    title('Learning')
end

parRegMat=mat(:,:,:,k);

STE=nanstd(parRegMat)./sqrt(size(parRegMat, 1));
tStat= squeeze(nanmean(parRegMat)./STE);

figure
imagesc(downSampTimes, [],tStat)
colorbar
set(gca, 'clim', [-5, 5])
if k==1
    title('Intercept t')
elseif k==2
    title('First Cue Location t')
elseif k==3
    title('Block Type t')
elseif k==4
    title('Surprise t')
elseif k==5
    title('Learning t')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   cluster based permutation test   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

connectThresh=40; 
thresh=.01;
numPerm=1000;  


% GET CONNECTION MATRIX specifying which electrodes are connected to which:
if exist('connectionMat.mat') & ~recreateConnectionMat
    load connectionMat.mat % if we've already got one made, load it.
else
    % otherwise, create one from scratch:
    
    % Then get the XYZ coordinates for each channel:
    allLocas=[[eegDat.EEG.chanlocs.X]; [eegDat.EEG.chanlocs.Y]; [eegDat.EEG.chanlocs.Z]]' ;
    % Loop through the channels and get the distance between that channel and
    % all other channels.
    for i = 1:length(allLocas)
        relDist(:,i)=sqrt(sum((allLocas-repmat(allLocas(i,:), length(allLocas), 1)).^2, 2));
    end
    
    % Set a threshold on distance... and mark channels that fall within that
    % threshold:
    connectionMat=relDist<connectThresh;
    connectionMat=connectionMat-eye(length(connectionMat));
    
    % CHECK OUT CONNECTIONS:
    %imagesc(connectionMat);
    %close all
    
    % check average number of connections:
    %numCons=nanmean(sum(connectionMat));
    
    % Look at connectivity weights across scalp, alla AC:
    % figure;
    % for e = 1:64
    % subplot(8,8,e)
    % topoplot(connectionMat(e,:),EEG.chanlocs)
    % end
    % saveas(gcf, 'connectionMatFig.eps', 'epsc2')
    % close all
    save connectionMat.mat connectionMat
end


%%
for i=1:size(mat,4)
%     [h,p,ci,stats]=ttest(mat(:,:,:,i),[],'Alpha',0.01,'Tail','both');
%     h=squeeze(h);
%     imagesc(downSampTimes,[],h)
    
    clusterInfo = getEEG_clusterSize(mat(:,:,:,i),connectionMat, thresh, 'both');
    
    matCoef=mat(:,:,:,i);
    
    for j=1:numPerm
        v=nan(size(matCoef));
        %v=nan(size(mat,1),1);
        for k=1:size(v,1)
            r=rand;
            if r<0.5
                v(k,:,:)=-1;
            else
                v(k,:,:)=1;
            end
        end
        
        permMat=matCoef.*v;
        permClusterInfo = getEEG_clusterSize(permMat,connectionMat, thresh, 'both');   
        
        maxClusterSize(i,j)=max(permClusterInfo.clustWtMap(:));
       
    end
    
    sig=prctile(maxClusterSize(i,:), 95);
    clusterInfo.clustWtMap(clusterInfo.clustWtMap<sig)=0;
    clusterW=unique(clusterInfo.clustWtMap);
    
    figure
    histogram(maxClusterSize(i,:))
    hold on
    ylim([0 700])
    if length(clusterW)>1
        for k=2:length(clusterW)
            xline(clusterW(k),'--b')
        end
    end
    xline(sig,'--r')
    hold off
    figDirHist=sprintf('~/Dropbox (Brown)/arousalLearningPerception/figure_eeg/CBPT/_hist_%s',num2str(i));
    saveas(gcf,figDirHist,'epsc')
    
    figNum
    figure(figNum)
    imagesc(downSampTimes, 1:64, clusterInfo.clustWtMap)
    colorbar
    figDirColor=sprintf('~/Dropbox (Brown)/EI_balance_and_EEG/figure_eeg/CBPT/_color_%s',num2str(i));
    saveas(gcf,figDirColor,'epsc')

end
%%

chanNames={eegDat.EEG.chanlocs.labels}
clusterInfo.tMap(64,:)

pause on
for i =1:length(downSampTimes)
    disp(i)
    figure
    hold on
    topoplot(clusterInfo.tMap(:,downSampTimes==downSampTimes(i)), eegDat.EEG.chanlocs)
    hold off
    pause(0.5)
    clf
end




%%%%%%%%%%%
%% Movie %%
%%%%%%%%%%%

vidObj=VideoWriter('topoplots400-1400.avi');
%currFrame=getframe;
%writeVideo(vidObj,currFrame);

open(vidObj);

for i=1:length(downSampTimes)
    
    currPlot=topoplot(clusterInfo.tMap(:,downSampTimes==downSampTimes(i)), eegDat.EEG.chanlocs);
    
    currFrame=getframe;
    
    writeVideo(vidObj,currFrame);
    clf
    delete(currPlot)
end

close(vidObj)

%%








