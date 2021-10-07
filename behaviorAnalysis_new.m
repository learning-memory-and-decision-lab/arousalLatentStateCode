
clear
close all

%%

whichComp=5;

if whichComp==1
    basePath='/Users/ttli/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/ttli/Dropbox (Brown)/sharedMatlabUtilities/';
    modelDataDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData/subCombined';
elseif whichComp==2
    basePath='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/';
elseif whichComp==4
    basePath='~/Dropbox/arousalLearningPerception/vwm_task';
    toolPath='~/Dropbox/sharedMatlabUtilities/';
elseif whichComp==5
    basePath='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/';
else
    basePath='/Users/CHANGE_HERE/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    toolPath='/Users/CHANGE_HERE/Dropbox (Brown)/sharedMatlabUtilities/';

end

dataDir=[basePath, '/behave_data'];

cd(basePath)
addpath(genpath(dataDir))
addpath(genpath(toolPath))

defaultPlotParameters

%%
subList=[2002:2007,2009:2010];

for s=1:length(subList)

sub=subList(s);
    
subNum=num2str(sub);

nBlock=4;

nStim=2;



clear allDataStruct

for block=3:nBlock
    
    fileName=sprintf('vwm_multiStim_color_%s%s.mat',subNum,num2str(block));
    fileDir=fullfile(dataDir,fileName);
    load(fileDir)
    
    nTrial=size(data.reportedValue,1);
    %keyboard
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
    
    estDeg=nan(nTrial,nStim);
    subPredErr=nan(nTrial,nStim);
    
    if ~data.prefs.doFeedback
        
        for r=1:nTrial
            for c=1:nStim
                %if data.chosenTarg(r,c)==1
                if data.chosenTargEst(r,c)==1
                    estDeg(r,c)=data.reportedValue(r,1);
                    est(r,c)=deg2rad(data.reportedValue(r,1));
                    estErr(r,c)=circ_dist(deg2rad(data.presentedColor(r,1)),deg2rad(data.reportedValue(r,1)));
                    trueEst(r,c)=deg2rad(data.presentedColor(r,1));
                    %generativeMeanRad(r,c)=deg2rad(data.generativeMean(r,1));
                else
                    estDeg(r,c)=data.reportedValue(r,2);
                    est(r,c)=deg2rad(data.reportedValue(r,2));
                    estErr(r,c)=circ_dist(deg2rad(data.presentedColor(r,2)),deg2rad(data.reportedValue(r,2)));
                    trueEst(r,c)=deg2rad(data.presentedColor(r,2));
                    %generativeMeanRad(r,c)=deg2rad(data.generativeMean(r,2));
                end
            end
        end
    end
    
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
                    
                    predictErr(r,c)=circ_dist(deg2rad(data.toPredict(r-1,1)),deg2rad(data.predictedValue(r-1,1)));
                    truePred(r,c)=deg2rad(data.toPredict(r-1,1));
                    predEstDiff(r,c)=circ_dist(deg2rad(data.reportedValue(r-1,1)),deg2rad(data.predictedValue(r-1,1)));
                    %predUpdate(r,c)=circ_dist(deg2rad(predVUpdated(r-1,1)),deg2rad(predV(r-1,1)));
                    surpriseTrial(r,c)=data.surprise(r-1,1);
                    %subPredErr(r,c)=circ_dist(deg2rad(data.reportedValue(r,1)),deg2rad(data.predictedValue(r-1,1)));
                else
                    predDeg(r,c)=data.predictedValue(r-1,2);
                    pred(r,c)=deg2rad(data.predictedValue(r-1,2));
                    predictErr(r,c)=circ_dist(deg2rad(data.toPredict(r-1,2)),deg2rad(data.predictedValue(r-1,2)));
                    truePred(r,c)=deg2rad(data.toPredict(r-1,2));
                    predEstDiff(r,c)=circ_dist(deg2rad(data.reportedValue(r-1,2)),deg2rad(data.predictedValue(r-1,2)));
                    %predUpdate(r,c)=circ_dist(deg2rad(predVUpdated(r-1,2)),deg2rad(predV(r-1,2)));
                    surpriseTrial(r,c)=data.surprise(r-1,2);
                    %subPredErr(r,c)=circ_dist(deg2rad(data.reportedValue(r,2)),deg2rad(data.predictedValue(r-1,2)));
                    
                end
            end
        end
        %keyboard
        predUpdate(2:end-1,:)=circ_dist(pred(3:end,:),pred(2:end-1,:));
        subPredErr=circ_dist(est,pred);
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
    blockData.predictErr=predictErr;
    blockData.predEstDiff=predEstDiff;
    %keyboard
    blockData.predUpdate=predUpdate;
    blockData.surpriseTrial=surpriseTrial;
    blockData.predDeg=predDeg;
    blockData.estDeg=estDeg;
    blockData.subPredErr=subPredErr;
    blockData.generativeMeanRad=deg2rad(blockData.generativeMean);
    if sub==2002 || sub==2007|| sub==2010
        blockData.condition=[1];
    else
        blockData.condition=[2];
    end
    
    blockData=straightStruct(blockData);
    
    if ~exist('allDataStruct')
        allDataStruct=blockData;
    else
        allDataStruct=catBehav(blockData, allDataStruct, true);
    end
    
    
end





alldata=allDataStruct;

%subNum='20000';

disp(subNum)

saveDir=[basePath, '/behaveData/subCombined'];

fn=fullfile(saveDir,[subNum,'_allBlockData.mat']);

save(fn,'alldata')
end



%% Data check:
 figure()
hold on
 plot(allDataStruct.generativeMeanRad(:,1), '--k')
% plot(allDataStruct.colorArray(:,2), 'or')

plot(allDataStruct.trueEst(:,1),'x')


%%

% sel=allDataStruct.block==2;
% blkRand_struct=selBehav(allDataStruct, sel);


    if allDataStruct.condition(1)==1
    allDataStruct.blockCond(1:length(allDataStruct.block)./2)=1;
    allDataStruct.blockCond(length(allDataStruct.block)./2 +1:length(allDataStruct.block))=-1;
    else
        allDataStruct.blockCond(1:length(allDataStruct.block)./2)=-1;
    allDataStruct.blockCond(length(allDataStruct.block)./2 +1:length(allDataStruct.block))=1;
    end
        
    
    
    %1=CP, -1 = oddball
if isfield(allDataStruct, 'condition')
allDataStruct=rmfield(allDataStruct, 'condition')
end

allDataStruct=straightStruct(allDataStruct);

sel=allDataStruct.block==4;
blkCP_struct=selBehav(allDataStruct, sel);
fileName=sprintf('vwm_multiStim_color_%s%s.mat',subNum,num2str(3));
fileDir=fullfile(dataDir,fileName);
load(fileDir)
blkCP_struct.generativeMean=data.generativeMean;

sel=allDataStruct.block==3;
blkOB_struct=selBehav(allDataStruct, sel);
fileName=sprintf('vwm_multiStim_color_%s%s.mat',subNum,num2str(4));
fileDir=fullfile(dataDir,fileName);
load(fileDir)
blkOB_struct.generativeMean=data.generativeMean;



%%

blkOB_struct.surpriseTrial(1,:)=0;
blkCP_struct.surpriseTrial(1,:)=0;

blkOB_struct.noSurpriseTrial=(blkOB_struct.surpriseTrial-1)*-1;
blkCP_struct.noSurpriseTrial=(blkCP_struct.surpriseTrial-1)*-1;

% muEstErrRand=nanmean(blkRand_struct.estErr,'all');
% sigmaEstErrRand=nanstd(blkRand_struct.estErr,0,'all');

muEstErrCP=nanmean(blkCP_struct.estErr,'all');
sigmaEstErrCP=nanstd(blkCP_struct.estErr,0,'all');

muEstErrOB=nanmean(blkOB_struct.estErr,'all');
sigmaEstErrOB=nanstd(blkOB_struct.estErr,0,'all');

muEstErrWPredict=nanmean(cat(1,blkOB_struct.estErr,blkCP_struct.estErr),'all');
sigmaEstErrWPredict=nanstd(cat(1,blkOB_struct.estErr,blkCP_struct.estErr),0,'all');

muEstErrWPredictSurprise=nanmean(cat(1,blkOB_struct.estErr.*blkOB_struct.surpriseTrial,blkCP_struct.estErr.*blkCP_struct.surpriseTrial),'all');
sigmaEstErrWPredictSurprise=nanstd(cat(1,blkOB_struct.estErr.*blkOB_struct.surpriseTrial,blkCP_struct.estErr.*blkCP_struct.surpriseTrial),0,'all');

muEstErrWPredictNoSurprise=nanmean(cat(1,blkOB_struct.estErr.*blkOB_struct.noSurpriseTrial,blkCP_struct.estErr.*blkCP_struct.surpriseTrial),'all');
sigmaEstErrWPredictNoSurprise=nanstd(cat(1,blkOB_struct.estErr.*blkOB_struct.noSurpriseTrial,blkCP_struct.estErr.*blkCP_struct.surpriseTrial),0,'all');

muPredErrCP=nanmean(blkCP_struct.predictErr,'all');
sigmaPredErrCP=nanstd(blkCP_struct.predictErr,0,'all');

muPredErrOB=nanmean(blkOB_struct.predictErr,'all');
sigmaPredErrOB=nanstd(blkOB_struct.predictErr,0,'all');

muPredErrWPredict=nanmean(cat(1,blkOB_struct.predictErr,blkCP_struct.predictErr),'all');
sigmaPredErrWPredict=nanstd(cat(1,blkOB_struct.predictErr,blkCP_struct.predictErr),0,'all');


%%


surpriseTrialCPstim1=find(blkCP_struct.surpriseTrial(:,1)==1);
surpriseTrialCPstim2=find(blkCP_struct.surpriseTrial(:,2)==1);

normalTrialCPstim1=find(blkCP_struct.surpriseTrial(:,1)==0);
normalTrialCPstim2=find(blkCP_struct.surpriseTrial(:,2)==0);

surpriseTrialOBstim1=find(blkCP_struct.surpriseTrial(:,1)==1);
surpriseTrialOBstim2=find(blkCP_struct.surpriseTrial(:,2)==1);

normalTrialOBstim1=find(blkCP_struct.surpriseTrial(:,1)==0);
normalTrialOBstim2=find(blkCP_struct.surpriseTrial(:,2)==0);










%%
figure()

subplot(4,2,1)
hold on
s=muEstErrRand-3*sigmaEstErrRand:0.01:muEstErrRand+3*sigmaEstErrRand;
plot(s,pdf('Normal',s,muEstErrRand,sigmaEstErrRand))

title('Randomized stim order Estimation Error')
xlim([-pi pi])
%ylim([0 1])
hold off

subplot(4,2,2)
hold on
histogram(blkRand_struct.estErr)
title('Randomized stim order Estimation Error')
xlim([-pi pi])
hold off

subplot(4,2,3)
hold on
s=muPredErrWPredict-3*sigmaPredErrWPredict:0.01:muPredErrWPredict+3*sigmaPredErrWPredict;
plot(s,pdf('Normal',s,muPredErrWPredict,sigmaPredErrWPredict))
title('Prediction Error')
xlim([-pi pi])
%ylim([0 1])
hold off

subplot(4,2,4)
hold on
histogram(cat(1,blkCP_struct.predictErr,blkOB_struct.predictErr))
title('Prediction Error')
xlim([-pi pi])
hold off

subplot(4,2,5)
hold on
s=muEstErrWPredictNoSurprise-3*sigmaEstErrWPredictNoSurprise:0.01:muEstErrWPredictNoSurprise+3*sigmaEstErrWPredictNoSurprise;
plot(s,pdf('Normal',s,muEstErrWPredictNoSurprise,sigmaEstErrWPredictNoSurprise))
title('Estimation Error with Prediction no Surprise')
xlim([-pi pi])
%ylim([0 1])
hold off

subplot(4,2,6)
hold on
histogram(nonzeros(cat(1,blkOB_struct.estErr.*blkOB_struct.noSurpriseTrial,blkCP_struct.estErr.*blkCP_struct.noSurpriseTrial)))
title('Estimation Error with Prediction no Surprise')
xlim([-pi pi])
hold off

subplot(4,2,7)
hold on
s=muEstErrWPredictSurprise-3*sigmaEstErrWPredictSurprise:0.01:muEstErrWPredictSurprise+3*sigmaEstErrWPredictSurprise;
plot(s,pdf('Normal',s,muEstErrWPredictSurprise,sigmaEstErrWPredictSurprise))
title('Estimation Error with Prediction with Surprise')
xlim([-pi pi])
%ylim([0 1])
hold off

subplot(4,2,8)
hold on
histogram(nonzeros(cat(1,blkOB_struct.estErr.*blkOB_struct.surpriseTrial,blkCP_struct.estErr.*blkCP_struct.surpriseTrial)))
title('Estimation Error with Prediction with Surprise')
xlim([-pi pi])
hold off

% figure()
% 
% subplot(3,2,1)
% hold on
% s=muEstErrRand-3*sigmaEstErrRand:0.01:muEstErrRand+3*sigmaEstErrRand;
% plot(s,pdf('Normal',s,muEstErrRand,sigmaEstErrRand))
% title('Randomized stim order Estimation Error')
% xlim([-pi pi])
% %ylim([0 1])
% hold off
% 
% subplot(3,2,2)
% hold on
% histogram(blkRand_struct.estErr)
% title('Randomized stim order Estimation Error')
% xlim([-pi pi])
% hold off
% 
% subplot(3,2,3)
% hold on
% s=muPredErrWPredict-3*sigmaPredErrWPredict:0.01:muPredErrWPredict+3*sigmaPredErrWPredict;
% plot(s,pdf('Normal',s,muPredErrWPredict,sigmaPredErrWPredict))
% title('Prediction Error')
% xlim([-pi pi])
% %ylim([0 1])
% hold off
% 
% subplot(3,2,4)
% hold on
% histogram(cat(1,blkCP_struct.predictErr,blkOB_struct.predictErr))
% title('Prediction Error')
% xlim([-pi pi])
% hold off
% 
% subplot(3,2,5)
% hold on
% s=muEstErrWPredict-3*sigmaEstErrWPredict:0.01:muEstErrWPredict+3*sigmaEstErrWPredict;
% plot(s,pdf('Normal',s,muEstErrWPredict,sigmaEstErrWPredict))
% title('Estimation Error with Prediction blocks (CP/OB)')
% xlim([-pi pi])
% %ylim([0 1])
% hold off
% 
% subplot(3,2,6)
% hold on
% histogram(cat(1,blkCP_struct.estErr,blkOB_struct.estErr))
% title('Estimation Error with Prediction blocks (CP/OB)')
% xlim([-pi pi])
% hold off


%%
figure()

subplot(3,1,1)
hold on
plot([-pi, pi], [0, 0], '--k')
plot(cat(1,blkCP_struct.predEstDiff,blkOB_struct.predEstDiff),cat(1,blkCP_struct.predUpdate,blkOB_struct.predUpdate),"ob")
xlabel('Subjective prediction error')
ylabel('Prediction update')
title('All with prediction')
diagline=refline(1,0);
diagline.Color='k';
hold off

subplot(3,1,2)
hold on
plot([-pi, pi], [0, 0], '--k')
plot(blkCP_struct.predEstDiff,blkCP_struct.predUpdate,"ob")
xlabel('Subjective prediction error')
ylabel('Prediction update')
title('Change point')
diagline=refline(1,0);
diagline.Color='k';
hold off

subplot(3,1,3)
hold on
plot([-pi, pi], [0, 0], '--k')
plot(blkOB_struct.predEstDiff,blkOB_struct.predUpdate,"ob")
xlabel('Subjective prediction error')
ylabel('Prediction update')
title('Odd ball')
diagline=refline(1,0);
diagline.Color='k';
hold off


%%

figure()

% subplot(2,2,1)
% %plot(blkCP_struct.predDeg(:,1),'-m')
% hold on
% plot(blkCP_struct.colorArray(:,1),'ok')
% %plot(blkCP_struct.estDeg(:,1),'-b')
% plot(blkCP_struct.generativeMean(:,1),'--r')
% title('stim 1 Change Point')
% xlabel('trial')
% ylabel('color (rad)')
% 
% subplot(2,2,3)
% %plot(blkCP_struct.predDeg(:,2),'-m')
% hold on
% plot(blkCP_struct.colorArray(:,2),'ok')
% %plot(blkCP_struct.estDeg(:,2),'-b')
% plot(blkCP_struct.generativeMean(:,2),'--r')
% title('stim 2 Change Point')
% xlabel('trial')
% ylabel('color (rad)')

subplot(2,2,2)
%plot(blkOB_struct.predDeg(:,1),'-m')
hold on
plot(blkOB_struct.colorArray(:,1),'ok')
%plot(blkOB_struct.estDeg(:,1),'-b')
plot(blkOB_struct.generativeMean(:,1),'--r')
title('stim 1 Odd Ball')
xlabel('trial')
xticks([10 20 30 40 50])
xlim([0 50])
ylabel('color (rad)')

% subplot(2,2,4)
% %plot(blkOB_struct.predDeg(:,2),'-m')
% hold on
% plot(blkOB_struct.colorArray(:,2),'ok')
% %plot(blkOB_struct.estDeg(:,2),'-b')
% plot(blkOB_struct.generativeMean(:,2),'--r')
% title('stim 1 Odd Ball')
% xlabel('trial')
% ylabel('color (rad)')


%%
blkOB_struct.predDegUpdate=circ_dist(deg2rad(blkOB_struct.predDeg(3:end,:)),deg2rad(blkOB_struct.predDeg(2:end-1,:)));




%%
figure()

subplot(2,2,1)
plot(blkCP_struct.predDeg(:,1),'-m')
hold on
plot(blkCP_struct.colorArray(:,1),'ok')
plot(blkCP_struct.estDeg(:,1),'-b')
plot(blkCP_struct.generativeMean(:,1),'--r')

x=[CPmysteryTstim1';CPmysteryTstim1'];
y=repmat([0;360],1,length(x));
plot(x,y,'-c')

title('stim 1 Change Point')
xlabel('trial')
ylabel('color (deg)')

subplot(2,2,3)
plot(blkCP_struct.predDeg(:,2),'-m')
hold on
plot(blkCP_struct.colorArray(:,2),'ok')
plot(blkCP_struct.estDeg(:,2),'-b')
plot(blkCP_struct.generativeMean(:,2),'--r')
x=[CPmysteryTstim2';CPmysteryTstim2'];
y=repmat([0;360],1,length(x));
plot(x,y,'-c')
xline(CPstrangeTStim2,'-g')
title('stim 2 Change Point')
xlabel('trial')
ylabel('color (deg)')

subplot(2,2,2)
plot(blkOB_struct.predDeg(:,1),'-m')
hold on
plot(blkOB_struct.colorArray(:,1),'ok')
plot(blkOB_struct.estDeg(:,1),'-b')
plot(blkOB_struct.generativeMean(:,1),'--r')
x=[OBmysteryTstim1';OBmysteryTstim1'];
y=repmat([0;360],1,length(x));
plot(x,y,'-c')
xline(OBstrangeTStim1,'-g')
title('stim 1 Odd Ball')
xlabel('trial')
ylabel('color (deg)')

subplot(2,2,4)
plot(blkOB_struct.predDeg(:,2),'-m')
hold on
plot(blkOB_struct.colorArray(:,2),'ok')
plot(blkOB_struct.estDeg(:,2),'-b')
plot(blkOB_struct.generativeMean(:,2),'--r')
x=[OBmysteryTstim2';OBmysteryTstim2'];
y=repmat([0;360],1,length(x));
plot(x,y,'-c')
xline(OBstrangeTStim2,'-g')
title('stim 2 Odd Ball')
xlabel('trial')
ylabel('color (deg)')



%% plot rad raw data

figure()

subplot(2,2,1)
plot(blkCP_struct.pred(:,1),'-m')
hold on
plot(blkCP_struct.trueEst(:,1),'ok')
plot(blkCP_struct.est(:,1),'-b')
plot(blkCP_struct.generativeMeanRad(:,1),'--r')
plot(blkCP_struct.predUpdate(:,1),'*')
%x=[CPmysteryTstim1';CPmysteryTstim1'];
% y=repmat([0;360],1,length(x));
% plot(x,y,'-c')

title('stim 1 Change Point')
xlabel('trial')
ylabel('color (rad)')

subplot(2,2,3)
plot(blkCP_struct.pred(:,2),'-m')
hold on
plot(blkCP_struct.trueEst(:,2),'ok')
plot(blkCP_struct.est(:,2),'-b')
plot(blkCP_struct.generativeMeanRad(:,2),'--r')
plot(blkCP_struct.predUpdate(:,2),'*')
%x=[CPmysteryTstim2';CPmysteryTstim2'];
% y=repmat([0;360],1,length(x));
% plot(x,y,'-c')
%xline(CPstrangeTStim2,'-g')
title('stim 2 Change Point')
xlabel('trial')
ylabel('color (rad)')

subplot(2,2,2)
plot(blkOB_struct.pred(:,1),'-m')
hold on
plot(blkOB_struct.trueEst(:,1),'ok')
plot(blkOB_struct.est(:,1),'-b')
plot(blkOB_struct.generativeMeanRad(:,1),'--r')
plot(blkOB_struct.predUpdate(:,1),'*')
%x=[OBmysteryTstim1';OBmysteryTstim1'];
% y=repmat([0;360],1,length(x));
% plot(x,y,'-c')
%xline(OBstrangeTStim1,'-g')
title('stim 1 Odd Ball')
xlabel('trial')
ylabel('color (rad)')

subplot(2,2,4)
plot(blkOB_struct.pred(:,2),'-m')
hold on
plot(blkOB_struct.trueEst(:,2),'ok')
plot(blkOB_struct.est(:,2),'-b')
plot(blkOB_struct.generativeMeanRad(:,2),'--r')
plot(blkOB_struct.predUpdate(:,2),'*')
%x=[OBmysteryTstim2';OBmysteryTstim2'];
% y=repmat([0;360],1,length(x));
% plot(x,y,'-c')
%xline(OBstrangeTStim2,'-g')
title('stim 2 Odd Ball')
xlabel('trial')
ylabel('color (rad)')








%%
figure()

subplot(3,3,1)
hold on
plot(blkRand_struct.trueEst,blkRand_struct.est,'ok')
xlabel('True color (rad)')
xlim([0 2*pi])
ylabel('Reported perception (rad)')
ylim([0 2*pi])
title('Randomized stim order')
diagline=refline(1,0);
diagline.Color='r';
hold off

subplot(3,3,2)
hold on
plot(blkCP_struct.trueEst,blkCP_struct.est,'ok')
xlabel('True color (rad)')
xlim([0 2*pi])
ylabel('Reported perception (rad)')
ylim([0 2*pi])
title('Change point')
diagline=refline(1,0);
diagline.Color='r';
hold off

subplot(3,3,3)
hold on
plot(blkOB_struct.trueEst,blkOB_struct.est,'ok')
xlabel('True color (rad)')
xlim([0 2*pi])
ylabel('Reported perception (rad)')
ylim([0 2*pi])
title('Odd ball')
diagline=refline(1,0);
diagline.Color='r';
hold off

subplot(3,3,5)
hold on
plot(blkCP_struct.truePred,blkCP_struct.pred,'ok')
xlabel('True color (rad)')
xlim([0 2*pi])
ylabel('Reported prediction (rad)')
ylim([0 2*pi])
title('Change point')
diagline=refline(1,0);
diagline.Color='r';
hold off

subplot(3,3,6)
hold on
plot(blkOB_struct.truePred,blkOB_struct.pred,'ok')
xlabel('True color (rad)')
xlim([0 2*pi])
ylabel('Reported prediction (rad)')
ylim([0 2*pi])
title('Odd ball')
diagline=refline(1,0);
diagline.Color='r';
hold off

subplot(3,3,8)
hold on
a=blkCP_struct.truePred(:,1);
b=blkCP_struct.truePred(:,2);
plotX=cat(1,a(surpriseTrialCPstim1),b(surpriseTrialCPstim2));
aa=blkCP_struct.pred(:,1);
bb=blkCP_struct.pred(:,2);
plotY=cat(1,aa(surpriseTrialCPstim1),bb(surpriseTrialCPstim2));
plot(plotX,plotY,'ok')
xlabel('True color (rad)')
xlim([0 2*pi])
ylabel('Reported prediction (rad)')
ylim([0 2*pi])
title('Change point surprise')
diagline=refline(1,0);
diagline.Color='r';
hold off

subplot(3,3,9)
hold on
a=blkOB_struct.truePred(:,1);
b=blkOB_struct.truePred(:,2);
plotX=cat(1,a(surpriseTrialOBstim1),b(surpriseTrialOBstim2));
aa=blkOB_struct.pred(:,1);
bb=blkOB_struct.pred(:,2);
plotY=cat(1,aa(surpriseTrialOBstim1),bb(surpriseTrialOBstim2));
plot(plotX,plotY,'ok')
xlabel('True color (rad)')
xlim([0 2*pi])
ylabel('Reported prediction (rad)')
ylim([0 2*pi])
title('Odd ball surprise')
diagline=refline(1,0);
diagline.Color='r';
hold off

%%
figure()

subplot(1,2,1)
hold on
a=plot(blkCP_struct.predictErr,blkCP_struct.estErr,'o', 'lineWidth', 1, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k')
set(a(2), 'markerFaceColor', cbColors(3,:))

title('Change point')
xlabel('Prediction error (rad)')
ylabel('Perception error (rad)')
diagline=refline(1,0);
diagline.Color='k';
hold off
set(gca, 'box', 'off')










subplot(1,2,2)
hold on
plot(blkOB_struct.predictErr,blkOB_struct.estErr,'o')
title('Odd ball')
xlabel('Prediction error (rad)')
ylabel('Perception error (rad)')
diagline=refline(1,0);
diagline.Color='k';
hold off

%% Perceptual bias fit CP

estErrWpredict=cat(1,blkCP_struct.estErr(:,1),blkCP_struct.estErr(:,2));
predErrAll=cat(1,blkCP_struct.predictErr(:,1),blkCP_struct.predictErr(:,2));

predErrSmall=(predErrAll<1&predErrAll>-1);

%trialType=cat(1,blkCP_struct.surpriseTrial(:,1),blkCP_struct.surpriseTrial(:,2),blkOB_struct.surpriseTrial(:,1),blkOB_struct.surpriseTrial(:,2));

%surpriseT=find(trialType==1);
%normalT=find(trialType==0);

x0=ones(length(estErrWpredict),1);
x1=predErrAll;
x1(predErrSmall)=0;
x2=(predErrAll);
x2(~predErrSmall)=0;
%x1(surpriseT)=0;
%x2=predErrAll;
%x2(normalT)=0;
%xMat_1=[x0,x1,x2];
xMat_1=[x0,x1,x2];
yVar_1=estErrWpredict;
isGood_1=isfinite(yVar_1)& all(isfinite(xMat_1), 2);

regData_1.Y=yVar_1(isGood_1);        % Get updates where things are finite!
regData_1.X=xMat_1(isGood_1,:);      % Get updates where things are finite!
regData_1.includeUniform = true;     % include the possibility that some trials are guesses?
% first param = concentration, 2-end-1 = predictors,
regData_1.whichParams = [true; true(size(xMat_1, 2),1);  false] ; % Last element is for frequency of uniform errors
regData_1.labels={'Concentration', 'Intercept', 'PredictionError normal','PredictionError surprise'};

regData_1.lb=[.01, -5, -5, -5]; % lb = lower bound
regData_1.ub=[ 50,  5,  5,  5]; % ub = upper bound

% Set values that you are NOT fitting here (other values will get replaced)
initStart=[10, zeros(1, size(xMat_1, 2)), .25];

bestNegLogLike=inf;
clear params
for ff = 1:50
    regData.startPoint(regData_1.whichParams==false)=initStart(regData_1.whichParams==false)
    regData.startPoint(regData_1.whichParams)= ...
        regData_1.lb(regData_1.whichParams) + ...
        (regData_1.ub(regData_1.whichParams)-regData_1.lb(regData_1.whichParams)).*rand(1, sum(regData_1.whichParams));
    
    
    [params_run, negLogLike_run]=fitLinearModWCircErrs(regData_1);
    negLogLike_run;
    if negLogLike_run<bestNegLogLike
        params=params_run;
        negLogLike=negLogLike_run;
        %allVM_likelihoods=allVM_likelihoods_run;
    end
    
    
    
end

figure()
hold on

%plot(predErrAll, estErrWpredict, 'ob')
plot(predErrAll, estErrWpredict, 'or')

exX=[-pi, pi];

bigPredErr=params(2)+params(3).*exX;
smallPredErr=params(2)+params(4).*exX;
%nonSurpPred=params(2)+params(3).*exX;
%surpPred=params(2)+params(4).*exX;

plot(exX, bigPredErr,'-r')
plot(exX, smallPredErr,'-b')

%plot(exX, nonSurpPred, '-b')
%plot(exX, surpPred, '-r')

plot(predErrAll(predErrSmall),estErrWpredict(predErrSmall),'ob')
%plot(predErrAll(surpriseT), estErrWpredict(surpriseT), 'or')
xlim(exX)
ylim(exX)

xlabel('prediction error')
ylabel('perceptual error')

%% Perceptual bias fit OB

estErrWpredict=cat(1,blkOB_struct.estErr(:,1),blkOB_struct.estErr(:,2));
predErrAll=cat(1,blkOB_struct.predictErr(:,1),blkOB_struct.predictErr(:,2));

predErrSmall=(predErrAll<1&predErrAll>-1);

%trialType=cat(1,blkCP_struct.surpriseTrial(:,1),blkCP_struct.surpriseTrial(:,2),blkOB_struct.surpriseTrial(:,1),blkOB_struct.surpriseTrial(:,2));

%surpriseT=find(trialType==1);
%normalT=find(trialType==0);

x0=ones(length(estErrWpredict),1);
x1=predErrAll;
x1(predErrSmall)=0;
x2=(predErrAll);
x2(~predErrSmall)=0;
%x1(surpriseT)=0;
%x2=predErrAll;
%x2(normalT)=0;
%xMat_1=[x0,x1,x2];
xMat_1=[x0,x1,x2];
yVar_1=estErrWpredict;
isGood_1=isfinite(yVar_1)& all(isfinite(xMat_1), 2);

regData_1.Y=yVar_1(isGood_1);        % Get updates where things are finite!
regData_1.X=xMat_1(isGood_1,:);      % Get updates where things are finite!
regData_1.includeUniform = true;     % include the possibility that some trials are guesses?
% first param = concentration, 2-end-1 = predictors,
regData_1.whichParams = [true; true(size(xMat_1, 2),1);  false] ; % Last element is for frequency of uniform errors
regData_1.labels={'Concentration', 'Intercept', 'PredictionError normal','PredictionError surprise'};

regData_1.lb=[.01, -5, -5, -5]; % lb = lower bound
regData_1.ub=[ 50,  5,  5,  5]; % ub = upper bound

% Set values that you are NOT fitting here (other values will get replaced)
initStart=[10, zeros(1, size(xMat_1, 2)), .25];

bestNegLogLike=inf;
clear params
for ff = 1:50
    regData.startPoint(regData_1.whichParams==false)=initStart(regData_1.whichParams==false)
    regData.startPoint(regData_1.whichParams)= ...
        regData_1.lb(regData_1.whichParams) + ...
        (regData_1.ub(regData_1.whichParams)-regData_1.lb(regData_1.whichParams)).*rand(1, sum(regData_1.whichParams));
    
    
    [params_run, negLogLike_run]=fitLinearModWCircErrs(regData_1);
    negLogLike_run;
    if negLogLike_run<bestNegLogLike
        params=params_run;
        negLogLike=negLogLike_run;
        %allVM_likelihoods=allVM_likelihoods_run;
    end
    
    
    
end


figure()
hold on

%plot(predErrAll, estErrWpredict, 'ob')
plot(predErrAll, estErrWpredict, 'or')

exX=[-pi, pi];

bigPredErr=params(2)+params(3).*exX;
smallPredErr=params(2)+params(4).*exX;
%nonSurpPred=params(2)+params(3).*exX;
%surpPred=params(2)+params(4).*exX;

plot(exX, bigPredErr,'-r')
plot(exX, smallPredErr,'-b')

%plot(exX, nonSurpPred, '-b')
%plot(exX, surpPred, '-r')

plot(predErrAll(predErrSmall),estErrWpredict(predErrSmall),'ob')
%plot(predErrAll(surpriseT), estErrWpredict(surpriseT), 'or')
xlim(exX)
ylim(exX)

xlabel('prediction error')
ylabel('perceptual error')

%% Make plot of update against subjective prediction error for CP
% Steps:
% 1) select data
% 2) identify big errors
% 3) 

% Analysis parameters:
crit=1.5; % CRITERION for separating big and small errors


predictUpdateAll=cat(1,blkCP_struct.predUpdate(:,1),blkCP_struct.predUpdate(:,2));
subjectivePredErr=cat(1,blkCP_struct.predEstDiff(:,1),blkCP_struct.predEstDiff(:,2));

%predErrAll=cat(1,blkCP_struct.predictErr(:,1),blkCP_struct.predictErr(:,2));

predErrSmall=(subjectivePredErr<crit&subjectivePredErr>-crit);


x0=ones(length(predictUpdateAll),1);
x1=subjectivePredErr;
x1(~predErrSmall)=0;
x2=subjectivePredErr;
x2(predErrSmall)=0;
xMat_1=[x0,x1,x2];
yVar_1=predictUpdateAll;
isGood_1=isfinite(yVar_1)& all(isfinite(xMat_1), 2);

regData_1.Y=yVar_1(isGood_1);        % Get updates where things are finite!
regData_1.X=xMat_1(isGood_1,:);      % Get updates where things are finite!
regData_1.includeUniform = true;     % include the possibility that some trials are guesses?
% first param = concentration, 2-end-1 = predictors,
regData_1.whichParams = [true; true(size(xMat_1, 2),1);  false] ; % Last element is for frequency of uniform errors
regData_1.labels={'Concentration', 'Intercept', 'PredictionError normal','PredictionError surprise'};

regData_1.lb=[.01, -5, -5, -5]; % lb = lower bound
regData_1.ub=[ 50,  5,  5,  5]; % ub = upper bound

% Set values that you are NOT fitting here (other values will get replaced)
initStart=[10, zeros(1, size(xMat_1, 2)), .25];

bestNegLogLike=inf;
clear params
for ff = 1:50
    regData.startPoint(regData_1.whichParams==false)=initStart(regData_1.whichParams==false)
    regData.startPoint(regData_1.whichParams)= ...
        regData_1.lb(regData_1.whichParams) + ...
        (regData_1.ub(regData_1.whichParams)-regData_1.lb(regData_1.whichParams)).*rand(1, sum(regData_1.whichParams));
    
    
    [params_run, negLogLike_run]=fitLinearModWCircErrs(regData_1);
    negLogLike_run;
    if negLogLike_run<bestNegLogLike
        params=params_run;
        negLogLike=negLogLike_run;
        %allVM_likelihoods=allVM_likelihoods_run;
    end
    
    
    
end




t=find(subjectivePredErr<crit&subjectivePredErr>-crit);
figure()
hold on

%plot(subjectivePredErr, predictUpdateAll, 'ob')

exX=[-pi, pi];

smallTrials=params(2)+params(3).*exX;
bigTrials=params(2)+params(4).*exX;

plot([-pi, pi], [0, 0], '--k')
plot([-pi, pi], [-pi,pi], '--k')



plot(exX, smallTrials,'-r')
plot(exX, bigTrials,'-b')

x=subjectivePredErr;
y=predictUpdateAll;

plot(x,y,'or', 'lineWidth', 1)
plot(x(t),y(t),'ob', 'lineWidth', 1)


xlim(exX)
ylim(exX)

xlabel('subjective prediction error')
ylabel('prediction update')

saveas(gcf, 'CP_updateVsSPE.eps','epsc2')

%% Make plot of update against subjective prediction error for OB

crit=1.5; % CRITERION for separating big and small errors

predictUpdateAll=cat(1,blkOB_struct.predUpdate(:,1),blkOB_struct.predUpdate(:,2));
subjectivePredErr=cat(1,blkOB_struct.predEstDiff(:,1),blkOB_struct.predEstDiff(:,2));

%predErrAll=cat(1,blkCP_struct.predictErr(:,1),blkCP_struct.predictErr(:,2));

predErrSmall=(subjectivePredErr<crit&subjectivePredErr>-crit);


x0=ones(length(predictUpdateAll),1);
x1=subjectivePredErr;
x1(~predErrSmall)=0;
x2=subjectivePredErr;
x2(predErrSmall)=0;
xMat_1=[x0,x1,x2];
yVar_1=predictUpdateAll;
isGood_1=isfinite(yVar_1)& all(isfinite(xMat_1), 2);

regData_1.Y=yVar_1(isGood_1);        % Get updates where things are finite!
regData_1.X=xMat_1(isGood_1,:);      % Get updates where things are finite!
regData_1.includeUniform = true;     % include the possibility that some trials are guesses?
% first param = concentration, 2-end-1 = predictors,
regData_1.whichParams = [true; true(size(xMat_1, 2),1);  false] ; % Last element is for frequency of uniform errors
regData_1.labels={'Concentration', 'Intercept', 'PredictionError normal','PredictionError surprise'};

regData_1.lb=[.01, -5, -5, -5]; % lb = lower bound
regData_1.ub=[ 50,  5,  5,  5]; % ub = upper bound

% Set values that you are NOT fitting here (other values will get replaced)
initStart=[10, zeros(1, size(xMat_1, 2)), .25];

bestNegLogLike=inf;
clear params
for ff = 1:50
    regData.startPoint(regData_1.whichParams==false)=initStart(regData_1.whichParams==false)
    regData.startPoint(regData_1.whichParams)= ...
        regData_1.lb(regData_1.whichParams) + ...
        (regData_1.ub(regData_1.whichParams)-regData_1.lb(regData_1.whichParams)).*rand(1, sum(regData_1.whichParams));
    
    
    [params_run, negLogLike_run]=fitLinearModWCircErrs(regData_1);
    negLogLike_run;
    if negLogLike_run<bestNegLogLike
        params=params_run;
        negLogLike=negLogLike_run;
        %allVM_likelihoods=allVM_likelihoods_run;
    end
    
    
    
end




t=find(subjectivePredErr<crit&subjectivePredErr>-crit);
figure()
hold on

%plot(subjectivePredErr, predictUpdateAll, 'ob')

exX=[-pi, pi];

smallTrials=params(2)+params(3).*exX;
bigTrials=params(2)+params(4).*exX;

plot([-pi, pi], [0, 0], '--k')
plot([-pi, pi], [-pi,pi], '--k')



plot(exX, smallTrials,'-r')
plot(exX, bigTrials,'-b')

x=subjectivePredErr;
y=predictUpdateAll;

plot(x,y,'or', 'lineWidth', 1)
plot(x(t),y(t),'ob', 'lineWidth', 1)


xlim(exX)
ylim(exX)

xlabel('subjective prediction error')
ylabel('prediction update')


%%

%Q:how to define range of small vs big errors 

%%
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

klAll1=cat(2,dataCP1.kl,dataOB1.kl);
klAll2=cat(2,dataCP2.kl,dataOB2.kl);

surpriseAll1=cat(2,dataCP1.surprise,dataOB1.surprise);
surpriseAll2=cat(2,dataCP2.surprise,dataOB2.surprise);

kl = [mean(klCP,1),mean(klOB,1)]';
surprise = [mean(surpriseCP,1),mean(surpriseOB,1)]';

x0=ones(length(alldata.block),1);
x1=klAll2';
x2=surpriseAll2';
xMat_1=[x0,x1,x2];
yVar_1=alldata.estErr(1,:);

regData_1.Y=yVar_1;        % Get updates where things are finite!
regData_1.X=xMat_1;      % Get updates where things are finite!
regData_1.includeUniform = true;     % include the possibility that some trials are guesses?
% first param = concentration, 2-end-1 = predictors,
regData_1.whichParams = [true; true(size(xMat_1, 2),1);  false] ; % Last element is for frequency of uniform errors
regData_1.labels={'Concentration', 'Intercept', 'PredictionError normal','PredictionError surprise'};

regData_1.lb=[.01, -5, -5, -5]; % lb = lower bound
regData_1.ub=[ 50,  5,  5,  5]; % ub = upper bound

% Set values that you are NOT fitting here (other values will get replaced)
initStart=[10, zeros(1, size(xMat_1, 2)), .25];

bestNegLogLike=inf;
clear params
for ff = 1:50
    regData.startPoint(regData_1.whichParams==false)=initStart(regData_1.whichParams==false)
    regData.startPoint(regData_1.whichParams)= ...
        regData_1.lb(regData_1.whichParams) + ...
        (regData_1.ub(regData_1.whichParams)-regData_1.lb(regData_1.whichParams)).*rand(1, sum(regData_1.whichParams));
    
    
    [params_run, negLogLike_run]=fitLinearModWCircErrs(regData_1);
    negLogLike_run;
    if negLogLike_run<bestNegLogLike
        params=params_run;
        negLogLike=negLogLike_run;
        %allVM_likelihoods=allVM_likelihoods_run;
    end
    
    
    
end



