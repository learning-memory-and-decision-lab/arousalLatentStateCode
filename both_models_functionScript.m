
whichComp=4;


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



%use this

%%


runLoop=1;

%list of subject numbers
subList=[2002:2007,2009:2010];

if runLoop
    behaveAll=[];
    
end


%%

for s=1:length(subList)  %loop through all subs
    
    subNum=subList(s);
    subStr=num2str(subNum);
    
    fileName=sprintf('behaveData/subCombined/%s_allBlockData.mat',subStr);
    
    DP = fullfile(baseDir, fileName);
    %DP='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/behaveData//subCombined/2005_allBlockData.mat';
    %DP="/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/Meera1.mat";
    vars = shared_variables(DP);
    vars.H = .15;
    numReps=40;
    dataPath=DP;
    
    data = load(vars.path);
    allDataStruct = data.alldata;
    
    %both = both_models_function(vars);
    
    
    %function params  = both_models_function(vars)
    %Change-point condition data
    % vars = shared_variables;
    %suffix = '_50';
    if vars.generate_data
        vars.sigma_x=30;
        dataCP = learningAndPerception(vars);
        
        vars.sigma_c=20;
        vars.sigma_x=5;
        %Oddball condition datap
        dataOB = oddBall_Inference(vars);
    else
        %     "specify path to data from both conditions"
        %     dataCP = importdata(vars.path_CP); % will replace this with path_CP
        %     dataOB = importdata(vars.path_OB);
        
        %     dataCP = learningAndPerception(vars);
        %
        %     %Oddball condition datap
        %     dataOB = oddBall_Inference(vars);
        
        % STeps:
        % 1 -- use colorArray in modeling script instead of presented color
        % 2 -- run modeling script lots of times, take "average" of everything
        %      you will use.
        % 3 -- create one structure to store data could be allDataStruct. Add
        %       fields for each variable you will extract from model. Plug in
        %        "averaged" model values to those fields for appropriate condition.
        % 4 -- rerun behavioral analyses with averaged model data. Make sure to
        %      do a check of whether new "surprise" measures make sense. For
        %      example plot update as a function of prediction error... plot
        %      trials with "high suprrise" in a different color. Are they the
        %      big prediction error trials? I hope so.
        
        
        % Once we've submitted -- goal: streamlining code so you can hit play,
        % load arbitrary number of datasets, perform all analyses and generate
        % all figures.
        
        dataCP = learningAndPerception(vars);
        
        allDataStruct.nTrialsCP = dataCP.nTrialsCP;
        allDataStruct.YCP = dataCP.YCP;
        allDataStruct.X = dataCP.X;
        allDataStruct.Mu = dataCP.Mu;
        allDataStruct.SCP = dataCP.SCP;
        allDataStruct.shannonCP = dataCP.shannonCP;
        allDataStruct.klCP = dataCP.klCP;
        allDataStruct.entropyCP = dataCP.entropyCP;
        allDataStruct.predictionErrorOnX = dataCP.predictionErrorOnX; %model obj PE
        allDataStruct.perceptualErrorOnX = dataCP.perceptualErrorOnX; 
        allDataStruct.predErrorCP = dataCP.predErrorCP; %model sub PE
        allDataStruct.muUpdate = dataCP.muUpdate;
        allDataStruct.posteriorPerceptCP = dataCP.posteriorPerceptCP;
        allDataStruct.posteriorCP = dataCP.posteriorCP;
        allDataStruct.expectationX = dataCP.expectationX;
        allDataStruct.maxLikePostMu = dataCP.maxLikePostMu;
        allDataStruct.maxLikePostX = dataCP.maxLikePostX;
        allDataStruct.surpriseCP = dataCP.surpriseCP;
        
        
        
        dataOB = oddBall_Inference(vars);
        
        allDataStruct.nTrialsOB = dataOB.nTrialsOB;
        allDataStruct.YOB = dataOB.YOB;
        allDataStruct.B = dataOB.B;
        allDataStruct.C = dataOB.C;
        allDataStruct.SOB = dataOB.SOB;
        allDataStruct.shannonOB = dataOB.shannonOB;
        allDataStruct.klOB = dataOB.klOB;
        allDataStruct.entropyOB = dataOB.entropyOB;
        allDataStruct.predictionErrorOnB = dataOB.predictionErrorOnB; %model obj PE
        allDataStruct.perceptualErrorOnB = dataOB.perceptualErrorOnB; 
        allDataStruct.predErrorOB = dataOB.predErrorOB; %model sub PE
        allDataStruct.cUpdate = dataOB.cUpdate;
        allDataStruct.posteriorPerceptOB = dataOB.posteriorPerceptOB;
        allDataStruct.posteriorOB = dataOB.posteriorOB;
        allDataStruct.expectationB = dataOB.expectationB;
        allDataStruct.maxLikePostC = dataOB.maxLikePostC;
        allDataStruct.maxLikePostB = dataOB.maxLikePostB;
        allDataStruct.surpriseOB = dataOB.surpriseOB;
        
        
        % create matrices to store surprise and entropy associated with
        % each observation over numReps "random" simulations of the model
        % -- each differing in the exact internal noise experienced on each
        % trial.
        allCPEntropy=[];
        allCPSurprise=[];
        allOBEntropy=[];
        allOBSurprise=[];
        allCPobjPE=[];
        allOBobjPE=[];
        allCPsubPE=[];
        allOBsubPE=[];
        
        
        if runLoop
            for i = 1:numReps
                % Get surprise and entropy for changepoint trials:
                dataCP = learningAndPerception(vars);
                allCPEntropy=[allCPEntropy; dataCP.entropyCP];
                allCPSurprise=[allCPSurprise; dataCP.surpriseCP];
                allCPobjPE=[allCPobjPE;dataCP.predictionErrorOnX];
                allCPsubPE=[allCPsubPE;dataCP.predErrorCP];
                
                % Get surprise and entropy for oddball trials:
                dataOB = oddBall_Inference(vars);
                allOBEntropy=[allOBEntropy; dataOB.entropyOB];
                allOBSurprise=[allOBSurprise; dataOB.surpriseOB];
                allOBobjPE=[allOBobjPE;dataOB.predictionErrorOnB];
                allOBsubPE=[allOBsubPE;dataOB.predErrorOB];
                
            end
        end
        allDataStruct.entropyCP=mean(allCPEntropy, 1);
        allDataStruct.surpriseCP=mean(allCPSurprise,1);
        allDataStruct.entropyOB=mean(allOBEntropy,1);
        allDataStruct.surpriseOB=mean(allOBSurprise,1);
        allDataStruct.subNum=nan(size(allDataStruct.isRand));
        allDataStruct.subNum(:)=subNum;
        allDataStruct.blockCond=nan(size(allDataStruct.isRand));
        allDataStruct.predictionErrorOnX=circ_mean(allCPobjPE);
        allDataStruct.predictionErrorOnB=circ_mean(allOBobjPE);
        allDataStruct.subPredErrorCP=circ_mean(allCPsubPE);
        allDataStruct.subPredErrorOB=circ_mean(allOBsubPE);
        
        t=length(allDataStruct.isRand);
        if allDataStruct.condition(1)==1
            allDataStruct.blockCond(1:t/2)=1;
            allDataStruct.blockCond(t/2+1:end)=-1;
        else
            allDataStruct.blockCond(1:t/2)=-1;
            allDataStruct.blockCond(t/2+1:end)=1;
        end
        
        %
        %     hold on
        %     plot(dataCP.Mu, '--k')
        %     plot(dataCP.X, 'or', 'markerFaceColor', 'r', 'lineWidth', 1, 'markerEdgeColor', 'k', 'markerSize', 6)
        %     hold on
        %     plot([dataCP.presentedColor(91:180,1); dataCP.presentedColor(91:180,2)], 'or', 'markerFaceColor', 'b', 'lineWidth', 1, 'markerEdgeColor', 'k', 'markerSize', 6)
        %     plot([dataCP.generativeMean(91:180,1); dataCP.generativeMean(91:180,2)], '--c')
        %
        %
        %Oddball condition datap
        
        
    end
    
    if isempty(behaveAll)
        behaveAll=allDataStruct;
    else
        behaveAll=catBehav(allDataStruct,behaveAll);
    end
    
end
%%
if vars.generate_data
    dataCP=straightStruct(dataCP);
    dataOB=straightStruct(dataOB);
else
    
    allDataStruct=straightStruct(allDataStruct);
    behaveAll=straightStruct(behaveAll);
    
    behaveAll.notMove=circ_dist(deg2rad(behaveAll.est),deg2rad(behaveAll.startPoint));
end



%%
paramsCircUpdateAll=[];
paramsCircBiasAll=[];


paramsUpdateAllCos=[];
paramsUpdateAllSin=[];

paramsBiasAllCos=[];
paramsBiasAllSin=[];

subPEall=[];
objPEall=[];

updateAll=[];
biasAll=[];

combineAllSubs=0;
makePlots=0;

%%
if ~combineAllSubs
    for sub=1:length(subList)
        
        %%
        subNum=subList(sub);
        subStr=num2str(subNum);
        % Get useful variables for descriptive analyses:
        
        decomposeSine=0;
        decomposeCosine=0;
        
        
        if vars.generate_data % use model data if generate_data is true:
            disp(8888)
            nTrials = length(dataCP.X);
            OBpredErrorSub = dataOB.predErrorOB;
            CPpredErrorSub = dataCP.predErrorCP;
            OBsurprise = dataOB.surpriseOB;
            CPsurprise = dataCP.surpriseCP;
            OBupdate = dataOB.cUpdate;
            CPupdate = dataCP.muUpdate;
            OBentropy = dataOB.entropyOB;
            CPentropy = dataCP.entropyCP;
            predictionErrorOnB = dataOB.predictionErrorOnB;
            predictionErrorOnX = dataCP.predictionErrorOnX;
            perceptualErrorOnX = dataCP.perceptualErrorOnX;
            perceptualErrorOnB = dataOB.perceptualErrorOnB;
            OBpredError=dataOB.objPredErrorOB;
            CPpredError=dataCP.objPredErrorCP;
            
            
            
            
            
        else % otherwise, use participant data
            nTrials = length(behaveAll.block(behaveAll.subNum==subNum))*.5;
            
            OBpredError = [behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
            %OBpredError(isnan(OBpredError))=0;
            CPpredError = [behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
            %CPpredError(isnan(CPpredError))=0;
            
            OBpredErrorSub = [behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
            %OBpredErrorSub(isnan(OBpredErrorSub))=0;
            CPpredErrorSub = [behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
            %CPpredErrorSub(isnan(CPpredErrorSub))=0;
            
            
            OBsurprise = behaveAll.surpriseOB(behaveAll.subNum==subNum);
            CPsurprise = behaveAll.surpriseCP(behaveAll.subNum==subNum);
            
            OBupdate = [behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
            %OBupdate(isnan(OBupdate))=0;
            CPupdate = [behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
            %CPupdate(isnan(CPupdate))=0;
            
            
            OBentropy = behaveAll.entropyOB(behaveAll.subNum==subNum);
            CPentropy = behaveAll.entropyCP(behaveAll.subNum==subNum);
            
            perceptualErrorOnX = [behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
            %perceptualErrorOnX(isnan(perceptualErrorOnX))=0;
            perceptualErrorOnB = [behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
            %perceptualErrorOnB(isnan(perceptualErrorOnB))=0;
            
            condition=behaveAll.condition(sub*2);
            
            conNum=behaveAll.blockCond(behaveAll.subNum==subNum);
            condNum=[conNum(1:nTrials);conNum(1:nTrials);conNum(nTrials+1:end);conNum(nTrials+1:end)];
            
            predictionErrorOnX=[behaveAll.predictionErrorOnX(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
            %predictionErrorOnX(isnan(predictionErrorOnX))=0;
            predictionErrorOnB=[behaveAll.predictionErrorOnB(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
            %predictionErrorOnB(isnan(predictionErrorOnB))=0;
            
            
            
            OBsurpriseTrialNum=sum(nansum(behaveAll.surpriseTrial(behaveAll.subNum==subNum & behaveAll.blockCond==-1,:)));
            CPsurpriseTrialNum=sum(nansum(behaveAll.surpriseTrial(behaveAll.subNum==subNum & behaveAll.blockCond==1,:)));
            
            OBmoveDist=[behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)];
            %OBmoveDist(isnan(OBmoveDist))=0;
            CPmoveDist=[behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.notMove(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)];
            %CPmoveDist(isnan(CPmoveDist))=0;
            
            
        end
        
        %     moveDistAll=[OBmoveDist;CPmoveDist];
        %     figure()
        %     histogram(moveDistAll,100)
        %     %keyboard
        
        %check how many surprise trials
        if decomposeSine
            CPpredError=sin(CPpredError);
            OBpredError=sin(OBpredError);
            perceptualErrorOnX=sin(perceptualErrorOnX);
            perceptualErrorOnB=sin(perceptualErrorOnB);
            OBupdate=sin(OBupdate);
            CPupdate=sin(CPupdate);
        elseif decomposeCosine
            CPpredError=cos(CPpredError);
            OBpredError=cos(OBpredError);
            perceptualErrorOnX=cos(perceptualErrorOnX);
            perceptualErrorOnB=cos(perceptualErrorOnB);
            OBupdate=cos(OBupdate);
            CPupdate=cos(CPupdate);
            
        end
        %     figure()
        %     plot([OBpredError;CPpredError],[predictionErrorOnB;predictionErrorOnX],'x')
        %     hold on
        %     plot([-4,4],[-4,4],'-')
        %     xlabel('subject')
        %     ylabel('model')
        %     hold off
        
        
        
        
        %     t=find(CPsurprise>0.95);
        %     for i=1:length(t)
        %         if t(i)>90
        %             t(i)=t(i)-90;
        %         end
        %     end
        
        subjectPEobj=[OBpredError;CPpredError];
        modelPEobj=[predictionErrorOnB;predictionErrorOnX];
        
        subModPEobjDiff=modelPEobj-subjectPEobj;
        
        %histogram(subModPEobjDiff)
        
        
        data.signedError=subModPEobjDiff;
        data.signedError_allTargs=[];
        
        %startPoint=[2, .3,  0, .5, 0, 0, 0, 0];
        whichParams=[1 1 0 0];
        data.doFit=true;
        
        
        
        % data.signedError  --> signed errors (or just angles if not for VWM task)
        % data.signedError_allTargs --> errors computed as if subject were
        % estimating all of the colors in the array (ie colorArray - subject
        % response).
        % data.xMat        --> all variables that could affect recall or precision
        % data.doFit  --> fit? if not, just evaluate at startPoint
        
        %params: maximum  likelihood model parameters,
        % 1= proportion gaussian,
        % 2= concentration (ie 1./sigma^2) of von mises,
        % 3= mean of gaussian (should be zero... but who knows!)
        % 4= proportion binding error (ie propr of gaussian that are evenly distributed across targets)
        
        % simplex order: gaussian, binding error, uniform
        
        
        
        
        [X, negLogLike, simplexP, uniformProb]=fit_VWM_mixtureModelTrial(data,whichParams);
        
        
        
        %%
        
%         if condition==1
%             PE=[CPpredErrorSub;OBpredErrorSub];
%             StateChangeProbUpdate=[CPsurprise;OBsurprise];
%             
%             data.Y=[CPupdate,OBupdate];
%             
%         else
%             PE=[OBpredErrorSub;CPpredErrorSub];
%             StateChangeProbUpdate=[OBsurprise; CPsurprise];
%             
%             data.Y=[OBupdate, CPupdate];
%         end
%         
%         figure()
%         hold on
%         a=plot(PE(condNum==1), data.Y(condNum==1), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 1, 'markerSize', 8 )
%         c=plot(PE(condNum==1&StateChangeProbUpdate>.99), data.Y(condNum==1&StateChangeProbUpdate>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'c', 'lineWidth', 1, 'markerSize', 12 )
%         b=plot(PE(condNum==-1), data.Y(condNum==-1), 'or','markerSize', 8, 'markerEdgeColor', 'k','markerFaceColor', 'r', 'lineWidth', 1, 'markerSize', 8  )
%         d=plot(PE(condNum==-1&StateChangeProbUpdate>.99), data.Y(condNum==-1&StateChangeProbUpdate>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'y', 'lineWidth', 1, 'markerSize', 12 )
%         %e=plot(PE(condNum==1&uniformProb>0.5), data.Y(condNum==1&uniformProb>0.5), 'og', 'markerEdgeColor', 'k', 'markerFaceColor', 'g', 'lineWidth', 1, 'markerSize', 8 )
%         %f=plot(PE(condNum==-1&uniformProb>0.5), data.Y(condNum==-1&uniformProb>0.5), 'og', 'markerEdgeColor', 'k', 'markerFaceColor', 'm', 'lineWidth', 1, 'markerSize', 8 )
%         
%         
%         ff=legend([a, b,c,d], 'changepoint', 'oddball','changepoint','oddball')
%         set(ff, 'location', 'northwest', 'box', 'off')
%         %ylim([-2.5 2.5])
%         xlabel('subjective PE')
%         ylabel('')
%         
%         normCoef=normalizedCoefUpdate(sub,3);
%         
%         titleForPlot=sprintf('subject %s normalized learning coef %s',subStr,num2str(normCoef));
%         title(titleForPlot)
%         
%         hold off
        
        
        
        %% Data check, make sure things were extracted properly
        
        
        %     %first plot out all stims based on OB/CP raw data in 4 quadrants
        %     figure()
        %
        %     subplot(2,2,1)
        %     plot(behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1),behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1),'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     %xlim([-1 1])
        %     title('oddball stim1')
        %
        %     subplot(2,2,2)
        %     plot(behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2),behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2),'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     %xlim([-1 1])
        %     title('oddball stim2')
        %
        %     subplot(2,2,3)
        %     plot(behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1),behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1),'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     %xlim([-1 1])
        %     title('changepoint stim1')
        %
        %     subplot(2,2,4)
        %     plot(behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2),behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2),'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     %xlim([-1 1])
        %     title('changepoint stim2')
        %
        %     % seems like we see positive relationship when oPE small in all
        %
        %     %then check if we're grabbing the right values
        %
        %     figure()
        %
        %     subplot(2,2,1)
        %     plot(OBpredError,perceptualErrorOnB,'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     title('OB extracted')
        %
        %     subplot(2,2,2)
        %     plot(CPpredError,perceptualErrorOnX,'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     title('CP extracted')
        %
        %
        %     subplot(2,2,3)
        %     plot([behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)],[behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2)],'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     title('OB combined raw data')
        %
        %     subplot(2,2,4)
        %     plot([behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.predictErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)],[behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);behaveAll.estErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2)],'x')
        %     xlabel('objective prediction error')
        %     ylabel('perceptual error')
        %     title('CP combined raw data')
        %
        %
        %
        %
        %prediction update vs sPE
%             figure()
%         
%         
%             subplot(2,2,1)
%         
%             hold on
%             plot(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1),behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1),'x')
%             [r,c]=find(abs(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1))>2)
%             x=behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);
%             y=behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,1);
%             plot(x(r),y(r),'rx')
%             lsline
%             xlabel('subjective prediction error')
%             ylabel('prediction update')
%             title('oddball stim1')
%             hold off
%         
%             subplot(2,2,2)
%         
%             hold on
%             plot(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2),behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2),'x')
%             [r,c]=find(abs(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2))>2)
%             x=behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2);
%             y=behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==-1,2);
%             plot(x(r),y(r),'rx')
%             lsline
%             xlabel('subjective prediction error')
%             ylabel('prediction update')
%             title('oddball stim2')
%             hold off
%         
%             subplot(2,2,3)
%         
%             hold on
%             plot(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1),behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,1),'x')
%             [r,c]=find(abs(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1))>2)
%             x=behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);
%             y=behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,1);
%             plot(x(r),y(r),'rx')
%             lsline
%             xlabel('subjective prediction error')
%             ylabel('prediction update')
%             title('changepoint stim1')
%             hold off
%         
%             subplot(2,2,4)
%         
%             hold on
%             plot(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2),behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,2),'x')
%             [r,c]=find(abs(behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2))>2)
%             x=behaveAll.subPredErr(behaveAll.subNum==subNum & behaveAll.blockCond==1,2);
%             y=behaveAll.predUpdate(behaveAll.subNum==subNum & behaveAll.blockCond==1,2);
%             plot(x(r),y(r),'rx')
%             lsline
%             xlabel('subjective prediction error')
%             ylabel('prediction update')
%             title('changepoint stim2')
%             hold off
        
        
        
        
        %%
        
        % ignore this section for now
%         if vars.average && vars.generate_data
%             for i = 1:vars.average_number-1
%                 dataCP = learningAndPerception(vars);
%                 dataOB = oddBall_Inference(vars);
%                 close all;
%                 OBpredError = OBpredError + dataOB.predError;
%                 CPpredError = CPpredError + dataCP.predError;
%                 OBkl = OBkl + dataOB.kl;
%                 CPkl = CPkl + dataCP.kl;
%                 OBsurprise = OBsurprise + dataOB.surprise;
%                 CPsurprise = CPsurprise + dataCP.surprise;
%                 OBupdate = OBupdate+  dataOB.cUpdate;
%                 CPupdate = CPupdate + dataCP.muUpdate;
%                 predictionErrorOnB = predictionErrorOnB + dataOB.predictionErrorOnB;
%                 predictionErrorOnX = predictionErrorOnX + dataCP.predictionErrorOnX;
%                 perceptualErrorOnX = perceptualErrorOnX + dataCP.perceptualErrorOnX;
%                 perceptualErrorOnB = perceptualErrorOnB + dataOB.perceptualErrorOnB;
%                 
%             end
%             OBpredError = OBpredError / vars.average_number;
%             CPpredError = CPpredError / vars.average_number;
%             OBkl = OBkl / vars.average_number;
%             CPkl = CPkl / vars.average_number;
%             OBsurprise = OBsurprise / vars.average_number;
%             CPsurprise = CPsurprise / vars.average_number;
%             OBupdate = OBupdate/ vars.average_number;
%             CPupdate = CPupdate / vars.average_number;
%             predictionErrorOnB = predictionErrorOnB / vars.average_number;
%             predictionErrorOnX = predictionErrorOnX / vars.average_number;
%             perceptualErrorOnX = perceptualErrorOnX / vars.average_number;
%             perceptualErrorOnB = perceptualErrorOnB / vars.average_number;
%         end
%         
%         
       
        %% Regression for the Update
        % medians = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_medians", suffix, '.mat'], ''));
        % medians_kl = medians.medians_kl;
        % params = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_params", suffix, '.mat'], ''));
        % params_kl = params.params_kl;
        
        % condNumCP=ones(length(CPsurprise),1);
        % condNumOB=ones(length(OBsurprise),1)*-1;
        %
        % condNum=[condNumOB;condNumCP];
        
        % xes = [ ones(4*nTrials,1), ...
        %   ([CPpredError;OBpredError]),...
        %   ([CPpredError;OBpredError] .* nanzscore([CPsurprise; OBsurprise].* condNum)),...
        %   ([CPpredError;OBpredError] .* nanzscore([CPsurprise; OBsurprise])),...
        %   ([CPpredError;OBpredError] .* condNum), ...
        %   ([CPpredError;OBpredError].* nanzscore([CPentropy;OBentropy]))];
        
        
        if vars.generate_data
            
            condNum=ones(2*nTrials,1);
            condNum(((length(condNum)/2)+1):end)=-1;
            
            xes = [ ones(2*nTrials,1), ...
                ([CPpredErrorSub;OBpredErrorSub]),...
                ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
                ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise])),...
                ([CPpredErrorSub;OBpredErrorSub] .* condNum), ...
                ([CPpredErrorSub;OBpredErrorSub].* nanzscore([CPentropy;OBentropy]))];%, ...
                %([CPpredErrorSub;OBpredErrorSub].* nanzscore(uniformProb))];

%              xes = [ ones(2*nTrials,1), ...
%                 ([CPpredErrorSub;OBpredErrorSub]),...
%                 ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
%                 ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise])),...
%                 ([CPpredErrorSub;OBpredErrorSub] .* condNum), ...
%                 ([CPpredErrorSub;OBpredErrorSub].* nanzscore([CPentropy;OBentropy])), ...
%                ];
%             
            Y=[ CPupdate;OBupdate];
            condition=1;
            [b_update, bint_update] = regress(Y, xes);
                
            
        else
            
           % keyboard
            if condition==1
                xes = [ ones(4*nTrials,1), ...
                    ([CPpredErrorSub;OBpredErrorSub]),...
                    ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
                    ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise])),...
                    ([CPpredErrorSub;OBpredErrorSub] .* condNum), ...
                    ([CPpredErrorSub;OBpredErrorSub].* nanzscore([CPentropy;OBentropy])), ...
                    ([CPpredErrorSub;OBpredErrorSub].* nanzscore(uniformProb))];

%                 xes = [ ones(4*nTrials,1), ...
%                 ([CPpredErrorSub;OBpredErrorSub]),...
%                 ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
%                 ([CPpredErrorSub;OBpredErrorSub] .* nanzscore([CPsurprise;OBsurprise])),...
%                 ([CPpredErrorSub;OBpredErrorSub] .* condNum), ...
%                 ([CPpredErrorSub;OBpredErrorSub].* nanzscore([CPentropy;OBentropy]))
%                ];
%                 
                Y=[ CPupdate; OBupdate];
                
            else
                xes = [ ones(4*nTrials,1), ...
                    ([OBpredErrorSub;CPpredErrorSub]),...
                    ([OBpredErrorSub;CPpredErrorSub] .* nanzscore([OBsurprise; CPsurprise].* condNum)),...
                    ([OBpredErrorSub;CPpredErrorSub] .* nanzscore([OBsurprise; CPsurprise])),...
                    ([OBpredErrorSub;CPpredErrorSub] .* condNum), ...
                    ([OBpredErrorSub;CPpredErrorSub].* nanzscore([OBentropy; CPentropy])), ...
                    ([OBpredErrorSub;CPpredErrorSub].* nanzscore(uniformProb))];

%                 xes = [ ones(4*nTrials,1), ...
%                 ([OBpredErrorSub;CPpredErrorSub]),...
%                 ([OBpredErrorSub;CPpredErrorSub] .* nanzscore([OBsurprise; CPsurprise].* condNum)),...
%                 ([OBpredErrorSub;CPpredErrorSub] .* nanzscore([OBsurprise; CPsurprise])),...
%                 ([OBpredErrorSub;CPpredErrorSub] .* condNum), ...
%                 ([OBpredErrorSub;CPpredErrorSub].* nanzscore([OBentropy; CPentropy])), ...
%                ];
%                 
                Y=[OBupdate, CPupdate];

            end
        end
        
        
        
        if ~decomposeSine && ~decomposeCosine
            
            regData.X=xes;
            regData.Y=Y;
          
            regData.startPoint=[5,zeros(1,size(regData.X,2)),0.05,0];
            %concentration, params, percent random
            regData.whichParams=logical([1,1,1,1,1,1,1,1,0,0]);
            %regData.whichParams=logical([1,1,1,1,1,1,1,0,0]);%no guess param
            regData.includeUniform=1;
            regData.doVarVar=1;
            
            regData.priorMean=[0,0,0,0,0,0,0];
            %regData.priorMean=[0,0,0,0,0,0];%no guess param
            
            %regData.priorWidth=[1,1,.01,.01,.01,.01]; %no guess param
            regData.priorWidth=[1,1,.01,.01,.01,.01,.01];
            
%             data.lb=[.0001, ones(1, size(data.X, 2)).*-.5,0,0];
%             data.ub=[100, ones(1, size(data.X, 2)).*1,1,1];
%             
            regData.lb=[.0001, ones(1, size(regData.X, 2)).*-.5,0,0];
            regData.ub=[100, ones(1, size(regData.X, 2)).*1,1,1];
            
            
            regData.Y(isnan(Y))=[];
            regData.X(isnan(Y),:)=[];
            
            %     data.X=xes;
            %     data.startPoint=[5,zeros(1,size(data.X,2)),0.05];
            %     %concentration, params, percent random
            %     data.whichParams=logical([1,1,1,1,1,1,1,0]);
            %     data.includeUniform=1;
            %
            %     data.priorMean=0;
            %     data.priorWidth=.01;
            %
            %     data.lb=[.0001, ones(1, size(data.X, 2)).*-.5];
            %     data.ub=[100, ones(1, size(data.X, 2)).*1];
            %
            %
            %     %plot([OBpredError;CPpredError], data.Y, '.')
            
            
            % Now lets fit circular regression in a loop -- trying different start
            % points:
%             if sub ==2
%                 keyboard
%             end
            
            nr=20;
            bestNegLogLikeUpdate=inf;
            storeParams=[];
            storeNegLogLike=[];
            for i = 1:nr
                
                
                [paramsUpdate, negLogLikeUpdate]=fitLinearModWCircErrsNewVar(regData);
                storeParams=[storeParams; paramsUpdate];
                storeNegLogLike=[storeNegLogLike;negLogLikeUpdate];
                
                if negLogLikeUpdate<bestNegLogLikeUpdate
                    bestParamsUpdate=paramsUpdate;
                    bestNegLogLikeUpdate=negLogLikeUpdate;
                end
                
                %  LRs= bestParams(3)+nanzscore([OBsurprise; CPsurprise].* condNum).*bestParams(4)+bestParams(5).*nanzscore([OBsurprise; CPsurprise])
                
                %data.startPoint(1:end-1) = data.lb' + rand(length(data.lb),1).*data.ub';
                regData.startPoint(regData.whichParams) = regData.lb(regData.whichParams)' + rand(length(regData.lb(regData.whichParams)),1).*regData.ub(regData.whichParams)';
                
            end
            
        end
        
        
        %%%%%%%%%%%%%
%         if ~decomposeSine&~decomposeCosine
            
            %             subPE=abs(data.X(:,2));
            %             Q=quantile(subPE,(1:9)/10);
            %
            %             Qs=nan(1,11);
            %             Qs(2:10)=Q;
            %             Qs(1)=min(subPE);
            %             Qs(end)=max(subPE);
            %
            %             for i=1:length(Qs)-1
            %
            %                 Qdata=abs(data.Y(subPE>Qs(i)&subPE<Qs(i+1)));
            %                 meanQs(i)=circ_mean(Qdata);
            %                 semQs(i)=std(Qdata)./sqrt(length(Qdata));
            %                 subPEbins(i)=mean([Qs(i) Qs(i+1)]);
            %
            %             end
            %
            %             figure()
            %             errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
            %             hold on
            %             xL=xlabel('subjective PE');
            %             xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
            %             ylabel('mean prediction update')
            %             xlim([0 3.2])
            %             xticks(Qs)
            %             xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
            %             xtickangle(45)
            %             hold off
            
            
            
            
            %%%%%%%%%%%%%%
            
            %             subPEall=cat(1,subPEall,data.X(:,2));
            %             updateAll=cat(1,updateAll,data.Y);
            %
            %
            %             subPE=data.X(:,2);
            %             Q=quantile(subPE,(1:9)/10);
            %
            %             Qs=nan(1,11);
            %             Qs(2:10)=Q;
            %             Qs(1)=min(subPE);
            %             Qs(end)=max(subPE);
            %
            %
            %             for i=1:length(Qs)-1
            %
            %                 Qdata=data.Y(subPE>Qs(i)&subPE<Qs(i+1));
            %                 meanQs(i)=circ_mean(Qdata);
            %                 semQs(i)=std(Qdata)./sqrt(length(Qdata));
            %                 subPEbins(i)=mean([Qs(i) Qs(i+1)]);
            %
            %             end
            %
            %             figure()
            %             errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
            %             hold on
            %             xL=xlabel('subjective PE');
            %             xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
            %             ylabel('mean prediction update')
            %             xlim([-3.2 3.2])
            %             xticks(Qs)
            %             xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
            %             xtickangle(45)
            %             hold off
            
            %storeParams
            %bestParamsBias
            
            
            
            
            
            
            %%%%%%%%%%%
            %OB only
            
%             subPE=OBpredErrorSub;
%             Q=quantile(subPE,(1:9)/10);
%             
%             Qs=nan(1,11);
%             Qs(2:10)=Q;
%             Qs(1)=min(subPE);
%             Qs(end)=max(subPE);
%             
%             
%             for i=1:length(Qs)-1
%                 
%                 Qdata=OBupdate(subPE>Qs(i)&subPE<Qs(i+1));
%                 meanQs(i)=circ_mean(Qdata);
%                 semQs(i)=std(Qdata)./sqrt(length(Qdata));
%                 subPEbins(i)=mean([Qs(i) Qs(i+1)]);
%                 
%             end
%             
%             figure()
%             errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
%             hold on
%             xL=xlabel('OB subjective PE');
%             xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
%             ylabel('OB mean prediction update')
%             yline(0)
%             refline(1,0)
%             xlim([-3.2 3.2])
%             xticks(Qs)
%             xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
%             xtickangle(45)
%             hold off
%             
%             %%%%%%%
%             %CP only
%             
%             subPE=CPpredErrorSub;
%             Q=quantile(subPE,(1:9)/10);
%             
%             Qs=nan(1,11);
%             Qs(2:10)=Q;
%             Qs(1)=min(subPE);
%             Qs(end)=max(subPE);
%             
%             
%             for i=1:length(Qs)-1
%                 
%                 Qdata=CPupdate(subPE>Qs(i)&subPE<Qs(i+1));
%                 meanQs(i)=circ_mean(Qdata);
%                 semQs(i)=std(Qdata)./sqrt(length(Qdata));
%                 subPEbins(i)=mean([Qs(i) Qs(i+1)]);
%                 
%             end
%             
%             figure()
%             errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
%             hold on
%             xL=xlabel('CP subjective PE');
%             xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
%             ylabel('CP mean prediction update')
%             yline(0)
%             refline(1,0)
%             xlim([-3.2 3.2])
%             xticks(Qs)
%             xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
%             xtickangle(45)
%             hold off
%             
%             
%         end
%         %%%%%%%%%%%%%%%
        
        % sub6 2007 has good update data?
        
        
        %monday 7/19
        %all sub data together create update and perceptual error binned plots
        %for CP and OB (2by2)
        
        %play with model data, change noise/drift rate in the task and model
        %parameters
        %OB higher drift, lower noise (will see larger learning when normal
        %trial)
        %CP higher noise (will see less learning normal trial)
        
        %run myself to see how behavior looks with conscious effort to create
        %the behavior
        
        %set up a variant of the task,including the changes above
        %change instruction to clarify structure,give instruction of OB right before OB and CP right before CP, add a training stage before
        %the OB/CP blocks
        %2 sessions, will help with EEG and eye data
        
        %3 sub 2 sessions each (6 sessions)
        
        
        
        % xes = [ ones(2*nTrials,1), ...
        %   ([CPpredError';OBpredError']),...
        %   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'].* condNum)),...
        %   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'])),...
        %   ([CPpredError';OBpredError'] .* condNum), ([CPpredError';OBpredError'].* nanzscore([CPentropy';OBentropy']))];
        %
        
        if decomposeSine
            [b_updateSin, bint_updateSin] = regress(Y, xes);
        elseif decomposeCosine
            [b_updateCos, bint_updateCos] = regress(Y, xes);
            
        end
        
        
        
        % change to other circular regression
        
%         if condition==1
%             PE=[CPpredErrorSub;OBpredErrorSub];
%             PE(isnan(Y))=[];
%             StateChangeProbUpdate=[ CPsurprise;OBsurprise];
%             StateChangeProbUpdate(isnan(Y))=[];
%             condnum=condNum;
%             condnum(isnan(Y))=[];
%         else
%             PE=[OBpredErrorSub;CPpredErrorSub];
%             PE(isnan(Y))=[];
%             StateChangeProbUpdate=[OBsurprise; CPsurprise];
%             StateChangeProbUpdate(isnan(Y))=[];
%             condnum=condNum;
%             condnum(isnan(Y))=[];
%         end
%         
%         OBpredErrorStim1=OBpredErrorSub(1:(length(OBpredErrorSub)./2));
%         OBpredErrorStim2=OBpredErrorSub(((length(OBpredErrorSub)./2)+1):end);
%         
%         CPpredErrorStim1=CPpredErrorSub(1:(length(CPpredErrorSub)./2));
%         CPpredErrorStim2=CPpredErrorSub(((length(CPpredErrorSub)./2)+1):end);
%         
%         
%         [OBmysteryTstim1,c]=find(abs(OBpredErrorStim1)>1.5);
%         [OBmysteryTstim2,c]=find(abs(OBpredErrorStim2)>1.5);
%         
%         [CPmysteryTstim1,c]=find(abs(CPpredErrorStim1)>1.5);
%         [CPmysteryTstim2,c]=find(abs(CPpredErrorStim2)>1.5);
%         
%         [CPstrangeTStim2,c]=find(round(CPpredErrorStim2,4)==2.3736);
%         [OBstrangeTStim1,c]=find(round(OBpredErrorStim1,4)==1.0297);
%         [OBstrangeTStim2,c]=find(round(OBpredErrorStim2,4)==1.0297);
%         
%         if makePlots
%             figure()
%             hold on
%             a=plot(PE(condnum==1), data.Y(condnum==1), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 1, 'markerSize', 8 )
%             c=plot(PE(condnum==1&StateChangeProbUpdate>.5), data.Y(condnum==1&StateChangeProbUpdate>.5), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'c', 'lineWidth', 1, 'markerSize', 12 )
%             b=plot(PE(condnum==-1), data.Y(condnum==-1), 'or','markerSize', 8, 'markerEdgeColor', 'k','markerFaceColor', 'r', 'lineWidth', 1, 'markerSize', 8  )
%             d=plot(PE(condnum==-1&StateChangeProbUpdate>.5), data.Y(condnum==-1&StateChangeProbUpdate>.5), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'y', 'lineWidth', 1, 'markerSize', 12 )
%             ff=legend([a, b], 'changepoint', 'oddball')
%             set(ff, 'location', 'northwest', 'box', 'off')
%             %ylim([-2.5 2.5])
%             xlabel('subjective PE')
%             ylabel('update')
%             
%             refline(1,0)
%             yline(0)
%             
%             hold off
%             
%         end
        
        %% Regression for Perceptual Error
        % medians = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_medians", suffix, '.mat'], ''));
        % medians_kl = medians.medians_kl;
        % params = load(join(["/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/neuro/lookups/regular_kl_params", suffix, '.mat'], ''));
        % params_kl = params.params_kl;
        
        % xes = [ ones(4*nTrials,1), ...
        %   ([CPpredErrorObj;OBpredErrorObj]),...
        %   ([CPpredErrorObj;OBpredErrorObj] .* nanzscore([CPsurprise; OBsurprise].* condNum)),...
        %   ([CPpredErrorObj;OBpredErrorObj] .* nanzscore([CPsurprise; OBsurprise])),...
        %   ([CPpredErrorObj;OBpredErrorObj] .* condNum), ([CPpredErrorObj;OBpredErrorObj].* nanzscore([CPentropy;OBentropy]))];
        
        %dropTrials
        
        
        %take the sine and cosine of PE and perceptual error
        
        
        %7/16
        
        %try to add new terms in the model regressions
        %add noise to model data, how much noise can we still see the
        %effect
        
        % do the bin plots for OB/CP for learning
        
        if vars.generate_data
            
            condNum=ones(2*nTrials,1);
            condNum(((length(condNum)/2)+1):end)=-1;
            
            xes = [ ones(2*nTrials,1), ...
                ([CPpredError;OBpredError]),...
                ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
                ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise])),...
                ([CPpredError;OBpredError] .* condNum), ...
                ([CPpredError;OBpredError].* nanzscore([CPentropy;OBentropy]))]%, ...
                %([CPpredError;OBpredError].* nanzscore(uniformProb))];

%               xes = [ ones(2*nTrials,1), ...
%                 ([CPpredError;OBpredError]),...
%                 ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
%                 ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise])),...
%                 ([CPpredError;OBpredError] .* condNum), ...
%                 ([CPpredError;OBpredError].* nanzscore([CPentropy;OBentropy])), ...
%                 ];
%             
            
            
            
            Y=[ perceptualErrorOnX; perceptualErrorOnB];
            condition=1
            [b_bias, bint_bias] = regress(Y, xes);
            
        else
            if condition==1
                xes = [ ones(4*nTrials,1), ...
                    ([CPpredError;OBpredError]),...
                    ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
                    ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise])),...
                    ([CPpredError;OBpredError] .* condNum), ...
                    ([CPpredError;OBpredError].* nanzscore([CPentropy;OBentropy])),...
                    ([CPpredError;OBpredError].* nanzscore(uniformProb))];
                

%                 xes = [ ones(4*nTrials,1), ...
%                 ([CPpredError;OBpredError]),...
%                 ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise].* condNum)),...
%                 ([CPpredError;OBpredError] .* nanzscore([CPsurprise;OBsurprise])),...
%                 ([CPpredError;OBpredError] .* condNum), ...
%                 ([CPpredError;OBpredError].* nanzscore([CPentropy;OBentropy])), ...
%                 ];
            
                
                
                Y=[perceptualErrorOnX; perceptualErrorOnB];
                
            else
                xes = [ ones(4*nTrials,1), ...
                    ([OBpredError;CPpredError]),...
                    ([OBpredError;CPpredError] .* nanzscore([OBsurprise; CPsurprise].* condNum)),...
                    ([OBpredError;CPpredError] .* nanzscore([OBsurprise; CPsurprise])),...
                    ([OBpredError;CPpredError] .* condNum), ...
                    ([OBpredError;CPpredError].* nanzscore([OBentropy; CPentropy])),...
                    ([OBpredError;CPpredError].* nanzscore(uniformProb))];

%                 xes = [ ones(4*nTrials,1), ...
%                 ([OBpredError;CPpredError]),...
%                 ([OBpredError;CPpredError] .* nanzscore([OBsurprise; CPsurprise].* condNum)),...
%                 ([OBpredError;CPpredError] .* nanzscore([OBsurprise; CPsurprise])),...
%                 ([OBpredError;CPpredError] .* condNum), ...
%                 ([OBpredError;CPpredError].* nanzscore([OBentropy; CPentropy])), ...
%                 ];
            
                
                
                Y=[perceptualErrorOnB; perceptualErrorOnX];
            end
        end
%         % figure()
%         % plot([CPpredError;OBpredError],nanzscore([CPsurprise;OBsurprise]),'x')
%         
%         
%         
%         
        if ~decomposeSine && ~decomposeCosine
            regData.X=xes;
            regData.Y=Y;
            regData.startPoint=[5,zeros(1,size(regData.X,2)),0.05,0];
            %regData.whichParams=logical([1,1,1,1,1,1,1,0,0]);%no guess param
            regData.whichParams=logical([1,1,1,1,1,1,1,1,0,0]);
            % 9: guess rate
            % 10: variance term
            regData.includeUniform=1;
            regData.doVarVar=1;
            
            regData.priorMean=[0,0,0,0,0,0,0];
            regData.priorWidth=[1,1,.01,.01,.01,.01,.01];
            
             
%             regData.priorMean=[0,0,0,0,0,0];%no guess param
%             regData.priorWidth=[1,1,.01,.01,.01,.01];%no guess param
% %             
            regData.lb=[.0001, ones(1, size(regData.X, 2)).*-.5,0,0];
            regData.ub=[100, ones(1, size(regData.X, 2)).*1,1,1];
            
            
            [nanR,c]=find(isnan(xes(:,2)));
            regData.X(nanR,:)=[];
            regData.Y(isnan(xes(:,2)))=[];
            
            
            %plot([OBpredError;CPpredError], data.Y, '.')
            
            
            % Now lets fit circular regression in a loop -- trying different start
            % points:
            
            %if ~decomposeSine&~decomposeCosine
            nr=20;
            bestNegLogLikeBias=inf;
            storeParams=[];
            for i = 1:nr
                %     data.Y=data.Y(abs(data.X(:,2))<1);
                %     data.X=data.X(abs(data.X(:,2))<1,:);
                [paramsBias, negLogLikeBias]=fitLinearModWCircErrsNewVar(regData);
                storeParams=[storeParams; paramsBias];
                
                if negLogLikeBias<bestNegLogLikeBias
                    bestParamsBias=paramsBias;
                    bestNegLogLikeBias=negLogLikeBias;
                end
                
                %  LRs= bestParams(3)+nanzscore([OBsurprise; CPsurprise].* condNum).*bestParams(4)+bestParams(5).*nanzscore([OBsurprise; CPsurprise])
                
                regData.startPoint(regData.whichParams) = regData.lb(regData.whichParams)' + rand(length(regData.lb(regData.whichParams)),1).*regData.ub(regData.whichParams)';
                
            end
        end
        
%         
%         %
%         %     Y_hat=bestParamsBias(2:end-1).*data.X;
%         %     err=abs(data.Y-Y_hat);
%         %     err(abs(data.X(:,2))>1)
%         %
%         %     % figure()
%         % plot(data.X(:,2),err,'x')
%         
%         
%         
%         %%%%%%%%%%%%%
%         if ~decomposeSine&~decomposeCosine
%             
%             %             objPEall=cat(1,objPEall,data.X(:,2));
%             %             biasAll=cat(1,biasAll,data.Y);
%             %
%             %             objPE=abs(data.X(:,2));
%             %             Q=quantile(objPE,(1:9)/10);
%             %
%             %             Qs=nan(1,11);
%             %             Qs(2:10)=Q;
%             %             Qs(1)=min(objPE);
%             %             Qs(end)=max(objPE);
%             %
%             %             for i=1:length(Qs)-1
%             %
%             %                 Qdata=abs(data.Y(objPE>Qs(i)&objPE<Qs(i+1)));
%             %                 meanQs(i)=circ_mean(Qdata);
%             %                 semQs(i)=std(Qdata)./sqrt(length(Qdata));
%             %                 objPEbins(i)=mean([Qs(i) Qs(i+1)]);
%             %
%             %             end
%             %
%             %             figure()
%             %             errorbar(objPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
%             %             hold on
%             %             xL=xlabel('objective PE');
%             %             xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
%             %             ylabel('mean perceptual error')
%             %             xlim([0 3.2])
%             %             xticks(Qs)
%             %             xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
%             %             xtickangle(45)
%             %             hold off
%             
%             
%             
%             
%             %%%%%%%%%%%%%%
%             
%             
%             objPE=data.X(:,2);
%             Q=quantile(objPE,(1:9)/10);
%             
%             Qs=nan(1,11);
%             Qs(2:10)=Q;
%             Qs(1)=min(objPE);
%             Qs(end)=max(objPE);
%             
%             
%             for i=1:length(Qs)-1
%                 
%                 Qdata=data.Y(objPE>Qs(i)&objPE<Qs(i+1));
%                 meanQs(i)=circ_mean(Qdata);
%                 semQs(i)=std(Qdata)./sqrt(length(Qdata));
%                 objPEbins(i)=mean([Qs(i) Qs(i+1)]);
%                 
%             end
%             
%             figure()
%             errorbar(objPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
%             hold on
%             xL=xlabel('objective PE');
%             xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
%             ylabel('mean perceptual error')
%             xlim([-3.2 3.2])
%             xticks(Qs)
%             xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
%             xtickangle(45)
%             hold off
%             
%             %storeParams
%             %bestParamsBias
%         end
%         %%%%%%%%%%%%%%%
%         
%         
%         %  xes = [ ones(2*nTrials,1), ...
%         %   ([CPpredError';OBpredError']),...
%         %   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'].* condNum)),...
%         %   ([CPpredError';OBpredError'] .* nanzscore([CPsurprise'; OBsurprise'])),...
%         %   ([CPpredError';OBpredError'] .* condNum), ([CPpredError';OBpredError'].* nanzscore([CPentropy';OBentropy']))];
%         
%         %         if decomposeSine||decomposeCosine
%         %
%         %             [b_perceptual, bint_perceptual] = regress([perceptualErrorOnB;perceptualErrorOnX], xes);
%         %             %b_perceptual
%         %         end
%         


        if decomposeSine
            [b_biasSin, bint_biasSin] = regress(Y, xes);
        elseif decomposeCosine
            [b_biasCos, bint_biasCos] = regress(Y, xes);
            
        end
        
%         
%         %         figure()
%         %         plot([CPpredError;OBpredError],[perceptualErrorOnB;perceptualErrorOnX],'x')
%         %         hold on
%         %         if decomposeSine
%         %             xlabel('sin(objPE)')
%         %             ylabel('sin(percaptual error)')
%         %         elseif decomposeCosine
%         %             xlabel('cos(objPE)')
%         %             ylabel('cos(percaptual error)')
%         %         else
%         %             xlabel('objPE')
%         %             ylabel('percaptual error')
%         %         end
%         
%         
%       
%         if condition==1
%             PE=[CPpredError;OBpredError];
%             PE(nanR)=[];
%             StateChangeProb=[ CPsurprise;OBsurprise];
%             StateChangeProb(nanR)=[];
%             condnum=condNum;
%             condnum(nanR)=[];
%         else
%             PE=[OBpredError;CPpredError];
%             PE(nanR)=[];
%             StateChangeProb=[OBsurprise; CPsurprise];
%             StateChangeProb(nanR)=[];
%             condnum=condNum;
%             condnum(nanR)=[];
%         end
% 
% 
%        
%         
%         if makePlots
%             figure()
%             hold on
%             a=plot(PE(condnum==1), data.Y(condnum==1), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 1, 'markerSize', 8 )
%             c=plot(PE(condnum==1&StateChangeProb>.99), data.Y(condnum==1&StateChangeProb>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'c', 'lineWidth', 1, 'markerSize', 12 )
%             b=plot(PE(condnum==-1), data.Y(condnum==-1), 'or','markerSize', 8, 'markerEdgeColor', 'k','markerFaceColor', 'r', 'lineWidth', 1, 'markerSize', 8  )
%             d=plot(PE(condnum==-1&StateChangeProb>.99), data.Y(condnum==-1&StateChangeProb>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'y', 'lineWidth', 1, 'markerSize', 12 )
%             ff=legend([a, b], 'changepoint', 'oddball')
%             set(ff, 'location', 'northwest', 'box', 'off')
%             %ylim([-2.5 2.5])
%             xlabel('objective PE')
%             ylabel('perceptual error')
%             
%             refline(1,0)
%             yline(0)
%             hold off
%         end
        
        %% Saving coefficient values
        
        if decomposeSine
            params.bint_biasSin = bint_biasSin;
            params.bint_updateSin = bint_updateSin;
            params.b_biasSin = b_biasSin;
            params.b_updateSin = b_updateSin;
        elseif decomposeCosine
            
            params.bint_biasCos = bint_biasCos;
            
            params.bint_updateCos = bint_updateCos;
            
            params.b_biasCos = b_biasCos;
            
            params.b_updateCos = b_updateCos;
        end
        %
        if ~decomposeSine && ~decomposeCosine
            params.circPerceptual=bestParamsBias;
            params.circUpdate=bestParamsUpdate;
            params.neglogLikePerceptual=bestNegLogLikeBias;
            params.neglogLikeUpdate=bestNegLogLikeUpdate;
            
        end
        
        
        filename=sprintf('%s_behaveParams',subStr);
        
        dir=[baseDir,'/behaveData','/behavCoef'];
        
        saveDir=fullfile(dir,filename);
        
        save(saveDir,'params')
        
        
        
        if ~decomposeSine && ~decomposeCosine
            paramsCircUpdateAll=cat(1,paramsCircUpdateAll,params.circUpdate);
            paramsCircBiasAll=cat(1,paramsCircBiasAll,params.circPerceptual);
        end
        
        %keyboard
        if decomposeSine
            paramsUpdateAllSin=cat(2,paramsUpdateAllSin,params.b_updateSin);
            paramsBiasAllSin=cat(2,paramsBiasAllSin,params.b_biasSin);
        elseif decomposeCosine
            paramsUpdateAllCos=cat(2,paramsUpdateAllCos,params.b_updateCos);
            paramsBiasAllCos=cat(2,paramsBiasAllCos,params.b_biasCos);
        end
        
    end
    close all
    
    
    aveSinCosUpdate=((paramsUpdateAllSin'+paramsUpdateAllCos')./2);
    aveSinCosBias=((paramsBiasAllSin'+paramsBiasAllCos')./2);
    
    figure()
    hold on
    plot(aveSinCosUpdate(:,3),paramsCircUpdateAll(:,4),'x')
    refline(1,0)
    yline(0)
    xline(0)
    xlim([-1 1])
    ylim([-1 1])
    
    figure()
    plot(aveSinCosBias(:,:),paramsCircBiasAll(:,4),'x')
    
    figure()
    plot(paramsCircUpdateAll(1:8,4),paramsCircUpdateAll(9:16,4),'o')
    
    [h,p]=ttest(paramsCircUpdateAll(1:8,4))
    
    %% plot 4 panel binned plot for model
    if ~decomposeSine&~decomposeCosine
        %%%%%%%%%%%
        %OB only update
        
        subPE=OBpredErrorSub;
        Q=quantile(subPE,(1:9)/10);
        
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(subPE);
        Qs(end)=max(subPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=OBupdate(subPE>Qs(i)&subPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            subPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        figure()
        a=subplot(2,2,3)
        errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        xL=xlabel('OB subjective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('OB mean prediction update')
        yline(0)
        refline(1,0)
        xlim(a,[-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
        %%%%%%%
        %CP only update
        
        subPE=CPpredErrorSub;
        Q=quantile(subPE,(1:9)/10);
        
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(subPE);
        Qs(end)=max(subPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=CPupdate(subPE>Qs(i)&subPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            subPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        b=subplot(2,2,1)
        errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        xL=xlabel('CP subjective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('CP mean prediction update')
        yline(0)
        refline(1,0)
        xlim(b,[-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
        
        
        %%%%%%%%%%%
        %OB only bias
        
        objPE=OBpredError;
        Q=quantile(objPE,(1:9)/10);
        
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(objPE);
        Qs(end)=max(objPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=perceptualErrorOnB(objPE>Qs(i)&objPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            objPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        c=subplot(2,2,4)
        errorbar(objPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        xL=xlabel('OB objective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('OB mean perceptual error')
        yline(0)
        ylim(c,[-0.5 0.5])
        refline(1,0)
        xlim(c,[-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
        %%%%%%%
        %CP only bias
        
        objPE=CPpredError;
        Q=quantile(objPE,(1:9)/10);
        
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(objPE);
        Qs(end)=max(objPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=perceptualErrorOnX(objPE>Qs(i)&objPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            objPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        d=subplot(2,2,2)
        errorbar(objPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        ylim(d,[-0.5 0.5])
        xL=xlabel('CP objective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('CP mean perceptual error')
        yline(0)
        
        refline(1,0)
        xlim(d,[-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
         ylim(c,[-0.5 0.5])
         ylim(d,[-0.5 0.5])
%         
        
        
    end
    
    
    %         %%%%%%
    
    
elseif combineAllSubs %if combine all subs data together
    %nTrials = length(behaveAll.block(behaveAll.subNum==subNum))*.5;
    
    OBpredErrorAll = [behaveAll.predictErr(behaveAll.blockCond==-1,1);behaveAll.predictErr( behaveAll.blockCond==-1,2)];
    
    CPpredErrorAll = [behaveAll.predictErr( behaveAll.blockCond==1,1);behaveAll.predictErr( behaveAll.blockCond==1,2)];
    
    
    OBpredErrorSubAll = [behaveAll.subPredErr( behaveAll.blockCond==-1,1);behaveAll.subPredErr( behaveAll.blockCond==-1,2)];
   
    CPpredErrorSubAll = [behaveAll.subPredErr( behaveAll.blockCond==1,1);behaveAll.subPredErr( behaveAll.blockCond==1,2)];
    
    
    
    OBsurpriseAll = behaveAll.surpriseOB;
    CPsurpriseAll = behaveAll.surpriseCP;
    
    OBupdateAll = [behaveAll.predUpdate( behaveAll.blockCond==-1,1);behaveAll.predUpdate( behaveAll.blockCond==-1,2)];
    
    CPupdateAll = [behaveAll.predUpdate( behaveAll.blockCond==1,1);behaveAll.predUpdate( behaveAll.blockCond==1,2)];
   
    
    
    OBentropyAll = behaveAll.entropyOB;
    CPentropyAll = behaveAll.entropyCP;
    
    perceptualErrorOnXAll = [behaveAll.estErr( behaveAll.blockCond==1,1);behaveAll.estErr( behaveAll.blockCond==1,2)];
    
    perceptualErrorOnBAll = [behaveAll.estErr( behaveAll.blockCond==-1,1);behaveAll.estErr( behaveAll.blockCond==-1,2)];
    
    
    conditionAll=behaveAll.condition(sub*2);
    
    %         conNum=behaveAll.blockCond;
    %         condNum=[conNum(1:nTrials);conNum(1:nTrials);conNum(nTrials+1:end);conNum(nTrials+1:end)];
    %
    predictionErrorOnXAll=[behaveAll.predictionErrorOnX( behaveAll.blockCond==1,1);behaveAll.estErr( behaveAll.blockCond==1,2)];
    
    predictionErrorOnBAll=[behaveAll.predictionErrorOnB( behaveAll.blockCond==-1,1);behaveAll.estErr( behaveAll.blockCond==-1,2)];
    
    
    OBsurpriseTrialNumAll=sum(nansum(behaveAll.surpriseTrial( behaveAll.blockCond==-1,:)));
    CPsurpriseTrialNumAll=sum(nansum(behaveAll.surpriseTrial( behaveAll.blockCond==1,:)));
    
    OBmoveDistAll=[behaveAll.notMove( behaveAll.blockCond==-1,1);behaveAll.notMove( behaveAll.blockCond==-1,2)];
    
    CPmoveDistAll=[behaveAll.notMove( behaveAll.blockCond==1,1);behaveAll.notMove( behaveAll.blockCond==1,2)];
    
    
    %
    
    %% plot the 4 pannel binned plots for combined sub data
    
    
    
    if ~decomposeSine&~decomposeCosine
        %%%%%%%%%%%
        %OB only update
        
        subPE=OBpredErrorSubAll;
        %subPE(isnan(subPE))=[];
        Q=quantile(subPE,(1:9)/10);
%         
%         update=OBupdateAll;
%         update(isnan(update))=[];
%         
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(subPE);
        Qs(end)=max(subPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=OBupdateAll(subPE>Qs(i)&subPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            subPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        figure()
        a=subplot(2,2,3)
        errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        xL=xlabel('OB subjective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('OB mean prediction update')
        yline(0)
        refline(1,0)
        xlim([-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
        %%%%%%%
        %CP only update
        
        subPE=CPpredErrorSubAll;
        Q=quantile(subPE,(1:9)/10);
        
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(subPE);
        Qs(end)=max(subPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=CPupdateAll(subPE>Qs(i)&subPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            subPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        b=subplot(2,2,1)
        errorbar(subPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        xL=xlabel('CP subjective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('CP mean prediction update')
        yline(0)
        refline(1,0)
        xlim(b,[-3.2 3.2])
        ylim(b,[-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
        
        
        %%%%%%%%%%%
        %OB only bias
        
        objPE=OBpredErrorAll;
        Q=quantile(objPE,(1:9)/10);
        
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(objPE);
        Qs(end)=max(objPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=perceptualErrorOnBAll(objPE>Qs(i)&objPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            objPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        c=subplot(2,2,4)
        errorbar(objPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        xL=xlabel('OB objective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('OB mean perceptual error')
        yline(0)
        ylim(c,[-0.5 0.5])
        refline(1,0)
        xlim(c,[-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
        %%%%%%%
        %CP only bias
        
        objPE=CPpredErrorAll;
        Q=quantile(objPE,(1:9)/10);
        
        Qs=nan(1,11);
        Qs(2:10)=Q;
        Qs(1)=min(objPE);
        Qs(end)=max(objPE);
        
        
        for i=1:length(Qs)-1
            
            Qdata=perceptualErrorOnXAll(objPE>Qs(i)&objPE<Qs(i+1));
            Qdata(isnan(Qdata))=[];
            meanQs(i)=circ_mean(Qdata);
            semQs(i)=std(Qdata)./sqrt(length(Qdata));
            objPEbins(i)=mean([Qs(i) Qs(i+1)]);
            
        end
        
        d=subplot(2,2,2)
        errorbar(objPEbins,meanQs,semQs,'^','markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 0.5, 'markerSize', 10 )
        hold on
        ylim(d,[-0.5 0.5])
        xL=xlabel('CP objective PE');
        xL.Position(2) = xL.Position(2) - abs(xL.Position(2) * 0.05);
        ylabel('CP mean perceptual error')
        yline(0)
        
        refline(1,0)
        xlim(d,[-3.2 3.2])
        xticks(Qs)
        xticklabels({'min','10th','20th','30th','40th','50th','60th','70th','80th','90th','max'})
        xtickangle(45)
        
        
%         ylim(c,[-0.3 0.3])
%         ylim(d,[-0.3 0.3])
%         
        
        
    end
    
end


%cbColors=[0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./256;
%%


defaultPlotParameters;

% num        = 1;
% wid        = 17; % total width
% hts        = [8, 6.5, 6.5 ]; % height of each row
% cols       = {[.50 .50], [.50 .50], [.50 .50]}; % width of columns
% [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1.5, 2.5, [], ''); % etc
% set(axs,'Units','normalized');
% set(gcf,'color','white')
% fig_.Position = [10 10 30 40];
% left_color = [0 0 255] / 256;
% right_color = [255 0 0] / 256;
% set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% % draw in each panel, one at a time
%
%
% lw=1.5
% lw2=3
% exSub=12

%plot behave coeffs for all subs

%Matt's 3 steps
%1.generate 5 random variables
%2.loop through columns
%3. in loop plot(rand+i,params(:,i)

semUpdate=std(paramsCircUpdateAll,1)./sqrt(size(paramsCircUpdateAll,2));

x=[0.1:0.1:0.8];

figure()

%plotting 5 coefficients

stdUpdate=std(paramsCircUpdateAll);


normalizedCoef=[]

for c=2:5
   
        %x+c
        
        
        normalizedCoef=paramsCircUpdateAll(:,c)./stdUpdate(c);
        normalizedCoefUpdate(:,c-1)=normalizedCoef;
        meanToPlot=mean(normalizedCoef);
        semToPlot=std(normalizedCoef)./sqrt(length(normalizedCoef));
        %keyboard
        plot(x+c,normalizedCoef,'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 10 )
        
        
        hold on
        errorbar(mean(x+c),meanToPlot,semToPlot,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 11 )
        
        xtickVal(c)=mean(x+c);
    
end

%xlim([0.1 0.6])
%ylim([-1 1])
title('Learning Regression','color',cbColors(4,:))
xticks(xtickVal)
xticklabels({'bleh','intercept','PE','PE x STP x CP/OB','PE x STP','PE x condition'})
xtickangle(45)
yline(0);
xlim([1.5 6])
ylim([-6 6])
hold off




semBias=std(paramsCircBiasAll,1)./sqrt(size(paramsCircBiasAll,2));

x=[0.1:0.1:0.8];

figure()

stdBias=std(paramsCircBiasAll);

for c=2:5
    
    
        normalizedCoef=paramsCircBiasAll(:,c)./stdBias(c);
        
        normalizedCoefBias(:,c-1)=normalizedCoef;
        
        meanToPlot=mean(normalizedCoef);
        semToPlot=std(normalizedCoef)./sqrt(length(normalizedCoef));   
    
    
    
        plot(x+c,normalizedCoef,'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 10 )
        hold on
        errorbar(mean(x+c),meanToPlot,semToPlot,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 11 )
        xtickVal(c)=mean(x+c);
    
end

%xlim([0.1 0.6])
%ylim([-0.1 0.1])
title('Bias Regression','color',cbColors(5,:))
xticks(xtickVal)
xticklabels({'bleh','intercept','PE','PE x STP x CP/OB','PE x STP','PE x condition'})
xtickangle(45)
yline(0);
xlim([1.5 6])
ylim([-4 4])

hold off

%%

defaultPlotParameters;

% num        = 1;
% wid        = 17; % total width
% hts        = [8, 6.5, 6.5 ]; % height of each row
% cols       = {[.50 .50], [.50 .50], [.50 .50]}; % width of columns
% [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1.5, 2.5, [], ''); % etc
% set(axs,'Units','normalized');
% set(gcf,'color','white')
% fig_.Position = [10 10 30 40];
% left_color = [0 0 255] / 256;
% right_color = [255 0 0] / 256;
% set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% % draw in each panel, one at a time
%
%
% lw=1.5
% lw2=3
% exSub=12

%plot behave coeffs for all subs

%Matt's 3 steps
%1.generate 5 random variables
%2.loop through columns
%3. in loop plot(rand+i,params(:,i)

semUpdate=std(paramsCircUpdateAll,1)./sqrt(size(paramsCircUpdateAll,2));

x=[0.1:0.1:0.8];

figure()

%plotting 5 coefficients
for c=2:6
    if c==6
        sem=semUpdate(end-1);
        plot(x+c,paramsCircUpdateAll(:,end-1),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth',1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircUpdateAll(:,end)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        
        xtickVal(c)=mean(x+c);
        
    else
        x+c
        sem=semUpdate(c);
        %keyboard
        plot(x+c,paramsCircUpdateAll(:,c),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircUpdateAll(:,c)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        
        xtickVal(c)=mean(x+c);
    end
end

%xlim([0.1 0.6])
ylim([-0.01 1])
title('Learning Regression','color',cbColors(4,:))
xticks([2:6])
xticklabels({'intercept','PE','PE x p(ST) x CP/OB','PE x p(ST)','PE x p(Guess)'})
xtickangle(45)
yline(0);
xlim([1 7])

hold off




semBias=std(paramsCircBiasAll,1)./sqrt(size(paramsCircBiasAll,2));

x=[0.1:0.1:0.8];

figure()

for c=2:6
    if c==6
        sem=semBias(end-1);
        plot(x+c,paramsCircBiasAll(:,end-1),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircBiasAll(:,end)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        xtickVal(c)=mean(x+c);
    else
        x+c
        sem=semBias(c);
        plot(x+c,paramsCircBiasAll(:,c),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircBiasAll(:,c)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        xtickVal(c)=mean(x+c);
    end
end

%xlim([0.1 0.6])
ylim([-0.05 0.08])
title('Bias Regression','color',cbColors(5,:))
xticks([2:6])
xticklabels({'intercept','PE','PE x p(ST) x CP/OB','PE x p(ST)','PE x p(Guess)'})
xtickangle(45)
yline(0);
xlim([1 7])

hold off


%% plotting beta for model


defaultPlotParameters;

% num        = 1;
% wid        = 17; % total width
% hts        = [8, 6.5, 6.5 ]; % height of each row
% cols       = {[.50 .50], [.50 .50], [.50 .50]}; % width of columns
% [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1.5, 2.5, [], ''); % etc
% set(axs,'Units','normalized');
% set(gcf,'color','white')
% fig_.Position = [10 10 30 40];
% left_color = [0 0 255] / 256;
% right_color = [255 0 0] / 256;
% set(fig_,'defaultAxesColorOrder',[left_color; right_color]);
% % draw in each panel, one at a time
%
%
% lw=1.5
% lw2=3
% exSub=12

%plot behave coeffs for all subs

%Matt's 3 steps
%1.generate 5 random variables
%2.loop through columns
%3. in loop plot(rand+i,params(:,i)

figure()
subplot(2,1,1)
plot([.1:.1:.6 ; .1:.1:.6], [bint_update(:,1)'; bint_update(:,2)'], '-', 'Color', 'black')
hold on
plot(.1:.1:.6, b_update, 'o', 'MarkerFaceColor', 'Cyan', 'MarkerEdgeColor', 'black', 'markerSize', 7)
title('Learning Regression','color',cbColors(4,:))
xticks([.1:.1:.6])
xticklabels({'intercept','PE','PE x p(ST) x CP/OB','PE x p(ST)','PE x p(Guess)'})
%set(gca,'XTick',[.1:.1:.6], 'xticklabel', labels);
xtickangle(45)
yline(0);
xlim([0 0.45])

subplot(2,1,2)
plot([.1:.1:.6 ; .1:.1:.6], [bint_bias(:,1)'; bint_bias(:,2)'], '-', 'Color', 'black')
hold on
plot(.1:.1:.6, b_bias, 'o', 'MarkerFaceColor', 'Cyan', 'MarkerEdgeColor', 'black','markerSize', 7)
title('Bias Regression','color',cbColors(5,:))
xticks([.1:.1:.6])
xticklabels({'intercept','PE','PE x p(ST) x CP/OB','PE x p(ST)','PE x p(Guess)'})
%set(gca,'XTick',[.1:.1:.6], 'xticklabel', labels);
xtickangle(45)
yline(0);
xlim([0 0.45])
ylim([-0.2 0.3])

% plot([.1:.1:.4 ; .1:.1:.4], [bint_update(:,1)'; bint_update(:,2)'], 's', 'markerSize', 8);
plot([0, .6], [0, 0], '--')
ylim([-0.3 0.9])
xlim([0 0.45])
set(gca,'XTick',[.1:.1:.6], 'xticklabel', labels);
ylabel("Value");
set(gca, 'box', 'off');


semUpdate=std(paramsCircUpdateAll,1)./sqrt(size(paramsCircUpdateAll,2));

x=[0.1:0.1:0.8];

figure()

%plotting 5 coefficients
for c=2:6
    if c==6
        sem=semUpdate(end-1);
        plot(x+c,paramsCircUpdateAll(:,end-1),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth',1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircUpdateAll(:,end)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        
        xtickVal(c)=mean(x+c);
        
    else
        x+c
        sem=semUpdate(c);
        %keyboard
        plot(x+c,paramsCircUpdateAll(:,c),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircUpdateAll(:,c)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        
        xtickVal(c)=mean(x+c);
    end
end

%xlim([0.1 0.6])
ylim([-0.01 1])
title('Learning Regression','color',cbColors(4,:))
xticks([2:6])
xticklabels({'intercept','PE','PE x p(ST) x CP/OB','PE x p(ST)','PE x p(Guess)'})
xtickangle(45)
yline(0);
xlim([1 7])

hold off




semBias=std(paramsCircBiasAll,1)./sqrt(size(paramsCircBiasAll,2));

x=[0.1:0.1:0.8];

figure()

for c=2:6
    if c==6
        sem=semBias(end-1);
        plot(x+c,paramsCircBiasAll(:,end-1),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircBiasAll(:,end)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        xtickVal(c)=mean(x+c);
    else
        x+c
        sem=semBias(c);
        plot(x+c,paramsCircBiasAll(:,c),'o','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1, 'markerSize', 11 )
        hold on
        errorbar(mean(x+c),mean(paramsCircBiasAll(:,c)),sem,'^','markerEdgeColor', 'k', 'markerFaceColor', cbColors(c,:), 'lineWidth', 1.5, 'markerSize', 13 )
        xtickVal(c)=mean(x+c);
    end
end

%xlim([0.1 0.6])
ylim([-0.05 0.08])
title('Bias Regression','color',cbColors(5,:))
xticks([2:6])
xticklabels({'intercept','PE','PE x p(ST) x CP/OB','PE x p(ST)','PE x p(Guess)'})
xtickangle(45)
yline(0);
xlim([1 7])

hold off




%%
% MRN added this code to look at the data to troubleshoot weird behavioral
% results.


PE=[OBpredError;CPpredError];
StateChangeProbUpdate=[OBsurprise; CPsurprise];

figure()
hold on
a=plot(PE(condNum==1), data.Y(condNum==1), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 1, 'markerSize', 8 )
c=plot(PE(condNum==1&StateChangeProbUpdate>.99), data.Y(condNum==1&StateChangeProbUpdate>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'c', 'lineWidth', 1, 'markerSize', 12 )
b=plot(PE(condNum==-1), data.Y(condNum==-1), 'or','markerSize', 8, 'markerEdgeColor', 'k','markerFaceColor', 'r', 'lineWidth', 1, 'markerSize', 8  )
d=plot(PE(condNum==-1&StateChangeProbUpdate>.99), data.Y(condNum==-1&StateChangeProbUpdate>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'y', 'lineWidth', 1, 'markerSize', 12 )
ff=legend([a, b], 'changepoint', 'oddball')
set(ff, 'location', 'northwest', 'box', 'off')
%ylim([-2.5 2.5])
xlabel('subjective PE')
ylabel('update')
hold off



PE=[OBpredErrorObj;CPpredErrorObj];
StateChangeProb=[OBsurprise; CPsurprise];

figure()
hold on
a=plot(PE(condNum==1), data.Y(condNum==1), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'b', 'lineWidth', 1, 'markerSize', 8 )
c=plot(PE(condNum==1&StateChangeProb>.99), data.Y(condNum==1&StateChangeProb>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'c', 'lineWidth', 1, 'markerSize', 12 )
b=plot(PE(condNum==-1), data.Y(condNum==-1), 'or','markerSize', 8, 'markerEdgeColor', 'k','markerFaceColor', 'r', 'lineWidth', 1, 'markerSize', 8  )
d=plot(PE(condNum==-1&StateChangeProb>.99), data.Y(condNum==-1&StateChangeProb>.99), 'ob', 'markerEdgeColor', 'k', 'markerFaceColor', 'y', 'lineWidth', 1, 'markerSize', 12 )
ff=legend([a, b], 'changepoint', 'oddball')
set(ff, 'location', 'northwest', 'box', 'off')
%ylim([-2.5 2.5])
xlabel('objective PE')
ylabel('Bias')
hold off























%
