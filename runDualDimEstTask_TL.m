


% run VSTM task:
% shell script to run VSTM task

%% CHANGES TO MAKE:

% Lets get orientation out of the picture -- we'll want a box centered on
% the stimulus location but colored gray to cue which stimulus to report.

% Lets make an "option" that controls whether rectangles/squares actually
% have a jittered orientation. For now we'll probably want to keep them all
% at the same orientation -- but down the line we may want to see how
% specific our findings are to the attended stimulus feature. --Done


% change number of items (1, 2, 3)? -- Done

% add generative dynamics to how colors are generated -- Done
% create two generative structures (changepoint & oddball) * in prefs
% struct --Done??

% force subject to predict color on each trial -- all colors, one at time.
% MAYBE have them report all colors too? last priority -- but could be
% useful. --Done





%% LAB MEMBER PILOT DATA TO COLLECT:

% Hoping to get data and feedback from a couple of lab members on this new
% visual working memory task. its not set up in any sort of fancy wrapper
% as I may want to change a lot of stuff when i get feedback. I'm looking
% for feedback on two versions of the task... one in which the subject is
% shown the color of the probed target and asked to report the orientation
% of that target (taskType=3) and another version where the subject is
% shown the orientation of a target and asked to report the color of the
% target that had that orientation (taskType=4)

% This version of the task also differs from my previous versions in that
% it has feedback (thumbs up for small error, down for large error).

% It also has an occasional post-decision wager... the post decision wager
% trials are designed to ask the subject whether they would like to bet an
% extra dollar on their memory from this trial... or whether they'd like to
% pass on this trial and not bet at all. OK, dollar is a little steep...
% but a point... which eventually should translate into a few cents.

% To run the color estimation version of the task the settings should be:
%       nTrials  = 100;
%       taskType = 4;

% The orientation version of the task should be similar
%       nTrials  = 100;
%       taskType = 3;


whichComp=3;

if whichComp==1
    basePath='/Users/ttli/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    colorPath='/Users/ttli/Dropbox (Brown)/sharedMatlabUtilities/';
    
elseif whichComp==2
    basePath='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    colorPath='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/';
    
elseif whichComp==3
    basePath='/Users/mattLab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    colorPath='/Users/mattLab/Dropbox (Brown)/sharedMatlabUtilities/';


end

% SN=input('Enter subject number:', 's');
% 
% ETname=input('Enter name for Eyetracking data file:','s');
% edfFile=[ETname '.edf'];
% 
% Cond=input('Enter condition for subject:','s');
% 
% CondN=str2double(Cond);
% 
% TN='vwm_multiStim_color_';

%cd('/Users/mattnassar/Dropbox/vwm_task/8_dualDimEstTask/')
cd(basePath)
addpath(genpath('./'))
addpath(genpath(colorPath))





%monitorData.gammaTable=
%load Gateway_monitor_normMatrix.mat



calibrationFile='labMonitorData_2019-10-30_426c';


load(fullfile(basePath,calibrationFile))
%load('/Users/mattnassar/Dropbox/vwm_task/code/labMonitorData_2014-10-10')


%prefs.taskType = 0;          %: Pick one of several task types, previously 3
%%
prefs.numBlocks=4;
%prefs.nTrials = 4;
prefs.nItems = 2;
prefs.stimulusDuration = .5;%3;%.5
prefs.retentionInterval = .900;
prefs.squareSize = 40;       % in pixels
prefs.radius = 160;
prefs.fixationSize = 3;
prefs.colorWheelRadius = 350;
prefs.propFixedSpace = 0;
prefs.fixedSpacings  = [360./prefs.nItems: 360./prefs.nItems:359];
prefs.fixedOSpacings = [360./prefs.nItems: 360./prefs.nItems:359];
prefs.useCieLab=true;
prefs.monitorData=monitorData;
prefs.rectSize1= 120;%80; %120;         %: length of stimulus rectangle
prefs.rectSize2= 120; %80; %120;         %: width of stimulus rectangle
prefs.debug = false;         %: Are you debugging the task? true = small screen
prefs.preEstPause=1;
%prefs.doFeedback=0;          %  Display + or - feedback after each trial?
prefs.fdbkProp=0;            %  Prop b feedback.
prefs.fdbk_oThresh=pi./6;    %  orientation error threshold
prefs.fdbk_cThresh=pi./6;    %  color threshold
prefs.fdbk_duration= 1;
prefs.ITI= .5;
prefs.pdw_prop=0;          %  Proportion of trials that will have a post decision wager
prefs.pdw_vals=[0 2];        %  Possible bet values [array of two integers]
prefs.openWindow=true;
prefs.closeWindow=false;
prefs.initRTthresh=15;
prefs.skipSyncTests=true;
prefs.jittered = false;
%prefs.doPredict =true;
prefs.showMask=true;
prefs.testGray = round(LabTo_calRGB([prefs.monitorData.lVal, 0,0], prefs.monitorData));
prefs.pointsPerDollar=50;
prefs.doEEG=false;
prefs.doET=false;
if whichComp==3
    prefs.fontSize=40;
    prefs.imageSizeX=1500;
    prefs.imageSizeY=1000;
else
    prefs.fontSize=20;
    prefs.imageSizeX=1200;
    prefs.imageSizeY=800;
end

SN=input('Enter subject number:', 's');


ETname=input('Enter name for Eyetracking data file:','s');
edfFile=[ETname '.edf'];


%1 is CP first, 2 is OB first; counter-balance across subjects
Cond=input('Enter condition for subject:','s');

CondN=str2double(Cond);

TN='vwm_multiStim_color_';



%%%%%%%%%%%%%%%%% Initialize EEG Calls %%%%%%%%%%%%%%%%%%%%%%%%%%
if prefs.doEEG
    config_parallel
    ioObject = io64;
    LPT1address = hex2dec('3FB8');
    %'E050'); %standard location of LPT1 port
    status = io64(ioObject);

    % Create structure of trigger numbers:
    triggerNum.start           =15;
    triggerNum.end             =14;
    
    triggerNum.practice           =16;
    %triggerNum.block           =17;
    
    triggerNum.block           =1;
    triggerNum.instructionsOn  =2;
    triggerNum.instructionsOff =3;
    triggerNum.stimOn          =4;
    triggerNum.stimOff         =5;
    triggerNum.cue1on          =6;
    triggerNum.responseMade1   =7;  %16;  %240;
    triggerNum.cue2on          =8;
    triggerNum.responseMade2   =9;  %32;  %224;
    triggerNum.predCue1on      =10;
    triggerNum.predResp1       =11;    %48;  %208;
    triggerNum.predCue2on      =12;
    triggerNum.predResp2       =13;   %64;  %192;
    
    % put structure in prefs:
    
    prefs.triggerNum=triggerNum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Test trigger
% % 
% while true
%     io64(ioObject,LPT1address,prefs.triggerNum.block);
%     pause(0.05)
%     disp('got here')
%     io64(ioObject,LPT1address,0);
% end

%io64(ioObject,LPT1address,0)

%%
% if prefs.doEEG
%          io64(ioObject,LPT1address,prefs.triggerNum.start);
%          WaitSecs(0.01)
%          io64(ioObject,LPT1address,0);
% end

%EyelinkSetup(1)
for bb=1:prefs.numBlocks
     % send trigger indicate the start of every block
%      if prefs.doEEG
%          io64(ioObject,LPT1address,prefs.triggerNum.block);
%          WaitSecs(0.01)
%          io64(ioObject,LPT1address,0);
%      end
%      
%      if prefs.doET
%          myMessage = 'block start';
%          Eyelink('SendMessage', myMessage);
%      end
     
    if bb>1
        prefs.openWindow=false;
        if bb==prefs.numBlocks
            prefs.closeWindow=true;
        end
    end
    
    fn=[TN SN num2str(bb) '.mat'];
    
    if CondN==2 && bb==3
        bb=4;
    elseif CondN==2 && bb==4
        bb=3;
    end
    
    [data, prefs]=colorOrient_estTask_TL_test(prefs, fn, edfFile,bb,CondN);
   
        
    
    blockScore(bb)=sum(nansum(data.score,2));
    %io64(ioObject,LPT1address,0);
end


disp(sprintf('You Earned %s Points for %s bonus dollars!!!', num2str(sum(blockScore)), ...
    num2str(round(sum(blockScore)./prefs.pointsPerDollar))))

% if prefs.doEEG
%          io64(ioObject,LPT1address,prefs.triggerNum.end);
%          WaitSecs(0.01)
%          io64(ioObject,LPT1address,0);
% end


