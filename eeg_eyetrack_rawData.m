%% EEG (Pupillometry) preprocessing script for R00

% Stolen from cannon preprocessing. Stolen shamelessly from Anne Collins who stole (or maybe politely borrowed)
% it from Jim Cavanagh. -----Matt

%2:37
%% 
clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                  set path and params                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%$$$$ pay attention to the path changes

% change whichComp accordingly!!
whichComp=3;

showPlots=true;

if whichComp==1 % iMac in lab manager office
    baseDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities';
elseif whichComp==2 % lab Macbook 
    baseDir='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities';
else % all other computers 
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities';
end


% change the local path to own username 
locpath=('/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');

behaveDir=fullfile(baseDir,'behave_data');
%datapath = '/Users/LabManager/Desktop/eeg_data';

if whichComp>2
    datapath=fullfile(baseDir,'eeg_data'); %eeg data folder 
else
    datapath=fullfile(baseDir,'eeg_data');
end
eyetrackDir=fullfile(baseDir,'ET_data');

addpath(genpath(baseDir))
addpath(genpath(eeglabDir))

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
cd(datapath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   load data                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% insert subject number below 
subnoStr = 'BV_demo_KaitlynPilot';

subno=str2num(subnoStr);

EEG = pop_loadbv(datapath,[subnoStr,'.vhdr']);
%EEG = pop_loadbv(datapath,['/BV_demo_' subnoStr,'.vhdr']); %EEG file names
EEG = pop_chanedit(EEG,  'lookup', locpath); % if this fails, use loc path as above but from local EEGLab
EEG = eeg_checkset( EEG );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     steps before first pass reject by inspection                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     triggerNum.start           =15;
%     triggerNum.end             =14;
%     
%     triggerNum.block           =1;
%     triggerNum.instructionsOn  =2;
%     triggerNum.instructionsOff =3;
%     triggerNum.stimOn          =4;
%     triggerNum.stimOff         =5;
%     triggerNum.cue1on          =6;
%     triggerNum.responseMade1   =7;  
%     triggerNum.cue2on          =8;
%     triggerNum.respnseMade2   =9;  
%     triggerNum.predCue1on      =10;
%     triggerNum.predResp1       =11;    
%     triggerNum.predCue2on      =12;
%     triggerNum.predResp2       =13;   

Triggers = {'S  4'}; % for stimulus locked 

EEG = pop_epoch( EEG, Triggers, [-8   8 ], 'BV_demo_KaitlynPilot', 'Epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );

EEG.data=EEG.data(1:63,:,:);
EEG.chanlocs = EEG.chanlocs(1,1:63);
EEG.nbchan=63;
% Add Cz  - but don't re-ref yet
% EEG = pop_chanedit(EEG, 'append',63,'changefield',{64 'CPz' 180 0.1266 -32.9279 [-4.0325e-15] 78.3630 -180 67.2080 85});
EEG = pop_chanedit(EEG, 'append',63,'changefield',{64, 'labels', 'Cz'});
%%180 0.1266 -32.9279 [-4.0325e-15] 78.3630 -180 67.2080 85});
EEG = pop_chanedit(EEG,  'lookup', locpath);

EEG = pop_rmbase(EEG,[],[]);% Remove mean

eeglab redraw
ALLEEG=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             First Pass Cleaning                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use GUI and click tools->reject data epochs->reject by inspection
% focus on period between S4 and S6              /
% set voltage range to 60

% Once you have highlighted bad epochs click "update marks" in lower right

if ~isempty(EEG.reject.rejmanual)
to_reject=find(EEG.reject.rejmanual);
disp(sprintf('These numbers should be rejected: \n %s', num2str(to_reject))); % copy this printout into the if statement below to save rejected channels
end

























