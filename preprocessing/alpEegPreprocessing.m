%% EEG preprocessing script for R00

% Stolen from cannon preprocessing. Stolen shamelessly from Anne Collins who stole (or maybe politely borrowed)
% it from Jim Cavanagh. -----Matt

%2:37

%% just hit "Run and Advance" for this cell 

clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                set path (directory) and params                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%$$$$ pay attention to the path changes

% change whichComp accordingly!!
whichComp=1;


showPlots=true;

if whichComp==1 % iMac in lab manager office
    baseDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    locpath=('/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
elseif whichComp==2 % lab Macbook 
    baseDir='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    locpath=('/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
elseif whichComp==3 % Meera's Macbook 
    baseDir='/Users/meerasingh/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/meerasingh/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/meerasingh/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    locpath=('/Users/meerasingh/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
elseif whichComp==4 %  Kaitlyn's Macbook
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    locpath=('/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
else % all other computers 
    % change the local path to own username 
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities';
    locpath=('/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
end


if whichComp>2
    datapath=fullfile(baseDir,'EEG_data'); %eeg data folder 
else
    datapath=fullfile(baseDir,'eeg_data');
end

%
eyetrackDir=fullfile(baseDir,'ET_data');

addpath(genpath(baseDir))
addpath(genpath(eeglabDir))
addpath(genpath(sharedFuncDir))  


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
cd(datapath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   load data                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% insert subject number below 
subnoStr = '2003';
subno=str2num(subnoStr);
% if subnoStr == '0311_022620'
%     subno=1999;
% 
% end



%EEG = pop_loadbv(datapath,[subnoStr,'.vhdr']);
%EEG = pop_loadbv(datapath,['/BV_demo_' subnoStr,'.vhdr']); %EEG file names
EEG = pop_loadbv(datapath,['/ALP_' subnoStr,'.vhdr']);


EEG = pop_chanedit(EEG,  'lookup', locpath); % if this fails, use loc path as above but from local EEGLab

EEG = eeg_checkset( EEG );

eeglab redraw

%after importing EEG data, use the gui to convert eye data from ascii to
%.mat

%then use the gui to import and synchronize ET data

%then move onto the filtering step

%% 
fn=fullfile(datapath, subnoStr,[subnoStr,'_RAW_w_EYE.mat']);
save(fn,'EEG','-v7.3');

fn=fullfile(datapath, subnoStr,[subnoStr,'_w_saccades.mat']);
save(fn,'EEG','-v7.3');

fn=fullfile(datapath, subnoStr,[subnoStr,'_newICAx.mat']);
save(fn,'EEG','-v7.3');



%% filter low freq before epoch

EEG=pop_eegfiltnew(EEG, [], 0.05,[],true,[],0);
%EEG=pop_eegfiltnew(EEG, [], 30,[],true,[],0);

mkdir(subnoStr)

fn=fullfile(datapath, subnoStr,[subnoStr,'_ALP_RawFilt.mat']);
save(fn,'EEG','-v7.3');
fprintf('filtered data saved for subject %s!\n', subnoStr);

beep

%% sanity check on raw data
chanLabs={EEG.chanlocs.labels};

fc3_ind=find(strcmp( 'FC3', chanLabs));
fc3=squeeze(EEG.data(fc3_ind,:,:));
fc3_plot=fc3(:,1:10:end);
figure
plot(EEG.times,fc3, 'r');
ylim([-1000 1000])
plot(EEG.times(:,1:10:end),fc3_plot, 'r');

EEG.data=EEG.data-nanmean(EEG.data,2);

eeglab redraw
ALLEEG=[];
pop_eegplot( EEG, 1, 1, 1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     steps before first pass reject by inspection                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just "Run and Advance" don't have to change stuff here

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
%     triggerNum.responseMade2   =9;  
%     triggerNum.predCue1on      =10;
%     triggerNum.predResp1       =11;    
%     triggerNum.predCue2on      =12;
%     triggerNum.predResp2       =13;   

Triggers = {'S  4'}; % for stimulus locked 

% epoching the continuous EEG data by trial, time 0 is S4 (stimulus onset) 
EEG = pop_epoch( EEG, Triggers, [-2  2], 'newname', 'Epochs', 'epochinfo', 'yes');
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

fn=fullfile(datapath, subnoStr,[subnoStr,'_ALP_rawEpoched.mat']);
save(fn,'EEG','-v7.3');
fprintf('epoched data saved for subject %s!\n', subnoStr);


%% sanity check
chanLabs={EEG.chanlocs.labels};

for i=1:length(EEG.chanlocs)
    chanLab=EEG.chanlocs(i).labels;
    chanData=squeeze(EEG.data(i,:,:));
    %bn_chanData=chanData-nanmean(chanData(EEG.times<0&EEG.times>-500,:));
    figure
    plot(EEG.times,nanmean(chanData, 2), 'b');
    title(chanLab)
end

o1_ind=find(strcmp( 'O1', chanLabs));
o1=squeeze(EEG.data(o1_ind,:,:));

bn_o1=o1-nanmean(o1(EEG.times<0&EEG.times>-500,:));

figure
plot(EEG.times,nanmean(bn_o1, 2), 'b');

pz_ind=find(strcmp( 'Pz', chanLabs));
pz=squeeze(EEG.data(pz_ind,:,:));

bn_pz=pz-nanmean(pz(EEG.times<0&EEG.times>-500,:));

figure
plot(EEG.times,nanmean(bn_pz, 2), 'r');

figure
plot(EEG.times,nanmean(EEG.data,3))

fc3_ind=find(strcmp( 'FC3', chanLabs));
fc3=squeeze(EEG.data(fc3_ind,:,:));
figure
plot(EEG.times,nanmean(fc3, 2), 'r');




data2plot=squeeze(nanmean(EEG.data(:,1,:)));
    
figure
topoplot(data2plot,EEG.chanlocs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             First Pass Cleaning                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do this for manual rejection

% use GUI and click tools->reject data epochs->reject by inspection
% focus on period between S4 and S6              /
% set voltage range to 60

% Once you have highlighted bad epochs click "update marks" in lower right

if ~isempty(EEG.reject.rejmanual)
    to_reject=find(EEG.reject.rejmanual);
    fprintf('These numbers should be rejected: \n %s\n', num2str(to_reject));
    % copy this printout into 'to_reject' below to save rejected channels
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             First Pass Cleaning                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do this for automatic rejection

[EEG,rmepochs]=pop_autorej(EEG,'threshold',1000,'startprob',5,'maxrej',5,'eegplot','on','nogui','on');

[EEG,indelec,measure,com]=pop_rejchan(EEG,'threshold',5,'norm','on','measure','kurt');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             Reject epoch/interpolate channels                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%if decide to interp chans, look up chan label from chan number by
    % ALLEEG.chanlocs.labels in workspace
    %to_reject var should have 2 visual rows, top row for epochs rejected 1st
    %pass,... then bottom row for epochs rejected in 2nd pass
    % example to_reject= [ #, #, #,...
    % #, #];

% HALT: IF THIS IS THE FIRST PASS... COPY THE DISPLAYED EPOCHS INTO THE
% APPROPRIATE to_reject definition!!!

% Otherwise (if this is the second pass and you've already appended 2nd
% stage rejections) proceed!


% add elseif statements for other subjects
if subno == 1999
    to_interp={'T8'};
    to_reject=[1    2    3    4    5    7    8   11   12   16   26   28   30   31   33   37   40   57   59   65   68   70   72   76   80   83   84   88   90   91  102  108  119  122];
    epochNumbers = 1:EEG.trials;
    epochNumbers(to_reject)=[];
elseif subno == 2001
    to_interp={'T7','P4','TP9','PO7','P6'};
    to_reject=[4    5    6   16   25   68   73   83  132];
    epochNumbers = 1:EEG.trials;
    epochNumbers(to_reject)=[];
else
    fprintf('Subject has already been cleaned... loading rejected epochs for subject %s\n', num2str(subno))
    fn=fullfile(datapath, subnoStr, [num2str(subno), '_rejectedEpochsAndChannels']);
    load(fn);
end

%Interpolate bad channels
if ~isempty(to_interp)
    for xi=1:size(to_interp,2)
        for ei=1:63
            if strmatch(EEG.chanlocs(ei).labels,to_interp{xi});
               badchans(xi)=ei;
            end
        end 
    end
    EEG.data=double(EEG.data);
    EEG = pop_interp(EEG,badchans,'spherical');
end

% Reject bad epochs
if ~isempty(to_reject)
    binarized=zeros(1,EEG.trials);
    binarized(to_reject)=1;
    EEG = pop_rejepoch(EEG,binarized,0);
end


% NOW re-ref - after interpolation (old ref: Cz)
EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'type',{[]},'theta',{0},'radius',{0},'X',{0},'Y',{0},'Z',{85},'sph_theta',{0},'sph_phi',{90},'sph_radius',{85},'urchan',{64},'ref',{''}));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        Save cleaned data & rejected channels and epochs                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mkdir(subnoStr)

fn=fullfile(datapath, subnoStr,[subnoStr,'_ALP_CLEANED.mat']);
save(fn,'EEG','epochNumbers','-v7.3');
fprintf('cleaned data saved for subject %s!\n', subnoStr);


% save rejected epochs and channels:
fn=fullfile(datapath,subnoStr, [subnoStr, '_rejectedEpochsAndChannels']);
save(fn,'to_interp','to_reject','epochNumbers', '-v7.3');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       run ICA                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ctrl Enter, Ignore popup box
% ICA -- Use separate script (batchICA) if running multiple subjects at once (recommended)

%this may take 10-15 minutes.30 Do not interrupt ICA.
EEG = pop_runica(EEG,'icatype','runica');

fn=fullfile(datapath, subnoStr,[subnoStr,'_ALP_ICA.mat']);
save(fn,'EEG','epochNumbers','-v7.3'); 
fprintf('ICA saved for subject %s!\n', subnoStr);

beep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 Reject eyeblinks                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subno = 603;
% subnoStr = '603';
clearvars -except locpath datapath subno subnoStr sharedFuncDir


cd([datapath,'/', subnoStr]);
%load([num2str(subno),'_CDA_chunk_ICA_Stim.mat']);
load([subnoStr,'_ALP_ICA.mat']);
pop_selectcomps(EEG, 1:20 ); % Displays topoplots
pop_eegplot( EEG, 0, 1, 1);   % Displays component scroll


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                         Reject ICA                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add elseif statements for other subjects
if subno == 1999 % Subjects 1999 = Meera (Pilot)
    to_ICAreject = [1];% LIST COMPONENTS FOR BLINKS AND EYE MOVEMENTS HERE!
elseif subno == 2001
    to_ICAreject = [1];
    
end


EEG = pop_subcomp(EEG, to_ICAreject, 0);

EEG = eeg_checkset( EEG );
eeglab redraw
ALLEEG=[];
pop_eegplot( EEG, 1, 1, 1);

fn=fullfile(datapath, subnoStr,[subnoStr,'_ALP_ICAx.mat']);
save(fn,'EEG','epochNumbers','-v7.3');
fprintf('eyeblinks saved for subject %s!\n', subnoStr);

fn=fullfile(datapath, subnoStr, [subnoStr, '_ICA_comp_to_remove']);
save(fn, 'to_ICAreject', '-v7.3');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                      Filtering                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ctrl Enter, Read comments below
% this Filters data:

addpath(genpath(sharedFuncDir))

cd([datapath,'/', subnoStr]);
%load([num2str(subno),'_CDA_chunk_ICA_Stim.mat']);
load([subnoStr,'_ALP_ICAx.mat']);

dims=size(EEG.data);
                   
for i=1:dims(3)  
    FEEDBACK(:,:,i)=EEG.data(:,:,i); % for feedback locked - JCH 7/30/14
end
                                           
% Filter
dims = size(FEEDBACK);
FILT=eegfilt_JFC(FEEDBACK,EEG.srate,[],30); % eegfilt % low pass
%FILT=eegfilt_JFC(FILT,EEG.srate,.1,[]); % eegfilt % high pass
FILT_FB=reshape(FILT,dims(1),dims(2),dims(3)); FILT=[];
EEG.data = FILT_FB;
filt_fn= [fn(1:end-19) '_FILT_STIM.mat'];
save(filt_fn,'EEG','epochNumbers','-v7.3');
fprintf('filtered data saved for subject %s!\n', num2str(subno));
beep

% IF YOU ARE ON THE SECOND PASS: Congrats. You are finished. almost.

% FINAL STEPS:

% 1) run scan through the 2nd pass cleaning quickly to make sure that you
%       didn't miss any gigantic artifacts


% 2) copy the rejected epochs and interpolated channels into the running
%       sheet on google drive. make sure you do this accurately

% 3) delete the elseif statement for this subject in the first pass
%    cleaning loop... so that next time someone runs this subject it just
%    loads your list of epochs and channels to reject directly.

% 4) copy the ICA components that were removed for this subject
%    (to_ICAreject) into the running sheet on dropbox

% 5) delete the elseif for the subject you just cleaned from the ICA
% component removal loop (~line 225)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   Second Pass Cleaning                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Pass Cleaning: Hit any key to move through channels

% CAREFUL: if the data is SUPER noisy this will generate a LOT of plots. 

subnoStr = '0311_022620';

whichComp=1;


% change the local path to own username 

if whichComp==1 % iMac in lab manager office
    baseDir='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    locpath=('/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
elseif whichComp==2 % lab Macbook 
    baseDir='/Users/mattlab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    eeglabDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eegAnalysis';
    locpath=('/Users/mattlab/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
else % all other computers 
    baseDir='/Users/kaitlynmi/Dropbox (Brown)/arousalLearningPerception/vwm_task/RAs';
    eeglabDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities/eeglab2019_1';
    sharedFuncDir='/Users/kaitlynmi/Dropbox (Brown)/sharedMatlabUtilities';
end


behaveDir=fullfile(baseDir,'behave_data');

if whichComp>2
    datapath=fullfile(baseDir,'EEG_data'); %eeg data folder 
else
    datapath=fullfile(baseDir,'eeg_data');
end
eyetrackDir=fullfile(baseDir,'ET_data');

addpath(genpath(baseDir))
addpath(genpath(eeglabDir))

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
cd(datapath);

clearvars -except locpath datapath EEG subnoStr
cd([datapath,'/',subnoStr]);

load([subnoStr,'_FILT_STIM.mat']);

thresh=50 % if subject has lots of alpha waves you may have to increase this 
          % threshold so that the loop below doesn't show you every trial!
 
t1=find(EEG.times>-1500, 1)
t2=find(EEG.times<4000, 1, 'last')

FILT_FB = EEG.data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   Generating plots                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ctrl Enter, press key to keep generating new plots, examine plot and note channel on notebook

% This stuff is unnecessary... you can do it to get a sense of what the data look like, but you don't need to
 
% % FOCUS MAINLY ON STUFF BETWEEN THE RED LINES!!!! %don't actually do
% anything but get an impression of the status of channels to pay attention
% to later in voltage fluctuations 
for electrode = 1:64
    toPlot_time = EEG.times(t1:t2);
    toPlot_chan = squeeze(FILT_FB(electrode,t1:t2,:));
    toPlot_chanMean = squeeze(mean(FILT_FB(electrode,t1:t2,:),3));

    Chanmax = max(toPlot_chan);
    Chanmin = min(toPlot_chan);

    [max_Chanmax,badep1] = max(Chanmax);
    [min_Chanmin,badep2] = min(Chanmin);

    if max_Chanmax > thresh || min_Chanmin < -thresh    
        figure
        hold on
        plot(toPlot_time,toPlot_chan,'k')
        plot([-1000, -1000], [-thresh, thresh], '--r')
        plot([1900, 1900], [-thresh, thresh], '--r')
        title([EEG.chanlocs(electrode).labels,' = channel ',num2str(EEG.chanlocs(electrode).urchan)]);
        pause
        close(gcf)
    else
        continue
    end
end
close(gcf)



 %% Ctrl + Enter
%SECOND PASS CLEANING!!!
eeglab redraw
% use GUI and click tools->reject data epochs->reject by inspection
% focus on p183  188  197  198  218  220  235  237  242  245  254  267  283  298  299  301  328  332  371  382  415  436  437  456eriod between S2 and S6/

% Once you have highlighted bad epochs click "update marks" in lower right
% Then run the code below. 
% paste the printout into the definition of "secRuntoReject" below.


% STEP 1) REJECT EPOCHS USING GUI AS DESCRIBED ABOVE
% STEP 2) GET LIST OF EPOCHS BY RUNNING 2 LINES OF CODE BELOW

% STEP 3) PASTE THAT LIST INTO LINE 400 secRuntoReject = [ PASTE HERE!]
% STEP 4) RUN THAT LINE OF CODE




%%  PART 1: reject epochs using GUI as described above... copy them into secRuntoReject

%the two lines of code below will print the new numbers in the command
%window to reject
to_reject=find(EEG.reject.rejmanual);
fprintf('copy these numbers: %s\n', num2str(to_reject)); 
% copy this printout into the secRuntoReject variable matrix below to save rejected channels

% clear var secRuntoReject, paste the printout into var secRuntoReject = [ PASTE HERE!], Ctrl Enter, epochs that contain extreme voltage fluctuations during critical epochs, keep pressing Enter to generate plots until no more appear
% if idenitied bad CHANNELS int he to_interp var above add the channel name
% if idenitied bad EPOCHS add them into the var newRejects in next cell
% to find the name for the bad CHANNEL open var ALLEEG > chanlocs > go to
% cell # you identified > click cell > read label

%%


secRuntoReject = [58];


%% Part 2: look at epochs flagged using a fixed threshold for outliers
% FOCUS MAINLY ON STUFF BETWEEN THE RED LINES!!!!

% One more search for potentially bad epochs:

alreadyDitched=false(size(squeeze(FILT_FB(1,t1:t2,:)),2),1);
alreadyDitched(secRuntoReject)=true;

np=0;
maxNp=15;
for channel = 1:64 % Enter channels to look at here: Can enter multiple at a time
    toPlot_chan = squeeze(FILT_FB(channel,t1:t2,:));
    Chanmax = max(toPlot_chan);
    Chanmin = min(toPlot_chan);
    chanMean= nanmean(toPlot_chan');
    
    badep1=find(Chanmax>thresh & ~alreadyDitched');
    badep2=find(Chanmin<-thresh & ~alreadyDitched');
 %   [max_Chanmax,badep1] = max(Chanmax);
 %   [min_Chanmin,badep2] = min(Chanmin);
    
    allBadEp=[badep1,badep2];
 
    for i=1:length(allBadEp)
          figure
          hold on
          a=plot(EEG.times(t1:t2), chanMean, 'g');
          b=plot(EEG.times(t1:t2), toPlot_chan(:,allBadEp(i)))
          title(sprintf('channel %d, epoch %s', channel, num2str(allBadEp(i))))
          legend([a, b], 'mean erp', 'trial erp')
          plot([0, 0], [-thresh, thresh], '--r')
          plot([2000, 2000], [-thresh, thresh], '--r')
           set(gcf, 'position', [400 278 860 420])  
          
          if np>maxNp
              input('Enter to close plots and continue!')
              close all
              np=0;
          end

          
          np=np+1;
    
    end   
    
end
%% Clear var newRejects, enter into newRejects epochs that contain extreme voltage fluctuations during critical epochs, Ctrl Enter, enter the numbers into the Second line of the to-reject matrix for that particular subject
% start script over and run all cells in script THROUGH the filter cell

% Write numbers of epochs containing extreme voltage fluctuations during
% critical epochs here:
newRejects=[];  % ADD epochs IDed through threshold method here! %the input for the matrix newRejects is the numbers of the epochs that were bad
secRuntoReject=unique([secRuntoReject, newRejects]);


% this code gets the epoch numbers from the original dataset (before rejecting). 
origEpochNumb = []; % DON'T ENTER ANYTHING HERE
for j = 1:size(secRuntoReject,2)
    origEpochNumb = cat(2,origEpochNumb,epochNumbers(secRuntoReject(j)));
end
disp('Enter these numbers back into top of script under to_reject (add to existing bad epochs)'); %"top of script" means in "Frist pass cleaning cell,...
%enter these into the Second line of the to-reject matrix for that particular subject
%numbers to be entered above (for the second time running through script): 60   86  126  188  197  218  242  245  299  324  332


% ADD THE NUMBERS THAT GET DISPLAYED IN COMMAND LINE TO THE EPOCH REJECTION
% LIST ABOVE!!!! (where you added epoch numbers from first pass)
disp(num2str(origEpochNumb));
% DON'T FORGOT TO ADD 2ND PASS EPOCHS TO TOP OF SCRIPT AND RE-RUN *through* FILTERING
% When you're all done, FILT_FB file should be most recently saved





















