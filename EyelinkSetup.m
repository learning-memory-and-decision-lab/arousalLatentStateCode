function [el] = EyelinkSetup(startOrStop,window)
%global DP;
%global SP;
%Local function to set up defaults for each Eyelink and start recording.
%Designed to be called before starting each block. Requires previous call
%to EyelinkInitDefaults to work appropriately

if startOrStop==1 %starting up the Eyelink
    el=EyelinkInitDefaults(window); %initialize the Eyelink default settings
    
    % Initialize Eyelink connection (real or dummy). The flag '1' requests
    % use of callback function and eye camera image display:
    if ~EyelinkInit([], 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;
        return;
    end
    
    % Send any additional setup commands to the tracker
    Eyelink('Command','calibration_type = HV13');
    Eyelink('Command','recording_parse_type = GAZE');
    %perform calls to the Eyelink host to get the correct queued eye data
    %Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');
    %Eyelink('command','inputword_is_window = ON');
    %Eyelink('Command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK');
    Eyelink('Command','sample_rate = 1000');
    %Eyelink('Command','screen_pixel_coords = 0 0 1151 863'); %need to fix the resolution input thing here
    Eyelink('Command','draw_line = x1, y1, x2, y2, 0');
    Eyelink('Command','file_event_data = GAZE,VELOCITY');
    Eyelink('Command','link_event_data = GAZE,VELOCITY');
    Eyelink('Command','saccade_velocity_threshold = 200');
    Eyelink('Command','saccade_motion_threshold = 0.15');
    result = Eyelink('StartSetup',1);
    Eyelink('StartRecording');

    disp('Eyelink up and running!');
    
elseif startOrStop==0 %done with Eyelink for now, shut it down
    
    %reset sampling rate for other experimenter's use and shutdown the Eyelink
%     Eyelink('Command','sample_rate = 250');
%     Eyelink('Command','calibration_type = HV9');
    disp('....now clearing eyelink! :)');
    Eyelink('Shutdown');
    el=-1;

end

