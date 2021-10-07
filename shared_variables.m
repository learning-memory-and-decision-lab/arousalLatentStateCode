function vars = shared_variables(datapath)
    vars.nTrials = 90;
    vars.H = .15;
    vars.replace = true;
    vars.sigma_c = 20;
    vars.sigma_x = 10;
    vars.sigma_y = 10;
    vars.nDegrees = 360;
    vars.path = datapath;
    vars.generate_data = false;
    vars.average = false;
    %"/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/2001_allBlockData.mat";
    save("shared_variables.mat");
end 


% function vars = shared_variables()
%     vars.nTrials = 50;
%     vars.H = .15;
%     vars.replace = true;
%     vars.sigma_c = 10;
%     vars.sigma_x = 10;
%     vars.sigma_y = 10;
%     vars.nDegrees = 360;
%     vars.generate_data = true;
%     vars.path = "~/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/2001_allBlockData.mat";
%     vars.average = true;
%     vars.average_number = 15;
% end 