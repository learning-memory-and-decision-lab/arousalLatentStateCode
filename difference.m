function [LR, UP, PE]=difference(outcomes, predictions, coordSys)
% Compute update (UP), prediction error (PE), and learning rate (LR=UP/PE)
% for either cartesian or polar coordinates
% Last modified by Ryan Thorpe on 2020.01.24

% Define default coordinate system
if nargin<4
    coordSys = 'cartesian';
end


UP=nan(length(outcomes),1);   %nans

switch coordSys
    case 'cartesian'
        UP(1:end-1) = predictions(2:end)-predictions(1:end-1);
        PE = outcomes-predictions;
    case 'polar'
        PE = circ_dist(outcomes,predictions); % Minimal distance between outcome and prediction orientations
        UP = circ_dist(predictions(2:end),predictions(1:end-1));
        UP_options(:,1) = circ_dist(predictions(2:end),predictions(1:end-1)); % Distance one way around the circle
        UP_options(:,2) = UP_options(:,1)- 2*pi; % Distance the other way around the circle
        UP_options(:,3) = UP_options(:,1)+ 2*pi; % Distance the other way around the circle
        
        
        [~,opt_ind] = min(abs(UP_options-PE(1:end-1)),[],2); % Index
        
        % Select the update that travels in the direction of PE
        for ti=1:length(opt_ind)
            UP(ti) = UP_options(ti,opt_ind(ti));
        end
    otherwise
        error("computeLR error: didn't recognize coordinate system!! Use coordSys='cartesian' or 'polar'.")
end

% Learing rate
LR=UP./PE;
end
