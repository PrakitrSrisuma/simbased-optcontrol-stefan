function [position,isterminal,direction] = Event_FinalIce(t,T,S_target)

position = T(end)-S_target;  % stop when ice volume is met
isterminal = 1; % Halt integration 
direction = 0; % The zero can be approached from either direction

end