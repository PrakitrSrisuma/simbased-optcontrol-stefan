function [position,isterminal,direction] = Event_CompleteMelting(t,T,input)

position = T(end)-input.tol;  % stop when ice volume is met
isterminal = 1; % Halt integration 
direction = 0; % The zero can be approached from either direction

end