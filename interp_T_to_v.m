function Ty = interp_T_to_v(Ny,yu,T)
% Interpolate Temp to the v-momementum cells

Ty = ( ( yu(3:Ny+2)-yu(2:Ny+1) )'.*T(3:Ny+2,:) + ...
         ( yu(2:Ny+1)-yu(1:Ny) )'.*T(2:Ny+1,:) ) ./ (yu(3:Ny+2)-yu(1:Ny))'; 