function src = get_finite_source_from_HyC(mag, stressdrop)
% Compute square-shaped finite source centered around hypocentre.

%stressdrop = 1e6;  % [1 MPa]
shearmod   = 30e9; % [30 GPa]

M0     = magnitude2moment(mag);
src.L  = ( 2*M0./(pi*stressdrop) ).^(1/3); % [m]
src.A  = src.L.^2; % [m^2]
src.D  = M0./(shearmod*src.L.^2);     % [m]