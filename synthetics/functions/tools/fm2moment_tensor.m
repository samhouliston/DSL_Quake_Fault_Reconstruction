function M = fm2moment_tensor(d, n, M0)
% Compute an equivalent moment tensor from focal mechanism vectors
% After Stein & Wysession, p. 242, equation (4)
% Men-Andrin Meier, 25/5/2023

% Inputs: 
%   d = slip vector
%   n = normal vector
%   M0 scalar seismic moment in Nm (default = )
%
% Output:
%   M = moment tensor

if nargin<3; M0 = 1e19;
end

M = M0*[2*n(1)*d(1),           n(1)*d(2)+n(2)*d(1), n(1)*d(3)+n(3)*d(1); ...
          n(2)*d(1)+n(1)*d(2), n(2)*d(2),           n(2)*d(3)+n(3)*d(2); ...
          n(3)*d(1)+n(1)*d(3), n(3)*d(2)+n(2)*d(3), 2*n(3)*d(3)];