function r = cauchyrnd(location, scale, m, n)
% cauchyrnd - Generate random samples from the Cauchy distribution
%
%   r = cauchyrnd(location, scale)           → returns a single sample
%   r = cauchyrnd(location, scale, m, n)     → returns m×n samples
%
%   Inputs:
%       location - Location parameter (x0)
%       scale    - Scale parameter (γ > 0)
%       m, n     - Output matrix size (optional)
%
%   Output:
%       r        - Random numbers sampled from Cauchy(location, scale)

    if nargin == 2
        r = location + scale * tan(pi * (rand - 0.5));
    else
        r = location + scale * tan(pi * (rand(m, n) - 0.5));
    end
end
