function coefM = TaylorCoefs(points,frac_output)
%   Input vector of distances of unit 'dx' for points under consideration
%   relative to the point at which the derivative is wanted, e.g.
%   points = [-2 -1 0] for a 3-point downwind scheme at u_i-2, u_i-1, u_i
%   Set symb = true for fractional output (default is decimal)

if nargin == 1
    frac_output = false; % Set default output to decimal format
end

N = length(points);
coefM = ones(N);
for j = 2:N
    coefM(j,:) = points.^(j-1)/factorial(j-1); 
end

if frac_output == true
    coefM = sym(coefM); % Converts output to pretty-printed fractions
end