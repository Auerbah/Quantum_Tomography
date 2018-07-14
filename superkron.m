function x = superkron(varargin)
x = varargin{1};
n = nargin;
for i=2:n
    x = kron(x, varargin{i});
end