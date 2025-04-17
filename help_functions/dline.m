function dline(col)
% dline(col)
if nargin < 1
    col = 'k';
end
range=axis;
hold on
plot([range(1) range(2)],[range(3) range(4)],col)
hold off
