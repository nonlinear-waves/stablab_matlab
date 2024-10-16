function [r,dr] = react(u)

if u > 0
    r = exp(-1/u);
    dr = exp(-1/u)/u^2;
else
    r = 0;
    dr = 0;
end
