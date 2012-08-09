%function value = myfun(h)

f = @(h) 0.04*acos((0.2-h)/0.2) - (0.2-h)*sqrt(0.4*h-h^2) - 0.3