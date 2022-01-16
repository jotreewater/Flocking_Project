function bump = Bump(z)
h = 0.2;
if (z < h && z ~= 0)
    bump = 1;
elseif (z > h && z < 1)
    bump = (1 / 2) * (1 + cos( pi * ((z - h)/(1 - h))));
else
    bump = 0;
end
end