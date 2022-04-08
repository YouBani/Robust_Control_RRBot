% Generation of a cubic trajectory
function a = CubicTraj(tspan, x0, xf)
    to = tspan(1);
    tf = tspan(2);
    
    A = [1 to to^2 to^3;
        0 1 2*to 3*(to^2);
        1 tf tf^2 tf^3;
        0 1 2*tf 3*(tf^2)];
    b = [x0; xf];
    
    a = A\b;
end