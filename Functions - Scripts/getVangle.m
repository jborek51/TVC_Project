function [theta_V,d] = getVangle(V)
if V(2) == 0
    if V(1) == 0
        theta_V = 0;
    elseif V(1) > 0
        theta_V = 90;
    else
        theta_V = 90;
    end
elseif V(1) == 0
    if V(2) > 0 
        theta_V = 0;
    else
        theta_V = 180;
    end
elseif V(1) > 0 && V(2) > 0
    theta_V = atan(V(1)/V(2))*180/pi;
elseif V(1) < 0 && V(2) > 0
    theta_V = -atan(V(1)/V(2))*180/pi;
elseif V(1) < 0 && V(2) < 0
    theta_V = 180-atan(V(1)/V(2))*180/pi;
else
    theta_V = 180+atan(V(1)/V(2))*180/pi;
end
theta_V = theta_V*pi/180;
if V(1) <= 0
    d = 1;
else
    d = -1;
end
end