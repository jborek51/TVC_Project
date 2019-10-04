function [v_app,alpha] = getVapp(v_z,v_y,v_w,theta)
V = [v_y v_z];
V_w = [v_w 0];
V_app = V_w-V;

[theta_V,d_V] = getVangle(V);
if norm(V) == 0
    theta_a = pi/2;
else
    theta_a = acos((-norm(V_w)^2+norm(V)^2+norm(V_app)^2)/(2*norm(V)*norm(V_app)));
end
if V_app(2) > 0
    theta_0 = theta_a-(theta_V*d_V);
else
    theta_0 = theta_a+theta_V*d_V;
end
if V_app(1) >= 0 && V_app(2) <= 0
    alpha = theta_0-theta;
elseif V_app(1) < 0 && V_app(2) <= 0
    alpha = theta_0-theta;
elseif V_app(1) < 0 && V_app(2) > 0
    alpha = -theta_0-theta;
% elseif V_app(1) >= 0 && V_app(2) > 0
else
    alpha = -theta_0-theta;
end
if alpha < -pi
    alpha = 2*pi+alpha;
end
v_app = norm(V_app);
end