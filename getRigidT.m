function A = getRigidT(param, type, DH)
    %Calculate rigid transformation from DH parameters by DH convention.
    if (strcmp(type, 'Revolute'))
        param(4) = param(4) + DH;
    elseif (strcmp(type, 'Prismatic'))
        param(3) = param(3) + DH;
    end
    
    %Retrieve the parameters
    a = param(1);
    alpha = param(2);
    d = param(3);
    theta = param(4);
    
    %Create frame i to frame j rigid tranformation.
    A = [cos(theta), -sin(theta), 0, a;
        cos(alpha)*sin(theta), cos(alpha)*cos(theta), -sin(alpha), -d*sin(alpha);
        sin(alpha)*sin(theta), sin(alpha)*cos(theta), cos(alpha), d*cos(alpha);
        0, 0, 0, 1];
end