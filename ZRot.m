function RM=ZRot(phi,ori)
% returning rotation matrix
if nargin~=2
    RM=[];
    errordlg('Invalid input ...');
else
    phi=phi/180*pi; % converted into radius
    switch ori
        case 'x'
            RM=[1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi);];
        case 'y'
            RM=[cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi);];
        case 'z'
            RM=[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1;];
        otherwise
            RM=ZRot(ori,'z')*ZRot(phi/pi*180,'x')*ZRot(-ori,'z');
    end
end
end
    