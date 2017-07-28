function [ w ] = induced( controlPoints, boundVortices )
% Commputes the induced velocities. The singularity is canceled neglecting
% the term induced by the bound vortex on himself

k = size(controlPoints,1);
w = cell(k,k);
u_inf = [1 0 0];

for i = 1:k
    for j = 1:k

%         R0 = boundVortices(j+1,:) - boundVortices(j,:);
        R1 = controlPoints(i,:) - boundVortices(j,:);
        R2 = controlPoints(i,:) - boundVortices(j+1,:);
        r1 = norm(R1);
        r2 = norm(R2);

        if i == j
            vAB = [0 0 0];
        else
            vAB = (r1+r2) * (cross(R1,R2)) / (r1*r2*(r1*r2 + dot(R1,R2)));
        end
        
        vINF = cross(u_inf,R2)/(r2*(r2-dot(u_inf,R2))) ...
             - cross(u_inf,R1)/(r1*(r1-dot(u_inf,R1)));
        
        w{i,j} = (vAB + vINF) / (4*pi);
        
    end
end

end
