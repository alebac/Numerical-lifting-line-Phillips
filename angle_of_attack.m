function [ alpha ] = angle_of_attack( Vtot, u_ni, u_ai )
% This function computes the angle of attack from the total velocity, the
% normal unit vector at control point i u_ni, and the spanwise unit vector 
% at control point i. Formula 9 of Phillips paper. The result is in
% radians.

k = size(Vtot,1);
alpha = zeros(k,1);

for i = 1:k
    alpha(i) = atan( dot(Vtot(i,:),u_ni) / dot(Vtot(i,:),u_ai) );
end

end

