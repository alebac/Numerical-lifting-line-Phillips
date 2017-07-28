function [ Vtot ] = total_velocity( gamma, inducedVelocities, Vrel )

% Computes the total velocity from the circulation distribution and the
% upstream velocity Vrel using formula 5 of Phillips paper

k = size(gamma,1);
Vtot = zeros(k,3);
for i = 1:k
    induction = [0 0 0];
    for t = 1:k
       induction  = induction + gamma(i)*inducedVelocities{i,t};
    end
    Vtot(i,:) = Vrel(i,:) + induction;
end

end

