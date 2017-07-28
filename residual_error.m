function [ residuals ] = residual_error( Vtot, gamma, dl, ...
    differentialArea, cl )

% Computes the residuals using formula 8 of the Hunsaker-Snyder paper

k = size(Vtot,1);
residuals = zeros(k,1);

for i = 1:k
    residuals(i) = 2 * norm( cross(Vtot(i,:),dl(i,:)) ) * gamma(i) ...
        - norm(Vtot(i,:))^2 * differentialArea(i) * cl(i);
end

end