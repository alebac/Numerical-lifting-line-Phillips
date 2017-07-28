clear all
close all
clc

%% Input data
k = 50;
iterations = 15;
relaxationFactor = 1;

sweepAngle = 0;
dihedralAngle = 0;

rootChord = 0.2;
taperRatio = 1;
wingspan = 0.8;
airspeed = [30 0 0];
density = 1.225;                % kg/m^3

alphaFly = 4; 
twist = 0;

alpha0lift = 0;
alphaRad = alphaFly * pi/180;
alpha0liftRad = alpha0lift * pi/180;
clalpha = 5.73; %2*pi;

%% Propellers
hubPositions = 0.4;
diameters = 0.15;
thrust = 5;
rpm = 25400;                                   % rounds per minute
propellerRotationSpeed = rpm/60 * 2*pi;         % radians per second

spinnerRadius = 0.1 * max(diameters);

%% PHILLIPS %%

% Wing twist
clear alpha;
alphaTwist = abs((0:k-1) / k - 1/2)*2 * twist;
alphaTwist = alphaTwist';
alpha = alphaFly + alphaTwist;

%% Propellers
propData = propellersPosition(hubPositions, diameters, k, wingspan, ...
    propellerRotationSpeed, airspeed, density, thrust);
inducedVelocity = induced_velocity_propeller(propData, airspeed, ...
    spinnerRadius);

Vrel = ones(k+1,1)*airspeed;
Vrel = Vrel + inducedVelocity;

%% 1. Mesh
%  The lifting line is specified by a matrix of the stations
%  Each line is the position of a station and it specifies the [x, y, z]
%  coordinates

liftingLineBorders = zeros(k+1,3);
liftingLine = zeros(k,3);

% Discretization of the lifting line in the x and y directions
station = 1:k;
dy = wingspan/k;

% Lifting line: control points at the center of the sections, 3/4 of chord
y = station'*dy - wingspan/2 - dy/2;
% Lifting line borders: nodes at the extremities of the sections, 1/4 of
% chord
station2 = 0:k;
y2 = station2'*dy - wingspan/2;

chord = ones(k,1) * rootChord; 
chord2 = ones(k+1,1) * rootChord;
% chord = rootChord - 2*abs(y4)/wingspan * (1-taperRatio)*rootChord;
% chord2 = rootChord - 2*abs(y5)/wingspan * (1-taperRatio)*rootChord;
% chord = rootChord * sqrt(1-(2*y4/wingspan).^2);
% chord2 = rootChord * sqrt(1-(2*y5/wingspan).^2);

surface = trapz(y2,chord2);

% Lifting line
sweepAngleRad = sweepAngle / 180 * pi();
x = abs(y) * tan(sweepAngleRad) + 0.25*rootChord;
dihedralAngleRad = dihedralAngle / 180 * pi();
z = abs(y) * tan(dihedralAngleRad);
liftingLine(:,1) = x;
liftingLine(:,2) = y;
liftingLine(:,3) = z;

% Lifting line borders
x2 = abs(y2) * tan(sweepAngleRad) + 0.25*rootChord;
z2 = abs(y2) * tan(dihedralAngleRad);
liftingLineBorders(:,1) = x2;
liftingLineBorders(:,2) = y2;
liftingLineBorders(:,3) = z2;

% dl and differential areas
dl = zeros(k,3);
differentialArea = zeros(k,1);
for i = 1:k
    dl(i,:) = liftingLineBorders(i+1,:) - liftingLineBorders(i,:);
    differentialArea(i) = dy * chord(i);
end

% % Surface
% %surface = rootChord*(1+taperRatio)*wingspan/2;
% % surface = pi*wingspan*rootChord/4;
% surface = trapz(y5,chord2);

%% 2. Induced vortices
inducedVelocities = induced(liftingLine, liftingLineBorders);

%% 3. Assume Gamma
gamma = 1 * sqrt(1-(2*y/wingspan).^2);

for z = 1:iterations
    tic
%% 4. Compute total velocity
Vtot = total_velocity(gamma, inducedVelocities, Vrel);

%% 5. Compute angle of attack and cl
u_ai = [1 0 0]; % chordwise unit vector at control point i 
u_ni = [0 0 1]; % normal unit vector at control point i
alphaiRad = angle_of_attack(Vtot, u_ni, u_ai);
alphaRad = alpha * pi/180;
alpha0liftRad = alpha0lift * pi/180;

cl = clalpha * (alphaRad - alpha0liftRad + alphaiRad);

%% 6. Compute residuals
residuals = residual_error(Vtot, gamma, dl, differentialArea, cl);

%% 7. Compute gradient
grad = gradient(Vtot, dl, u_ni, u_ai, clalpha, differentialArea, ...
    inducedVelocities, gamma, cl);

%% 8. Solve system
deltaGamma = grad\-residuals;
%% 9. Update gamma
gamma = gamma + relaxationFactor * deltaGamma;

toc
end

%% 10. Last computations
CL = 2 / (norm(airspeed)*surface) * trapz(y,gamma);
CDi = 2 / (norm(airspeed)*surface) * trapz(y,gamma.*(-alphaiRad));
AR = wingspan^2 / surface;
e = CL^2/(pi*AR*CDi);

alphai = alphaiRad*180/pi;   %% this makes the rotation on the induced angle of attack different from the other 2 methods

%% Plots
plots
