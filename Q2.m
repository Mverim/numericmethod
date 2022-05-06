clc,clear
m1 = input('Enter first mass(kg): ');           % Getting mass inputs
m2 = input('Enter second mass(kg): ');
m3 = input('Enter third mass(kg): ');
k = input('Enter spring coefficient(kg/s2): '); % Getting friction coef. input
g = 9.81;           % Gravity constant (m/s2)
Mass = [m1 m2 m3];  % Mass array
MgMass = Mass*g;    % Mass array fixed for 'g' constant

%% Matrices for calculating the displacements of the masses
% K_Matrix is the stiffness matrix
K_Matrix = [3*k -2*k 0; -2*k 3*k -1*k; 0 -1*k 1*k];
% MATLAB 'inv' function performs LU decomposition
% Check mathworks documentation on 'inv' for more info
K_Inversed = inv(K_Matrix); % K_Inversed is inversed of K_Matrix
X = zeros(3,1);             % 'X' array for holding displacement data

%% Calculating the displacements between masses
for i = 1:3
    for j=1:3
        X(i) = X(i) + K_Inversed(j,i) * MgMass(j);
    end
end

fprintf('X1=%f\nX2=%f\nX3=%f\n',X(1,1),X(2,1),X(3,1));

