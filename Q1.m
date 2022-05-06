clc, clear, clearvars

%% Initializing the choosen methods 
method = input(['Choose desired method:\nType 1 for bracketing' ...
    '  method\nType 2 for open method\n']);
% Using Swamee-Jain eq. to find approximate root value
RootApproximate = vpa(1.325/log((1.5E-6/3.7*5E-3)+(5.74/13743^0.9))^2,32);

if method == 1
finalroot = BisectionFunction();
end
if method == 2
finalroot = NewtonRaphsonFunction();
end

%% Calculating Root Using Bisection Method
function [finalroot] = BisectionFunction()   

clc
% Setting x as symbolic variable
syms x;

% Inputs and initializing variables
trueval = 0.0289678; % Got the truevalue from the book
y = 1/sqrt(x) + 2*log10( 1.50e-06/(3.7*0.0050 ) + 2.51/ ...
    ( 1.374301675977654e+04*sqrt(x) ) ); % Friction function
a = input('Enter lower bracket guess: ');
b = input('Enter upper bracket guess: ');
e = input('Enter desired Max error(epsilon): ');

% Evaluating and finding values at inputs
fa = eval(subs(y,x,a));
fb = eval(subs(y,x,b));

% Implementing Bisection Method

if fa*fb > 0                % Checking if the values bracket the actual root
    disp('Root is not inside of the given values.');
else
    c = (a+b)/2;            % Setting root as (Xupper+Xlower)/2
    fc = eval(subs(y,x,c)); % Evaluating the function in the supposed root
    fprintf('\n\nLower\t\tUpper\t\tRoot\t\t\tError\t\t\tf(c)\n'); % Printing 'a b c er fc' tables

    while abs(fc)>e
        er=abs(((c-trueval)/trueval)*100);              % Percent relative error
        fprintf('%f\t%f\t%f\t%f\t%f\n',a,b,c,er,fc);    % Printing values of current 'a b c er fc'

        if fa*fc< 0
            b =c; % Setting new Xupper as the middle value
        else
            a =c; % Setting new Xlower as the middle value
        end
        c = (a+b)/2;            % Getting the middle value
        fc = eval(subs(y,x,c)); % Substituting 'c' instead of 'x' in function and evaluating
    end
    finalroot = vpa(c,32);      % Setting precision to 32 decimals and printing
    fprintf('\nRoot is: %f\n', finalroot); 
end

end

%% Calculating Root Using Newton-Raphson Method
function [finalroot] = NewtonRaphsonFunction()
clc
% Setting x as symbolic variable
syms x;

% Getting inputs, setting constants and core function
trueval = 0.0289678; % Got the truevalue from the book
y = 1/sqrt(x) + 2*log10( 1.50e-06/(3.7*0.0050 ) + 2.51/ ...
    ( 1.374301675977654e+04*sqrt(x) ) ); % Friction function
a = input('Enter initial guess: ');
e = input('Enter desired Max error(epsilon): ');
N = input('Enter maximum number of steps: ');
% Initializing step counter
step = 1;
 
% Finding derivate of given function
g = diff(y,x);

% Finding Functional Value
fa = eval(subs(y,x,a));

while abs(fa)> e
    fa = eval(subs(y,x,a)); % Function evaluated at point 'a'
    ga = eval(subs(g,x,a)); % Function's derivative evaluated at 'a'

    if ga == 0 % If the evaluated derivative value is 0, function errors
        disp('Division by zero.');
        break; % Since the algorithm has fa/ga part cannot be divided by 0
    end

    b = a - fa/ga;                     % Actual Newton-R. method formula
    er = abs((b-trueval)/trueval*100); % Percent relative error
    fprintf('Step=%d\tRoot=%f\tError=%f\tf(a)=%f\n',step,a,er,fa);
    a = b; 
    
    if step>N % If the iteration diverges, this shows the error
       disp('Not convergent'); 
       break;
    end
    step = step + 1;
end
finalroot = vpa(a,32); % Setting precision to 32 decimals and printing
fprintf('Root is %f\n', finalroot);
end
