
function res = lossFunction(param)
% For other test-functions, see https://en.wikipedia.org/wiki/Test_functions_for_optimization

x = param(1);
y = param(2);

disp([x, y]);

%% Mishra's Bird function - constrained
f = @(x, y) sin(y)*exp((1-cos(x))^2) + cos(x)*exp((1-sin(y))^2) + (x-y)^2;

if (x+5)^2 + (y+5)^2 >= 25 ||...
        x < -10 || x > 0 ||...
        y < -6.5 || y > 0
    inBounds = false;
else 
    inBounds = true;
end

%% Rosenbrock 
% f = @(x, y) (1-x)^2 + 100*(y-x^2)^2;
% if x^2 + y^2 > 2 ||...
%         abs(x) > 1.5 || ...
%         abs(y) > 1.5
%     inBounds = false;
% else 
%     inBounds = true;
% end    

%%
if inBounds == true
    res = f(x, y);
else
    % If not in bounds, give a *large* value
    res = 1E10;
end
