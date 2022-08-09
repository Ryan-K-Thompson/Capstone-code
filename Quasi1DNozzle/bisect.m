function [x] = bisect(f,a,b,n)
% From Intro to Numerical Methods Ohio University
%Performs n iterations of the bisection method to find root for function f
%Inputs: f -- an anonymous function 
%        a,b -- left and right edges of interval
%        n -- number of bisections performed
%Outputs: x -- the estimated solution of f(x)=0
%         e -- an upper bound on the error
% Commented out the portion that displays the coordinate of the root.
format long; %Use as many decimal places as possible
c = f(a); %Set c equal to function at a
d = f(b); %Set d equal to function at b
if c*d > 0.0 %If c times d is positive, bisection method will not work.
    fprintf('\nError\n');
    error('Function has the same sign at both endpoints.Bisection method will not work');
end
%disp('xy');%Displace xy to show coordinate of root
for k = 1:n %Perform n iterations of bisection method
    x = (a+b)/2; %new x is half the interval
    y = f(x); %Evaluate f(x)
    %disp([x y]); %display value of (x,y) coordinate
    if y == 0.0 %If y equals zero, have found root 
        e = 0; %Set initial error equal to zero
        break %Exit for loop
    end 
    if c*y<0 %Opposite sign
        b = x;%Set new right limit to x ((a+b)/2)
    else a = x;%Otherwise, set left limit to x (a=x);
    end
end
end