%%
clear;                   %clearing all variables
x=0:.01:10;              %initilizing x and then introducing the functions.
f=exp(-x).*sin(5*x);
g=exp(x).*sin(5*x);
t=0:.01:5;

h=sin(pi*t)./(pi*t);
l=sin(t.^2);
%opening a new figure for each function with command "figure"
%and adding title and labels.
figure(1);plot(x,f);title('f(x)=e^(-x)*sin(5x)');xlabel('x');ylabel('y');    
figure(2);plot(x,g);title('g(x)=e^x*sin(5x)');xlabel('x');ylabel('y');
figure(3);plot(t,h);title('h(x)=sin(pi*x)/(pi*x)');xlabel('x');ylabel('y');
figure(4);plot(t,l);title('l(x)=sin(x^2)');xlabel('x');ylabel('y');


%%
%ploting all 4 functions in a single figure.
figure(5);
subplot(2,2,1);plot(x,f);title('f(x)=e^(-x)*sin(5x)');xlabel('x');ylabel('y');
subplot(2,2,2);plot(x,g);title('g(x)=e^x*sin(5x)');xlabel('x');ylabel('y');
subplot(2,2,3);plot(t,h);title('h(x)=sin(pi*x)/(pi*x)');xlabel('x');ylabel('y');
subplot(2,2,4);plot(t,l);title('l(x)=sin(x^2)');xlabel('x');ylabel('y');
%%
figure;
%ploting functions f and g ,which have considerable difference in amount, in a single figure.
%with the scale setting to show both functions completely
%since this system is intended to show both functions completely,it is not so useful when 
%working with functions that have high magnitudes of difference,just as in this case,
%where g has the amount of a thousand and f only oscilates around zero.
%the result will be that the function that has the small magnitude will not be shown clearly.
plot(x,f); 
hold on
plot(x,g);
title('plotting two functions with extreme difference in just one scale system');
 %plotting them in the same figure but with two different axes and scale.
figure(2);[h,a,b]=plotyy(x,f,x,g);title('the same functions using different scale for each');   
ylabel(h(1),'f(x)=e^(-x)*sin(5x)');
ylabel(h(2),'f(x)=e^(x)*sin(5x)');
%%
%this section checks whether or not a matrix is magic.
%this code is run three times and each time the output is shown.
%first run with magic(4) as input matrix.
clear
a=input('\ninput matrix to determine whether is magic or not:');

%checking whether or not the matrix is square
asize=size(a);
if mod(sum(asize),2)~=0 fprintf('\ninput matrix is not a square matrix\n');  
else 
    diagsum=sum(diag(a));%calculating sum of each row and column and the diagnol.
    columnsum=sum(a);
    rowsum=sum(a');
      %checking whether or not the mentioned sums are equal or not
    if diagsum==columnsum & diagsum==rowsum ,fprintf('YES\n');
    else fprintf('NO\n'); end
end
%second run with an arbitrary matrix:[1 2;3 4]
clear
a=input('\ninput matrix to determine whether is magic or not:');
asize=size(a);
if mod(sum(asize),2)~=0 fprintf('\ninput matrix is not a square matrix\n'); 
else 
    diagsum=sum(diag(a));
    columnsum=sum(a);
    rowsum=sum(a');
    if diagsum==columnsum & diagsum==rowsum ,fprintf('YES\n');
    else fprintf('NO\n'); end
end
%third run with a nonsquare matrix:[1,2]
clear
a=input('\ninput matrix to determine whether is magic or not:');
asize=size(a);
if mod(sum(asize),2)~=0 fprintf('\ninput matrix is not a square matrix\n'); 
else 
    diagsum=sum(diag(a));
    columnsum=sum(a);
    rowsum=sum(a');
    if diagsum==columnsum & diagsum==rowsum ,fprintf('YES\n');
    else fprintf('NO\n'); end
end
%%
%this is the answer to the system AX=B
clear
A=[2 9 3 0 6 2 9 4 5 7;4 9 7 5 3 2 8 0 6 4;
    0 0 6 8 1 6 8 5 4 0;9 7 5 6 4 2 2 4 8 2;
    1 2 6 1 4 8 5 6 7 1;1 4 6 3 1 9 0 6 9 2;
    3 5 1 4 5 7 4 6 5 4;1 9 1 9 2 3 3 0 3 5;
    4 4 9 1 3 5 1 0 1 4;3 9 1 8 5 1 1 3 6 8];
B=[217;279;227;236;156;285;225;262;164;168];
X=A\B
%We could also calculate X using inv(A)*B
%%
fprintf('The inverse is:\n');disp(inv(A))
fprintf('sum of each row is:\n');z=sum(A');disp(z');
fprintf('average of each column is:');format bank;disp(mean(A))
fprintf('The variance of the diagnol of A is:');format;disp(var(diag(A)))
%%
clear
[x,y]=meshgrid(-1:.01:1);    %initilizing x and y and then introducing z
z=5*x.^2+3*y.^2;
 %opening a figure for each plot and adding titles and labels.
figure(1);contour(x,y,z);title('contourplot of z=5x^2 + 3y^2');xlabel('x');ylabel('y');
figure(2);surf(x,y,z);title('surfplot of z=5x^2 + 3y^2');xlabel('x');ylabel('y');zlabel('z');
figure(3);plot3(x,y,z);title('plot3 of z=5x^2 + 3y^2');xlabel('x');ylabel('y');zlabel('z');

%%
clear
 %generating a random matrix of dimentions 10 by 10 with only integers 0 to 10
A=floor(11*rand(10,10));                
t=input('to display the matrix,enter 1.otherwise enter 0:');
if t==1 ,disp(A);end
%the function histogram takes as input, in this case particularly, first the data set matrix and second the number of rectangles
histogram(A,10);

