% file name = BetaLogHyperGeom
% date of creation = 17/03/2024
%
% Study of Ln(X) and of Ln(-X) when X~B(a,b)
%
clear all
close all
pkg load statistics
pkg load gsl
%
% parameters
%
a1=5;
b1=3;
a2=3;
b2=7;
n=2000; % number o randomly generated numbers
select = 1 % 1 for using the hypergeometric function, 2 for an approximate summation
%
% random generation and histagram
%
X1=betarnd(a1,b1,1,n);
Y=log(X1);
M=mean(Y);
num_bins = 30;
[counts, bin_centers] = hist(Y, num_bins);
bin_width = bin_centers(2) - bin_centers(1);
%
% analytic density of Y=Ln(X1)
%
len=50;
y=linspace(5*M,0,len);
for i=1:len
  fa(i)=exp(y(i)*a1)*(1-exp(y(i)))^(b1-1);
endfor
fa=fa/beta(a1,b1);
%
% plotting
%
figure
subplot(1,3,1)
bar(bin_centers,counts/(n*bin_width),'FaceColor','yellow','EdgeColor','black','BarWidth',1)
hold on
plot(y,fa,'-k','Linewidth',1)
hold off
title(strcat('density of Y=Ln(X_{1}), X_{1}~Beta(',num2str(a1),',',num2str(b1),')'),'fontsize',15)
xlabel('Ln(X_{1})','fontsize',15)
ylabel('Probability Density','fontsize',15)
%
% histogram and density of Y=Ln(1/X2)
%
X2=betarnd(a2,b2,1,n);
Y=-log(X2);
M=mean(Y);
[counts,bin_centers]=hist(Y,num_bins);
bin_width = bin_centers(2) - bin_centers(1);
%
y=linspace(0,5*M,len);
for i=1:len
  fa(i)=exp(-y(i)*a2)*(1-exp(-y(i)))^(b2-1);
endfor
fa=fa/beta(a2,b2);
%
% plotting
%
subplot(1,3,2)
bar(bin_centers,counts/(n*bin_width),'FaceColor','green','EdgeColor','black','BarWidth',1)
hold on
plot(y,fa,'-k','Linewidth',1)
hold off
title(strcat('density of Y=Ln(1/X_{2}),X_{2}~Beta(',num2str(a2),',',num2str(b2),')'),'fontsize',15)
xlabel('Ln(1/X_{2})','fontsize',15)
ylabel('Probability Density','fontsize',15)
%
% histogram of Y=Ln(X1/X2)
%
Z=log(X1./X2);
M=mean(Z);
[counts,bin_centers]=hist(Z,num_bins);
bin_width = bin_centers(2) - bin_centers(1);
%
% analytic density of Y=Ln(X1/X2)
%
minimo=min(Z);
massimo=max(Z);
len=50;
z=linspace(minimo,massimo,len);
clear fa
if (select==1)
  for i=1:len
    fa(i)=BetaLogDensityHG(a1,b1,a2,b2,z(i));
  endfor
endif
if (select==2)
  for i=1:len
    fa(i)=BetaLogDensity(a1,b1,a2,b2,z(i));
  endfor
endif
%
% plotting
%
subplot(1,3,3)
bar(bin_centers,counts/(n*bin_width),'FaceColor','blue','EdgeColor','black','BarWidth',1)
hold on
plot(z,fa,'-k','Linewidth',1)
hold off
title(strcat('density of Y=Ln(X_{1}/X_{2})','X_{1}~Beta(',num2str(a1),',',num2str(b1),')',...
',X_{2}~Beta(',num2str(a2),',',num2str(b2),')'),'fontsize',15)
xlabel('Ln(X_{1}/X_{2})','fontsize',15)
ylabel('Probability Density','fontsize',15)
