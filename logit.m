%% Matching Logistic Curve to Data
clc;
clear all;
close all;

%data = xlsread('pop_colorado.xls','Sheet1'); p0=[0.1,160.,1.];
%data(:,2) = data(:,2)/1e6;
%ttl='Colorado Population Growth'; xlbl='Time (year)'; ylbl='Population (millions)';

%data = xlsread('pop_connecticut.xls','Sheet1'); p0=[0.1,16.,1.];
%data(:,2) = data(:,2)/1e6;
%ttl='Connecticut Population Growth'; xlbl='Time (year)'; ylbl='Population (millions)';

%data = xlsread('sales_iphone.xls','Sheet1'); p0=[1.,1.,1.];
%ttl='iPhone Sales World-Wide'; xlbl='Time (year)'; ylbl='Units Sold (millions)';

data = xlsread('sales_pcs.xls','Sheet1'); p0=[1.,1.,1.];
ttl='Personal Computer Sales World-Wide'; xlbl='Time (year)'; ylbl='Units Sold (millions)';

xdata = data(:,1);
ydata = data(:,2);

% Scale data
xdata = xdata-xdata(1);
ydata = ydata./ydata(1);


% MATLAB Levenberg-Marquardt Method
[p1 err] = lsqcurvefit(@logistic,p0,xdata,ydata)
p1 = p1.*[1,data(1,2),data(1,2)]
err = err*data(1,2)

% Own implementation of LM-method
%[p2 err] = lm_opt(p0,@fnc,@grad,xdata,ydata)

figure('position', [0, 1000, 800, 400]) 
t=linspace(0,xdata(length(xdata)),400);
%xlim=([xdata_orig(1), xdata_orig(length(xdata_orig))]);
y1=logistic(p1,t);
%y2=ydata_orig(1)*logistic(p2,t-xdata_orig(1));
%plot(xdata,ydata,'o',t,y1,'r')
plot(data(:,1),data(:,2),'o',t+data(1,1),y1,'r');
title(ttl,'FontSize',16);
xlabel(xlbl,'FontSize',14);
ylabel(ylbl,'FontSize',14);

function yprime = logistic(param,t)
    r=param(1);
    k=param(2);
    y0=param(3);
    yprime = y0.*k./(y0+(k-y0).*exp(-r.*t));
end
function f = error(x,xdata,ydata)
    N = length(xdata);
    r = x(1);
    k = x(2);
    y0 = x(3);
    sum = 0;
    for i = 1:N
        sum = sum + (ydata(i)-fnc_logistic(x,xdata(i))).^2;
    end
    f = sum;
end
function f = fnc(p,t)
	r=p(1);
    k=p(2);
    y0=p(3);
    f = y0.*k./(y0+(k-y0).*exp(-r.*t));
end
function g = grad(x,t)
    N = length(x);
    r = x(1);
    k = x(2);
    y0 = x(3);
    e = exp(r*t);
    denom = (k+y0*(e-1.)).^2;
    dr = - (k*t*y0*(k-y0)*e)/denom;
    dk = - (y0*y0*e*(e-1.))/denom;
    dy = - (k*k*e)/denom;
    g = [dr dk dy];
end