function [F]=gaussian2D(x,coords)
%% x = [Amp, x0, sigmax, y0, sigmay, theta]

X=coords(:,:,1);
Y=coords(:,:,2);

a= (cos(x(6))^2)/(2*x(3)^2) + (sin(x(6))^2)/(2*x(5)^2);
b= - sin(2*x(6))/(4*x(3)^2) + sin(2*x(6))/(4*x(5)^2);
c= (sin(x(6))^2)/(2*x(3)^2) + (cos(x(6))^2)/(2*x(5)^2);

F=x(1)*exp(-(a.*(X-x(2)).^2 - 2.*b.*(X-x(2)).*(Y-x(4)) + c.*(Y-x(4)).^2));

end