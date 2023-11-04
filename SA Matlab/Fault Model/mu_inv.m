function x = mu_inv(y,mu)
%inverse of mu-law
x = (((1+mu).^abs(y)- 1)/mu).*sign(y);
