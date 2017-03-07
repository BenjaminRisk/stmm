function y = biweight(x)
    y=(15/16).*(1-x.^2).^2.*(abs(x)<=1);
end