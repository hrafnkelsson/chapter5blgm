function [cov0]=cov12(a,b)
 %m=length(a);
 %cov0=((mean((a-mean(a)).*(b-mean(b))))*m)/(m-1);
 cov01=cov(a,b);
 cov0=cov01(1,2);


