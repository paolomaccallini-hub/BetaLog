% file name = BetaLogHyperGeom
% date of creation = 17/03/2024
%
% It calculates the density of Y=ln(X1/X2) in z, with
% X1~Beta(a1,b1), X2~Beta(a2,b2). It uses the hypergeometric function.
%
function result = BetaLogDensityHG(a1,b1,a2,b2,z)
  result=0;
  if (z>0)
    result=gsl_sf_hyperg_2F1(1-b2,a1+a2,a1+a2+b1,exp(-z));
    result=exp(-z*a2)*result;
    result=result/beta(a1,b1);
    result=result/beta(a2,b2);
    result=result*beta(a1+a2,b1);
  endif
  if (z<=0)
    result=gsl_sf_hyperg_2F1(1-b1,a1+a2,a1+a2+b2,exp(z));
    result=exp(z*a1)*result;
    result=result/beta(a1,b1);
    result=result/beta(a2,b2);
    result=result*beta(a1+a2,b2);
  endif
end
