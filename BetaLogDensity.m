% file name = BetaLogHyperGeom
% date of creation = 17/03/2024
%
% It calculates the density of Y=ln(X1/X2) in z, with
% X1~Beta(a1,b1), X2~Beta(a2,b2). It uses a finite summation.
%
function result = BetaLogDensity(a1,b1,a2,b2,z)
    a1=round(a1);
    b1=round(b1);
    a2=round(a2);
    b2=round(b2);
    result=0;
    if (z>0)
      for k=0:(b2-1)
        result = result + nchoosek(b2-1,k)*exp(-k*z)*beta(a1+a2+k,b1)*(-1)^k;
      endfor
      result=exp(-z*a2)*result;
      result=result/beta(a1,b1);
      result=result/beta(a2,b2);
    endif
    if (z<=0)
      for k=0:(b1-1)
        result = result + nchoosek(b1-1,k)*exp(k*z)*beta(a1+a2+k,b2)*(-1)^k;
      endfor
      result=exp(z*a1)*result;
      result=result/beta(a1,b1);
      result=result/beta(a2,b2);
    endif
end

