function p=hyge2cdf(k,n,k1,n1)
  p=0;
  for ki=0:k
      p=p+hyge2pdf(ki,n,k1,n1); % See hyge2pdf code in repository
  end
  p = max(0,min(1,p));
end