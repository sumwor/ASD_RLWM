function p = x2testOverlap(var1,var2,N)
        n1 = sum(var1&var2); n2 = sum(~var1&var2);
       n3 = sum(~var2&var1); n4 = sum(~var1 & ~var2); 
       % Pooled estimate of proportion
       p1 = (n1+n3)/N; p2 = (n1+n2)/N;
       e1 = p1*p2*N; e2 = (1-p1)*p2*N;
       e3 = (1-p2)*p1*N; e4 = (1-p1)*(1-p2)*N;
       
       % Chi-square test, by hand
       observed = [n1 n2 n3 n4];
       expected = [e1 e2 e3 e4];
       chi2stat = sum((observed-expected).^2 ./ expected);
       p = (1 - chi2cdf(chi2stat,1))/2;
end
