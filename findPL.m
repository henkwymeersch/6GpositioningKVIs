function PL=findPL(sigma_pos,TIR)
  r_min=0;
  r_max=10*sigma_pos;
  tol=1e-8;
  maxIter=100;
  if (sigma_pos==+Inf)
      PL=+Inf; 
      return;
  end    
  for iter=1:maxIter
      r_mid=(r_min+r_max)/2;
      
    
      C_min=2*qfunc(r_min/sigma_pos)-TIR;
      C_max=2*qfunc(r_max/sigma_pos)-TIR;
      C_mid=2*qfunc(r_mid/sigma_pos)-TIR;
    
      if abs(C_mid) < tol || (r_max - r_min) / 2 < tol
        PL= r_mid;
        return 
      end
    
      if (C_min*C_mid <0)
          r_max=r_mid;
      else
          r_min=r_mid;
      end
  end
  
  PL=r_mid;