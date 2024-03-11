function [stk2, dip2, rak2, fDipDir2] = ComputeSecondPlane(stk1, dip1, rak1)
      
% Do not modify strike (Dip and rake need to be modified in order to
% minimize computational errors)
dip1 = dip1 + 0.000001;
rak1 = mod((rak1 + 0.000001) + 360, 360);
[stk2, dip2, rak2, fDipDir2, nError] = focal_pl2pl(stk1, dip1, rak1);