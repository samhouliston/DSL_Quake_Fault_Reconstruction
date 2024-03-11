function fm = get_all_FM_vectors(stk,dip,rak)


% First plane .............................................................
fm.stk = stk;
fm.dip = dip;
fm.rak = rak;

[n1,n2,n3,s1,s2,s3] = fm2nsb_vectors(stk,dip,rak);
fm.s1 = [s1,s2,s3];
fm.n1 = [n1,n2,n3];
fm.b1 = cross(fm.n1, fm.s1);
fm.t1 = fm.n1+fm.s1;
fm.p1 = fm.n1-fm.s1;

[paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(stk,dip,rak);
fm.p1_azi = paz;
fm.p1_dip = pdip;
fm.t1_azi = taz;
fm.t1_dip = tdip;
fm.b1_azi = baz;
fm.b1_dip = bdip;



% Second plane ............................................................
[stk2,dip2,rak2,~,~] = focal_pl2pl(stk,dip,rak);
fm.stk2 = stk2;
fm.dip2 = dip2;
fm.rak2 = rak2;

[n1,n2,n3,s1,s2,s3] = fm2nsb_vectors(stk2,dip2,rak2);
fm.s2 = [s1,s2,s3];
fm.n2 = [n1,n2,n3];
fm.b2 = cross(fm.n2, fm.s2);
fm.t2 = fm.n2+fm.s2;
fm.p2 = fm.n2-fm.s2;


[paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(stk2,dip2,rak2);
fm.p2_azi = paz;
fm.p2_dip = pdip;
fm.t2_azi = taz;
fm.t2_dip = tdip;
fm.b2_azi = baz;
fm.b2_dip = bdip;