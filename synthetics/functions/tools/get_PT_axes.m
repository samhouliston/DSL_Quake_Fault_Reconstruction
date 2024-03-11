function cat = get_PT_axes(cat)

neq = numel(cat.lat);
cat.p_azi = zeros(neq,1);
cat.p_dip = zeros(neq,1);
cat.t_azi = zeros(neq,1);
cat.t_dip = zeros(neq,1);
cat.b_azi = zeros(neq,1);
cat.b_dip = zeros(neq,1);

for ieq = 1:neq
    
    [paz,pdip,taz,tdip,baz,bdip] = sdr2ptb(cat.stk(ieq),cat.dip(ieq),cat.rak(ieq));
    
    cat.p_azi(ieq) = paz;
    cat.p_dip(ieq) = pdip;
    cat.t_azi(ieq) = taz;
    cat.p_dip(ieq) = tdip;
    cat.b_azi(ieq) = baz;
    cat.b_dip(ieq) = bdip;
end