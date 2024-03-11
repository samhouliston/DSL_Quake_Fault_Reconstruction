function fm = get_mode_FM(cat,plotme)

if nargin<2; plotme=0;
end

fprintf(1,'Selecting mode of FM angles, after selecting nodal plane with s1(2)>0.\n')
% clf; histogram(cat.stk); mode(cat.stk)
% clf; histogram(cat.dip); mode(cat.dip)
% clf; histogram(cat.rak); mode(cat.rak)

stk = cat.stk;
dip = cat.dip;
rak = cat.rak;

% flipme = cat.fm_s1(:,2)>0;
% 
% fm.comment = 'Flipping FMs with u(2)>0';
% stk(flipme) = cat.stk2(flipme);
% dip(flipme) = cat.dip2(flipme);
% rak(flipme) = cat.rak2(flipme);

fm.stk = mode(stk);
fm.dip = mode(dip);
fm.rak = mode(rak);

fm.string = sprintf('FM_{ref} = %i / %i / %i', fm.stk, fm.dip, fm.rak);

if plotme; beachball(fm.stk,fm.dip,fm.rak)
end
