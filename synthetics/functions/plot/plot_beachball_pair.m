function plot_beachball_pair(stk1,dip1,rak1,stk2,dip2,rak2)

hf = figure(9001); clf; hold on; grid off; box off; axis equal
txt1 = beachball_mod(stk1,dip1,rak1,1,1,1,hf);
txt2 = beachball_mod(stk2,dip2,rak2,4,1,1,hf);

title({['left: ',txt1]; ['right: ',txt2]},'fontWeight','normal')

col = get(gcf, 'color');
set(gca, 'xtick',[], ...
         'ytick',[], ...
         'xcolor', col, ...
         'ycolor', col, ...
         'color', 'none');