stk = 20;
dip = 40;
rak = 40;

rota = 1:2:180;

r1 = 0;
r2 = 0;
r3 = 1;

xlm = [0 max(rota)];
ylm = [0 max(rota)];

% Rotation around fixed  axis
% Rotation around random axis
% With    FP ambiguity
% Without FP ambiguity
nrot  = numel(rota);
a_fix       = zeros(nrot,1);
a_rnd       = zeros(nrot,1);
a_fix_woFPA = zeros(nrot,1);
a_rnd_woFPA = zeros(nrot,1);
for irot = 1:nrot
    
    [stk_fix,dip_fix,rak_fix] = FMrot_known_axis (stk,dip,rak,rota(irot),r1,r2,r3);
    [stk_rnd,dip_rnd,rak_rnd] = FMrot_random_axis(stk,dip,rak,rota(irot));

    [a_fix(irot),~,~] = kagan([stk,dip,rak],[stk_fix,dip_fix,rak_fix]);
    [a_rnd(irot),~,~] = kagan([stk,dip,rak],[stk_rnd,dip_rnd,rak_rnd]);
    
    a_fix_woFPA(irot) = get_fault_plane_rot(stk,dip,rak,stk_fix,dip_fix,rak_fix);
    a_rnd_woFPA(irot) = get_fault_plane_rot(stk,dip,rak,stk_rnd,dip_rnd,rak_rnd);
end


figure(3010); clf; 

s1 = subplot(2,1,1); hold on; grid on; box on;
plot([0 180],[0 180],'-','color',[.7 .7 .7])
p1 = plot(rota,a_rnd      ,'s','markerEdgeColor',[.1 .1 .1], 'markerFaceColor','r');
p2 = plot(rota,a_rnd_woFPA,'o','markerEdgeColor',[.1 .1 .1], 'markerFaceColor',[.4 .4 .4]);
tstring = sprintf('Using random rotation axis   -   FM: %i / %i / %i',stk,dip,rak);
title(tstring,'fontWeight','normal')
ylabel('Computed rotation angle [deg]')
set(gca,'xlim',xlm,'ylim',ylm)

s2 = subplot(2,1,2); hold on; grid on; box on;
plot([0 180],[0 180],'-','color',[.7 .7 .7])
p1 = plot(rota,a_fix      ,'s','markerEdgeColor',[.1 .1 .1], 'markerFaceColor','r');
p2 = plot(rota,a_fix_woFPA,'o','markerEdgeColor',[.1 .1 .1], 'markerFaceColor',[.4 .4 .4]);
tstring = sprintf('Using fixed rotation axis: r_i = %i/%i/%i',r1,r2,r3);
title(tstring,'fontWeight','normal')
xlabel('Imposed rotation angle [deg]')
ylabel('Computed rotation angle [deg]')
set(gca,'xlim',xlm,'ylim',ylm)


subplot(s1)
hl = legend([p1;p2],'Using FP ambiguity','Not using FP ambiguity');
hl.Location = 'northWest';


if false
    figName = '~/fig/proj/fm/new/FM_plane_ambiguity';
    set(gcf,'PaperPosition',[0 0 6 8])
    print('-dpng',figName)
end