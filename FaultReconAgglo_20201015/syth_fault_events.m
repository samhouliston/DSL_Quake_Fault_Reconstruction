function [big_mat,clustIDX_true]=syth_fault_events(NUM_FAULTS,rseed,PER_NO,DENS_KM2,SPC_FAC, MODE_ELIP, DO_PLOT)
% Generates random points sampling synthetic faults in 3D
% NUM_FAULTS    : Number of faults;
% rseed         : Random seed generator for repeatable synthetics
% PER_NO        : Background noise percentage
% DENS_FAC      : Event density per km2
% SPC_FAC       : Spacing factor
% MODE_ELIP     : 0: Rectangular planes, 1: Gaussian ellipses 
% DO_PLOT       : 0: no plotting 1: marked plotting 2: unmarked plotting

% Y.Kamer 20201015

if(nargin<7)
    DO_PLOT     = 0;
end
rng(rseed); 
v_ST=(rand(NUM_FAULTS,1)-0.5)*180;  % STRIKE -90:+90;
v_DP= 45 + rand(NUM_FAULTS,1)*90;   % DIP 45-90;
v_L = 20 + rand(NUM_FAULTS,1)*20;   % LENGTH 20-40
v_W = 5  + rand(NUM_FAULTS,1)*10;   % WIDTH 5-15
v_D = rand(NUM_FAULTS,1)*10;        % DEPTH 0-10;

v_DN    = ones(NUM_FAULTS,1)*DENS_KM2; % Event density per km2;
LOC_ERR = 1; %thickness*2

OFV = rand(2,1)*SPC_FAC;
v_X = OFV(1)* rand(NUM_FAULTS,1)*NUM_FAULTS*mean(v_L);%1
v_Y = OFV(2)* rand(NUM_FAULTS,1)*NUM_FAULTS*mean(v_W);%1
v_Z = rand(NUM_FAULTS,1)*10;

cmap = jet(NUM_FAULTS+1);
cmap = cmap(randperm(NUM_FAULTS),:);
big_mat =   [];
clustIDX_true  =   [];
if(DO_PLOT) 
    figure;
    PTS_SIZE    = 10;
end
for i=1:NUM_FAULTS
    pnt_mat     = fault_plane(v_ST(i),v_DP(i),v_L(i),v_W(i),v_D(i),v_DN(i),rseed,LOC_ERR,MODE_ELIP,0);
    pnt_mat     = pnt_mat + ones(1,size(pnt_mat,1))'*[v_X(i) v_Y(i) v_Z(i)];
    big_mat     = [big_mat; pnt_mat];
    clustIDX_true      = [clustIDX_true; ones(size(pnt_mat,1),1)*i];
    
    if(DO_PLOT)
        if(DO_PLOT==1) 
            f_clr       = cmap(i,:);
        elseif(DO_PLOT==2)
            f_clr       = [1 1 1]*0.2;
        end
        scatter3(pnt_mat(:,1),pnt_mat(:,2),pnt_mat(:,3),PTS_SIZE,'ko','MarkerFaceColor',f_clr);
        hold on;
    end
end
NUM_NO  = round((PER_NO/100)*size(big_mat,1)/(1-(PER_NO/100)));
if(NUM_NO)
    clustIDX_true      = [clustIDX_true; ones(NUM_NO,1)*(i+1)];
    noise_mat   = rand(NUM_NO,3).*(ones(1,NUM_NO)'*range(big_mat)) + ones(1,NUM_NO)'*min(big_mat);
    
    if(DO_PLOT==1) 
        scatter3(noise_mat(:,1),noise_mat(:,2),noise_mat(:,3),PTS_SIZE,'.k');
    elseif(DO_PLOT==2)
        scatter3(noise_mat(:,1),noise_mat(:,2),noise_mat(:,3),PTS_SIZE,'ko','MarkerFaceColor',f_clr);
    end
    big_mat=[big_mat; noise_mat];
end

if(DO_PLOT) 
    daspect([1 1 1]);
    set(gca,'color',[1 1 1]*1)
    axis tight;
    TTL     = get(get(gca,'title'),'string');
    if(strcmp(TTL,''))
        TTL ={[num2str(NUM_FAULTS) ' Faults, ' num2str(numel(clustIDX_true)) ' Points ' ...
            ' Background Noise: ' num2str(PER_NO) '%']};
    else
        TTL{numel(TTL)+1,1} = [num2str(NUM_FAULTS) ' Faults, ' num2str(numel(clustIDX_true)) ' Points ' ...
            ' Background Noise: ' num2str(PER_NO) '%' ];
    end
    title(TTL,'fontsize',14);
end
rng('shuffle')
end