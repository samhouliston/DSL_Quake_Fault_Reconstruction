function plot_mix_model(param,Pts_Train,Pts_Val,prf,BKG_ON,BBOX_ON)
if(nargin<6)
    BBOX_ON=0;
    if(nargin<5)
        BKG_ON=1;
    end
end
figure;
nKernels    = size(param.m,2);
clr         = parula(nKernels);
%clr         = gray(nKernels+1);
%clr         = flipud(clr(1:nKernels,:));
padC        = ceil(log10(nKernels));
vec_dens    = log10(param.w./prod(param.sig(:,1:nKernels),1));
%{
tmp             = prod(param.sig(:,1:nKernels),1);
vec_dens        = log10(tmp);
%}

if(all(param.bkg))
    norm_w      = (param.w-min(param.w))./range(param.w);   
    %norm_w      = (vec_dens-min(vec_dens))./range(vec_dens);   
else
    %norm_w      = (param.w-min(param.w(~param.bkg)))./range(param.w(~param.bkg));
    norm_w      = (vec_dens-min(vec_dens(~param.bkg)))./range(vec_dens(~param.bkg));
end
%}
%

%ID          = [1067,1181,1149,2093,2590,2392,2717,2704,2142,454,2621,1217,2699,1200;];
%ID          = [1217];
%norm_w      = (param.w-min(param.w(ID)))./range(param.w(ID));
norm_w(isnan(norm_w))=0;
facs = [1 2 3 4
        5 6 7 8
        4 3 6 5
        3 2 7 6
        2 1 8 7
        1 4 5 8];
PLOT_alpha  =  0.50;%25
    for j=1:nKernels
        strCINFO= sprintf(['Cluster %0' num2str(padC) 'd : w= %6.4f'],j,param.w(j));
        if(~BKG_ON && param.bkg(j))
            continue;
        end
        if(param.bkg(j))
            DrawCuboidYK(param.m(:,j),param.covar(:,:,j));
            hold on;
            %patch('Faces',facs,'Vertices',param.bbox(:,:,j),'FaceColor',[1 1 1]*0.5,'FaceAlpha',0.25);
        else
            clr_p   = clr(floor(norm_w(j)*(nKernels-1))+1,:);
            %
            h       = plot_gaussian_ellipsoid(param.m(:,j), param.covar(:,:,j),sqrt(12)/2);%
            set(h,'edgecolor','none',...
                'facealpha',PLOT_alpha,'FaceColor',clr_p,...
                'displayname',strCINFO);
            %}
            %DrawCuboidYK(param.m(:,j),param.covar(:,:,j)*12,clr_p,PLOT_alpha);
        end
        if(BBOX_ON && ~param.bkg(j))
            DrawCuboidYK(param.m(:,j),param.covar(:,:,j)*12);
            patch('Faces',facs,'Vertices',param.bbox(:,:,j),'FaceColor',[1 1 1]*0.5,'FaceAlpha',0.25);
        end
        
        %norm_w(j)*0.4+0.1
    end

h = camlight('right');
lighting gouraud
material dull
daspect([1 1 1])
axis tight;
view(3);
grid on;
%
h   = colorbar('southoutside');
mnmx = minmax(vec_dens);
tc  = ceil(mnmx(1)):1:floor(mnmx(2));
set(h,'ticks',tc);
caxis(mnmx);
%}
%{
for i = 1:20;
 camorbit(10,0)
 camlight(h,'right')
 drawnow;
 pause(0.25)
end
%}
%{
colorbar;
if(nKernels>1)
    caxis([min(param.w) max(param.w)]);
end
%}
if(~isempty(Pts_Train))
    plot3(Pts_Train(1,:),Pts_Train(2,:),Pts_Train(3,:),'.k');
end
%scatter3(Pts_Val(1,:),Pts_Val(2,:),Pts_Val(3,:),10,'k');
hold on;
%{
NUM_PER     = 1000;
per_pdf     = zeros(size(Pts_Val,2),NUM_PER);
for p=1:NUM_PER
    Pts_per=randn(size(Pts_Val'))*0.05;
    for j=1:nKernels
        switch param.mix_type
            case 'gaus'
                per_pdf(:,p)     = per_pdf(:,p) + mvnpdf(Pts_Val'+Pts_per,param.m(:,j)',param.covar(:,:,j))*param.w(j);
            case 'stud'
                per_pdf(:,p)     = per_pdf(:,p) + mvtpdf(((Pts_Val+Pts_per)-param.m(:,j)*ones(1,nPnt))',param.covar(:,:,j),param.v(j))*param.w(j);
        end
    end
end
s_per_pdf=sum(per_pdf);
lb=min(s_per_pdf);
ub=max(s_per_pdf);
%}
%{
ViewAZ=-35;%-75;%-90;%-105
ViewEL=55;%80;%90;%35
setGCA=['set(gca,''YDir'',''reverse'',''ZDir'',''normal'',''XDir'',''normal'',''Color'',ones(1,3)*1);'...
        'grid on;'];
view(ViewAZ,ViewEL);
eval(setGCA);
ylim([0 60]);
xlim([0 110]);
%}
%if(isfield(param,'BIC'))
    title(['Clusters: ' num2str(numel(param.w)) ', BKG = ' num2str(sum(param.w(param.bkg)),'%.2f')],'fontsize',14);
%end
axis tight;
set(gca, 'color', [1 1 1]*0.5)
view(2);
end

