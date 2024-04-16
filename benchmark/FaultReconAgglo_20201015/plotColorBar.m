cellFiles = {'SoCal_Faults\LOCAL_CH_5_FLTS_1407_150828_2102.mat',...
            'SoCal_Faults\GLOBAL_CH_5_FLTS_237_150913_0619.mat',...
            'SoCal_Faults\LOCAL_CH_30_FLTS_1482_150828_1806.mat',...
            'SoCal_Faults\GLOBAL_CH_30_FLTS_385_150912_0339.mat'};
figure;
for i=1:numel(cellFiles)
    subplot(2,4,i);
    ld = load(cellFiles{i}); 
    nKernel= numel(ld.param_All_FIN_LOC.w);
    valVec = log10(size(ld.chunk_mat_All,1)*ld.param_All_FIN_LOC.w./...
                prod(ld.param_All_FIN_LOC.sig(:,1:nKernel),1));
    colormap(jet(50));
    mnmx = minmax(valVec);
    caxis(mnmx);
    tc = ceil(mnmx(1)):1:floor(mnmx(2));
    h=colorbar;
    set(h,'ticks',tc);
    set(h,'ticklabels',num2str(10.^tc'));
    
end