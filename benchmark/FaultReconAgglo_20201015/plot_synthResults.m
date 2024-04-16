%load ws_local_global;
vecNoisePrc     = [5 10 20];
vecElips        = [0 1];
numN        = numel(vecNoisePrc);
clrN    = [0 0 0; 0.25 0.25 0.25; 0.5 0.5 0.5; ];
figure;
for e=[0 1]
subplot(1,2,e+1);
    for p=[5 1]
        for i=1:numN 
        idSel = matNDE(:,1) == vecNoisePrc(i) & ...
            matNDE(:,3) == e ;
            switch p
                case 1
                    pChar = 'o--';
                    critN = 'Local';
                case 5
                    pChar = '^-';
                    critN = 'Global';
            end
            plot(matNDE(idSel,2),...
                 outRandLG_elps(idSel,p),...
                 pChar,'color',clrN(i,:),...
                 'linewidth',1,...
                 'markersize',5,...
                 'displayName',[critN ', noise: ' num2str(vecNoisePrc(i)) ' %']);
            hold on;
        end
    end
    grid on;
    switch e
        case 0
            ttlStr = 'Fault type: Rectangle';
            ylabel('Rand index');
        case 1
            ttlStr = 'Fault type: 2D Gaussian';
            legend('Location','southeast')
    end
    title(ttlStr);
    ylim([0.4 1]);
    pbaspect([1.25 1 1])
    xlabel('Density (points/km^2)');
end
set(findall(gcf,'-property','FontSize'),'FontSize',8);
%print(gcf,'-dpng','-r300','synth_res.png');