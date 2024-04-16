function hard_All=hardJoin(hard_All,hard_Inc)
% The ID vector size should match the new accumulated matrix, padd with
% zeros
    numALL      = numel(hard_All);
    numINC      = numel(hard_Inc);
    sizALL      = numel(hard_All{1}.ID);
    sizINC      = numel(hard_Inc{1}.ID);
    hard_All    = [hard_All hard_Inc];
    for i=1:numel(hard_All)
        if(i<=numALL)
            hard_All{i}.ID = [hard_All{i}.ID; false(sizINC,1)];
        else
            hard_All{i}.ID = [false(sizALL,1); hard_All{i}.ID];
        end
        
    end
end