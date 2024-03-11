function FMs = get_FMs_from_rotation_around_fixed_axis(...
    nfm, ...
    stk0, ...
    dip0, ...
    rak0, ...
    drot, ...
    nrot)


% Sample rotation angles from a Von Mises distribution
oneSigma = drot*pi/180;
kappa    = 1/oneSigma^2;             % kappa =! 1/sigma^2 [1/rad]
rota_rad = circ_vmrnd(0, ...
                      kappa, ...
                      nfm);
rota = rota_rad*180/pi;
%clf; hold on; histogram(rota)

for ifm = 1:nfm

     [stk, dip, rak] = FMrot_known_axis(...
         stk0, ...
         dip0, ...
         rak0, ...
         rota(ifm), ...
         nrot(1), ...
         nrot(2), ...
         nrot(3) );

    FMs(ifm) = get_all_FM_vectors(stk, dip, rak);
end