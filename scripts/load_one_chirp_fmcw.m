vdat = fmcw_load('../data/sample/thwaites/DATA2023-02-12-0437.DAT'); %defaults to just the first burst
[rc,~,spec_cor,spec] = fmcw_range(vdat,2,2500,@blackman);
plot(rc,20*log10(abs(spec_cor(1,:))))
%save('../data/sample/thwaites/DATA2023-02-12-0437_p4.mat')
figure
plot(vdat.vif(1,:))