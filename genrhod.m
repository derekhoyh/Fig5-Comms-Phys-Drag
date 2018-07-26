function genrhod

% maxNumCompThreads(1);

Tvec=[70 130 190];

load('scaleddataTeq130K.mat');
load('rhodhomogdata-T130.mat','rhod');
factor=scaleddataTeq130K(176,2)/rhod(601);
clear scaleddataTeq130K;
clear rhod;


alpha=0.568231496731503;
echarge=1.6021766208e-19; 
hbar=1.054571800e-34;   

cons=echarge.^2./hbar;

for k=1:length(Tvec)
    load(['sigmadragvals-EMT-T' num2str(Tvec(k)) '.mat']);
    rhodfun=@(nq) interp1(ng1,-sigmadeff./(sigma1eff.*sigma2eff)./cons.*4.*alpha.^2.*pi,nq,'pchip',NaN);
    nplt=(-60:0.2:60);
    rhod=rhodfun(nplt).*factor;
%     plot(nplt,rhod,'LineStyle','-','LineWidth',3)
    save(['rhodemtdata-T' num2str(T) '.mat'],'nplt','rhod');
end


end