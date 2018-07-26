function sigmadnrhodEMT

maxNumCompThreads(1);

load('scaleddataTeq130K.mat');
load('rhodhomogdata-T130.mat','rhod');
factor=scaleddataTeq130K(176,2)/rhod(601);
clear scaleddataTeq130K;
clear rhod;

Tvec=[70 130 190];
nrms1=6;
nrms2=6;
eta=@(x,T) (-0.99+0.005.*(T-70).*theta(T-70)).*exp(-0.5.*x).*theta(-(-0.99+0.005.*(T-70)));

ng1=(-60:1:60);
ng2=zeros(1,length(ng1));
figure;hold on;box on;

alpha=0.568231496731503;
echarge=1.6021766208e-19; 
hbar=1.054571800e-34;   

cons=echarge.^2./hbar;

for k=1:length(Tvec)
    T=Tvec(k);
load(['monolayersigmas_T' num2str(T) 'K.mat'],'n', 'sigmamono');
nplus=n(2:length(n));
nfull1=[-fliplr(nplus) n];
sigmaplus=sigmamono(2:length(sigmamono));
sigmafull1=[fliplr(sigmaplus) sigmamono];
monocond=@(x) interp1(nfull1,sigmafull1,x,'pchip',NaN);

load(['draggrid-T' num2str(T) '.mat'])
sigmaDinterp=@(n1,n2) interp2(nA,nP,sigmaDgrid,n1,n2,'spline',NaN);%change nimp from 5x15^10 to 10x10^10

% Now do the EMT.
fun=@(n,m,s) (n-m).^2./s.^2; 
P=@(n1,n2,ng1,ng2) exp(- 1./(2.*(1-( eta(sqrt(abs(ng1.*ng2)),T) ).^2)) .* ...
    (fun(n1,ng1,nrms1) + fun(n2,ng2,nrms2) -2.*eta(sqrt(abs(ng1.*ng2)),T).*(n1-ng1).*(n2-ng2)./(nrms1.*nrms2) ) ); 

sigmadeff=zeros(1,length(ng1));
sigma1eff=zeros(1,length(ng1));
sigma2eff=zeros(1,length(ng2));
save(['sigmadragvals-EMT-T' num2str(T) '.mat'])

for j=1:length(ng1)
    diary('sigmaDemt.txt')
    tic
    sigma1eff(j)=EMTmono(ng1(j),nrms1,monocond);
%     sigma2eff(j)=EMTmono(ng2(j),nrms2,T);
    sigma2eff(j)=EMTmono(ng2(j),nrms2,monocond);
    
    numerator= integral2(@(n1,n2) P(n1,n2,ng1(j),ng2(j)) .* sigmaDinterp(n1,n2) .* sigma1eff(j) ./ ... 
        ( (sigma1eff(j)+monocond(n1)).*(sigma2eff(j)+monocond(n2)) ),ng1(j)-5.*nrms1,ng1(j)+5.*nrms1,ng2(j)-5.*nrms2,ng2(j)+5.*nrms2,'Method','iterated','RelTol',1e-4);
    denominator=integral2(@(n1,n2) P(n1,n2,ng1(j),ng2(j)) .* monocond(n1) ./ ... 
        ( (sigma1eff(j)+monocond(n1)).*(sigma2eff(j)+monocond(n2)) ),ng1(j)-5.*nrms1,ng1(j)+5.*nrms1,ng2(j)-5.*nrms2,ng2(j)+5.*nrms2,'Method','iterated','RelTol',1e-4);
    sigmadeff(j)=numerator./denominator;
    toc
    save(['sigmadragvals-EMT-T' num2str(T) '.mat'],'-append','sigmadeff','sigma1eff','sigma2eff')
    
    
    diary off;
end

% rhodfun=@(nq) interp1(ng1,-sigmadeff./(sigma1eff.*sigma2eff-(sigmadeff.*4.*alpha.^2.*pi).^2)./cons.*4.*alpha.^2.*pi,nq,'pchip',NaN);
%     nplt=(-60:1:60);
%     rhod=rhodfun(nplt).*factor;
% %     plot(nplt,rhod,'LineStyle','-','LineWidth',3)
%     save(['rhodemtdata-T' num2str(T) '.mat'],'nplt','rhod');

end

% legend({'$T=40K$' '$T=70K$' '$T=100K$' '$T=130K$' '$T=160K$' '$T=190K$' '$T=220K$' '$T=240K$' '$T=300K$'}, 'Interpreter', 'latex','FontSize',20, 'Location','NorthEast','Orientation','Vertical')

% print('-dpdf','Fig4c.pdf')


end