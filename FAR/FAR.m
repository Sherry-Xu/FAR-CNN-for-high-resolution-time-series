function y=FAR(data,str,str_res,str_yhat,mn,traint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: this function is to perform FAR for out of sample FAR
%forecast. The function will use supporting functions in the folder and
%subfolders in fdaM
%
%Input: data   : data matrix of size 24*days
%       str    : a text string indicating the name of the file that the results will be saved, e.g. 'FARout.mat'
%       mn     : sieve truncation, e.g. mn=11, s.t. num_of_basis=2*mn+1=23
%       traint : number of days for the rolling window, e.g. 28 days (4 weeks)
%
%Output: FARout.mat : save all the forecast errors
%        y       : FAR forecasts
%
%Author: Wee Song CHUA, 22 Jan 2018

addpath(genpath('fdaM'))

gasflow=data;
gasflow0=gasflow; % can subtract first hour gas flow here for change in gasflow from 6am, but need add back at line 84

% total number of days
T=size(gasflow,2); 
 
%Create the fd objects by Fourier transformation
CAnbasis=2*mn+1;
CAfourierb=create_fourier_basis([0,1],CAnbasis,1);
 
% x-axis
time=0:23;
timescaled=zeros(24,1);
for i=1:24
    timescaled(i)=time(i)/23;
end

% matrix to store Fourier coefficients for y
yt_CWS=zeros(CAnbasis,T);
 
% extract Fourier coefficients for y
for t=1:T
    % y
    CAfourierfdy=smooth_basis(timescaled,gasflow(:,t),CAfourierb);
    yt_CWS(:,t)=getcoef(CAfourierfdy);
end
% 1-step ahead forecast
%traint=28;
newfT=T-traint;
xaxis=timescaled;
gasflowHat=zeros(24,newfT);
for ft=1:newfT
    % training set
    trainy=yt_CWS(:,ft:ft+traint-1);
     
    %estimation using FAR estimator
    hc0y= ( ((sum(trainy(1,1:end-1)))*(sum(trainy(1,2:end)))-(traint-1)*(sum(trainy(1,1:end-1).*trainy(1,2:end))))) / ((sum(trainy(1,1:end-1)))^2-(traint-1)*(sum(trainy(1,1:end-1).^2)));
    hcky=zeros(mn,1);
    for k=1:mn
        hcky(k)= ( (sqrt(2)*(traint-1)*sum(trainy(2*k+1,2:end).*trainy(2*k+1,1:end-1)+trainy(2*k,2:end).*trainy(2*k,1:end-1))-sqrt(2)*((sum(trainy(2*k,2:end)))*(sum(trainy(2*k,1:end-1)))+(sum(trainy(2*k+1,2:end)))*(sum(trainy(2*k+1,1:end-1))))) ) / ((traint-1)*sum(trainy(2*k,1:end-1).^2+trainy(2*k+1,1:end-1).^2)-((sum(trainy(2*k,1:end-1)))^2+(sum(trainy(2*k+1,1:end-1)))^2)) ;
    end
 
    hp0=(sum(trainy(1,2:end))-hc0y*sum(trainy(1,1:end-1))) / (traint-1);
    hs0=(sum((trainy(1,2:end)-repmat(hp0,1,traint-1)-hc0y*trainy(1,1:end-1)).^2)) / (traint-1);
    hpk=zeros(mn,1);  
    hqk=zeros(mn,1);
    hsk=zeros(mn,1);    
    for k=1:mn
       hpk(k)=( sum(trainy(2*k+1,2:end))-(1/sqrt(2))*hcky(k)*sum(trainy(2*k+1,1:end-1)) ) / (traint-1);
       hqk(k)=( sum(trainy(2*k,2:end))-(1/sqrt(2))*hcky(k)*sum(trainy(2*k,1:end-1)) ) / (traint-1);
       hsk(k)=(sum( (trainy(2*k,2:end)-repmat(hqk(k),1,traint-1)-(1/sqrt(2))*hcky(k)*trainy(2*k,1:end-1)).^2 + (trainy(2*k+1,2:end)-repmat(hpk(k),1,traint-1)-(1/sqrt(2))*hcky(k)*trainy(2*k+1,1:end-1)).^2 )) / (2*(traint-1));
    end
     
    % forecast
    tempYHat=zeros(size(xaxis));
    for tau=1:24
        est=0;
        for k=1:mn
            est=est+(hqk(k)+(1/sqrt(2))*hcky(k)*trainy(2*k,traint))*(sqrt(2)*sin(2*pi*k*xaxis(tau)))+(hpk(k)+(1/sqrt(2))*hcky(k)*trainy(2*k+1,traint))*(sqrt(2)*cos(2*pi*k*xaxis(tau)));
        end
        est=est+hc0y*trainy(1,traint)+hp0; %est is the FAR estimate of Y_t(tau)
        tempYHat(tau)=est;
    end
     
    % FAR estimated gasflow
    gasflowHat(:,ft)= tempYHat;
end
 
%RMSE
resd=gasflowHat-gasflow0(:,(traint+1:end));
RMSE=sqrt(mean2(resd.^2))
 
%MAPE (remove zero gasflow in the calculation)
gasflow2=gasflow0(:,traint+1:end);
gasflowVec=gasflow2(:);
zeropos=find(gasflowVec==0);
resdVec=resd(:);
resdVec(zeropos)=[];
gasflowVec(zeropos)=[];
MAPE=mean(abs(resdVec./gasflowVec))
 
%R2
gasflowbar=mean2(gasflow0);
gasflow2=gasflow0(:,traint+1:end);
gasflowVec=gasflow2(:);
SST=sum((gasflowVec-gasflowbar).^2);
gasflowHatVec=gasflowHat(:);
SSE=sum((gasflowHatVec-gasflowVec).^2);
R2=1-SSE/SST %

% by hr
RMSEhr=sqrt(mean(resd.^2,2));
MAPEhr=zeros(24,1);
gasflow2=gasflow0(:,traint+1:end);
for i=1:24
    gasflowVec=gasflow2(i,:);
    zeropos=find(gasflowVec==0);
    resdVec=resd(i,:);
    resdVec(zeropos)=[];
    gasflowVec(zeropos)=[];
    MAPEhr(i)=mean(abs(resdVec./gasflowVec));
end
RMSEhr
MAPEhr
y=gasflowHat;
save(str)
csvwrite(str_res,resd)
csvwrite(str_yhat,y)
end