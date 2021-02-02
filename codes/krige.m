function [sim]=krige(tsdata,trdata, mja,mra,rmajor,rminor,eangle,nuggetj,sillj)
%This is a function for kriging estimate
%data file should have x-cor, ycor, variable to be kriged
tsxcor=tsdata(1:end,1);
tsycor=tsdata(1:end,2);
trdl=size(trdata);
trxcor=trdata(1:end,1);
trycor=trdata(1:end,2);
trsal=trdata(1:end,3);
%Guassian covariance function is used for this study
coj=nuggetj;
cj=sillj-nuggetj;
amajor=mja;
aminor=mra;
%search ellipsoid used
rmajor=rmajor;
rminor=rminor;
eangle=eangle;
aniratio=1;
% Transformation of co-ordinates along major axis
tsnxcor=tsxcor.*cos(eangle)+tsycor.*sin(eangle);
tsnycor=-tsxcor.*sin(eangle)+tsycor.*cos(eangle);
trnxcor=trxcor.*cos(eangle)+trycor.*sin(eangle);
trnycor=-trxcor.*sin(eangle)+trycor.*cos(eangle);
sim=[];
%ap=length(tsnxcor);
for i=1:size(tsdata,1)
   neighbour=[];
   for j=1:trdl(1,1)    
      dx1=tsnxcor(i)-trnxcor(j);
      dy1=tsnycor(i)-trnycor(j);
      dis=(dx1*dx1./rmajor^2)+ (dy1*dy1./rminor^2);
      if(dis<=1)
         neighbour=[neighbour; trnxcor(j) trnycor(j) trsal(j)];
      end      
   end
   if(size(neighbour,1)<=4)
      simuval=randn(1,1);
   else
       dnl=size(neighbour,1);
       xcorn=neighbour(1:end,1);
       ycorn=neighbour(1:end,2);
       saln=neighbour(1:end,3);
                %xcorn=fneighbour(1:end,1);
                %ycorn=fneighbour(1:end,2);
                %saln=fneighbour(1:end,3);
      %sample to sample covariance matrix
       sscov=[];
       h1=[];
       h2=[];
       for k=1:dnl
         for l=1:dnl
            h1(k,l)=sqrt((xcorn(k)-xcorn(l)).^2 +((ycorn(k)-ycorn(l))/aniratio).^2);
            if(h1(k,l)>amajor)
                yg=sillj;
            else
             yg=coj+(cj).*((1.5.*h1(k,l)./amajor)-(0.5.*h1(k,l).^3./amajor.^3));
            end;
            sscovm(k,l)=sillj-yg;
            sscov(k,l)=sscovm(k,l);
         end;
       end;
       %sscovs=[];
       spcov=[];
       lsp=[];
       for m=1:dnl
         hs1(m)=sqrt((tsnxcor(i)-xcorn(m)).^2 + ((tsnycor(i)-ycorn(m))/aniratio).^2);
         %hs2(m)=sqrt(((tsnxcor(i)-xcorn(m))/aniratio).^2 + (tsnycor(i)-ycorn(m)).^2);
         %spcov1(m)=cj*(1-(1-exp(-3*hs1(m)*hs1(m)/amajor^2)));
         if(hs1(m)>amajor)
             yg=sillj;
         else
             yg=coj+(cj).*((1.5.*hs1(m)./amajor)-(0.5.*(hs1(m)).^3./(amajor.^3)));
         end;
         %yg=coj+(cj).*(1-exp((-3.*hs1(m)^2)/amajor^2)); 
         spcov1(m)=sillj-yg;
         %spcov2(m)=cr*(1-(1-exp(-3*hs2(m)*hs2(m)/aminor^2)));
         %spcov(m)=spcov1(m)+spcov2(m);
         spcov(m)=spcov1(m);
       end
       spcov=spcov';
      % For Ordinary kriging
      %lsp=length(spcov);
      %spcov(lsp+1)=1;
      %spcovl=[];              
% Starting of Ordinary kriging
       %w=[];
       %we=[];
      % weight is a function for calculating weight matrix 
       sscov=[sscov,ones(length(sscov),1)];
       sscov=[sscov;[ones(1,(length(sscov)-1)),0]];
       spcov=[spcov;1];
      %lencov=length(sscov);
      %lambda=.01
      %imatrix=.1*eye(lencov);
      %eig(sscov);
       %spcov;
      
       w=pinv(sscov)*spcov;
      %w=weight(sscov,spcov);
      %we=w(1:(end-1));
      %sum(we);
      %wl=length(we);
      %flag=[];
      %for l=1:wl
       %   if(we(l)<0)
       %       flag=1;
       %   end;
       %end;
       %  if(length(flag)>0)
       %           minw=[]; 
       %           minw=abs(min(we));
       %           wn=we+minw;          
       %           wnew=wn./sum(wn);
       %           we=[];
       %           we=wnew;
       %       end;    
       w(end)=[];
       %       we
       estok=w'*saln;
      %spcovr=spcov(1:end-1); 
      %krigvar=sillj-we'*spcovr-w(end)
      % krigvar=sillj-spcov'*inv(sscov)*spcov;
     % if(krigvar<0)
     %     break;
     % end;
     %  krigstd=sqrt(krigvar);
       simuval=estok;%+krigstd*randn(1,1);      
      % if(krigvar<0)
        % simuval=randn(1,1);
     %  end;
     %clear dx1
     %clear dy1
    % clear dis
    % clear neighbour
    % clear h1
    % clear hs1
    % clear sscov
    % clear spcov
   %  clear w
    % clear neighbour
    % clear
   end
   sim=[sim;simuval];
end

sim
% End of ordinary kriging
  
      
     
      


