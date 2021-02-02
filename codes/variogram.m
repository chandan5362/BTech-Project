function [llag,gamlag]=variogram(xcor,ycor,data,nlag,minlag,laginv,azm,atol,maxbandw)
%This function is for variogram calculation
dl=length(data);
lowangle=azm-atol;
upangle=azm+atol;
gamlag=0;
n=0;
for k=1:nlag
    gamlag(k)=0;
    n(k)=0;
end;
for i=1:dl
    for j=i:dl
        if (i~=j)
            d=sqrt((xcor(i)-xcor(j))^2 + (ycor(i)-ycor(j))^2);
            if((xcor(j)-xcor(i))==0)
                theta=90;
            else
              theta=180.*atan((ycor(j)-ycor(i))./(xcor(j)-xcor(i)))./pi;
            end
            for k=1:nlag
                %lag=minlag+k*xlag;
                llag(k)=minlag+(k-1)*laginv;
                hlag(k)=llag(k)+laginv;
                 b=abs(d.*sin(theta-azm));
                if((llag(k)<=d & d<hlag(k))&(lowangle<=theta &theta<upangle)&(b<=maxbandw))
                        gamlag(k)=gamlag(k)+(data(i)-data(j)).^2;
                        n(k)=n(k)+1;
                    end;
                    end;
                end;    
            end;        
    end;
    for k=1:nlag
           gamlag(k)=gamlag(k)./(2.*n(k));
        end;
gamlag;
llag;
