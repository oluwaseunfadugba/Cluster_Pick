function [c_all,icontour_all, numstart_all, numfin_all] = ...
    calc_P_nodal_planes(strike,dip,rake,vp,vs,dens, z_displ, r_displ, t_displ)
% This function determines the P-nodal planes for the strike, dip and rake.
% The strike, dip and rake can be a vector but of the same size.
% 
% Outputs: 
%            c will contain no of roles = 2 x no of FM
%            The other outputs will have no roles = no of FM
%
N = length(strike);
%
c_all = zeros(2*N,5000);
icontour_all = zeros(N,1);
numstart_all = zeros(N,5);
numfin_all = zeros(N,5);
%
for no_of_FM = 1:N
        eps=+1.0;
        %   calculate the moment tensor for a point dislocation
        m=dismom(strike(no_of_FM),dip(no_of_FM),rake(no_of_FM));
        %
        %  Compute P wave radiation pattern
        az=linspace(0,360,1000);
        inc=linspace(0,90,300);
        pamp=prad(az,inc,vp,vs,dens,eps,m, z_displ);
        
        %
        %  find the zero amplitude contor, if it exists
        c=contourc(az,inc,pamp,[0 0]);
        %
        %  find all of the contours
        [ix,iy]=size(c);
        icontour=0;
        num=0;
        numstart=0;
        numfin=0;
        if iy > 0
        num1=c(2,1);
        if (num1 > 1)
            icontour=1;
            num(1)=num1;
            numstart(1)=2;
            numfin(1)=numstart(1)+num1-1;
        else
            icontour=0;
        end
        size(numfin);
        istop=0;
        ik=0;
        if (icontour > 0)
            test=numfin(1);
            while (test < iy)
                ik=ik+1;
                icontour=icontour+1;
                num(ik+1)=c(2,numfin(ik)+1);
                numstart(ik+1)=numfin(ik)+2;
                numfin(ik+1)=numstart(ik+1)+num(ik+1)-1;
                test=numfin(ik+1);
            end
        end
        end
 
        [nsr, nsc] = size(numstart);
        [nfr, nfc] = size(numfin);
        icontour_all(no_of_FM,:) = icontour;
        numstart_all(no_of_FM,1:nsc) = numstart;
        numfin_all(no_of_FM,1:nfc) = numfin;
        
        [mmm,nnn]=size(c);
        a = 2*no_of_FM - 1;        b = 2*no_of_FM;
    
        c_all(a:b,1:nnn) = c;
        
end
