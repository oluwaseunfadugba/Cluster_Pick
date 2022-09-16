function preferred_FM = calc_preferred_FM(best_strikes,best_dips,best_rakes,plot_FM)
    
ang_thresh = 30;
    
N_fm = length(best_strikes);

for i = 1:2
    if i == 2 
        for j =1:N_fm
            [best_strikes(j),best_dips(j),best_rakes(j)] = ...
                calc_AuxPlane(best_strikes(j),best_dips(j),best_rakes(j));
        end     
    end
    
    [strike_avg, dip_avg, rake_avg,ang_rms,perc_win_angthresh] = ...
        average_FMs(best_strikes,best_dips,best_rakes,ang_thresh);

    strike_avgn(i) = round(strike_avg,0);
    dip_avgn(i) = round(dip_avg,0);
    rake_avgn(i) = round(rake_avg,0);
    ang_rmsn(i) = round(ang_rms,2);
    perc_win_angthreshn(i) = round(perc_win_angthresh,2);
    
end

[ang_rms,I] = min(ang_rmsn);

strike1 = strike_avgn(I);
dip1 = dip_avgn(I);
rake1 = rake_avgn(I);
perc_win_30deg = perc_win_angthreshn(I);

[strike2,dip2,rake2] = calc_AuxPlane(strike1,dip1,rake1);

strike2 = round(strike2,0);
dip2 = round(dip2,0);
rake2 = wrapTo360(round(rake2,0));

if strike1 <= strike2
    preferred_FM = [strike1 dip1 rake1 strike2 dip2 rake2 ang_rms perc_win_30deg];
else
    preferred_FM = [strike2 dip2 rake2 strike1 dip1 rake1 ang_rms perc_win_30deg];
end

if plot_FM ==1
    fig = figure ('Position', get(0, 'Screensize'),'visible', 'on', 'color','w'); clf 

    plot_P_nodal_planes(preferred_FM(1),preferred_FM(2),preferred_FM(3),...
        best_strikes,best_dips,best_rakes)

    print(fig, 'fig_avg_FMs.png', '-dpng'); 

end

end


%% Helper functions
function [strike_avg, dip_avg, rake_avg,ang_rms,perc_win_angthresh] = ...
    average_FMs(best_strikes,best_dips,best_rakes,ang_thresh)

best_strikes_orig = real(best_strikes);
best_dips_orig = real(best_dips);
best_rakes_orig = real(best_rakes);

best_strikes_tmp = real(best_strikes);
best_dips_tmp = real(best_dips);
best_rakes_tmp = real(best_rakes);

max_ang = ang_thresh;

while max_ang>=ang_thresh
    N = length(best_strikes_tmp);
    
    [s_avg,d_avg,r_avg] = MECH_AVG(best_strikes_tmp,best_dips_tmp,best_rakes_tmp);

    ang_btw_sdr = zeros(N,1);

    for i = 1:N
        sdr_avg = [s_avg d_avg r_avg];
        sdr_i = [best_strikes_tmp(i) best_dips_tmp(i) best_rakes_tmp(i)];

        ang_tmp = kagan(sdr_avg,sdr_i);
        ang_btw_sdr(i,1) = abs(ang_tmp);
    end

    max_ang = max(ang_btw_sdr);
    
    if max_ang >= ang_thresh
        N;
        indexx=find(ang_btw_sdr == max_ang);
        best_strikes_tmp(indexx(1))=[];
        best_dips_tmp(indexx(1))=[];
        best_rakes_tmp(indexx(1))=[];
    end

end

[strike_avg,dip_avg,rake_avg] = MECH_AVG(best_strikes_tmp,best_dips_tmp,best_rakes_tmp);

% determining rms and perc_win_angthresh
N_orig = length(best_strikes_orig);
ang_btw_sdr_rms = zeros(N_orig,1);

for i = 1:N_orig
    sdr_avg = [strike_avg dip_avg rake_avg];
    sdr_i = [best_strikes_orig(i) best_dips_orig(i) best_rakes_orig(i)];
    
    ang_tmp = kagan(sdr_avg,sdr_i);
    ang_btw_sdr_rms(i,1) = abs(ang_tmp);

end

ang_rms = rms(ang_btw_sdr_rms);
perc_win_angthresh = (length(ang_btw_sdr_rms(ang_btw_sdr_rms<ang_thresh))/N_orig)*100;

end

function [strike_avg,dip_avg,rake_avg] = MECH_AVG(best_strikes,best_dips,best_rakes)

% c ------------------------------------------------------------ 
% c
% c subroutine MECH_AVG determines the average focal mechanism of a set
% c   of mechanisms
% c
% c  Inputs:  nf     =  number of fault planes
% c           norm1(3,nf) = normal to fault plane
% c           norm2(3,nf) = slip vector
% c  Output:  norm1_avg(3)    = normal to avg plane 1
% c           norm2_avg(3)    = normal to avg plane 2
% c
% c    Written  10/4/2000 by Jeanne Hardebeck                              
% c    Modified 5/14/2001 by Jeanne Hardebeck                              
% c    Modified 4/29/2020 by Oluwaseun fadugba (Convert to MATLAB and used kagan.m to determine ROTA)
%
%       subroutine MECH_AVG(nf,norm1,norm2,norm1_avg,norm2_avg)
%            
%       real dot1,fract1
%       real misf,maxmisf,avang1,avang2
%       real norm1(3,nf),norm2(3,nf),temp1(3),temp2(3)
%       real norm1_avg(3),norm2_avg(3),ln_norm1,ln_norm2
%       real theta1,theta2,ref1(3),ref2(3)
%       integer nf
% 

nf = length(best_strikes);

norm1(:,1) = -sin(best_dips*pi/180).*sin(best_strikes*pi/180);
norm1(:,2) =  sin(best_dips*pi/180).*cos(best_strikes*pi/180);
norm1(:,3) = -cos(best_dips*pi/180);

norm2(:,1) =  cos(best_rakes*pi/180).*cos(best_strikes*pi/180) + ...
    cos(best_dips*pi/180).*sin(best_rakes*pi/180).*sin(best_strikes*pi/180);
norm2(:,2) =  cos(best_rakes*pi/180).*sin(best_strikes*pi/180) - ...
    cos(best_dips*pi/180).*sin(best_rakes*pi/180).*cos(best_strikes*pi/180);
norm2(:,3) = -sin(best_rakes*pi/180).*sin(best_dips*pi/180);   

degrad=180./3.1415927;
%       
% c if there is only one mechanism, return that mechanism
%    
if (nf==1) 

    norm1_avg=norm1;
    norm2_avg=norm2;

else
    % c find the average normal vector for each plane - determine which
    % c nodal plane of each event corresponds to which by finding the
    % c minimum focal mechanism rotation
    %
    norm1_avg = sum(norm1);
    norm2_avg = sum(norm2);
    
    ln_norm1=sqrt(norm1_avg(1).^2 + norm1_avg(2).^2 + norm1_avg(3).^2);
    ln_norm2=sqrt(norm2_avg(1).^2 + norm2_avg(2).^2 + norm2_avg(3).^2);

    norm1_avg=norm1_avg/ln_norm1;
    norm2_avg=norm2_avg/ln_norm2;

    % c determine the RMS observed angular difference between the average 
    % c normal vectors and the normal vectors of each mechanism
    % 
    avang1=0.;
    avang2=0.;

    for i=1:nf

        temp1=norm1(i,:);
        temp2=norm2(i,:);

        %rota = MECH_ROT(norm1_avg,temp1,norm2_avg,temp2);
        d11=temp1(1)*norm1_avg(1)+temp1(2)*norm1_avg(2)+...
                            temp1(3)*norm1_avg(3);
        d22=temp2(1)*norm2_avg(1)+temp2(2)*norm2_avg(2)+...
                            temp2(3)*norm2_avg(3);
                        
        if (d11>=1.); d11=1.; end
        if (d11<=-1.); d11=-1.; end
        if (d22>=1.); d22=1.; end
        if (d22<=-1.); d22=-1.; end

        a11=acos(d11);
        a22=acos(d22);
        avang1=avang1+a11*a11;
        avang2=avang2+a22*a22;
    end

    avang1=sqrt(avang1/nf);
    avang2=sqrt(avang2/nf);

    % c the average normal vectors may not be exactly orthogonal (although
    % c usually they are very close) - find the misfit from orthogonal and 
    % c adjust the vectors to make them orthogonal - adjust the more poorly 
    % c constrained plane more
 
    if ((avang1+avang2)> 0.0001)

        maxmisf=0.01;
        fract1=avang1/(avang1+avang2);

        for icount=1:100  
            dot1=norm1_avg(1)*norm2_avg(1)+norm1_avg(2)...
                 *norm2_avg(2)+norm1_avg(3)*norm2_avg(3);

            misf=90.-acos(dot1)*degrad;

            if (abs(misf)<=maxmisf); break; end

            theta1=misf*fract1/degrad;
            theta2=misf*(1.-fract1)/degrad;

            for j=1:3
              temp=norm1_avg(j);
              norm1_avg(j)=norm1_avg(j)-norm2_avg(j)*sin(theta1);
              norm2_avg(j)=norm2_avg(j)-temp*sin(theta2);
            end

            ln_norm1=sqrt(norm1_avg(1).^2 + norm1_avg(2).^2 + norm1_avg(3).^2);
            ln_norm2=sqrt(norm2_avg(1).^2 + norm2_avg(2).^2 + norm2_avg(3).^2);

            norm1_avg=norm1_avg/ln_norm1;
            norm2_avg=norm2_avg/ln_norm2;

        end

    end

end

% Convert the normal and slip vectors to strike, dip and rake values
[strike_avg,dip_avg,rake_avg,~,~,~] = strike_dip_rake(norm1_avg,norm2_avg);

end

function [strike1,dip1,rake1,strike2,dip2,rake2] = strike_dip_rake(n,u)
% *************************************************************************%
%   program STRESSINVERSE                                                 %
%                                                                         %
%   joint iterative inversion for stress and faults from focal mechanisms %
%                                                                         %
%   Vavrycuk, V., 2014. Iterative joint inversion for stress and fault    %
%   orientations from focal mechanisms, Geophys. J. Int., 199, 69-77,     %
%   doi: 10.1093/gji/ggu224                                               %
%                                                                         %
%                                                                         %
%   function STRIKE_DIP_RAKE                                              %
%                                                                         %
%   calculation of strike, dip and rake from the fault normals and slip   %
%   directions                                                            %
%                                                                         %
%   input:  fault normal n                                                %
%           slip direction u                                              %
%                                                                         %
%   output: strike, dip and rake                                          %
%                                                                         %

%*************************************************************************%

n1 = n;
u1 = u;
	
if (n1(3)>0) n1 = -n1; u1 = -u1; end % vertical component is always negative!
%if (n1(3)>0) n1 = -n1; end; % vertical component is always negative!
        
n2 = u;
u2 = n;
    
if (n2(3)>0) n2 = -n2; u2 = -u2; end  % vertical component is always negative!

% ------------------------------------------------------------------------
% 1st solution
%--------------------------------------------------------------------------
dip    = acos(-n1(3))*180/pi;
strike = asin(-n1(1)/sqrt(n1(1)^2+n1(2)^2))*180/pi;

% determination of a quadrant
if (n1(2)<0) strike=180-strike; end

rake = asin(-u1(3)/sin(dip*pi/180))*180/pi;

% determination of a quadrant
cos_rake = u1(1)*cos(strike*pi/180)+u1(2)*sin(strike*pi/180);
if (cos_rake<0) rake=180-rake; end

if (strike<0   ) strike = strike+360; end 
if (rake  < 0) rake   = rake  +360; end %if (rake  < -180) rake   = rake  +360; end;
%if (rake  > 180) rake   = rake  -360; end;  % rake is in the interval -180<rake<180
    
strike1 = real(strike); dip1 = real(dip); rake1 = real(rake);

% ------------------------------------------------------------------------
% 2nd solution
%--------------------------------------------------------------------------
dip    = acos(-n2(3))*180/pi;
strike = asin(-n2(1)/sqrt(n2(1)^2+n2(2)^2))*180/pi;

% determination of a quadrant
if (n2(2)<0) strike=180-strike; end

rake = asin(-u2(3)/sin(dip*pi/180))*180/pi;

% determination of a quadrant
cos_rake = u2(1)*cos(strike*pi/180)+u2(2)*sin(strike*pi/180);
if (cos_rake<0) rake=180-rake; end

if (strike<0   ) strike = strike+360; end
if (rake  < 0) rake   = rake  +360; end %if (rake  <-180) rake   = rake  +360; end;
%if (rake  > 180) rake   = rake  -360; end;  
    
strike2 = real(strike); dip2 = real(dip); rake2 = real(rake);

end

function plot_P_nodal_planes(strike1,dip1,rake1,best_strikes,best_dips,best_rakes)

p_pol= [ 0	  190.00	   0	    5.00	    1.00	   -1.00];
vp=6.57; vs=3.71; dens=2.79; z_displ=1; r_displ=1; t_displ=1;simul_tag='syndata';

line_color=[0.5,0.5,0.5; 1, 0, 1; 0, 0, 1; 1, 0, 0; 0, 1, 0];

N = length(strike1)+1;

    for ii = 1:N 

        if ii == 1
        [c_all,icontour_all, numstart_all, numfin_all] = ...
            calc_P_nodal_planes(best_strikes,best_dips,best_rakes,vp,vs,dens, z_displ, r_displ, t_displ);
        else
        [c_all,icontour_all, numstart_all, numfin_all] = ...
            calc_P_nodal_planes(strike1(ii-1),dip1(ii-1),rake1(ii-1),vp,vs,dens, z_displ, r_displ, t_displ);
        end
        
        hold on; %p_pol=[1];
        foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,...
            p_pol,1,line_color(ii,:)); hold on;
    end
    
    con=pi/180.;
    radius=1.0;
    a=linspace(0,360,1000)*con;
    yc=radius*cos(a);
    xc=radius*sin(a);

    plot(xc,yc,'linewidth',2,'color','r');


    title('P Wave Nodal Surfaces', 'FontSize', 20);
    
    markersize = 8;
    textfontsize = 14;
    linewidth = 2;
    
    text(0.8-0.3-0.7,1.4,simul_tag, 'FontSize', textfontsize,'FontWeight','bold');
        
    %  Put polarity info onto plot
    plot(0.9,1.2,'MarkerSize',markersize,'marker','o','markeredgecolor','k','markerfacecolor','k');
    text(1.0,1.2,'Compr.', 'FontSize', textfontsize);
    plot(0.9,1.05,'MarkerSize',markersize,'marker','x','markeredgecolor','k','markerfacecolor','k');
    text(1.0,1.05,'Nodal', 'FontSize', textfontsize);
    plot(0.9,0.90,'MarkerSize',markersize,'marker','o','markeredgecolor','k','markerfacecolor','w');
    text(1.0,0.90,'Dilat.', 'FontSize', textfontsize);
    
    dx = -0.8; dy = -1.3;
    text(0.6-dx,-0.85-dy,'Legend', 'FontSize', textfontsize);
    
    plot([0.6 0.8]-dx,-[1.05 1.05]-dy,'linewidth',linewidth,'color',line_color(1,:));
    text(0.85-dx,-1.05-dy,'Acceptable FMs', 'FontSize', textfontsize);
    
    text(0.8-dx,-1.25-dy,'Preferred FMs', 'FontSize', textfontsize);
    
    plot([0.6 0.8]-dx,-[1.4 1.4]-dy,'linewidth',linewidth,'color',line_color(2,:));
    text(0.85-dx,-1.4-dy,'Nodal planes', 'FontSize', textfontsize); 
    
    axis off;
     
    %axis off; 
    set(gca,'xtick',[]); set(gca,'xticklabel',[])
    set(gca,'ytick',[]); set(gca,'yticklabel',[])
    
end