function FM_seeded_planes(n0,FAULT_FLAG)
%  Construct n0 faults using nearby focal mechanisms when we are plitting
%  the thickest fault in OADC_3D
%  just need the vertice locations

%  FAULT_FLAG = kthick (from splitfault.m),  Split the thickest fault into
%               two random planes of length L/2

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt yt zt Nt xb yb zb lambda3
global Strike Dip FM_file plot_avg_FM cat_Tmatrix

% addpath('/Users/oluwaseunfadugba/Documents/Waveform_Modeling_Now/focal_2/focal_for_CSZ_project/src')

% Check if to use FM constraints and if FM_file exists.
% If not run the random fault generator function

if ~isempty(FM_file)
    
    if isfile(FM_file) 
        
        % File exists. Load FM data set
        fid = fopen(FM_file) ; %'All_write_mechanisms.txt'
        data = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',2) ;
        data = cell2mat(data); fclose(fid);

        data = data(data(:,13) >= 4,:);  % no of sta threshold
        data = data(data(:,5) >= 2,:);   % magnitude threshold

        [N_fm,~]=size(data);

        xsf = (data(:,3) + 70.6)*75.778;
        ysf = (data(:,2) - 47.2)*111.1743;
        zsf = -data(:,4);

        strikes=data(1:N_fm,6);
        dips = data(1:N_fm,7);
        rakes = data(1:N_fm,8);

        if FAULT_FLAG == 0
            NtF = length(xs);
            xthick = xs; 
            ythick = ys;
            zthick = zs;

        else
            NtF = Nt(FAULT_FLAG);
            xthick = xt(FAULT_FLAG,1:Nt(FAULT_FLAG)); 
            ythick = yt(FAULT_FLAG,1:Nt(FAULT_FLAG));
            zthick = zt(FAULT_FLAG,1:Nt(FAULT_FLAG));

        end

        % Find the closet FM to the earthquake in the thick cluster
        kk= 0; dst = zeros(1,N_fm); 

        for i = 1:NtF % per hypocenter

            for m=1:N_fm  %per FM
                dst(m)= sqrt((xsf(m)-xthick(i))^2 + (ysf(m)-ythick(i))^2 + (zsf(m)-zthick(i))^2);           
            end

            %  find the closest fault plane
            [mindist, index] = min(dst);

            dist2FM_threshold=2;

            if mindist < dist2FM_threshold 
                kk = kk + 1;

                xs_fm(kk) = xthick(i);
                ys_fm(kk) = ythick(i);
                zs_fm(kk) = zthick(i);

                strikes_fm(kk) = strikes(index);
                dips_fm(kk) = dips(index);
                rakes_fm(kk) = rakes(index);     

            end 
        end

        % At this point, we know if there are FMs close to the thick cluster.
        % If no FM is present, OADC_3D will use random-seeded faults.
      
        if kk==0 % No FM near the thick fault
            fprintf('No focal mechanisms near the ''thick'' fault\n');
            fprintf('Using randomly-seeded planes\n');
            randfaults(n0,FAULT_FLAG)   

        else
            fprintf('Averaging %i focal mechanisms\n',kk);
            fprintf('Using focal-mechanism-seeded planes\n');

            preferred_FM = calc_preferred_FM(strikes_fm',dips_fm',rakes_fm',plot_avg_FM);

            if preferred_FM(1) <= preferred_FM(4)
                strike_pr = preferred_FM(1); 
                dip_pr = preferred_FM(2);
                %rake_pr = preferred_FM(3);
            else
                strike_pr = preferred_FM(4);
                dip_pr = preferred_FM(5);
                %rake_pr = preferred_FM(6);
            end

            generate_plane_using_FM(n0,FAULT_FLAG,strike_pr,dip_pr,kk,xs_fm,ys_fm,zs_fm);
        end
    else
        
        error('Error: Focal mechanism file does not exit.')
    end
    
else
     % File does not exist.
     fprintf('Focal mechanisms doesn''t exist\n');
     fprintf('Using randomly-seeded planes\n');
     randfaults(n0,FAULT_FLAG)
end

end
