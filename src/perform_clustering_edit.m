function perform_clustering_edit()

global nkh cdf_rnd vr
global nks vs cdf_cat
global xs_decl ys_decl zs_decl 
global xs_sort ys_sort zs_sort vol_sort
global h1_decl h2_decl h3_decl h4_decl h_perc_edit
global V05 NV05_cat

global xs_clu ys_clu zs_clu
global xs_diffuse ys_diffuse zs_diffuse

PROB = str2double(get(h_perc_edit,'String'))/100;

%  Calculate the volume of tetrahedra at the 5% probability level (0.05) 
%  from the random catalog
for kh=1:nkh-1
    if cdf_rnd(kh) < PROB && cdf_rnd(kh+1) >= PROB
        ncalc=kh;
    end
end
ncalc;
V05=vr(ncalc);
fprintf('ncalc = %i, Volume at 5 percent probability (V05, km cubed) = %g\n',ncalc,V05);

%V05=0.5;
%**************************************************************************
%  Calculate N(V(0.05)) from original catalog

for ks=1:nks-1
    if vs(ks) < V05 && vs(ks+1) >= V05
        ncat=ks;
    end
end
ncat;
NV05_cat=cdf_cat(ncat);
fprintf('ncat = %i, Probability at target volume (NV05) = %g\n',ncalc,NV05_cat);

semilogx(h3_decl,vr,cdf_rnd,'-k','LineWidth',1.5); hold on;
set(h3_decl,'NextPlot','add');
semilogx(h3_decl,vs,cdf_cat,'-b','LineWidth',1.5); hold off;

set(h3_decl,'NextPlot','add');
semilogx(h3_decl,[V05 V05],[0 1],'-r','LineWidth',1.5);
set(h3_decl,'NextPlot','add');
semilogx(h3_decl,[min([vs; vr]) max([vs; vr])],[PROB PROB],'--r','LineWidth',1.5);
grid(h3_decl,'MINOR')
legend(h3_decl,'Random','Observed','V05',num2str(round(PROB,2)),'Location','best');
xlim(h3_decl,[min([vs; vr]) max([vs; vr])]);
set(h3_decl,'NextPlot','replace');
% 
xs_clu = xs_sort(vol_sort <= V05);
ys_clu = ys_sort(vol_sort <= V05);
zs_clu = zs_sort(vol_sort <= V05);

xs_diffuse = xs_sort(vol_sort > V05);
ys_diffuse = ys_sort(vol_sort > V05);
zs_diffuse = zs_sort(vol_sort > V05);

scatter3(h2_decl,xs_clu,ys_clu,zs_clu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 0 1])
axis(h2_decl,'equal');
grid(h2_decl,'on');
xlabel(h2_decl,'X km');
ylabel(h2_decl,'Y km');
zlabel(h2_decl,'Z km');

xx = round((max(xs_decl)-min(xs_decl))/5);
yy = round((max(ys_decl)-min(ys_decl))/5);
zz = round((max(zs_decl)-min(zs_decl))/5);
    
xlim(h2_decl,[min(xs_decl)-xx max(xs_decl)+xx])
ylim(h2_decl,[min(ys_decl)-yy max(ys_decl)+yy])
zlim(h2_decl,[min(zs_decl)-zz max(zs_decl)+zz])

viewpt=h1_decl.View;
view(h2_decl,viewpt(1),viewpt(2));
hold on;
set(h2_decl,'NextPlot','replace');

% ------------------------------
scatter3(h4_decl,xs_diffuse,ys_diffuse,zs_diffuse,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 0.75 0])   
set(h4_decl,'NextPlot','add');
scatter3(h4_decl,xs_clu,ys_clu,zs_clu,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 0 1])

axis(h4_decl,'equal');
grid(h4_decl,'on');

xlabel(h4_decl,'X km');
ylabel(h4_decl,'Y km');
zlabel(h4_decl,'Z km');

xx = round((max(xs_decl)-min(xs_decl))/5);
yy = round((max(ys_decl)-min(ys_decl))/5);
zz = round((max(zs_decl)-min(zs_decl))/5);
    
xlim(h4_decl,[min(xs_decl)-xx max(xs_decl)+xx])
ylim(h4_decl,[min(ys_decl)-yy max(ys_decl)+yy])
zlim(h4_decl,[min(zs_decl)-zz max(zs_decl)+zz])

viewpt=h1_decl.View;
view(h4_decl,viewpt(1),viewpt(2));
hold on;
set(h4_decl,'NextPlot','replace');

% -----------------------
Link = linkprop([h1_decl, h2_decl h4_decl], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);
