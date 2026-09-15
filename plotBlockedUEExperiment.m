function figs = plotBlockedUEExperiment(e)
%PLOTBLOCKEDUEEXPERIMENT Five study figures; roads are display-only patches.
c=e.config; s=e.scenario; k=e.targetBlockedUE.UEIndex;
visibility='off'; if c.showFigures, visibility='on'; end
figs=gobjects(5,1);
figs(1)=figure('Name','Urban scenario','Visible',visibility);
drawCity(s,c); hold on;
regular=~e.classification.allAPBlockedMask;
blocked=e.classification.serviceableMask;
other=e.classification.allAPBlockedMask & ~blocked;
scatterGroup(s.APpos,70,'^',[0.1 0.35 0.85],'AP');
scatterGroup(s.UEpos(regular,:),18,'o',[0.25 0.65 0.35],'Ordinary UE');
scatterGroup(s.UEpos(blocked,:),38,'o',[0.9 0.3 0.1],'Serviceable blocked UE');
scatterGroup(s.UEpos(other,:),25,'x',[0.5 0.3 0.5],'Other all-AP blocked UE');
scatterGroup(s.UEpos(k,:),160,'p',[1 0.75 0],'Target blocked UE');
title(sprintf('Urban: %d APs, %d UEs, %d buildings; target UE %d',c.M,c.K,numel(s.obstacles),k));
legend('Location','eastoutside');
figs(2)=figure('Name','Target NLoS paths','Visible',visibility);
drawCity(s,c); hold on; scatterGroup(s.UEpos(k,:),160,'p',[1 0.75 0],'Target UE');
links=e.channelRankingTable; links=links(links.FeasibleNLoSLink,:);
zoomPoints=s.UEpos(k,:);
seen=strings(0);
for i=1:height(links)
    m=links.APIndex(i); ap=s.APpos(m,:); ue=s.UEpos(k,:);
    zoomPoints(end+1,:)=ap; %#ok<AGROW>
    scatter3(ap(1),ap(2),ap(3),70,'^','filled','MarkerFaceColor',[0.1 0.35 0.85],'HandleVisibility','off');
    text(ap(1),ap(2),ap(3)+2,sprintf('AP%d (#%d)',m,links.Rank(i)));
    paths=e.rayTracingResult.pathInfo(m,k).paths;
    for l=1:numel(paths)
        p=paths(l);
        if p.type=="reflection"
            via=p.reflectionPoint; color=[0.75 0.05 0.65];
        elseif p.type=="diffraction"
            via=p.diffractionPoint; color=[0 0.65 0.8];
        else
            continue;
        end
        xyz=[ap;via;ue]; show='off';
        zoomPoints(end+1,:)=via; %#ok<AGROW>
        if ~any(seen==p.type), show='on'; seen(end+1)=p.type; end %#ok<AGROW>
        plot3(xyz(:,1),xyz(:,2),xyz(:,3),'-','Color',color,'LineWidth',1, ...
            'DisplayName',char(p.type),'HandleVisibility',show);
    end
end
title(sprintf('Retained reflection / diffraction paths to UE %d (no direct LoS)',k));
xlim([max(c.region(1,1),min(zoomPoints(:,1))-20),min(c.region(1,2),max(zoomPoints(:,1))+20)]);
ylim([max(c.region(2,1),min(zoomPoints(:,2))-20),min(c.region(2,2),max(zoomPoints(:,2))+20)]);
legend('Location','eastoutside');
links=links(links.EffectiveGainLinear>0,:);
figs(3)=figure('Name','Channel ranking','Visible',visibility);
floorDB=floor(min(links.ChannelGain_dB)/10)*10-5;
bar(links.Rank,[links.ChannelGain_dB,links.EffectiveGain_dB],'BaseValue',floorDB); grid on;
ylim([floorDB,ceil(max(links.ChannelGain_dB)/10)*10+5]);
xticks(links.Rank); xticklabels(compose('AP%d',links.APIndex));
xlabel('AP ordered by effective gain'); ylabel('Gain [dB]');
legend('Channel gain','MRT effective gain','Location','best');
title('NLoS AP channel ranking (unit-norm MRT)');
figs(4)=figure('Name','Gain and required power','Visible',visibility);
t=e.singleAPMinimumPowerTable;
t=t(isfinite(t.EffectiveGain_dB) & isfinite(t.MinimumPowerFor10Mbps_W),:);
scatter(t.EffectiveGain_dB,t.MinimumPowerFor10Mbps_W,50,'filled'); grid on;
set(gca,'YScale','log');
text(t.EffectiveGain_dB,t.MinimumPowerFor10Mbps_W*1.2,compose(' AP%d',t.APIndex));
xlabel('Effective Channel Gain [dB]'); ylabel('Minimum Required Power for 10 Mbps [W, log scale]');
yline(1,'--','1 W');
title(sprintf('N=%d; Pearson=%.4f; Spearman=%.4f', ...
    e.correlationResults.N,e.correlationResults.Pearson,e.correlationResults.Spearman));
figs(5)=figure('Name','Feasible power combinations','Visible',visibility);
t=e.feasiblePowerCombinationTable;
scatter(t.TotalPower_W,t.Rate_Mbps,6,t.NumActiveAPs,'filled'); grid on; hold on;
if e.minimumPowerSolution.Found
    t0=e.minimumPowerSolution.Combination;
    scatter(t0.TotalPower_W,t0.Rate_Mbps,100,'rp','filled');
end
yline(c.targetRate/1e6,'--','QoS'); xlabel('Total Power [W]'); ylabel('Achieved Rate [Mbps]');
cb=colorbar; cb.Label.String='Active AP count';
cb.Ticks=1:numel(e.minimumPowerSolution.CandidateAPIndices);
title(sprintf('QoS feasible combinations: %d',height(t)));
if c.saveFigures
    if ~exist(c.resultsDir,'dir'), mkdir(c.resultsDir); end
    names={'urban_scenario','target_nlos_paths','channel_ranking','gain_vs_required_power','power_vs_rate'};
    for i=1:numel(figs)
        exportgraphics(figs(i),fullfile(c.resultsDir,[names{i},'.png']),'Resolution',160);
        savefig(figs(i),fullfile(c.resultsDir,[names{i},'.fig']));
    end
end
if ~c.showFigures, close(figs); end
end
function drawCity(s,c)
hold on;
if isfield(s,'corridors')
    for i=1:size(s.corridors,1)
        b=s.corridors(i,:);
        patch([b(1) b(2) b(2) b(1)],[b(3) b(3) b(4) b(4)],zeros(1,4), ...
            [0.8 0.8 0.8],'EdgeColor','none','HandleVisibility','off');
    end
end
for i=1:numel(s.obstacles)
    o=s.obstacles(i);
    v=[-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1].* ...
        [o.width o.depth o.height]/2;
    R=[cos(o.yaw),-sin(o.yaw),0;sin(o.yaw),cos(o.yaw),0;0 0 1];
    v=v*R.'+o.center;
    patch('Vertices',v,'Faces',[1 2 3 4;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8], ...
        'FaceColor',[0.55 0.6 0.7],'FaceAlpha',0.35,'EdgeColor',[0.4 0.4 0.4], ...
        'HandleVisibility','off');
end
axis equal; grid on; view(35,45); xlim(c.region(1,:)); ylim(c.region(2,:)); zlim(c.region(3,:));
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
set(gcf,'Position',[50 50 1100 750]);
end
function scatterGroup(p,sizePoint,marker,color,name)
if isempty(p), return; end
scatter3(p(:,1),p(:,2),p(:,3),sizePoint,marker,'MarkerEdgeColor',color, ...
    'MarkerFaceColor',color,'DisplayName',name);
end
