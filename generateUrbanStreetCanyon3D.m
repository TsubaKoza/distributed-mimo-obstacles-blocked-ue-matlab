function scenario = generateUrbanStreetCanyon3D(cfg)
%GENERATEURBANSTREETCANYON3D Buildings first, then outdoor APs and UEs.
% Corridors and open gaps are geometry metadata, never ray-tracing obstacles.
assert(cfg.enableObstacles,'Urban study requires buildings.');
validateattributes(cfg.urbanCorridorCount,{'numeric'},{'vector','numel',2,'integer','positive'});
validateattributes(cfg.urbanCorridorWidth,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.urbanBuildingGap,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.minimumAPSeparation,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(cfg.urbanYawJitterDeg,{'numeric'},{'scalar','nonnegative','<=',15});
validateattributes(cfg.urbanNLoSProneUEFraction,{'numeric'},{'scalar','>=',0,'<=',1});
ranges={cfg.urbanBuildingWidthRange,cfg.urbanBuildingDepthRange,cfg.urbanBuildingHeightRange};
for i=1:3
    validateattributes(ranges{i},{'numeric'},{'vector','numel',2,'positive','finite'});
    assert(ranges{i}(2)>=ranges{i}(1));
end
centers=cell(1,2); blocks=cell(1,2); corridors=zeros(0,4);
for d=1:2
    edges=linspace(cfg.region(d,1),cfg.region(d,2),cfg.urbanCorridorCount(d)+2);
    centers{d}=edges(2:end-1);
    low=centers{d}-cfg.urbanCorridorWidth/2;
    high=centers{d}+cfg.urbanCorridorWidth/2;
    blocks{d}=[[cfg.region(d,1),high].',[low,cfg.region(d,2)].'];
    assert(all(diff(blocks{d},1,2)>cfg.urbanBuildingGap),'Corridors leave no viable blocks.');
    for j=1:numel(low)
        if d==1
            corridors(end+1,:)=[low(j),high(j),cfg.region(2,:)]; %#ok<AGROW>
        else
            corridors(end+1,:)=[cfg.region(1,:),low(j),high(j)]; %#ok<AGROW>
        end
    end
end
template=struct('center',[],'width',[],'depth',[],'height',[], ...
    'yaw',[],'reflectionCoefficient',[],'active',true);
obstacles=repmat(template,0,1); envelopes=zeros(0,4);
for ix=1:size(blocks{1},1)
    for iy=1:size(blocks{2},1)
        bx=blocks{1}(ix,:); by=blocks{2}(iy,:);
        nx=max(1,ceil(diff(bx)/(cfg.urbanBuildingWidthRange(2)+cfg.urbanBuildingGap)));
        ny=max(1,ceil(diff(by)/(cfg.urbanBuildingDepthRange(2)+cfg.urbanBuildingGap)));
        xe=linspace(bx(1),bx(2),nx+1); ye=linspace(by(1),by(2),ny+1);
        for i=1:nx
            for j=1:ny
                slot=[xe(i)+cfg.urbanBuildingGap/2,xe(i+1)-cfg.urbanBuildingGap/2, ...
                    ye(j)+cfg.urbanBuildingGap/2,ye(j+1)-cfg.urbanBuildingGap/2];
                accepted=false;
                for trial=1:cfg.maxPlacementAttempts
                    o=template;
                    o.width=sampleRange(cfg.urbanBuildingWidthRange);
                    o.depth=sampleRange(cfg.urbanBuildingDepthRange);
                    o.height=sampleRange(cfg.urbanBuildingHeightRange);
                    o.yaw=deg2rad((2*rand-1)*cfg.urbanYawJitterDeg);
                    extent=[abs(cos(o.yaw))*o.width+abs(sin(o.yaw))*o.depth, ...
                        abs(sin(o.yaw))*o.width+abs(cos(o.yaw))*o.depth];
                    if extent(1)>diff(slot(1:2)) || extent(2)>diff(slot(3:4)), continue; end
                    % Small placement variation inside its own disjoint slot.
                    xy=[slot(1),slot(3)]+extent/2+ ...
                        ([diff(slot(1:2)),diff(slot(3:4))]-extent).*rand(1,2);
                    o.center=[xy,o.height/2];
                    o.reflectionCoefficient=sampleReflectionCoefficient(cfg);
                    obstacles(end+1,1)=o; %#ok<AGROW>
                    envelopes(end+1,:)=[xy(1)-extent(1)/2,xy(1)+extent(1)/2, ...
                        xy(2)-extent(2)/2,xy(2)+extent(2)/2]; %#ok<AGROW>
                    accepted=true; break;
                end
                assert(accepted,'Building ranges do not fit the street-block slots. Adjust ranges/gap/corridors.');
            end
        end
    end
end
APpos=zeros(cfg.M,3);
for m=1:cfg.M
    accepted=false;
    for trial=1:cfg.maxPlacementAttempts
        xy=sampleCorridor(centers,cfg,true);
        q=[xy,cfg.zAP];
        if isOutdoor(q,obstacles,cfg.obstacleClearance) && ...
                all(vecnorm(APpos(1:m-1,:)-q,2,2)>=cfg.minimumAPSeparation)
            APpos(m,:)=q; accepted=true; break;
        end
    end
    assert(accepted,'Cannot satisfy AP spacing on the configured corridors.');
end
UEpos=zeros(cfg.K,3); placementClass=strings(cfg.K,1);
for k=1:cfg.K
    prone=rand<cfg.urbanNLoSProneUEFraction;
    accepted=false;
    for trial=1:cfg.maxPlacementAttempts
        if prone
            % All facades are eligible, independent of AP/LoS information.
            o=obstacles(randi(numel(obstacles))); side=randi(4);
            margin=cfg.obstacleClearance+1+3*rand;
            if side<=2
                xy=[(2*side-3)*(o.width/2+margin),(rand-0.5)*(o.depth+4)];
            else
                xy=[(rand-0.5)*(o.width+4),(2*(side-2)-3)*(o.depth/2+margin)];
            end
            c=cos(o.yaw); s=sin(o.yaw);
            xy=xy*[c s;-s c]+o.center(1:2);
        else
            xy=sampleCorridor(centers,cfg,false);
        end
        q=[xy,cfg.zUE];
        if all(xy>=cfg.region(1:2,1).' & xy<=cfg.region(1:2,2).') && ...
                isOutdoor(q,obstacles,cfg.obstacleClearance) && ...
                all(vecnorm(APpos-q,2,2)>=cfg.minAPUEDistance) && ...
                all(vecnorm(UEpos(1:k-1,:)-q,2,2)>=cfg.minAPUEDistance)
            UEpos(k,:)=q; accepted=true; break;
        end
    end
    assert(accepted,'Could not place an outdoor UE within the attempt limit.');
    if prone, placementClass(k)="building perimeter"; else, placementClass(k)="street"; end
end
scenario=struct('APpos',APpos,'UEpos',UEpos,'obstacles',obstacles, ...
    'antennaPos',generateAPArrayPositions(APpos,cfg),'mode',cfg.scenarioMode, ...
    'region',cfg.region,'corridors',corridors,'buildingEnvelopes',envelopes, ...
    'UEPlacementClass',placementClass,'seed',cfg.seed);
validateUrbanScenario3D(scenario,cfg);
end

function xy=sampleCorridor(centers,cfg,nearSide)
d=randi(2); other=3-d;
xy=zeros(1,2); xy(other)=sampleRange(cfg.region(other,:));
if nearSide
    offset=(2*(rand>0.5)-1)*(cfg.urbanCorridorWidth/2-1);
else
    offset=(rand-0.5)*(cfg.urbanCorridorWidth-2);
end
xy(d)=centers{d}(randi(numel(centers{d})))+offset;
end
function value=sampleRange(bounds)
value=bounds(1)+diff(bounds)*rand;
end
function yes=isOutdoor(q,obstacles,margin)
yes=true;
for i=1:numel(obstacles)
    o=obstacles(i); local=transformToObstacleLocal(q,o);
    if all(abs(local(1:2))<=[o.width,o.depth]/2+margin)
        yes=false; return;
    end
end
end
