function Offspring = Exploration(Problem,Population,Polpu,nND,N)

    PCObj  = Population.objs;
    NPCObj = Polpu.objs;
    
    %% Normalization
    fmax   = max(PCObj,[],1);
    fmin   = min(PCObj,[],1);
    PCObj  = (PCObj-repmat(fmin,size(PCObj,1),1))./repmat(fmax-fmin,size(PCObj,1),1);
    NPCObj = (NPCObj-repmat(fmin,size(NPCObj,1),1))./repmat(fmax-fmin,size(NPCObj,1),1);

    %% Determine the size of the niche
    d  = pdist2(PCObj,PCObj);
    d(logical(eye(length(d)))) = inf;
    d  = sort(d,2);
    r0 = mean(d(:,min(3,size(d,2))));
    r  = nND/N*r0;
    
    %% Detect the solutions in PC to be explored
    d = pdist2(PCObj,NPCObj);
    S = find(sum(d<=r,2)<=1);
    
    %% Generate new solutions
    if ~isempty(S)
        MatingPool = randi(length(Population),1,length(S));
        Offspring  = OperatorDEhalf(Problem,Population([S',MatingPool]));
    else
        Offspring = [];
    end
end