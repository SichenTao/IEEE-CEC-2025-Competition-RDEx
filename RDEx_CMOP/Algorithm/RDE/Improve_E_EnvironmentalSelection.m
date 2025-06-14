function [return_pop,return_Fitness] = Improve_E_EnvironmentalSelection(Population,N,VAR)

input_cons = Population.cons;
input_cons(input_cons<0) = 0;
input_cons = sum(input_cons,2);
findex = find(input_cons<=VAR);
ifindex = find(input_cons>VAR);
fPopulation = Population(findex);
ifPopulation = Population(ifindex);
if isempty(fPopulation)
    ifFitness = CalFitness(ifPopulation.objs,ifPopulation.cons);
    Next2 = ifFitness < 1;
    if sum(Next2) <= N
        [~,Rank] = sort(ifFitness);
        Next2(Rank(1:N )) = true;
    elseif sum(Next2) > N
        Del  = Truncation(ifPopulation(Next2).objs, sum(Next2)-N );
        Temp = find(Next2);
        Next2(Temp(Del)) = false;
    end
    ifPopulation = ifPopulation(Next2);
    ifFitness    = ifFitness(Next2);
    [ifFitness,rank] = sort(ifFitness);
    ifPopulation = ifPopulation(rank);
    fPopulation = [];
    fFitness = [];
elseif length(fPopulation) <= N
    cons = fPopulation.cons;
    cons(cons<0)=0;
    cons = sum(cons,2);
    fFitness = CalFitness([fPopulation.objs,cons]);
    Next = fFitness < 1;
    [~,Rank] = sort(fFitness);
    Next(Rank(1:length(fPopulation) )) = true;
    fPopulation = fPopulation(Next);
    fFitness    = fFitness(Next);
    [fFitness,rank] = sort(fFitness);
    fPopulation = fPopulation(rank);
    ifFitness = CalFitness(ifPopulation.objs,ifPopulation.cons);
    Next2 = ifFitness < 1;
    if sum(Next2) <= N - length(fPopulation)
        [~,Rank] = sort(ifFitness);
        Next2(Rank(1:N - length(fPopulation) )) = true;
    elseif sum(Next2) > N - length(fPopulation)
        Del  = Truncation(ifPopulation(Next2).objs, sum(Next2)-(N - length(fPopulation)) );
        Temp = find(Next2);
        Next2(Temp(Del)) = false;
    end
    ifPopulation = ifPopulation(Next2);
    ifFitness    = ifFitness(Next2) + max(fFitness);
    [ifFitness,rank] = sort(ifFitness);
    ifPopulation = ifPopulation(rank);
elseif length(fPopulation) > N
    cons = fPopulation.cons;
    cons(cons<0)=0;
    cons = sum(cons,2);
    fFitness = CalFitness([fPopulation.objs,cons]);
    Next = fFitness < 1;
    if sum(Next) <= N
        [~,Rank] = sort(fFitness);
        Next(Rank(1:N )) = true;
    elseif sum(Next) > N
        Del  = Truncation(fPopulation(Next).objs, sum(Next)-N );
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    fPopulation = fPopulation(Next);
    fFitness   = fFitness(Next);
    [fFitness,rank] = sort(fFitness);
    fPopulation = fPopulation(rank);
    ifPopulation = [];
    ifFitness = [];
end

return_pop = [fPopulation,ifPopulation];
return_Fitness = [fFitness,ifFitness];
end

function Del = Truncation(PopObj,K)
Distance = pdist2(PopObj,PopObj);
Distance(logical(eye(length(Distance)))) = inf;
Del = false(1,size(PopObj,1));
while sum(Del) < K
    Remain   = find(~Del);
    if isempty(Remain)
        keyboard
    end
    Temp     = sort(Distance(Remain,Remain),2);
    [~,Rank] = sortrows(Temp);
    Del(Remain(Rank(1))) = true;
end
end