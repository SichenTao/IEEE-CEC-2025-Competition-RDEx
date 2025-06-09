function [ Offspring ] = OperatorDE(Population, popsize, Problem, Fitness)

[~, sorted_index] = sort(Fitness, 'ascend');
r0 = 1 : popsize;
[r1, r2] = gnR1R2(popsize, popsize, r0);
fit_r1 = Fitness(r1);
fit_r2 = Fitness(r2);
pNP = max(floor(popsize * (1-0.99*Problem.FE/Problem.maxFE)), 2);
randindex = ceil(rand(1, popsize) .* pNP);
randindex = max(1, randindex);
pbest = Population(sorted_index(randindex));

NewDec = OperatorDEx(Problem,Population(r0).decs,pbest.decs,Population(r1).decs,Population(r2).decs,fit_r1,fit_r2);
Offspring = Problem.Evaluation(NewDec);
end

function [r1, r2] = gnR1R2(NP1, NP2, r0)
NP0 = length(r0);
r1 = zeros(1, NP0);
r2 = zeros(1, NP0);

for i = 1:NP0
    pool = setdiff(1:NP1, r0(i));
    r1(i) = pool(randi(length(pool)));
end

for i = 1:NP0
    pool = setdiff(1:NP2, [r0(i), r1(i)]);
    r2(i) = pool(randi(length(pool)));
end
end