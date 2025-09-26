function Offspring = OperatorDE(Problem, Population, fitness)

    if isa(Population(1), 'SOLUTION')
        evaluated = true;
        Popu = Population.decs;
    else
        evaluated = false;
    end
    [~, index] = sort(fitness, 'ascend');
    Offspring = DE(Popu, Problem, index);
    if evaluated
        Offspring = Problem.Evaluation(Offspring);
    end
end

function Offspring = DE(Popu, Problem, index)
    [N, D] = size(Popu);
    Fm = [0.6,0.8,1.0];
    CRm=[0.1,0.2,1.0];
    Index = randi([1,length(Fm)],Problem.N,1);
    F = Fm(Index);
    F = F';
    F = F(:, ones(1,Problem.D));
    Index = randi([1,length(CRm)],Problem.N,1);
    Cr = CRm(Index);
    Cr = Cr';
    P = max(round(0.17 * N * (1 - 0.9 * Problem.FE/Problem.maxFE)), 2);
    pbest = Popu(index(randi(P, N, 1)), :);
    [r1, r2] = gnR1R2(N, 1:N);
    V = Popu + F .* ...
        (pbest - Popu + Popu(r1,:) - Popu(r2,:));
    jrand = randi(D, N, 1);
    CrossMask = rand(N,D) < Cr;
    for i = 1:N
        CrossMask(i, jrand(i)) = true;
    end
    Offspring = Popu;
    Offspring(CrossMask) = V(CrossMask);
    PerturbMask = rand(N,D) < 0.2;
    Offspring(PerturbMask) = cauchyrnd(Offspring(PerturbMask), 0.2);
    Lower = repmat(Problem.lower, N, 1);
    Upper = repmat(Problem.upper, N, 1);
    Offspring = min(max(Offspring, Lower), Upper);
end

function [r1, r2] = gnR1R2(N, r0)
    NP = length(r0);
    maxTrials = 1000;

    %% Generate r1 ≠ r0
    for trial = 1:maxTrials
        r1 = randi(N, 1, NP);
        if all(r1 ~= r0)
            break;
        end
        % fix conflicts only
        conflict = r1 == r0;
        r1(conflict) = randi(N, 1, sum(conflict));
        if trial == maxTrials
            error('Failed to generate r1 after %d trials', maxTrials);
        end
    end

    %% Generate r2 ≠ r0 and ≠ r1
    for trial = 1:maxTrials
        r2 = randi(N, 1, NP);
        if all((r2 ~= r0) & (r2 ~= r1))
            break;
        end
        % fix conflicts only
        conflict = (r2 == r0) | (r2 == r1);
        r2(conflict) = randi(N, 1, sum(conflict));
        if trial == maxTrials
            error('Failed to generate r2 after %d trials', maxTrials);
        end
    end
end
