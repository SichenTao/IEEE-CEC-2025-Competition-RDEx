function Offspring = OperatorDEhalf(Problem,Population)
    if isa(Population(1),'SOLUTION')
        evaluated = true;
        Population    = Population.decs;
    else
        evaluated = false;
    end
    Popu   = Population(1:floor(end/2),:);
    PopuG   = Population(floor(end/2)+1:floor(end/2)*2,:);
    Offspring = DE(Popu,PopuG,Problem);
    if evaluated
        Offspring = Problem.Evaluation(Offspring);
    end
end

function Offspring = DE(Popu,PopuG,Problem)
    [N, D] = size(Popu);
    F = zeros(N,1);
    Cr = zeros(N,1);
    r1 = randperm(N,N);
    for i = 1:N
        F(i)  = max(0, min(cauchyrnd(0.7, 0.2), 1));
        Cr(i) = max(0, min(normrnd(0.5, 0.1), 1));
    end
    V = Popu + repmat(F,1,D) .* (Popu(r1,:) - PopuG);
    jrand = randi(D, N, 1);
    CrossMask = rand(N, D) < repmat(Cr,1,D);
    for i = 1:N
        CrossMask(i, jrand(i)) = true;
    end
    Offspring = Popu;  
    Offspring(CrossMask) = V(CrossMask);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Offspring = min(max(Offspring, Lower), Upper);
end