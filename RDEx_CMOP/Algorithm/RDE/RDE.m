classdef RDE < ALGORITHM

    methods
        function main(Algorithm,Problem)

            %% Generate random population
            Population = Problem.Initialization();
            cons = Population.cons;
            cons(cons <= 0)=0;
            conss = sum(cons,2);
            epsilon0 = max(conss);
            if epsilon0 == 0
                epsilon0 = 1;
            end
            Fitness = CalFitness_E(Population.objs,Population.cons,epsilon0);
            %% Optimization
            while Algorithm.NotTerminated(Population)

                cp=(-log(epsilon0)-6)/log(1-0.5);
                epsilon = epsilon0*(1-Problem.FE/Problem.maxFE)^cp;
                
                Offspring = OperatorDE(Population, Problem.N, Problem, Fitness);

                [Population,Fitness] = Improve_E_EnvironmentalSelection([Population,Offspring],Problem.N,epsilon);
                

            end
        end
    end
end

