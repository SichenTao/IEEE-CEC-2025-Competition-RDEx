classdef RDE < ALGORITHM
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            eplion = Algorithm.ParameterSet(0.05);
            %% Generate random population
            Popul      = Problem.Initialization();
            [Population,N] = Selection(Popul,Problem.N);
            Population2=Population;

            %% Optimization
            while Algorithm.NotTerminated(Population)

                if Problem.FE>= Problem.maxFE*0.5 && rand()<0.5
                    NewPopulation = Exploration(Problem,Population,Popul,N,Problem.N);
                    Population2 = EnvironmentalSelection([Population2,NewPopulation],Problem.N,eplion);
                    MatingPool = TournamentSelection(2,Problem.N,-CalFitness(Population2.objs,eplion));
                    NewPopul     = OperatorDE(Problem,Population2(MatingPool),-CalFitness(Popul(MatingPool).objs,eplion));
                    Population2        = EnvironmentalSelection([Population2,NewPopul],Problem.N,eplion);
                else
                    NewPopulation = Exploration(Problem,Population,Popul,N,Problem.N);
                    Popul = EnvironmentalSelection([Popul,NewPopulation],Problem.N,eplion);
                    MatingPool = TournamentSelection(2,Problem.N,-CalFitness(Popul.objs,eplion));
                    NewPopul     = OperatorDE(Problem,Popul(MatingPool),-CalFitness(Popul(MatingPool).objs,eplion));
                    Popul        = EnvironmentalSelection([Popul,NewPopul],Problem.N,eplion);
                end
                [Population,N] = Selection([Population,NewPopul,NewPopulation,Population2],Problem.N);
                Population2=Population;
            end
        end
    end
end