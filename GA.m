clear
clf
close all
%% Read Excel File
% There are three different test cases we are running
% 'European_Cities_Big.xlsx'     = 41 cities
% 'European_Cities_Medium.xlsx'  = 20 cities
% 'European_Cities_Short.xlsx'   = 10 cities

filename = 'European_Cities_Medium.xlsx';
cities = table2struct(readtable(filename));
price_matx = xlsread(filename, 'Price');
%% ---- Paramters to be defined ---- %
% Number of different routes to search
popsize = 500;

% Number of Generations to search along
generation = 200;

% Percentage of population that will be choosen for crossover
crossover_perc = 30;

% Type of Crossover to use
% Key 1: Order 1 Crossover
% Key 2: Sequential Constructive Crossover
% Key otherwise: Do Nothing
key = 2;

% Percentage of population that will be choosen for mutation
mutation_perc = 100;

% After how many iterations the mutation should occur
mutation_itr_limit = 10;

% To define what cost funtion to use
% For Distance = 1
% For Price    = 0
alpha = 1;






%% ---- Variable for GA ---- %
% DO NOT CHANGE %

% Number of cities for GA
numcities = size(cities, 1);
% Initializing fitness for population
max_fit = Inf;
% Fitness of each chromosome in a population
genration_fit = zeros(popsize, 1);
% Initialize the population
pop = zeros(popsize, numcities);
% Record of the best fitness in across all generations
fit_hist = zeros(generation, 2);
% Initializing the variable to check for mutation
mutation_iter = 0;
% Variable for scaling
scale = 1;
% Number of of total chromosome to take for crossover
crossover_sz = (crossover_perc/100) * popsize;
% Number of of total chromosome to take for mutation
mutation_sz = (mutation_perc/100) * popsize;

% % Figures Handles
map = figure('Name','Best Route in current Generation', 'Numbertitle','off');
hist = figure('Name', 'Generation History', 'Numbertitle','off');


%% COST MATRIX

% COST MATRIX: PRICE
% Scale price between (0,1)
max_price = max(price_matx(:));
price_matx = price_matx / max_price;

% COST MATRIX: DISTACNE
distance_matx = dmatrix(cities);
% Scale distance between (0,1)
max_dist = max(distance_matx(:));
distance_matx = distance_matx / max_dist;



% Selecting scale for different cost matrix,
% The variable is only needed for printing actual cost
% Does not take part in actual GA
if(alpha == 1)
    scale = max_dist;
elseif (alpha == 0)
    scale = max_price;
end


% ---- The final cost matrix ---- %
% If alpha = 1, then only distance and if alpha = 0, then only price
cost_matx = alpha*(distance_matx) + (1-alpha)*price_matx;

%% ------------------------------------- GENETIC ALGORITHM ---------------------------- %%

% Initialize the population with random permutation
for t = 1:popsize
    pop(t,:) = randperm(numcities);
end

%% Plot World Map
figure(map)
createmap(cities)

%% Generation Iterations
% Sanity check
if(crossover_sz >=2)
    for itr = 1:generation
        %% Evaluate Fitness
        tic
        for i = 1:popsize
            genration_fit(i) = fitness(pop(i,:), cost_matx);
        end
        [best_fit, best_chr_id] = min(genration_fit(:,1));
        [worst_fit, ~] = max(genration_fit(:,1));
        median_fit = median(genration_fit(:,1));
        
        % Update the fitness history
        fit_hist(itr,:) = [median_fit, best_chr_id];
        
        %% Plotting figures
        
        % Plotting World Map with routes
        best_chr = pop(best_chr_id, :);
        lat = zeros(length(best_chr) ,2);
        lon = zeros(length(best_chr) ,2);
        figure(map)
        for i = 1:(length(best_chr) -1)
            lat(i,:) = [cities(best_chr(i)).Lat, cities(best_chr(i+1)).Lat];
            lon(i,:) = [cities(best_chr(i)).Lon, cities(best_chr(i+1)).Lon];
        end
        lat(length(best_chr),:) = [cities(best_chr(length(best_chr))).Lat, cities(best_chr(1)).Lat];
        lon(length(best_chr),:) = [cities(best_chr(length(best_chr))).Lon, cities(best_chr(1)).Lon];
        pt = plotm(lat, lon, 'r-', 'linewidth',1);
        title(sprintf('Generation: %d / %d | Best Value: %g', itr, generation, best_fit*scale));
        
        
        % Plotting the Fitness History
        figure(hist)
        subplot(1,2,1)
        fit = [worst_fit, median_fit, best_fit];
        fit_name = {'Worst'; 'Median'; 'Best'};
        barh(fit, 0.9);
        title(sprintf('Generation: %d / %d', itr, generation));
        xlabel('Fitness')
        set(gca,'yticklabel',fit_name);
        axis([0  15 -inf inf]);
        
        subplot(1,2,2)
        hold on
        plot(itr, fit_hist(itr), 'ro');
        ylabel('Fitness')
        xlabel('Generation')
        axis([0 generation 0 inf]);
        drawnow
        
        % Time delay to see the result
        pause(0.01)
        
        
        % To delete the last plot lines to draw new ones
        if (itr~=generation)
            delete(pt);
        end
        
        %% Crossover Operators
        
        % Percentage of population choosen for cross over
        [~, best_pos] = sort(genration_fit);
        crossover_chr = pop(best_pos(1:crossover_sz), :);
        
        % ---- Generating new population ---- %
        pop = zeros(popsize, numcities);
        for i = 1:popsize
            switch key
                case 1 % Order 1 Crossover Operator
                    pop(i,:) = order1crossover(crossover_chr);
                case 2 % Sequential Constructive Crossover Operator (SCX)
                    pop(i,:) = scx(crossover_chr, cost_matx);
                otherwise
                    % do nothing
            end
        end
        
        %% Mutation: Reciprocal Exchange Mutation
        
        % Do mutation only when the fitness remains same for mutation_itr_limit consecutive iterations
        if(itr ~=1 && fit_hist(itr) == fit_hist(itr-1))
            mutation_iter = mutation_iter + 1;
        end
                
        if(mutation_iter == mutation_itr_limit)
            mutation_iter = 0;
            mutation_chr = pop(1:mutation_sz, :);
            
            % Mutate the chromosomes
            mutation_chr = mutation(mutation_chr);
            
            % Add it to the population
            pop(1:mutation_sz, :) = mutation_chr;
        end
    end
toc
disp('Final Solution:')
disp(int32(best_fit*scale))



else
    disp('Number of chromosome choosen for crossover')
    disp(crossover_chr)
    disp('Wrong number of chromosome for crossover should be more than 2')
end