%% cleaning
clf;
clc;
clear;
close all;

%% parameters (different settings for test)
Temp_init = 50;
Temp_min = 1;
alpha = 0.999;
n_iter = 100;
                                            
%% import cities
% There are three different test cases we are running
% 'European_Cities_Big.xlsx'     = 41 cities
% 'European_Cities_Medium.xlsx'  = 20 cities
% 'European_Cities_Short.xlsx'   = 10 cities

filename = 'European_Cities_Medium.xlsx';
sol = readtable(filename);
cities = table2struct(readtable(filename));
sol = [sol; sol(1,:)];
[num_cities, ~] = size(sol);




% Generate Figures
map = figure('Name','Best Route in current Generation', 'Numbertitle','off');
hist = figure('Name', 'Distance Evolution', 'Numbertitle','off');

% Create World Map
figure(map)
createmap(cities)

%% loops
Temp = Temp_init;
count = 0;
convergence = 0;
while Temp > Temp_min
    if Temp == Temp_init
        [prev_cost, dist] = cost(sol, num_cities);
    end   
    iter = 1;
    while iter < n_iter
        new_sol = neigh(sol, num_cities);
        [new_cost, new_dist] = cost(new_sol, num_cities);
        ap = acc_prob(Temp, prev_cost, new_cost);
        if ap > rand
            sol = new_sol;
            prev_cost = new_cost;
            dist = new_dist;
            convergence = 0;
        end
        iter = iter + 1;
	convergence = convergence + 1;
       lat = zeros(size(sol ,1) - 1, 2);                                    % All lat in the best population
    lon = zeros(size(sol ,1) - 1, 2);                                    % All lon in best population

    figure(map)
    for i = 1:(size(sol,1) - 1)
        lat(i,:) = [sol.Lat(i), sol.Lat(i+1)]; 
        lon(i,:) = [sol.Lon(i), sol.Lon(i+1)];
    end
   
    
    % To close the loop, and adding the first and last cities lat lon
    pt = plotm(lat, lon, 'r-', 'linewidth',1);
    title(sprintf('Temperature: %g | Best Distance: %g', Temp, dist));
 

    
        figure(hist)
        hold on
        plot(count, dist/1000, 'ro');
        ylabel('Distance')
        xlabel('Iterations')
        axis([0 inf 0 inf]);
        title('Cost Evolution');
        drawnow
             pause(0.01)
    if (Temp > Temp_min)
        delete(pt);                                                         % To delete the last plot lines drawn on the map
    end
        
    end
    
    count = count + 1;
   
    %-------------------------------- Plotting -------------------------------%
    % First plotting the lines on the map
    
    %plan evolution 


    Temp = alpha*Temp;
 
       
end



