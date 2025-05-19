% DESCRIPTION
% This Project implements a Genetic Algorithm to plan a path for a robot, 
% across an environment, with the intention of avoiding obstacles in minimal
% distance.

% QUESTIONS
%The following lines of Code pose three questions to the user, asking to
%specify their desired selection algorithm, crossover algorithm and mutation
%algorithm, for use in the Genetic Algorithm.
selectionType = input("Enter Selection Algorithm 0)Roulette Wheel, 1)Tournament, 2)Linear Rank?");
crossoverType = input("Enter Crossover Algorithm 0)1-Point, 1)Uniform?");
mutationType = input("Enter Mutation Algorithm 0)Value Encode, 1)Bit Flip?");

% PARAMETERS
%The following 7 lines of code establish variables that will largely remain
%static throughout the development of this genetic algorithm. The first
%three concern the confines of the environment that needs to be traversed,
%establishing the maximum and minimum coordinates of the image in the map
%variable.
map=im2bw(imread('random_map.bmp'));
start = [1 1];
finish=[500 500];
%The 3 lines below establish the maximum amount of turning points for the
%robot, as well as amount of times the genetic algorithm will loop to get
%an optimum candidate solution, and the size of the population each time
%the genetic algorithm is processed.
noOfPointsInSolution = 10;
iterations = 20;
populationSize = 100;

% 1 INITIAL POPULATION
%Our initial random population is generated using the createPopulation
%function. This population has 100 chromosomes of randomly generated x and
%y coordinates.
population = createPopulation(populationSize,map,noOfPointsInSolution);
population = [population zeros(populationSize,1)];
%The fittest variable, will be used to store the fittest chromosome of each
%generation, with an extra row, to store the final fittest chromosome to be
%plotted upon our map.
fittest = zeros(iterations+1,(size(population,2)));

% GENETIC ALGORITHM
%The following for loop repeats for the amount of iterations stated above.
for k=1:iterations
    
    % 2 EVALUATE FITNESS
    %In order to define the fittest chromosomes to be passed onto the next
    %generation, we must first define the value of each chromosome. This
    %entails measuring their suitability as a candidate solution to our
    %problem. 
    %The for loop below, cycles through our population and returns a set of
    %fitness values, adjacent to each chromosome.
    for i=1:populationSize
        fitness = findFitness(population, map, noOfPointsInSolution);
    end
    %The following three lines of code then take the generated fitness
    %values, add them to the 21st column of our population, in
    %correspondence to their measured chromosome, and then sorts them
    %descendingly, with the highest fitness on the bottom of our
    %population. The Fittest (At the bottom of our population) is then
    %stored as this generation's fittest variable in the fittest matrix.
    population(:,21) = fitness;
    population = sortrows(population,21);
    fittest(k,1:21) = population(end,1:21);

    %ELITISM
    %In order to generate our new population as part of the iterative
    %nature of our genetic algorithm, me must set aside the top 30% of our
    %prior population, to ensure that we generate a new population with a
    %generally higher fitness than the prior. Each generations population
    %should in turn get fitter.
    newPopulation = zeros(populationSize,20);    
    newPopulation(1:(0.3*populationSize),:) = population(populationSize-(0.3*populationSize-1):populationSize,1:20);
    newPopulationSize = (0.3*populationSize);
    
    %As we have filled places in our new population with the fittest 30
    %chromosomes, we should seek to fill the remaining 70 spaces with
    %relatively fit chromosomes. In order to ensure fitness and diversity,
    %we use selection and mutation algorithms to produce the rest of our
    %chromosomes.
    while(newPopulationSize<populationSize)
                
        % 3 SELECTION
        if (selectionType==0)
            %In order to perform the Roulette Wheel Selection, we must
            %first accumulate the weight of each chromosome, by calculating
            %the proportion of their fitness to the rest of the population.
            %This is displayed in the line below.
            weights=population(:,21)/sum(population(:,21));
            choice1=RouletteWheelSelection(weights);
            choice2=RouletteWheelSelection(weights);
        end
        if (selectionType==1)
            %Unlike the Roulette Wheel, the Tournament Selection algorithm
            %hanldes only the fitness of each chromosome of the population.
            %Hence why the only parameter is the population, with its 21st
            %column containing each chromosome's fitness.
            choice1=TournamentSelection(population);
            choice2=TournamentSelection(population);
        end
        if (selectionType==2)
            %The final type of selection is linear rank selection, which
            %functions similarly to the roulette wheel, but doesnt use the
            %same weight calculated previously in the program. It instead
            %calculates a different one using the fitness, as will be seen
            %explained in the function.
            choice1=LinearRankSeletion(population);
            choice2=LinearRankSeletion(population);
        end
        %After two new chromosomes have been nominated from our population,
        %we are ready to cross them over.
        temp_chromosome_1 = population(choice1, 1:20);
        temp_chromosome_2 = population(choice2, 1:20);
        
        % 4 CROSSOVER
        if (rand<0.6)
            if(crossoverType==0)
                %Given the nature of the crossover, the algorithm should
                %take both nominated chromosomes to crossover and generate
                %two new unique chromosomes to be added to the population.
                [offspring]=OnePointCrossover(temp_chromosome_1,temp_chromosome_2);
                temp_chromosome_1 = offspring(1,1:20);
                temp_chromosome_2 = offspring(2,1:20);

            end
            if (crossoverType==1)
                %The uniform crossover also takes both chromosomes to
                %substitute Alleles with one another, hence the following
                %three lines below.
                [offspring]=UniformCrossover(temp_chromosome_1,temp_chromosome_2);
                temp_chromosome_1 = offspring(1,1:20);
                temp_chromosome_2 = offspring(2,1:20);
            end
        end       
        
        % 5 MUTATION
        %Given the natural limitations of just using crossover functions,
        %we employ mutation to add diversity to our population. I've used a
        %relatively low mutation rate of 0.2 for both chromosomes, to
        %ensure that the best candidate solutions are not mutated and lost.
        if (mutationType==0)
            if(rand<0.3)
                %The Value Encoding mutation functions by taking the
                %chromosome and adding or subtracting a random value from
                %one of its Alleles, as can be seen explained in the
                %function's code.
                temp_chromosome_1= ValueEncodingMutation(temp_chromosome_1);
            end
            if (rand < 0.3)
                temp_chromosome_2 = ValueEncodingMutation(temp_chromosome_2);
            end
        end
        if(mutationType==1)
            if(rand<0.3)
                %Flip mutation functions differently to value encode
                %mutation as it manipulates the order of the alleles to
                %achieve mutation. An explaination of this can be seen in
                %the functions code.
                temp_chromosome_1= FlipMutation(temp_chromosome_1);
            end
            if (rand < 0.3)
                temp_chromosome_2 = FlipMutation(temp_chromosome_2);
            end
        end
        %Now that we have selected, crossed-over and mutated our
        %chromosomes, we are able to add them back to our new population,
        %already containing the top 30% of previously generated chromosomes
        %that had the highest fitness of the last generation.       
         
        newPopulationSize=newPopulationSize+1;
        newPopulation(newPopulationSize,1:20)=temp_chromosome_1;
        %The following if statement ensures that we do not exceed the
        %confides of the set population size each time k loops. Exceeding
        %would result in an error as we would generate populations of
        %different sizes each generation.
        if(newPopulationSize<populationSize)
            newPopulationSize=newPopulationSize+1;
            newPopulation(newPopulationSize,:)=temp_chromosome_2;
        end
    end
    %The following line of code replaces last generations population with
    %our newly generated new population. This is now ready to undergo the
    %same algorithms as the previous generation.
    population(:,1:20)=newPopulation;
end
%Now that we have left the for loop, we have ourselves our final
%population, containing our optimum candidate solution. This is what we will
%use to plot our points on the map.
%Firstly we must measure the fitness of each chromosome once again, to
%determine the fittest.
for i=1:populationSize
    fitness = findFitness(population, map, noOfPointsInSolution);
end
%Once we have measured the fitness, we can select the final generation's
%fittest chromosome to be added to our fittest variable. The final row
%containing our optimum solution.
population(:,21) = fitness;
population = sortrows(population,21);
fittest(end,1:21) = population(end,1:21);
fittest = sortrows(fittest,21);
solutionOnMap = displayMapAndLine(fittest,map,start,finish);

%% 1 Population Creation
function population = createPopulation(populationSize, map, noOfPointsInSolution)
  %Initially when creating our population, we should ensure that each value
  %adheres to the predetermined boundaries of our map, those being what we
  %have below.
  lowerBounds = ones(1,noOfPointsInSolution);
  upperBounds = size(map,1)*ones(1,noOfPointsInSolution); 
  %Given that each of the ten turning points must have both a x coordinate
  %and a y coordinate, I create a population containing 20 alleles for each
  %chromosome.
  population = zeros(populationSize,noOfPointsInSolution*2);
  %The following for loop then generates two sets of ten values, the first
  %ten being the x coordinates and the second ten being the y coordinates.
  for i=1:populationSize
      xCoordinate = randi([min(lowerBounds),max(upperBounds)],1,noOfPointsInSolution);
      yCoordinate = randi([min(lowerBounds),max(upperBounds)],1,noOfPointsInSolution);
      path=[sort(xCoordinate(:));sort(yCoordinate(:))]';      
      population(i,:)=path; 
      %The above coordinates are then sorted ascendingly. Given the nature
      %of our problem and the characteristics of the environment I assume
      %that it would be best to ensure that both the x and y coordinates
      %are increasing at each turning point, as we are traversing a path
      %starting at 1,1 and ending at 500,500. 
      %The sorted coordinates are then added in two parts to each
      %chromosome to our population, and then added to our path to be
      %parsed to the main method.
  end
end

%% 2 MEASURE FITNESS AND EVALUATION
%Within the context of our project, I decided to define a chromosomes
%fitness along its ability to avoid obstacles in a minimum distance. This
%would mean accumulating data that would connote whether the point lies in
%an object, whether the line between it and the previous point intersect
%with an object, and finally whether the length of said line is short
%enough to compete with other chromosomes.
function fitness = findFitness(population, map, noOfPointsInSolution)
  %Initially, we must predefine the number of chromosomes we are measuring
  %for, hence the two lines below.
  indPopulation = size(population,1);
  fitness=zeros(indPopulation,1);
  %The following for loop places each chromosome through a function that
  %calculates the chromosomes overall line fitness, and the fitness of its
  %individual points.
  for k=1:indPopulation
      xCoordinate = population(k,1:noOfPointsInSolution);
      yCoordinate = population(k, noOfPointsInSolution+1:end);
      fitness(k) = 1/(calculatePointFitness(xCoordinate,yCoordinate,map)+calculateGradientFitness(xCoordinate,yCoordinate,map));
      %In order to make the fitness scores clear, I divide the value by
      %one, to ensure the "fittest" chromosomes have the largest fitness
      %value assigned to them (in oppose to the lowers)
  end
end

function pathFitness = calculatePointFitness(xCoordinate,yCoordinate,map)
  %This function uses the binary map we had previously read into the
  %program and sees if each of our x and y chromosomes lie on a point that
  %equals 0. On the binary map, objects are represented with the number 0.
  numPoints=length(xCoordinate);
  pathFitness=0;
  for l=1:numPoints
      if map(yCoordinate(l),xCoordinate(l))==0
          pathFitness=pathFitness+100;
          %Each time it lies on an object, a penalty score is assigned to
          %the chromosome. This will be further computed and combined with
          %other penalties to give us the fitness of the chromosome.
      end
  end
end

function gradientFitness = calculateGradientFitness(xCoordinate,yCoordinate,map)
  %This function is called to find the "Gradient fitness" between two
  %points. This means it takes all the values that lie on the path between
  %our coordinate the previous coordinate, and measures whether or not it
  %lies upon a value 0 (An object).
  gradientFitness=0;
  for l=1:11
      if l==1
          x1=1;
          y1=1;
          x2=xCoordinate(l);
          y2=yCoordinate(l);
          %An exception is made for our first coordinate as its previous
          %point (for which it will be connected to) remains static
          %throughout the program and will not be randomly generated.
      end
      if(l==11)          
          x1=xCoordinate(l-1);
          y1=yCoordinate(l-1);
          x2=500;
          y2=500;
          %Similarly, our final turning point is also given exceptional
          %circumstances as the turning point following it will also be
          %static throughout the program.
      end
      if (l>1&&l<11)
          x1=xCoordinate(l-1);
          y1=yCoordinate(l-1);
          x2=xCoordinate(l);
          y2=yCoordinate(l);
          %For the rest of our coordinates, as seen in the code above, they
          %are stored in two x and two y coordinates for measurement of the
          %line connecting them.
      end
      %The following four lines of code calculate the points between each x
      %and y coordinates and store them in the variable inLineCoords. 
      %The first of the four takes the amount of points that lie on a line
      %between the two coordinates using Distance Formula, and stores it in
      %the variable called numPoints.
      numPoints=max(abs(x2-x1), abs(y2-y1))+1;
      %The next two lines then use the amount of points, alongside the
      %function "linspace" which is able to precisely calculate what values
      %lie between our x and y points. These values are stored in the
      %variables xpoints and ypoints.
      xpoints=round(linspace(x1,x2,numPoints));
      ypoints=round(linspace(y1,y2, numPoints));
      %Finally they are reassembled in the vector below, on two different
      %rows connoting their axis.
      inLineCoords=[xpoints', ypoints'];
      %The following for loop then relays through the line coordinates that
      %lie on our line and utilises the same concept as previously states
      %in our calculatePointFitness function to assign a penalty score to
      %any point that lies on a map value equal to 0.
      for i=1:size(inLineCoords,1)
          x=inLineCoords(i,1);
          y=inLineCoords(i,2);
          if map(y,x)==0
              gradientFitness=gradientFitness+100;
          end
          %Each time the for loop loops, a further penalty is assigned to
          %each chromosome as this is in direct correlation to the length
          %of the line. This is displayed below It is only a small penalty
          % as avoiding objects is still considered to be the priority of 
          % our genetic algorithm.
          gradientFitness=gradientFitness+1;
      end
  end
end

%% 3.1 Selection 1: Roulette Wheel
function choice = RouletteWheelSelection(weights)
  %The first step to creating a Roulette Wheel Selection is by totalling
  %the sum of all of the fitnesses of each chromosome in our population, in
  %order to find the proportion of each chromosome's fitness to the overall
  %fitness of all of them.
  accumulation = cumsum(weights);
  %To find this we use the function cumsum as stated above, to find the
  %cumulative sum of all of the weights passed from our main method.
  p = rand();
  chosen_index = -1;
  %Ultimately the nature of the roulette wheel selection is to be random
  %within the confides of our weighted environment, so we generate a random
  %value to nominate our chromosome.
  for index = 1 : length(accumulation)
    if (accumulation(index) > p)
        %The line below will choose our chromosome and be sent back to the
        %main method for crossover and mutation.
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end

%% 3.2 Selection 2: Tournament Selection
function choice = TournamentSelection(population)
  %Tournament selection operates differently to Roulette Wheel Selection by
  %nominating two chromosomes completely randomly. It then compares each
  %chromosomes fitness with the other and returns the chromosome with the
  %highest fitness.
  populationSize = size(population,1);
  randomTournamentIndex = randperm(populationSize, 2);
  chosen_index=0;
  %As the fitness scores are already stored in column 21 of our population
  %matrix, we compare each of our random chromosomes by said value.
  if(population(randomTournamentIndex(1,1),21)>population(randomTournamentIndex(1,2),21))
      chosen_index = randomTournamentIndex(1,1);
  end
  if(population(randomTournamentIndex(1,2),21)>=population(randomTournamentIndex(1,1),21))
      chosen_index = randomTournamentIndex(1,2);
  end
  %The fitter chromosome is returned to the main method for crossover and
  %mutation.
  choice=chosen_index;
end

%% 3.3 Selection 3: Linear Rank Selection
function choice = LinearRankSeletion(population)
  %Linear Rank Selection functions similar to Roulette Wheel Selection,
  %but instead of ranking purely on each chromosome's fitness, it measures
  %the fitness of each chromosome and assigns an incremental value to each
  %chromosome in our ordered list.
  populationSize = size(population,1);
  %The below formula works by reassigning a value of incremental fitness
  %scores to each chromosome, to make their ranking linear in nature.
  seletionProbabilities =(1-2/populationSize)+(2*(1:populationSize)*(2-(1+2/populationSize)))/populationSize^2;
  %This information is then parsed to the RouletteWheelSelection function
  %as it will randomly select a value based on the new ranked fitness, in
  %oppose to its proportional weighting.
  choice=RouletteWheelSelection(seletionProbabilities);
end

%% 4.1 Crossover 1: 1-Point Crossover
function offspring = OnePointCrossover(parent1,parent2)
  %The one point crossover acts by taking our two selected chromosomes and
  %selecting a random point on each of the chromosome. It then selects all
  %the alleles following this point on the first one and swaps it with
  %those following the point on the second chromosome.
  offspring = zeros(2,20);
  for i=1:2      
      crossoverPoint = randperm(20,1);
      offspring(i,1:crossoverPoint)=parent1(1:crossoverPoint);
      offspring(i,crossoverPoint+1:end)=parent2(crossoverPoint+1:end);
      %The three lines above display how a random allele is selected in our
      %chromosomes, and those values to the right of it are swapped with
      %its adjacent points in the other chromosome.
  end
end

%% 4.2 Crossover 2: Uniform Crossover
function offspring = UniformCrossover(parent1,parent2)
  %Uniform crossover entails generating a crossover probability and
  %swapping alleles on an independent basis, unlike K Point Crossovers
  %which tend to swap multiple alleles at a time.
  child1=zeros(1,20);
  child2=zeros(1,20);
  %As the new offspring will not have fitness scores yet, the two lines
  %above connote that there will be columns in each matrix, exclusively
  %for the x and y coordinates.
  temp1=parent1;
  temp2=parent2;
  for i=1:20
      if rand>=0.5
          %Each time the for loop initiates, a new crossover probability is
          %generated and if the selected chromosome lands on position i, it
          %is swapped
          temp1(i)=parent2(i);
          temp2(i)=parent1(i);
      end
      child1(i)=temp1(i);
      child2(i)=temp2(i);
      %Each time the for loop loops, each of the offspring is gradually
      %constructed, to be returned to the main method.
  end
  offspring(1,1:20)=child1;
  offspring(2,1:20)=child2;
  %The two lines above show each new child being placed into the
  %offspring matrix.
end

%% 5.1 Mutation 1: Value Encode Mutation
function mutatedChromosome = ValueEncodingMutation(chromosome)
  %Permanent encode mutation, mutates a chromosome by randomly selecting a
  %point in the chromosome and adding or subtracting a small amount from
  %it.
  point = randi([1, length(chromosome)]);
  mutation = randi([1,10]);
  x=rand;
  %The three lines above connote how a random point is selected in our
  %chromosome to be mutated, and by how much of a value. The variable x
  %will be the randomly generated decider in whether the mutation value is
  %subtracted or added to the allele selected randomly in the variable
  %point.
  if x>0.5 && chromosome(point)+mutation<500
      chromosome(point)=chromosome(point)+mutation;
  end
  if (x<=0.5)&&chromosome(point)-mutation>0  
      chromosome(point)=chromosome(point)-mutation;
  end
  mutatedChromosome=chromosome;
  %The new mutated chromosome is then returned to the main method.
end

%% 5.2 Mutation 2: Flip Mutation
function flippedChromosome = FlipMutation(chromosome)
  %Flip mutation functions in selecting two random locations in our
  %chromosome and swapping the values, as well as the alleles that lie in
  %between them, thus reversing the order for part of our chromosome.
  point1 = randi([1, length(chromosome)]);
  point2 = randi([1, length(chromosome)]);
  %The above two lines show the points in our chromosome being randomly
  %generated.
  while(point2 == point1)
      point2 = randi([1, length(chromosome)]);
      %This while loop ensures that if the two randomly generated points
      %are the same, a new point will be generated to be used.
  end
  if(point2 < point1)
      temp=point1;
      point1=point2;
      point2=temp;
      %If the circumstances are met for the mutation to be possible, the
      %two points are swapped, thus creating mutated 
  end
  tobeflipped = chromosome(point1:point2);
  chromosome(point1:point2) = fliplr(tobeflipped);
  %The above two lines, then take our values and flip the values between
  %our point1 and point2.
  flippedChromosome = chromosome;
  %The flipped chromosome is then returned to the main method.
end

%% Display Map and Line
function solution = displayMapAndLine(fittest,map,start,finish)
  %The final function of the project will construct the path of our optimum
  %candidate solution, that should hopefully traverse the environment,
  %avoiding obstacles in minimal distance.
  solution1=fittest(end,1:20);
  finalXCoordinates=solution1(1,1:10);
  finalYCoordinates=solution1(1,11:20);
  solution=zeros(12,2);
  solution(2:11,1)=finalXCoordinates(1,1:10);
  solution(2:11,2)=finalYCoordinates(1,1:10);
  solution(1,:)=start(:);
  solution(12,:)=finish(:);
  %The above 8 lines of code create a two columned matrix with the x axis
  %coordinates on the left and y axis coordinates on the right. It also
  %adds our beginning coordinate 1,1 at the beginning and the final
  %coordinate 500,500 at the end, making 12 coordinates overall.
  clf;
  imshow(map);
  rectangle('position',[1 1 size(map)-1],'edgecolor','k');
  line(solution(:,1),solution(:,2));
  %The above four lines construct our image containing our final fittest
  %turning points, with a line that connected them, thus displaying our
  %final path.
end
