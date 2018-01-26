%% 213 Final Project: Analysis of Galaxies
% Alexander Caravaggio
% 10187672
% Due Date: April 12 2017
% All data retrieved from http://leda.univ-lyon1.fr/leda/fullsql.html
% Using the search parameters 'SELECT mod0, vgsr WHERE mod0 IS NOT NULL'
% The data obtained is from the Centre de Recherche Astrophysique de Lyon
% which is associated with the Universty of Lyon

%% Age of the Universe

gdata = csvread('galaxies1.csv'); % Retrieves the galaxy data from the csv file
modx = c2r(gdata(:,1)); % Retrieves the distance data in units of modular distance and turns the column vector into a row vector
v = c2r(gdata(:,2)); % Retrieves the distance data in units of m/s and turns the column vector into a row vector
xpars = 10.^((modx./5)+1); % Converts the distance from modular distance to parsecs 
xm = xpars.*((3.086.*10^13)); % Converts the distance from parsecs to meters

figure(1) % Opens the first figure
plot(xm,v,'.') % Plots the velocity of the galxies vs. the distance they are from the earth
hold on;

[r1,m,b1] = regression(v,xm); % preforms a regression on the data
regv = (1./m).*xm; % Uses the regression data to get new velocity points so that a regression line can be plotted
plot(xm,regv) % Plots the regression line
xlabel('Distance [km]')
ylabel('Velocity [km/s]')
title('Velocity of Galaxies vs. Position')
legend('Galaxies','Regression Line')
hold off

t = m*3.17098E-8; % Calculates the age of the universe in years
dt = t-(r1*t); % Calculates the uncertainty of the age using the regression value

tbill = t*10^-9; % Formats the age of the universe in billions of years for readability
dtbill = dt*10^-9; % Formats the uncertainty in the age of the universe in billions of years for readability

fprintf('The age of the universe is %.1f +/- %.1f billion years\n',tbill,dtbill)

%% Future of the Universe Part I
% Plots the galaxies at their proper distances in a random direction away
% from the milky way. Model does not take into account that gravity causes
% the expansion to not be uniform. 

% Initializes all arrays
test = zeros(size(v)); 
theta = zeros(size(v));
x = zeros(size(v));
y = zeros(size(v));

j = 1; % Initiates the counter

while j<= length(v)
    test(j) = rand; % Creates a random number for each galaxy
    theta(j) = test(j) * 2*pi; % Multiplies that number by 2*pi to get an angle
    x(j) = xm(j)*cos(theta(j)); % Retrieves the x-coordinate for a given angle
    y(j) = xm(j)*sin(theta(j)); % Retrieves the y-coordinate for a given angle
    j = j+1;
end

H0 =(67.15./(3.086.*10^22)); % The hubble expansion constant in units of 1/s
ly = 9.461.*10.^15; % For conversion of light years into meters;

% Finds the radius and x and y coordinates for the observavle universe
eoou = ly.*t; 
phi = linspace(0,2.*pi);
xo = eoou.*cos(phi);
yo = eoou.*sin(phi);

% The earth's coordinates remain the same and is plotted at the origin
ex = 0;
ey = 0;

figure(2)
whitebg('black') % Changes the background to black so that it looks like space
q = plot(x,y,'y.'); % Plots the position of the galaxies
title('The Expansion of the Universe')

% Holds the axis constant
xlim([-10^23,10^23]) 
ylim([-10^23,10^23])
hold on

plot(ex,ey,'b*') % Plots the Milky Way Galaxy
hold on

H0 =((67.15.*10^3)./(3.086.*10^22)); % The hubble expansion constant in units of 1/y
c = 3.0*10^8;

N = 500; % Determines the speed of the animation

r0 = xm; %Creates a new variables for initial position

tg = ones(length(xm),N); % Creates a matrix that will contain all of the time indices
count = 1;

phold = (linspace(t,100*10^9,N)).*(3.154.*10^7); % Array that contains all the time entries for which the calculations will be made

while count <= length(xm)
    
    tg(count,:) = tg(count,:).*phold; % Applies each element of the array to a column
    count = count+1;

end
tg = transpose(tg); % Flips the matrix so that all entries in the same row have the same value

% Intializes matrices in which 
coeff = ones(N,length(xm));
rfin=zeros(N,length(xm));
vfin=ones(N,length(xm));
a = ones(N,length(xm));
xc = ones(N,length(xm));
yc = ones(N,length(xm));

i = 2; 
coeff(1,:) = 1./exp(H0*tg(1,:)); % Calculates the first row of coefficients
rfin(1,:) = xm; % Sets the first distance value to be the one provided from the data file

while i< N % Will loop for all 500 different time increments
    
   a(i-1,:) = (coeff(i-1,:).*(exp(H0.*tg(i,:))-exp(H0.*tg(i-1,:)))); % Sets the expansion coefficients
   rfin(i,:)= r0 +(a(i-1,:).*r0); % Calculates the distance from the milky way given an expansion coefficient
   vfin(i,:) = (rfin(i,:)-r0)./(tg(i,:)-tg(i-1,:)); % Calculates the velocity after expansion
   
   coeff(i,:) = 1./(exp(H0*tg(i,:))); %Calculates the new coeficcient that is used to calculate the expansion coefficient
   
   % Calculates the new x and y coordiante for each galaxy
   xc = rfin(i,:).*cos(theta); 
   yc = rfin(i,:).*sin(theta);
   
   % Updates the position of the galaxies
   drawnow 
   q.XData = xc;
   q.YData = yc;
   
   r0 = rfin(i,:); % Sets the old final distance as the new initial distance
   i = i+1; % Adds to the counter
   
   pause(10^-10) % Stops calculations for a bit after each iteration so that the animation runs at the right pace
end

%% Future of the Universe Part II
% Calculates the point at which no other galaxies will be visible from
% earth with the naked eye as well as the radius of the hubble sphere, the
% sphere beyonnd which all objects moving away faster than the speed of 
% light and therefore can no longer be casually linked with our universe.
%  All galaxies are assumed to only output yellow light

ylambda = 570*10^-9; % The wavelength at which yellow light is found
rlambda = 750*10^-9; % The wavelenght at which red light is found

vcrs = ((((rlambda/ylambda)^2)-1)/(((rlambda/ylambda)^2)+1))*c; % Calculates the speed something must be travelling so that yellow light 
% is redshifted out of the visual spectrum

rrs = vcrs.* (tg(2,1)-tg(1,1)); % Calculates the distance at which that redshift occurs

rhs = c/H0; % Calculates the radius of the hubble sphere

trs = (log((rrs-xm(1,1))/(coeff(1,1)*xm(1,1)))/H0)*3.17098E-7; % Calcluates the time at which the nearest galaxy will be redshifted out 
% of the visible spectrum in years
trstril = trs/(10^12); % Converts the time to trillions of years for readability

tsol = (log((rhs-xm(1,1))/(coeff(1,1)*xm(1,1)))/H0)*3.17098E-7; % Calculates the time at which the nearest galaxy is outside the hubble
% sphere in years
tsoltril = tsol/(10^12); % Converts the time to trillions of years for readability

fprintf('In %.1f trillion years no other galaxies will be visible with the naked eye in the night sky due to redshift\n',trstril)
fprintf('In %.1f trillion years no other galaxies will be casually linked to our own\n',tsoltril)
