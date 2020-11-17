%Generate images of the test functions

%House cleaning
close all; clear all; clc



%% Givens

%Function to use
% fxn = 'Ackley';
% fxn = 'Griewank';
% fxn = 'Rosenbrock';
fxn = 'Schwefel';



%% Generate Data

%---Get the info---
switch fxn
    
    case {'Ackley'}
        dim_1 = [-5:0.05:5];     %Dimension 1 - Limits
        dim_2 = dim_1;
        min = [0,0,0];
        dim_3 = zeros(length(dim_1),length(dim_2)); %Preallocate Z values

        for i1 = 1:length(dim_1)
            for i2 = 1:length(dim_2)
                dim_3(i1,i2) = ackley([dim_1(i1);dim_2(i2)]);
            end
        end
 
    case {'Griewank'}
        dim_1 = [-10:0.05:10];     %Dimension 1 - Limits
        dim_2 = dim_1;
        min = [0,0,0];
        dim_3 = zeros(length(dim_1),length(dim_2)); %Preallocate Z values

        for i1 = 1:length(dim_1)
            for i2 = 1:length(dim_2)
                dim_3(i1,i2) = griewank([dim_1(i1);dim_2(i2)]);
            end
        end
 
    case {'Rosenbrock'}
        dim_1 = [-10:0.05:10];     %Dimension 1 - Limits
        dim_2 = dim_1;
        min = [0,0,0];
        dim_3 = zeros(length(dim_1),length(dim_2)); %Preallocate Z values

        for i1 = 1:length(dim_1)
            for i2 = 1:length(dim_2)
                dim_3(i1,i2) = rosenbrock([dim_1(i1);dim_2(i2)]);
            end
        end
 
    case {'Schwefel'}
        dim_1 = [-500:1:500];       %Dimension 1 - Limits
        dim_2 = dim_1;
        min = [420.9687,420.9687,0];
        dim_3 = zeros(length(dim_1),length(dim_2)); %Preallocate Z values

        for i1 = 1:length(dim_1)
            for i2 = 1:length(dim_2)
                dim_3(i1,i2) = schwefel([dim_1(i1);dim_2(i2)]);
            end
        end
        
    otherwise
        disp('Incorrect function selected.')
        disp('Valid options are "Ackley", "Griewank", "Rosenbrock", or "Schwefel"')
        return
        
end



%% Generate Plots

%---Surface Plot---
figure('units','normalized','outerposition',[0 0 1 1])
mesh(dim_1,dim_2,dim_3','edgecolor','none','linewidth',1); hold on; grid on; grid minor
shading interp  %Makes the shading of the plot colors smooth
alpha 1       %Makes the plot somewhat transparent
plot3(min(1),min(2),min(3),'k.','markersize',50)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis (Cost)')
% title(fxn + "'s Function")
set(gca,'fontsize',48)
set(gca,'fontname','times')  % Set it to times


%---Contour Plot---
figure('units','normalized','outerposition',[0 0 1 1])
contour(dim_1,dim_2,dim_3','linewidth',3); hold on; grid on; grid minor
plot(min(1),min(2),'k.','markersize',50)
xlabel('X Axis')
ylabel('Y Axis')
% title(fxn + "'s Function")
set(gca,'fontsize',48)
set(gca,'fontname','times')  % Set it to times


