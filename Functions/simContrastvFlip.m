function [OptimumTheta MaxContrast]=simContrastvFlip(T1pair,TR,varargin)
%%
% function [OptimumTheta MaxContrast]=simContrastvFlip(T1pair,TR)
% T1pair should be a vector with 2 T1 values
% TR should be the value of the repetition time used
%
% this function plots the signal of a GRE sequence run with repetition time
% TR as a function of the flip angle
% additionally it calculates and plots the contrast that can be generated between
% these two sample and its optimum flip angle
%
% the ouptput from the function are OptimumTheta and MaxContrast


if ~exist('T1pair','var') || ~exist('TR','var')
    disp('ERROR: T1pair and TR must be specified')
    return
end

if length(T1pair)~=2
    disp('ERROR: T1pair must contain 2 values')
    return
end

if ~exist('tColours','var')
    tColours = {[66 122 223]/255,[108 158 80]/255,[223 76 76]/255};
end

if nargin>2
    ShowFigure=varargin{1};
else
    ShowFigure=1;
end;

%% doing all the calculations

theta = linspace(0,180,1000);
sig = zeros(length(theta),2);
for iT1 = 1:2
    sig(:,iT1) = (1-exp(-TR/T1pair(iT1))).*sin(theta*pi/180)./(1-cos(theta*pi/180)*exp(-TR/T1pair(iT1)));
end
contrast=abs(sig(:,1)-sig(:,2));
[maximum position]=max(contrast);
OptimumTheta=theta(position);
MaxContrast=maximum;


%% plotting the results
if ShowFigure
    figure
    % set(gcf,'Position',[    50   679   740   409])
    set(gcf,'Position',[   36   719   517   400])
    set(gcf,'Color',[1 1 1]);
    plot(theta,sig(:,1),'linewidth',2,'color',tColours{1}); hold all
    plot(theta,sig(:,2),'linewidth',2,'color',tColours{2});
    h1 = plot(theta,contrast);
    set(h1,'linewidth',2,'color','k')
    xlabel('Flip Angle (degrees)')
    ylabel('Relative signal')
    grid on
    title(['TR = ' num2str(TR) 'ms'])
    legend(['T1 = ' num2str(T1pair(1)) 'ms'], ['T1 = ' num2str(T1pair(2)) 'ms'], 'Contrast');
    xlim([0 180])
    text(45 ,maximum,['Optimum Flip=',num2str(round(OptimumTheta)),'degrees'])
    text(45 ,maximum*2,['Max Contrast=',num2str(round(maximum*1000)/1000),' '])
    fontScale(1.2)
end