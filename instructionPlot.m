whichComp=1;

if whichComp==1
    basePath='/Users/ttli/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    colorPath='/Users/ttli/Dropbox (Brown)/sharedMatlabUtilities/';
    
elseif whichComp==2
    basePath='/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    colorPath='/Users/LabManager/Dropbox (Brown)/sharedMatlabUtilities/';
    
elseif whichComp==3
    basePath='/Users/mattLab/Dropbox (Brown)/arousalLearningPerception/vwm_task';
    colorPath='/Users/mattLab/Dropbox (Brown)/sharedMatlabUtilities/';
end


addpath(genpath('./'))
addpath(genpath(colorPath))


calibrationFile='labMonitorData_2019-10-30_426c';


load(fullfile(basePath,calibrationFile))

%%

length = 30;
red = [1, 0, 0];
red2= [215, 96, 20]/255;
pink = [240, 227, 90]/255;

colors_r = [linspace(red(1),red2(1),length)', linspace(red(2),red2(2),length)', linspace(red(3),red2(3),length)'];
colors_p = [linspace(red(1),pink(1),length)', linspace(red(2),pink(2),length)', linspace(red(3),pink(3),length)'];

% plot random markers on the map and assign them the colors created
S=9;   % marker size
figure()
hold on
% xlim([1 30])
% xticks(1:1:30)
%axes('YColor','none')
set(gca,'ytick',[],'ycolor','none')
set(gca,'xtick',[],'xcolor','none')

%set(gca,'ytick',[])
colorAll=randperm(30);

for i=1:30
    odd=rand(1);
    if odd>0.1 && i<=15
        plot(i,10,'o','MarkerEdgeColor','k','MarkerFaceColor',colors_r(colorAll(i),:),'markersize',S)
    elseif odd>0.1 && i>15
        plot(i,10,'o','MarkerEdgeColor','k','MarkerFaceColor',colors_p(colorAll(i),:),'markersize',S)
    elseif odd<0.1
        c=rand(1);
        if c<=1/3
            plot(i,10,'o','MarkerEdgeColor','k','MarkerFaceColor','g','markersize',S)
        elseif c>1/3 && c<2/3
            plot(i,10,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',S)
        else
            plot(i,10,'o','MarkerEdgeColor','k','MarkerFaceColor','c','markersize',S)
        end
    end
end



%%

length = 10;

red = [1, 0, 0];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),length)', linspace(red(2),pink(2),length)', linspace(red(3),pink(3),length)'];

blue = [0, 0, 1];
lightblue = [189, 240, 255]/255;
colors_b = [linspace(blue(1),lightblue(1),length)', linspace(blue(2),lightblue(2),length)', linspace(blue(3),lightblue(3),length)'];

green= [0, 1, 0];
lightgreen = [224, 255, 189]/255;
colors_g = [linspace(green(1),lightgreen(1),length)', linspace(green(2),lightgreen(2),length)', linspace(green(3),lightgreen(3),length)'];

yellow=[1,1,0];
lightyellow=[255, 253, 189]/255;
colors_y = [linspace(yellow(1),lightyellow(1),length)', linspace(yellow(2),lightyellow(2),length)', linspace(yellow(3),lightyellow(3),length)'];

colorAll=cat(4,colors_p,colors_g,colors_b,colors_y);

S=8;   
figure()
hold on
% xlim([1 30])
% xticks(1:1:30)
set(gca,'ytick',[],'ycolor','none')
set(gca,'xtick',[],'xcolor','none')



for i=1:30
    change=rand(1);
    a=Randi(10);
    if i==1 || change<0.1
        mu=Randi(4);
        C(i)=mu;
        plot(i,10,'o','MarkerEdgeColor','k','MarkerFaceColor',colorAll(a,:,mu),'markersize',S)
    else
        mu=C(i-1);
        C(i)=mu;
        plot(i,10,'o','MarkerEdgeColor','k','MarkerFaceColor',colorAll(a,:,mu),'markersize',S)
    end
   
    
end












