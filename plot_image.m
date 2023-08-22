function plot_image(A, fig, name, location)

%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

screen_size = get(0,'ScreenSize');
screen_size(1:2)=[];

horizontal_step =floor(0.89*screen_size(1)/2);
vertical_step = floor(0.95*screen_size(2)/2);

cur_fig=figure(fig);
clf
imagesc(A);
hold on
set(fig,'NumberTitle','off')
set(fig,'Name',name)
bar_handle=colorbar;
colormap gray
axis tight
axis image
%title(name)
set(cur_fig,'position',[5+(location(1)-1)*(8+horizontal_step) -30+(location(2)-1)*vertical_step horizontal_step vertical_step],'Color','k','MenuBar','none')
set(get(gca,'Title'),'Color','w')
set(bar_handle,'YColor','w','FontSize',8,'Box','off')
drawnow