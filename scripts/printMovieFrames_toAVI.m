% make fish movie from fish movie frames

load('/Users/audreysederberg/Dropbox/postdoctoral projects/theory project 2014 to 2015/code for respository/fishmovie20110225_frames.mat')


v = VideoWriter('fishmovie.avi');
v.FrameRate = 60;
v.Quality = 100;
open(v);

[nFrames, fHeight, fWidth] = size(frames);

hFig = makeMyFigure(fWidth/10, fHeight/10);
hFig.Color = [1 1 1];

colormap(gray)
imagesc(squeeze(frames(1, :, :)))
axis image
xyoff
set(gca, 'nextplot', 'replacechildren');


for k = 1:nFrames
    imagesc(squeeze(frames(k, :, :)))
    axis image
    xyoff
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);