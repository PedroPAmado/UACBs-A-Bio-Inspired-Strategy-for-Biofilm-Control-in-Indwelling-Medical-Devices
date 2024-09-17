function [ ] = mat2avi( FrameMat, frameRate, VideoName, PlayVideo)
%MAT2AVI Summary converts the frames contained in the matrix FrameMat
%to a video in format .avi with name VideoName (set e.g 'Video_processed') and with frame rate
%specified by frameRate.
%If PlayVideo is set to true ( PlayVideo = 1) Matlab will play the video after the
%conversion.
%   Author: ARTORG Center, CVE-Group
%       - Damiano Schirinzi

%Write Video
video = VideoWriter(VideoName);
video.FrameRate = frameRate;
open(video);
for i = 1:size(FrameMat,3)
    writeVideo(video,FrameMat(:,:,i));
end
close(video)

%Play Video if PlayVideo is true
if(PlayVideo)
    i=1;
    videoAvi = VideoReader(VideoName);
    while hasFrame(videoAvi)
        mov(i)= im2frame(readFrame(videoAvi));
        i=i+1;
    end
    
    figure
    imshow(mov(1).cdata,'Border','tight')
    movie(mov,1,videoAvi.FrameRate)
    close(gcf)
end
end

