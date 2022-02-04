clear;


% vid = VideoWriter("vid");
% open(vid);
% y = zeros(length(S), 1);
% for i = 1:length(t)
%     for j = 1:length(y)
%         y(j) = S{j}(3, i) - S{j}(2, i);
%     end
%     bar(x, y, 1);
%     title(sprintf("Time: %8.2f year", t(i)/3600/24/365.25*2));
%     frame = getframe(gcf());
%     writeVideo(vid, frame);
% end
% close(vid);
