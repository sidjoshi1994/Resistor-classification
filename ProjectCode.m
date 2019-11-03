clc
clear
%% Read Image
I = imread('Test_resistor1.jpg');
%I = imread('Test_resistor3.jpg');
%I = imread('Test_resistor2.jpg');  % not accurate segmentation
%I = imread('Test_resistor4.jpg');
%I = imread('Test_resistor5.jpg');
%I = imread('Test_resistor6.jpg');  % not accurate segmentation
%I = imread('Test_resistor7.jpg');
%I = imread('Test_resistor8.jpg'); % not accurate segmentation
figure
imshow(I);
Gray = rgb2gray(I);
%% Noise Filtering
M = medfilt2(Gray); % To filter out potential salt and pepper noise
G = imgaussfilt(M);
se = strel('square',5);
F1 = imopen(G,se);
%imshowpair(F1,Gray,'montage');
%% Line Detection
F2 = imerode(F1,se);
E = edge(F2,'canny',[0.15 0.6]);
%imshow(E);
%E1 = edge(F1,'canny');
%imshowpair(E1,E,'montage');
[H,theta,rho] = hough(E);

figure
imshow(imadjust(rescale(H)),[],...
       'XData',theta,...
       'YData',rho,...
       'InitialMagnification','fit');
xlabel('\theta (degrees)')
ylabel('\rho')
axis on
axis normal 
hold on
colormap(gca,hot)

P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = theta(P(:,2));
y = rho(P(:,1));
plot(x,y,'s','color','black');

lines = houghlines(E,theta,rho,P,'FillGap',5,'MinLength',7);

figure, imshow(Gray), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
%imshow(Gray);

% Get the strongest line
% First column is rho index, second is theta index
R = houghpeaks(H,1);
% Find out how much it needs to rotate
thetaPeak = theta(R(1,2))-90;
% Fix the image
I3 = imrotate(I,thetaPeak);
%figure
%imshow(I3);

%% Localisation
sigma_0 = 70;
d_0 =  20;
%count = 0;

I3_filtered = medfilt3(I3);
I3_filtered_lab = rgb2lab(I3_filtered);

%Coloumn standard deviation 
C_deviation = std(I3_filtered_lab(:,:,1)) + 10 * (std(I3_filtered_lab(:,:,2))...
    +std(I3_filtered_lab(:,:,3)));
r_region = zeros(size(C_deviation));

for i = 1:length(C_deviation)
    if (C_deviation(i) > sigma_0)
        r_region(i) = 1;
        %count = count + 1;
    end
end
k = find(r_region);
xmin = k(1);
width = k(length(k))-k(1);

%Row standard deviation 
C_deviation_v = std(I3_filtered_lab(:,:,1),0,2) + 10 * (std(I3_filtered_lab(:,:,2),0,2)...
    +std(I3_filtered_lab(:,:,3),0,2));
r_region_v = zeros(size(C_deviation_v));

for i = 1:length(C_deviation_v)
    if (C_deviation_v(i) > 60)
        r_region_v(i) = 1;
    end
end

k_v = find(r_region_v);
ymin = length(r_region_v)-k_v(length(k_v));
height = k_v(length(k_v))-k_v(1)+10;

I_final = imcrop(I3_filtered,[xmin ymin width height]);
figure
I4 = imrotate(I_final,90);
imshow(I4);
%%