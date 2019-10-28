%This script reads a set of input images of  Hemotoxylin & Eosin stained 
%transverse/ maximal lung alveolar sections that are present in a given folder. 
%
% image criterion: the images must be taken under brightfield microscope with a 20x objective.
% image labeling criterion: the image names of images in a given folder should containing a running count at its end 
%
%the program is capable of reading multiple images , perform
%binarization using automated thresholding based on Otsu's method, perform erosion, dilation to fill the holes in the
%alveolar walls and invert the image such that the alveolar spaces are black and the walls are
%white. Subsequently, the program uses Laplacian operator to extract the boundaries of the
%alveolar spaces, and using the bwboundaries function, detects alveolar spaces as distinct
%objects. The object detection enables computation of the minimum distance between each
%alveolar space in an iterative manner where distances between each alveolar space with the
%neighboring spaces are compared and minimum distances corresponding to the narrowest
%segments is extracted. About 100- 200 measurements per field are taken and the average value
%representing the mean interseptal wall thickness, standard deviation and the corresponding
%dispersion index is returned to the user. In addition, the module computes parameters such as
%alveolar area, average alveolar diameters can be extracted depending on user’s interest.
% 
% Script execution rules: 
%   1. the file handling of this script is simplistic. Hence, this code
%   should be copied to the folder containing the images
%   2.the first line of the code shoul have the updated directory (same folder name as images)
%   3. the filename template should be specified. 
% flow chart:
% read input image one at a time ---> binarization by thresholding based on otsu's method
%  ---> image preprcessing:closing operation and area filtering --->
%  calculate laplacian to detect edges/ segmentation ---> detect boundaries
%  of alveoli---> filter boundaries to remove capillaries --> for every
%  alveoli in the image , compute distances with every other alveoli (one at a time), calculate
%  the minimal distances and the overall minimum distances---> store overall
%  minimum distance as the septal thickness of that alveoli---> repeat this
%  for every alveoli---> extract septal wall thickness and calculate mean
%  septal wall thickness ---> store it in global variables corresponding to
%  that image.
%  

% 
% Input
%     sourceImage: 8 bits grayscale image;
%     threshold: threshold [0 1] for binarization.
%     
% Output:
%     output: 3 array variables :
%     alveolarseptalthickness,alveolarsurfacearea,alveolardiameter. copy
%     and paste results in excel sheet.
%     
%     
% If you have any ideas or improvements about this method, welcome to contact me 'aravind245@gmail.com'.




clc; % clear command window
clear all;% delete all variables
close all;

dir('D:/Users/zeiss/Documents/MATLAB/52 wo Female WT Non-Breeders/'); % specify the directory in which the images to be analyzed are kept. Past this code in that folder and execute it from that folder.
n_fields = 8;% number of images to be analyzed in the folder
filename =('52 wo F WT nonbreed m19-126a 4466 '); % filename template. check spelling and spaces. 

% initializing global variables
alveolarseptalthickness = [];
alveolarsurfacearea = [];
meanalveolardiameter = [];
alveolardiameter = [];
scale = 8.24;% conversion factor for converting pixel measurements to microns(can change depending on imaging scope)
allareas = [];
alldiameters = [];
 

for o=1:n_fields
    % read image one at a time
    frame = o;
    I = imread([filename num2str(o) '.tif']);
    figure;imshow(I);

    level = graythresh(I);
    %I1 = im2bw(I,level);% Thresholding is done according to Otsu's method.
    %change thresholding parameter manually if needed. Have a threshold value greater than level
    %if more regions of alveoli has to be included and vice versa. This has
    %to be determined emperically based on the image lighting quality.
    %decreasing threshold value increases stringency.
    I1 = im2bw(I, 0.8);
    I1 = bwareaopen(I1,1000);% area filter for removing small unwanted remnants of binarization.
    I1 = imcomplement(I1);
    sc = strel('disk',10);
    I1 = imclose(I1,sc);% performs closing operation to close the gaps between the alveolar walls.
    I1 = bwareaopen(I1,1000);% area filter for removing small unwanted remnants of closing operation
    figure; imshow(I1);

    I2 = I1(:,:,1);% converting I1 to a gray scale 2D image. 
   
    I3 = del2(double(I2));% calculating the laplacian of the image for edge detection.
    figure; imshow (I3);

    [d1,d2] = size(I3);%extracting dimensions of the image 
    % initializing variables
    A= [];
    AA = [];
    IA =[];
    OD = [];
    ID = [];
    thickness = {};
    wallthickness = {};
    count = 1;
    area2 = [];
    OD2 = [];

    [B1,L1] = bwboundaries(I3,'noholes');% no holes parameter enables extraction of outer boundaries
    %imshow(label2rgb(L1, @jet, [.5 .5 .5]))% displaying the labels in color
    hold on

    stats = regionprops(L1,'Area','Centroid');
    for k = 1:length(B1)% for calculating measurements based on outer boundary
      boundary = B1{k};%extracting boundaries of objects one at a time)
      plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)% plotting boundary on top of image.
      % alveoalr perimeter extraction.
      delta_sq = diff(boundary).^2;
      circumference = sum(sqrt(sum(delta_sq,2))); 
      area = stats(k).Area;%outer area calculation
      outerdiameter = 2 * sqrt(area/pi);% calculation of the outer diameter
      OD = vertcat(OD,outerdiameter);% storing outerdiameter values
      AA = vertcat(AA, area);% storing area values.
      sk = sprintf('%2.2d', k);
      text(boundary(1,2)-10, boundary(1,1)-10,sk,'color','m');%numbering and labelling the boundaries
    end


    lenB = length(B1);% extracting number of objects
     B2 = B1;% initializing another variable to hav esame value as B1
  
% the following while loop is to drop the boundary values of small capilaries
% that are embedded in the walls of alveoli for more accurate detection of septal wall thickness.
    while count <= lenB-1
    boundary1 = B1{count};
    area1 = stats(count).Area;
     if area1 <= 700% from experience and manual measurements in imageJ, 
         %the area of capillaries was determined to be less than 700
      B2{count} = [];
     end
     if area1 >=700
         area2{count} = area1;
     end
    count = count+1;
    end 

    % removing the empty values in B2 and storing it in B3
    B3 = cell(lenB,1);
    m = 1;
    n =1;
    area3 = [];
    ncount =1; 

   while m <= length(B2)-1
        tempboundary1 = B2{m};        
        if isempty(tempboundary1) ==0
            B3{n} = tempboundary1;
            n= n+1;
        end
        m = m+1;
   end
   
   % removing corresponding measured values of area, diameter of alveoli.
   cc= 1;

   while cc <= length(area2)-1
       temparea = area2{cc};
   if isempty(temparea) ==0
            area3{ncount} = temparea;
            diameter2 = 2 *sqrt(temparea/pi);
         OD2{ncount} = diameter2;
         alldiameters = vertcat(alldiameters,diameter2);
            ncount = ncount+1; 
   end
        cc = cc+1;
   end
    A = vertcat(A,area3);
    
    % display the new corrected boundary list not including the capillaries.
    figure; imshow (I3);
    hold on

    for k2 = 1:length(B3)
      boundaryb = B3{k2};%extracting boundaries of objects one at a time
      if isempty(boundaryb) ==0
      plot(boundaryb(:,2), boundaryb(:,1), 'r', 'LineWidth', 2)    
      sk2 = sprintf('%2.2d', k2);
      text(boundaryb(1,2)-10, boundaryb(1,1)-10,sk2,'color','m');%numbering the boundaries  
      end
    end
    
    % this part of the program calculates the minimal distance between
    % alveoli after examining all possible distances between object pairs of a single alveoli
    % at a given time
    X1 = {};
    X2 = {};
    Y1 = {};
    Y2 = {};  

    for i = 1: length(B3)
       boundary1 = B3{i};%extracting boundaries of objects one at a time
        if isempty(boundary1) ==0
        boundary1x = boundary1(:,2);%storing the X coordinate of the boundary of ith object
         boundary1y = boundary1(:,1);%storing the Y coordinate of the boundary of ithe object
        x1=1;
        y1=1;
        x2=1;
        y2=1;
        overallMinDistance = inf;
         for j = i+1:length(B3) %j is a counter variable that comapares ith alveoli with the rest of alveoli
        boundary2 = B3{j};
        if isempty(boundary2) ==0
        for p = 1: size(boundary2,1)
         boundary2x = boundary2(p,2);% storing the X coordinate of the boundary of jth object
         boundary2y = boundary2(p,1);% storing the Y coordinate of the boundary of jth object
         allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2); % computing all pairwise distance between boundaries
        [minDistance(p), indexOfMin] = min(allDistances); % extracting the minimum distance between two alveoli
      
        %stores the coordinates and overall minimum distance for alveoli i when compared with the remaining alveoli
        if minDistance(p) < overallMinDistance
            x1 = boundary1x(indexOfMin); 
            y1 = boundary1y(indexOfMin);
            x2 = boundary2x;
            y2 = boundary2y;
            overallMinDistance = minDistance(p);
       end
       end
         minDistance = min(minDistance);
        end
         end 

 thickness {i} = overallMinDistance;%extracting th eseptal wall thickness
    X1{i} = x1;
    X2{i} = x2;
    Y1{i} = y1;
    Y2{i} = y2;
    end
    end
    
    % converting the cell variables to double variables
thickness = cell2mat(thickness');
X1 = cell2mat(X1');
X2 = cell2mat(X2');
Y1 = cell2mat(Y1');
Y2 = cell2mat(Y2');
count =1;
dummy =1;
k=1;
max = d1/40;

for count = 1:length(thickness)
% filtering the thickness vector to prevent artifact measurements
    t = thickness(count);
    if t < max && t ~=Inf
        line([X1(count), X2(count)], [Y1(count), Y2(count)], 'Color', 'y', 'LineWidth', 2);
        si = sprintf('%2.2d', count);
        text(X1(count)-5,Y1(count)-5,si,'color','g');
        count = count+1;
        wallthickness {k} = t;
        k= k+1;
    end
    dummy = dummy+1;
end

 
wallthickness = cell2mat(wallthickness');%converting cell variable to double
meanwallthickness = mean(wallthickness);% calculating mean wall thickness for the image
stdwallthickness = std(wallthickness);% calculating standard deviation of wall thickness of image
A = cell2mat(A');
dia = cell2mat(OD2);
meanalveolararea = mean(A);% calculating mean alveolar area for the image
stdalveolararea = std(A);% calculating standard deviation of alveolar areas for the image
meanalveolardia = mean(dia);% calculating mean alveolar diameter for image


    alveolarseptalthickness = horzcat(alveolarseptalthickness,meanwallthickness/scale);% converting pixel mesaurements to microns 
    alveolarsurfacearea = horzcat(alveolarsurfacearea,meanalveolararea);% converting pixel mesaurements to microns 
    meanalveolardiameter = horzcat(meanalveolardiameter, meanalveolardia);% converting pixel mesaurements to microns 
    allareas = vertcat(allareas,A);

end

 alveolarseptalthickness = alveolarseptalthickness';% output 1 
 alveolarsurfacearea = alveolarsurfacearea';% output 2 
 alveolardiameter = meanalveolardiameter';% output 3