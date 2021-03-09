function [pointcloud] = retrieve_data(directory)

%This function retrieve_data checks the file extension of the point cloud 
%file in the respective directory, which is given as the input

% input:  - directory [path]

% output: - pointcloud: the point cloud matrxi [n x 3] where X, Y, Z are
%                       the first, second and third column entries
           
%           -> Remark: Please note that only four file extensions are
%           supported, which are, .asc, .txt, .ply, and .pcd

% =========================================================================

%% retrieving data

[~,~,filetype] = fileparts(directory);

switch filetype 
    
     case '.asc'
        
        data = readmatrix(directory,'FileType', 'text');
        pointcloud=data(:,1:3); %XYZ Cartesian coordinates of all the points in the point cloud

         
     case '.txt'
        
        data = readmatrix(directory,'FileType', 'text');
        pointcloud=data(:,1:3); %XYZ Cartesian coordinates of all the points in the point cloud

            
     case '.pcd'
        
        data = pcread(directory);
        pointcloud = double(data.Location); %XYZ Cartesian coordinates of all the points in the point cloud
        
        
     case '.ply'
        
        data = pcread(directory);
        pointcloud = double(data.Location); %XYZ Cartesian coordinates of all the points in the point cloud
        
    otherwise 
        warning('The file format not supported')
        
end


end

