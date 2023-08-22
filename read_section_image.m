function [A, resolution, slice_z, filename] = read_section_image

%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

slice_z=0;
A=0;
resolution=0;
%% UI Dialog to open scan file
[filename, pathname] = uigetfile('*.dcm; *.M0*', 'Select scan file');
if isequal(filename,0) %case where nothing valid is selected
   disp('User selected Cancel')
   return
else
    full_name = fullfile(pathname, filename); %full_name is the full name of the file

%% DICOM format 
    if findstr(full_name,'.dcm')>1 %case where we have DICOM input
        info = dicominfo(full_name);
        A = (dicomread(info));        
        A=double(A);
        resolution = info.PixelSpacing';
    end
    
%% STRATEC pQCT format
    if findstr(full_name,'.M')>1   %case where we have pQCT input
        [A, resolution, slice_z] = parse_pQCT (full_name);
        A=wiener2(A,[8 8]);     %apply adaptive filtering  
    end
    

%%    
end
