function wmpeak=find_WM_peak(input_image,modality,brainmask)

% function wmpeak=find_WM_peak(input_image,modality,brainmask)
% input_image       Either stripped or unstripped T1/T2/PD/FL image. If the
%                   image is unstripped, a brain mask is required for 
%                   correct estimation of white matter peak
% modality          T1/T2/PD/FL
% brainmask         Not required if the input image is stripped. For images
%                   with skull, it is required.
% 
if nargin==2
    brainmask=[];
end
if ~strcmpi(modality,'T1') && ~strcmpi(modality,'T2') ...
        && ~strcmpi(modality,'PD') && ~strcmpi(modality,'FL') 
    fprintf('Modality must be either of T1, T2, PD, or FL (FLAIR). \n');
    fprintf('Returning -1. Exiting.\n');
    wmpeak=-1;
    return;
end
if ischar(input_image)
    subj=load_untouch_nii(input_image);
elseif isnumeric(input_image)
    subj.img=double(input_image);
end
if ~isempty(brainmask)
    if ischar(brainmask)
        mask=load_untouch_nii(brainmask);
    elseif isnumeric(brainmask)
        mask.img=double(brainmask);
    end
    U=unique(mask.img);
    if length(U)~=2
        fprintf('Brainmask volume (%s) must be binary 0 and 1 volume.\n',brainmask);
        fprintf('Returning -1. Exiting.\n');
        wmpeak=-1;
        return;
    end
    subj.img=double(subj.img).*double(mask.img);
end


[IntensityLimit, ~ , ~, ~]=OutlierReduction(subj.img,0,0,[1 2]);
[temp , ~]=Crop3D(subj.img,0,0,0);
temp(temp<IntensityLimit(1) & temp>0)=0;
temp(temp>IntensityLimit(2))=0;

delta=(IntensityLimit(2)-IntensityLimit(1))/80;
[pA,ff]=find_peaks(single(temp),[],[],0,2,delta);
if strcmpi(modality,'T1')
    wmpeak=pA(end);
else    
    indx=find(ff==max(ff));
    wmpeak=pA(indx(1));    
end
fprintf('%.4f\n',wmpeak);

