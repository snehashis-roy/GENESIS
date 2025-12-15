
function success = run_hallucination(atlas_vol_list,...
    subj_vol_list,isskull,numsubjpatch,numatlaspatch,numNN,psize,numcpu,outl,...
    modal1,modal2,atlasstartslice,subjectstartslice,atlaswmpeaks,subjectwmpeaks,...
    outputvol)
% function run_hallucination(atlas_vol_list,...
%     subj_vol_list,isskull,numsubjpatch,numatlaspatch,numNN,psize,numcpu,outl,...
%     inputmodal,outputmodal,atlasstartslice,subjectstartslice,atlaswmpeaks,subjectwmpeaks,...
%     outputvol)
% 
% ATLAS_VOL_LIST : A comma separated list of atlas volumes. Example:
%          /home/user/ute1.nii.gz,/home/user/ute2.nii.gz,/home/user/ct.nii.gz
%          Paths should not have any special character like space. Number of atlas 
%          volumes must be 1 greater than number subject volumes.
% 
% SUBJ_VOL_LOST : A comma separated list of atlas volumes. Example:
%          /home/user/ute1.nii.gz,/home/user/ute2.nii.gz,/home/user/ct.nii.gz
%          Paths should not have any special character like space. Number of atlas 
%          volumes must be 1 greater than number subject volumes.
%
% ISSKULL :  Does the subject contain skull? Yes or no. If ISSKULL is yes, then
%          ATLASWMPEAKS and SUBJECTWMPEAKS must be provided (see below). Otherwise,
%          the peaks will be estimated automatically from skullstripped image
%          based on modality (see below)
%
% NUMSUBJECTPATCH :  Number of subject patches used to generate EM parameters, 
%          reasonable values are 1E4-1E5.
% 
% NUMATLASPATCH : Number of atlas patches used to generate EM parameters, a
%        reasonable value is between 1E4-1E5.
% 
% NUMNN : Number of nearest neighbors to be used. Very small number 
%         leads to noisy result (FAST), very high number leads to better
%         result (SLOW). Reasonable values are 10-40. Default is 15. 
%         Usually the complexity of the originally proposed algorithm is O(N M^2) 
%         where M=NumAtlasPatches and N=NumSubjectPatches. However, M and N
%         are usually very large. So the reduced complexity is O(MN D^2) where
%         D=NUMNN. The dictionary size is MD^2, so choose M and D such that
%         MD^2 <1E6. Otherwise, the code can be very slow.
%
% PATCHSIZE : Size of patches to be used. A reasonable value is [3 3 3] or [3 3 1]
%
% NUMCPU : Number of parallel processing cores
% 
% OUTLIER : A flag to identify how much outlier is present in the image.
%         Default is 'none'. There are 4 options, 'low','moderate','high', and
%         'none'.
% 
% INPUTMODALITY : Options are T1,T2,FL or CT. For multichannel inputs, use a
%         comma separated list, e.g. t1,t2,fl. In this case with 3 input channels,
%         the ATLAS_VOL_LIST must have 4 files (3 inputs and 1 output), and
%         SUBJECT_VOL_LIST has 3 files
%         
% 
% OUTPUTMODALITY : Options are T1, T2, FL or CT.
% 
% ATLASSTARTSLICE :  Atlas slices where the patches are collected from. 0-->Default.
%         Choose slices from the middle of the atlas image. The range is
%         usually [ATLASSTARTSLICE, ATLASSTARTSLICE+15]
% 
% SUBJECTSTARTSLICE :  Subject slices where the patches are collected from. 0-->Default.
%          Choose slices from the middle of the subject image. The range is
%          usually [SUBJECTSTARTSLICE, SUBJECTSTARTSLICE+15]
% 
% ATLASWMPEAKS : For images with skull, enter WM peaks in comma separated format,
%          e.g. 199.0,299.0. If ISSKULL is false, the peaks will automatically
%          be calculated based on modality. If 0, it will be computed automatically.
% 
% SUBJECTWMPEAKS :  For images with skull, enter WM peaks in comma separated format,
%          e.g. 199.0,299.0. If 0, it will be computed automatically.
% 
% OUTPUTVOL : Full path and name of the output volume, e.g.
% /home/user/synthct.nii


% (For reference, version=12)



fprintf('Setting Params\n');
%atlas.vol1=atlas_vol1;
%atlas.vol2=atlas_vol2;
%atlas.mask=atlas_mask;
%subj.vol=subj_vol;
%subj.mask=subj_mask;
fprintf('Atlas volumes = %s.\n',atlas_vol_list);
temp=strsplit(atlas_vol_list,',');
T=length(temp);
fprintf('%d atlas inputs detected.\n',T);
atlas.vol=cell(T,1);
subj.vol=cell(T-1,1);
for t=1:T
    atlas.vol{t}=temp{t};
end
temp=strsplit(subj_vol_list,',');
if length(temp)~=T-1 || T==1
    fprintf('ERROR: Number of atlas inputs (%d) must be 1 greater than number \n',T);
    fprintf('ERROR: of subject inputs (%d). Exiting.\n',length(temp));
    success=0;
    return;
end
T=length(temp);

fprintf('%d subject inputs detected.\n',T);
for t=1:T
    subj.vol{t}=temp{t};
end

if strcmp(psize,'3x3x1')
	param.PatchSize=[3 3 1];   
elseif strcmp(psize,'3x3x3')
	param.PatchSize=[3 3 3];   
end
fprintf('Using %s sized patches.\n',psize);
param.MinNumBasis = str2num(numNN);
fprintf('Using %d nearest neighbors.\n',param.MinNumBasis);
param.NumAtlasPatches = str2num(numatlaspatch);
fprintf('Using %d atlas patches.\n',param.NumAtlasPatches);

param.NumSubjectPatches = str2num(numsubjpatch);
fprintf('Using %d subject patches.\n',param.NumSubjectPatches);

param.atlasoutlierratio = [0.01 1.5];
param.subjectoutlierratio = [0.01 1.5];

try
    x=atlasstartslice;
catch e
    fprintf('Using default atlas slices.\n');
    atlasstartslice=0;
end
if str2num(atlasstartslice)>1   
    param.AtlasSliceRange=[str2num(atlasstartslice) str2num(atlasstartslice)+15];
    fprintf('Atlas slice range =[%d, %d].\n',param.AtlasSliceRange(1),...
       param.AtlasSliceRange(2) );
end
try
    x=subjectstartslice;
catch e
    fprintf('Using default atlas slices.\n');
    subjectstartslice=0;
end
if str2num(subjectstartslice)>1   
    param.SubjectSliceRange=[str2num(subjectstartslice) str2num(subjectstartslice)+15];
    fprintf('Subject slice range =[%d, %d].\n',param.SubjectSliceRange(1),...
       param.SubjectSliceRange(2) );
end


if strcmpi(isskull,'no')
    param.isskull=0;
elseif strcmpi(isskull,'yes')
    param.isskull=1;
end
if param.isskull
    fprintf('Subject and atlas images contain skull. WM peaks must be provided.\n');
else
    fprintf('Subject and atlas images do not contain skull. WM peaks are optional.\n');
end


temp=strsplit(atlaswmpeaks,',');
if str2num(temp{1})==0
    fprintf('Atlas WM peaks not mentioned. Using automatic ones.\n');
else
    if length(temp)~=T
        fprintf('Number of atlas WM peaks must be %d. Exiting.\n',T);
        success=0;
        return;
    else
        for t=1:T
            param.atlaswmpeak(t)=str2num(temp{t});
        end
    end
      fprintf('Atlas WM peaks are: %.2f',param.atlaswmpeak(1));
    for t=2:T
        fprintf(', %.2f',param.atlaswmpeak(t));
    end
    fprintf('\n');
end
temp=strsplit(subjectwmpeaks,',');
if str2num(temp{1})==0
    fprintf('Subject WM peaks not mentioned. Using automatic ones.\n');
else
    if length(temp)~=T
        fprintf('Number of subject WM peaks must be %d. Exiting.\n',T);
        success=0;
        return;
    else
        for t=1:T
            param.subjectwmpeak(t)=str2num(temp{t});
        end
    end
    fprintf('Subject WM peaks are: %.2f',param.subjectwmpeak(1));
    for t=2:T
        fprintf(', %.2f',param.subjectwmpeak(t));
    end
    fprintf('\n');
end


param.outlier=outl;
fprintf('Using %s outlier removal.\n',outl);

param.NumCPU=str2num(numcpu);
fprintf('Using %d cpu.\n',param.NumCPU);

%Default values
% param.isinputT2=0;
% param.isT2=0;
% param.isCT=0;
% 
% if strcmpi(modal1,'T2') || strcmpi(modal1,'CT')
%     param.isinputT2=1;
% end
% if strcmpi(modal2,'T2')
%     param.isT2=1;
% elseif strcmpi(modal2,'CT')
%     param.isCT=1;
% end
% if param.isT2
% fprintf('Atlas 2nd modality is T2.\n');
% else
% fprintf('Atlas 2nd modality is not T2.\n');
% end
% if param.isCT
% fprintf('Atlas 2nd modality is CT.\n');
% else
% fprintf('Atlas 2nd modality is not CT.\n');
% end

temp=strsplit(modal1,',');
if length(temp)~=T
    fprintf('ERROR: Number of input modalities must %d (i.e. same as number of input images).\n',T);
    success=0;
    return;
end
param.inputmodal=modal1;
param.outputmodal=modal2;
fprintf('Input modalities: %s\n',modal1);
fprintf('Output modality :   %s\n',modal2);
if strcmpi(modal2,'CT')  % For output modality of T2 and CT, extra precaution needed 
    param.isCT=1;
else
    param.isCT=0;
end
if strcmpi(modal2,'T2')
    param.isT2=1;
else
    param.isT2=0;
end

% if param.isinputT2
% fprintf('Atlas 1st modality is T2.\n');
% else
% fprintf('Atlas 1st modality is %s.\n',modal1);
% end


param.outname=outputvol;
% fprintf('Output will be written in \n');
% fprintf('%s\n',param.outname);


image_hallucination(subj,atlas,param);


success=1;

