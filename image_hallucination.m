%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outvol=image_hallucination(SUBJ,ATLAS,PARAM)
%
% 
% This is for image normalization as well as synthesis using probabilistic
% model and 2-sparsity (GENESIS)
%
%
% SUBJ.vol{1},..SUBJ.vol{T}, ATLAS.vol{1}..ATLAS.vol{T+1} contain the 
% subject and atlas volumes (full path), respectively. Number of ATLAS
% volumes must be 1 greater than number SUBJ volumes.
%
% PARAM is a structure, containing the following fields:
%
% param.MinNumBasis = Number of nearest neighbors to be used. Very small 
% number leads to noisy result (FAST), very high number leads to better
% result (SLOW). Reasonable values are 10-40. Default is 15.
%
% param.NumAtlasPatches = Number of atlas patches used to generate EM
% parameters, reasonable values are 1E4-1E5.
%
% param.NumCPU = Number of parallel processing cores
% 
% param.NumSubjectPatches = Number of subject patches used to generate EM
% parameters, reasonable values are 1E4-1E5.
%
% Usually the complexity of the originally proposed algorithm is O(N M^2) where
% M=param.NumAtlasPatches and N=param.NumSubjectPatches. However, M and N
% are usually very large. So the reduced complexity is O(MN D^2) where
% D=param.MinNumBasis. The dictionary size is MD^2, so choose M and D such that
% MD^2 <1E6. Otherwise, the code might be very slow.
%
% param.PatchSize = size of patches to be used. A reasonable value is
% [3 3 3] or [3 3 1]
%
% param.outname = output volume name, must be .xml or .nii. The output can
% be nii if and only if the input volumes are nifti, otherwise the output
% will be written as raw FLOAT32 files.
%
% param.subjectwmpeak or param.atlaswmpeak = Subject or Atlas white matter
% peak values. If not mentioned, then peaks are calculated based on KDE of
% the histograms. THIS should be made if atlas or subject is a BRAINWEB
% PHANTOM, because peak detection usually FAILS on BRAINWEB PHANTOMS.
%
% param.atlasoutlierratio and param.subjectoutlierratio are two 1x2 vectors
% containing the outlierratios. A typical value is [0.01 2], meaning
% 2% maximal values are going to be reduced from computation. It affects 
% the peak finding algorithm, so use these values carefully.
%
% param.AtlasSliceRange = Range of atlas slices to be chosen to construct
% atlas patch cloud. Default (=0) will choose slices from the middle of the
% image.
%
% param.SubjectSliceRange = Range of subject slices to be chosen to construct
% subject patch cloud. By default (=0), they are chosen to be middle 11 slices.
%
% param.isT2 = a flag to identify if the C2 contrast is T2-w image. Default
% is 0 (false). Always use this flag as 1 is the C2 contrast is T2-w to get
% better results.
%
% param.isinputT2 = a flag to identify if the C1 contrast is T2
%
% param.isCT = a flag to identify if the output contrast is CT image. Default
% is 0 (false). Always use this flag as 1 is the output contrast is CT to get
% better results.
%
% param.outlier = a flag to identify how much outlier is present in the image.
% Default is 'none'. There are 4 options, 'low', 'moderate','high','none';

% (For reference, version=12)

function ret=image_hallucination(subj,atlas,param)

% fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
% fprintf('------------------------WARNING WARNING--------------------\n');
% fprintf('If any of your input or atlas volume contains a Brainweb phantom \n');
% fprintf('check the WM peak calculated inside the code (will be displayed next).\n');
% fprintf('Most likely it would be wrong. So use the options param.subjectwmpeak\n');
% fprintf('or param.atlaswmpeak to feed the WM peak info.\n');
% fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
% fprintf('----------------------------------------------------------------\n');

disp(param);

if ~isfield(param,'iswrite')
    param.iswrite=1;
end
if ~isfield(param,'outname')
    s=(subj.vol{1});
    fprintf('param.outname is not mentioned. Using default.\n');
    param.outname=strcat(s(1:end-4),'_synth.nii');
end
fprintf('Output will be written in\n');
fprintf('%s\n',param.outname);
% try
NumCPU=param.NumCPU;
%     if NumCPU>feature('numCores')
%         NumCPU=feature('numCores');
%     end
% catch e
%     NumCPU=feature('numCores');
%     fprintf('Using default %d parallel processes. \n',NumCPU);
%     fprintf('If there are LESS than %d physical processors available, decrease the option param.NumCPU.\n',NumCPU);
% end
% fprintf('Using %d parallel processes. \n',NumCPU);
if isempty(gcp('nocreate'))
    pooldirname=tempname(fullfile(getenv('HOME'),'.matlab','local_cluster_jobs','R2022a'));
    mkdir(pooldirname);
    cluster=parallel.cluster.Local();
    cluster.NumWorkers=NumCPU;
    cluster.JobStorageLocation=pooldirname;
    fprintf('Temp Job directory = %s\n',pooldirname);
    pl=parpool(cluster,NumCPU);
    flag=1;
else
    flag=0;
end
% try
%     if matlabpool('size') < NumCPU
%         matlabpool close;
%         c = parcluster('local');
%         c.NumWorkers = NumCPU;
%         matlabpool open;
%         flag=1;
%     else
%         fprintf('Using existing %d matlabpools.\n',matlabpool('size'));
%         flag=0;
%     end
% catch e
%     try
%         matlabpool close;
%     catch e
%     end
%     c = parcluster('local');
%     c.NumWorkers = NumCPU;
%     matlabpool open;
% %     matlabpool(NumCPU);
%     flag=1;
% end

if ~isfield(param,'PatchSize')
    param.PatchSize=[3 3 3];
end
if ~isfield(param,'atlasoutlierratio')
    param.atlasoutlierratio = [0 2];
end
if ~isfield(param,'subjectoutlierratio')
    param.subjectoutlierratio = [0 2];
end

if ~isfield(param,'isskull')
    param.isskull=0;
end
if ~isfield(param,'NumCPU')
    param.NumCPU=8;
end
if ~isfield(param,'AtlasSliceRange')
    param.AtlasSliceRange=[];
end
if length(param.AtlasSliceRange)==1
        param.AtlasSliceRange=[param.AtlasSliceRange  param.AtlasSliceRange+20];
end
if ~isfield(param,'SubjectSliceRange')
    param.SubjectSliceRange=[];
end
if length(param.AtlasSliceRange)==1
        param.AtlasSliceRange=[param.AtlasSliceRange  param.AtlasSliceRange+20];
end
% if ~isfield(param,'isCT')   
%     param.isCT=0;
% end
% try
%     x=param.isinputT2;
% catch e
%     param.isinputT2=0;
% end
try
    x=param.isskull;
catch e
    param.isskull=0;
end
try
    x=param.outlier;
catch e
    param.outlier='moderate';
end
% try
%     x=param.isT2;
% catch e
%     param.isT2=0;
% end
if param.isCT
    imagethreshold=-1024;
    fprintf('Minimum intensity is set to -1024 for CT modality.\n');
else
    imagethreshold=0;
    fprintf('Minimum intensity is set to 0 for non-CT modality.\n');
end


if ~isfield(param,'NumAtlasPatches')
    param.NumAtlasPatches=5E4;
end
if ~isfield(param,'NumSubjectPatches')
    param.NumAtlasPatches=5E4;
end
if param.isskull
    if ~isfield(param,'subjectwmpeak')
        fprintf('Subject with skull must be associated with WM peaks.\n');
        fprintf('Enter WM peaks for subject in param.subjectwmpek as 1x%d vector.\n',T);
        ret=0;
        return;
    end
    if ~isfield(param,'atlaswmpeak')
        fprintf('Atlas with skull must be associated with WM peaks.\n');
        fprintf('Enter WM peaks for atlas in param.atlaswmpek as 1x%d vector.\n',T);
        ret=0;
        return;
    end
end
Ta=param.MinNumBasis;
if ~isfield(param,'outname')
    s=(subj.vol{1});
    fprintf('param.outname is not mentioned. Using default.\n');
    param.outname=strcat(s(1:end-4),'_synth.nii');
end
T=length(subj.vol);
if length(atlas.vol)~=T+1
    fprintf('Atlas has %d files and subject has %d files. \n',length(atlas.vol),length(subj.vol));
    fprintf('Atlas must have one more contrast than subject volumes.\n');
    ret=0;
    return;
end

atlas_M=cell(T+1,1);
for t=1:T+1
    % This is not a good way to check nii or nii.gz extension, what if file is
    % called 1.nii
    if strcmpi(atlas.vol{t}(end-2:end),'nii')  || strcmpi(atlas.vol{t}(end-5:end),'nii.gz')
        temp=load_untouch_nii(atlas.vol{t});
        atlas_M{t}=double(temp.img);
    elseif strcmpi(atlas.vol{t}(end-2:end),'xml')
        atlas_M{t}=ReadXml(atlas.vol{t});
    elseif isnumeric(atlas.vol{t})
        atlas_M{t}=atlas.vol{t};
    else
        fprintf('ERROR: Only NIFTI and XML files are supported for Atlas volumes file.\n');
        ret=0;
        return;
    end
end


if ~param.isCT
    [I, ~ ,~ ,~]=OutlierReduction(atlas_M{T+1},0,0,param.atlasoutlierratio);
    atlas_M{T+1}(atlas_M{T+1}>I(2))=I(2);
    atlas_M{T+1}(atlas_M{T+1}<I(1) & atlas_M{T+1}>0)=I(1);
end

temp=strsplit(param.inputmodal,',');
for t=1:T
    try
        pA=param.atlaswmpeak(t);
    catch e
        pA=find_WM_peak(atlas_M{t},temp{t});
    end
    [I, ~ ,~ ,~]=OutlierReduction(atlas_M{t},0,0,param.atlasoutlierratio);
%     temp=atlas_M{t};
%     temp(temp>I(2))=0;
%     temp(temp<I(1) & temp>0)=0;
%     
%         if pA(1)==0
%             delta=(I(2)-I(1))/80;
%             [pA,ff]=find_peaks(temp,[],[],0,2,delta);
%             if t>1
%                 indx=find(ff==max(ff));
%                 pA(end)=pA(indx(1));
%             end
%         end
%     catch e
%         pA=find_peaks(temp,[],[],0,2,(I(2)-I(1))/80);
%     end
    
    %     Outlier removal
    atlas_M{t}(atlas_M{t}>I(2))=I(2);
    atlas_M{t}(atlas_M{t}<I(1) & atlas_M{t}>0)=I(1);
    atlas_M{t}=atlas_M{t}/pA(1);
    fprintf('Normalizing the atlas by %.3f\n',pA(1));
    
%     if ~param.isinputT2
%         atlas_M{t}=atlas_M{t}/pA(end);
%         fprintf('Normalizing the atlas by %.3f\n',pA(end));
%     else
%         atlas_M{t}=atlas_M{t}/pA(1);
%         fprintf('Normalizing the atlas by %.3f\n',pA(1));
%     end
    
end
clear temp;
%---------------------------------------------------------------------------
% Low intensity removal
%---------------------------------------------------------------------------
if strcmpi(param.outlier,'low')
    outthresh=0.05;
elseif strcmpi(param.outlier,'moderate')
    outthresh=0.10;
elseif strcmpi(param.outlier,'high')
    outthresh=0.20;
elseif strcmpi(param.outlier,'none')
    outthresh=0;
end
fprintf('Using %s outlier removal.\n',param.outlier);
if ~strcmpi(param.outlier,'none')
    if param.isT2 && ~param.isskull
        len=0;
        for t=1:T
            for k=1:8
                q(1)=quantile(atlas_M{t}(atlas_M{t}>0),0.6/k);
                q(2)=quantile(atlas_M{T+1}(atlas_M{T+1}>0),0.12*k);
                indx=find(atlas_M{t}<q(1) & atlas_M{T+1}<q(2) & atlas_M{t}>0);

                if (len+length(indx))/sum(atlas_M{t}(:)>0)<=outthresh;
                    len=len+length(indx);
                    atlas_M{T+1}(indx)=0;
                    atlas_M{t}(indx)=0;
                else
                    break;
                end
            end
            fprintf('%.3f percent outliers removed for %d th atlas image.\n',100*len/sum(atlas_M{t}(:)>0),t);
        end
    else
        fprintf('Using MPRAGE parameters.\n');
    end

    if param.isCT && ~param.isskull
        for t=1:T
            q(1)=quantile(atlas_M{t}(atlas_M{t}>0),0.1);
            q(2)=quantile(atlas_M{T+1}(atlas_M{T+1}>0),0.15);
            indx=find(atlas_M{t}<q(1) & atlas_M{T+1}>q(2) & atlas_M{t}>0);
            atlas_M{T+1}(indx)=0;
            atlas_M{t}(indx)=0;
            len=length(indx);
            q(1)=quantile(atlas_M{t}(atlas_M{t}>0),0.75);
            q(2)=120;  % Max CT number of brain tissues
            indx=find(atlas_M{t}<q(1) & atlas_M{T+1}>q(2) & atlas_M{t}>0);
            atlas_M{T+1}(indx)=0;
            atlas_M{t}(indx)=0;
            len=len+length(indx);
            fprintf('%.3f percent outliers removed for CT image.\n',100*len/sum(atlas_M{t}(:)>0));
        end
    elseif param.isskull && param.isCT
        for t=1:T
            qq(1)=quantile(atlas_M{t}(atlas_M{t}>0),0.4); % 0.4 = Max possible fraction of CSF+GM in a normal subject
            qq(2)=120;  % Max CT number of brain tissues (CSF,WM,GM)

            indx=find(atlas_M{T+1}>qq(2) & atlas_M{t}>qq(1));
            atlas_M{T+1}(indx)=0;
            atlas_M{t}(indx)=0;
            len=length(indx);

            qq(1)=quantile(atlas_M{t}(atlas_M{t}>0),0.4); %
            qq(2)=120;  % Max CT number of soft tissue = 120
            indx=find(atlas_M{T+1}<qq(2) & atlas_M{t}<qq(1) & atlas_M{t}>0);
            atlas_M{T+1}(indx)=0;
            atlas_M{t}(indx)=0;
            len=len+length(indx);

            qq(2)=-350;  % Min CT number of FAT = -150
            indx=find(atlas_M{T+1}<qq(2) & atlas_M{t}>qq(1));
            atlas_M{T+1}(indx)=0;
            atlas_M{t}(indx)=0;
            len=len+length(indx);

            fprintf('%.3f percent outliers removed for CT image.\n',100*len/sum(atlas_M{t}(:)>0));
        end
    else
        fprintf('Using MPRAGE parameters.\n');
    end
end
%---------------------------------------------------------------------------

subj_M=cell(T,1);
for t=1:T
    % This is not a good way to check nii or nii.gz extension, what if file is
    % called 1.nii
    if strcmp(subj.vol{t}(end-2:end),'nii') || strcmpi(subj.vol{t}(end-5:end),'nii.gz')
        temp=load_untouch_nii(subj.vol{t});
        subj_M{t}=double(temp.img);
        pixd=temp.hdr.dime.pixdim(2:4);
    elseif strcmp(subj.vol{t}(end-2:end),'xml')
        [subj_M{t}, hd]=ReadXml(subj.vol{t});
        pixd=hd.res;
    elseif isnumeric(subj.vol{t})
        subj_M{t}=subj.vol{t};
        pixd=[1 1 1];
    else
        fprintf('ERROR: Only NIFTI and XML files are supported for Subject file.\n');
        ret=0;
        return;
    end
end

temp=strsplit(param.inputmodal,',');
for t=1:T
    try
        pS=param.subjectwmpeak(t);
    catch e
        pS=find_WM_peak(subj_M{t},temp{t});
    end
    [I, ~ ,~ ,~]=OutlierReduction(subj_M{t},0,0,param.subjectoutlierratio);
    fprintf('Using effective subject intensity range = [%.2f,%.2f]\n',I(1),I(2));
%     temp=subj_M{t};
%     temp(temp>I(2))=0;
%     temp(temp<I(1) & temp>0)=0;
%     try
%         pS=param.subjectwmpeak(t);
%         if pS(1)==0
%             delta=(I(2)-I(1))/80;
%             [pS,ff]=find_peaks(temp,[],[],0,2,delta);
%             if t>1
%                 indx=find(ff==max(ff));
%                 pS(end)=pS(indx(1));
%             end
%         end
%     catch e
%         pS=find_peaks(temp,[],[],0,2,(I(2)-I(1))/80);
%     end
    % Outlier removal
    subj_M{t}(subj_M{t}>I(2))=I(2);
    subj_M{t}(subj_M{t}<I(1) & subj_M{t}>0)=I(1);
    subj_M{t}=subj_M{t}/pS(end);
    fprintf('Normalizing the subject by %.3f\n',pS(1));
%     if ~param.isinputT2
%         subj_M{t}=subj_M{t}/pS(end);
%         fprintf('Normalizing the subject by %.3f\n',pS(end));
%     else
%         subj_M{t}=subj_M{t}/pS(1);
%         fprintf('Normalizing the subject by %.3f\n',pS(1));
%     end
end
clear temp;



% number of atlas patches to be considered

M=param.NumAtlasPatches;

psize=param.PatchSize;
% psize=[3 3 1];
dsize=floor(psize/2);
L=psize(1)*psize(2)*psize(3);
P=[floor(L/2)+1 floor(L/2)+1+L];

mask=ones(size(atlas_M{1}));
for t=2:T
    mask=mask.*double(atlas_M{t}>0);
end

for t=1:T
    atlas_M{t}=atlas_M{t}.*mask;
end

clear mask;

dim=size(atlas_M{1});
Atlas=zeros(L*T,M);
AtlasM2=zeros(L,M);
% AtlasIndxList=zeros(3,M);
try
    if param.AtlasSliceRange(1)==0 || param.SubjectSliceRange(2)==0
        param.AtlasSliceRange=[];
    end
end

fprintf('Generating atlas dictionary..\n');

atlasslice=param.AtlasSliceRange;
if length(atlasslice)==1
        atlasslice=[atlasslice atlasslice];
end


if isempty(atlasslice) || atlasslice(1)==0
    [~,cropParams]= Crop3D(atlas_M{1},0,0,0);
    a=cropParams.boundingBox;
    atlasslice=a(5)+floor(0.65*(a(6)-a(5)));
    atlasslice=[atlasslice-10  atlasslice+10];
    atlasslice(1)=max(1,atlasslice(1));
    atlasslice(2)=min(size(atlas_M{1},3),atlasslice(2));
end
atlasslice=[max(1,atlasslice(1)) min(atlasslice(2),size(atlas_M{1},3))];
fprintf('Atlas slices %d to %d will be used.\n', atlasslice(1),atlasslice(2));
count=1;

for k=atlasslice(1):atlasslice(2)
    for i=dsize(1)+1:dim(1)-dsize(1)
        for j=dsize(2)+1:dim(2)-dsize(2)
            
            m=cell(T,1);
            if atlas_M{1}(i,j,k)>0
                for t=1:T            
                    m{t}=reshape(atlas_M{t}(i-dsize(1):i+dsize(1),...
                        j-dsize(2):j+dsize(2),k-dsize(3):k+dsize(3)),[L 1]);
                end
%                 m1=reshape(atlas_M2(i-dsize(1):i+dsize(1),j-dsize(2):j+dsize(2),k-dsize(3):k+dsize(3)),[L 1]);
                
                for t=1:T
                    Atlas((t-1)*L+1:t*L,count)=m{t};
                end
                AtlasM2(:,count)=reshape(atlas_M{T+1}(i-dsize(1):i+dsize(1),...
                    j-dsize(2):j+dsize(2),k-dsize(3):k+dsize(3)),[L 1]);
%                 AtlasM2(1:L,count)=m1;
                count=count+1;
            end
%             Atlas(L+1:2*L,count)=T;
            if count>M
                break;
            end
        end
        if count>M
            break;
        end
    end
    if count>M
        break;
    end
end


if count<M
    Atlas=Atlas(:,1:count-1);
    AtlasM2=AtlasM2(:,1:count-1);
end
fprintf('Atlas slices [%d,%d] have been used to generate atlas patch dictionary.\n',...
    atlasslice(1),k);


N=param.NumSubjectPatches;

% SubjIndxList=ones(3,N);
Subj=zeros(L*T,N);
fprintf('Generating subject dictionary..\n');
subslice=param.SubjectSliceRange;


if isempty(subslice) || subslice(1)==0
    [~,cropParams]= Crop3D(subj_M{1},0,0,0);
    a=cropParams.boundingBox;
    subslice=a(5)+floor(0.65*(a(6)-a(5)));  % heuristics to get bounding box
    subslice=[subslice-10  subslice+10];
    subslice(1)=max(1,subslice(1));
    subslice(2)=min(size(subj_M{1},3),subslice(2));
end
count=1;
subslice=[max(1,subslice(1)) min(subslice(2),size(subj_M{1},3))];
fprintf('Subject slices %d to %d will be used to generate subject patch cloud.\n',...
    subslice(1),subslice(2));
for k=subslice(1):subslice(2)
    for i=dsize(1)+1:dim(1)-dsize(1)
        for j=dsize(2)+1:dim(2)-dsize(2)
            m=cell(T,1);
            if subj_M{1}(i,j,k)>0
                for t=1:T
                    m{t}=reshape(subj_M{t}(i-dsize(1):i+dsize(1),...
                        j-dsize(2):j+dsize(2),k-dsize(3):k+dsize(3)),[L 1]);
                end
                for t=1:T
                    Subj((t-1)*L+1:t*L,count)=m{t};                    
                end
                count=count+1;
            end
            if count>N
                break;
            end
        end
        if count>N
            break;
        end
    end
    if count>N
        break;
    end
end
if count<N
    Subj=Subj(:,1:count-1);
    %     SubjIndxList=SubjIndxList(:,1:count-1);
end
fprintf('Subject slices [%d,%d] have been used to generate subject patch dictionary.\n',subslice(1),k);

% initial estimates of W and S1, S2
N=size(Subj,2);
M=size(Atlas,2);
Nt=3*N;
% Nt is the new size of the atlas dictionary, containing T closest patch-pairs
% for each subject patch
dim=size(subj_M{1});

tempSubjM2=zeros(L,N);
SubjM2=zeros(L,N);   % synthetic M2 contrast of Subj


fprintf('Generating atlas kdtree..\n');
% tic
treeA = createns(Atlas','distance', 'cityblock');
% toc

fprintf('Finding nearest neighbors on the atlas patch cloud..\n');
tic
tempIIB=cell(NumCPU,1);
INDXLISTS1=zeros(N,3);
Del=(N-mod(N,NumCPU))/NumCPU;
parfor count=1:NumCPU
    tempIIB{count}=knnsearch(treeA,Subj(:,Del*(count-1)+1:Del*count)',...
        'distance', 'cityblock','K',3);
end

count=1;
for i=1:NumCPU
    for j=1:length(tempIIB{i})
        INDXLISTS1(count,:)=tempIIB{i}(j,:);
        count=count+1;
    end
end
for i=N-mod(N,NumCPU)+1:N
    INDXLISTS1(i,:)=knnsearch(treeA,Subj(:,i)','distance','cityblock','K',3);
end
clear tempIIB treeA;
% INDXLISTS1 = knnsearch(treeA,Subj','dist','euclidean','K',3);
toc
fprintf('Generating a new pairwise Atlas..\n');

newAtlasM1_1=zeros(L*T,Nt);  % M1 contrast of 1st part of the atlas
newAtlasM1_2=zeros(L*T,Nt);  % M1 contrast of 2nd part of the atlas
newAtlasM2_1=zeros(L,Nt);
newAtlasM2_2=zeros(L,Nt);

for i=1:N
    %     newAtlasM1_1(:,2*i-1)=Atlas(:,INDXLISTS1(i,1));
    %     newAtlasM1_2(:,2*i-1)=Atlas(:,INDXLISTS1(i,2));
    %     newAtlasM1_1(:,2*i)=Atlas(:,INDXLISTS1(i,3));
    %     newAtlasM1_2(:,2*i)=Atlas(:,INDXLISTS1(i,4));
    
    newAtlasM1_1(:,3*i-2)=Atlas(:,INDXLISTS1(i,1));
    newAtlasM1_1(:,3*i-1)=Atlas(:,INDXLISTS1(i,2));
    newAtlasM1_1(:,3*i)=Atlas(:,INDXLISTS1(i,3));
    newAtlasM1_2(:,3*i-2)=Atlas(:,INDXLISTS1(i,2));
    newAtlasM1_2(:,3*i-1)=Atlas(:,INDXLISTS1(i,3));
    newAtlasM1_2(:,3*i)=Atlas(:,INDXLISTS1(i,1));
    
    %     newAtlasM2_1(:,2*i-1)=AtlasM2(:,INDXLISTS1(i,1));
    %     newAtlasM2_2(:,2*i-1)=AtlasM2(:,INDXLISTS1(i,2));
    %     newAtlasM2_1(:,2*i)=AtlasM2(:,INDXLISTS1(i,3));
    %     newAtlasM2_2(:,2*i)=AtlasM2(:,INDXLISTS1(i,4));
    
    newAtlasM2_1(:,3*i-2)=AtlasM2(:,INDXLISTS1(i,1));
    newAtlasM2_1(:,3*i-1)=AtlasM2(:,INDXLISTS1(i,2));
    newAtlasM2_1(:,3*i)=AtlasM2(:,INDXLISTS1(i,3));
    newAtlasM2_2(:,3*i-2)=AtlasM2(:,INDXLISTS1(i,2));
    newAtlasM2_2(:,3*i-1)=AtlasM2(:,INDXLISTS1(i,3));
    newAtlasM2_2(:,3*i)=AtlasM2(:,INDXLISTS1(i,1));
end
% Now checking which one of these pairs are relevant for i^th subject patch

fprintf('Generating modified atlas kdtree..\n');
% tic
tree=KDTreeSearcher(newAtlasM1_1','dist','cityblock');

fprintf('Finding nearest neighbors on the modified atlas patch cloud..\n');
tic
tempIIB=cell(NumCPU,1);
INDXLISTS=zeros(N,Ta);
Del=floor((N-mod(N,NumCPU))/NumCPU);
% fprintf('Parallelizing processes...\n');
parfor count=1:NumCPU
    tempIIB{count}=knnsearch(tree,Subj(:,Del*(count-1)+1:Del*count)',...
        'dist','cityblock','K',Ta);
end
% fprintf('Back to single cpu...\n');
count=1;
for i=1:NumCPU
    for j=1:length(tempIIB{i})
        INDXLISTS(count,:)=tempIIB{i}(j,:);
        count=count+1;
    end
end
for i=N-mod(N,NumCPU)+1:N
    INDXLISTS(i,:)=knnsearch(tree,Subj(:,i)','distance','cityblock','K',Ta);
end
clear tempIIB;
toc
% INDXLISTS is of size N x Ta, containing numbers from [1,Nt], Nt=3N

% Inititalization of Sigmas
S1=zeros(1,Nt);
% S2=zeros(1,Nt);
R=cell(Nt,1);
C=cell(Nt,1);
fprintf('Generating subject-pair kdtree..\n');
% tic
ns1 = KDTreeSearcher(INDXLISTS(:));
% ns2 = KDTreeSearcher(tmpINDXLISTS2(:));
% toc
fprintf('Finding nearest neighbor pairs on the subject patch clouds..\n');
pp=zeros(1,Nt);
indx1=rangesearch(ns1,[1:Nt]',0);
% indx2=rangesearch(ns2,[1:Nt]',0);
fprintf('Parallelizing processes...\n');
parfor j=1:Nt-mod(Nt,NumCPU)
    [R{j}, C{j}]=ind2sub([N Ta],indx1{j});
    pp(j)=length(R{j});
end
fprintf('Back to single cpu...\n');
for j=Nt-mod(Nt,NumCPU)+1:Nt
    [R{j}, C{j}]=ind2sub([N Ta],indx1{j});
    pp(j)=length(R{j});
end
fprintf('WARNING: %.2f percent of the Sigma''s are going to be null.\n',100*length(find(pp==0))/Nt);
if length(find(pp==0))>Nt/10
    fprintf('!!!!Try increasing param.MinNumBasis or use more subject patches!!!!\n');
    fprintf('If the amount of null sigma is still large, the initial normalization\n');
    fprintf('between the atlas and subject C1 contrast volumes could be bad. \n');    
    fprintf('MAKE SURE THAT THE WM PEAKS OF THE SUBJECT AND THE ATLAS ARE CORRECT.\n');
end

% remove redundant stuff from R and C
fprintf('Removing some redundant stuff.\n');
fprintf('Number of redundant sigmas = %d.\n',length(find(pp>2*param.MinNumBasis)));
indx=find(pp>2*param.MinNumBasis);
for i=1:length(indx)
    j=indx(i);
    R{j}=R{j}(1:param.MinNumBasis);
    C{j}=C{j}(1:param.MinNumBasis);
end
clear pp ns1;

% Initializing the mixing coefficient Alpha and posteriors W
Alpha=0.5*ones(N,Ta);
W=ones(N,Ta)/(Ta);
clear temp* indx1 ;
% Initializing Sigma1
% matlabpool(NumCPU);
fprintf('Initializing Sigma1..\n');
% fprintf('Parallelizing processes...\n');
tic
parfor j=1:Nt-mod(Nt,NumCPU)
    temp=zeros(L*T,length(R{j}));
    for i=1:length(R{j})
        temp(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM1_1(:,j)+...
            (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM1_2(:,j);
    end
    temp=Subj(:,R{j})-temp;
    S1(j)=sum(sum(temp.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
    
end
% fprintf('Back to single cpu...\n');
for j=Nt-mod(Nt,NumCPU)+1:Nt
    temp=zeros(L*T,length(R{j}));
    for i=1:length(R{j})
        temp(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM1_1(:,j)+...
            (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM1_2(:,j);
    end
    temp=Subj(:,R{j})-temp;
    S1(j)=sum(sum(temp.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
    
end

S1=sqrt(S1);
S1(S1==0)=1E-4;
S1(isnan(S1)==1)=1E-4;


% Intialization of W
% fprintf('Initializing W..\n');
% fprintf('Parallelizing processes...\n');
parfor i=1:N-mod(N,NumCPU)
    indx=INDXLISTS(i,:);
    temp1=bsxfun(@times,Alpha(i,:),newAtlasM1_1(:,indx));
    temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM1_2(:,indx));
    temp=bsxfun(@minus,temp1+temp2,Subj(:,i));
    temp=sum(temp.^2,1);
    temp=(1./(S1(indx))).*exp(-temp./(2*S1(indx).^2));
    temp(isnan(temp)==1)=0;
    W(i,:)=temp/sum(temp);
end
% fprintf('Back to single cpu...\n');
for i=N-mod(N,NumCPU)+1:N
    indx=INDXLISTS(i,:);
    temp1=bsxfun(@times,Alpha(i,:),newAtlasM1_1(:,indx));
    temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM1_2(:,indx));
    temp=bsxfun(@minus,temp1+temp2,Subj(:,i));
    temp=sum(temp.^2,1);
    temp=(1./(S1(indx))).*exp(-temp./(2*S1(indx).^2));
    temp(isnan(temp)==1)=0;
    W(i,:)=temp/sum(temp);
end

% Initialization of SubjM2
% fprintf('Parallelizing processes...\n');
S2=ones(1,Nt);
parfor i=1:N-mod(N,NumCPU)
    indx=INDXLISTS(i,:);
    temp1=bsxfun(@times,Alpha(i,:),newAtlasM2_1(:,indx));
    temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM2_2(:,indx));
    temp=temp1+temp2;
    temp=bsxfun(@times,temp,W(i,:)./S2(indx).^2);
    SubjM2(:,i)=sum(temp,2)/sum(W(i,:)./S2(indx).^2,2);
end
% fprintf('Back to single cpu...\n');
for i=N-mod(N,NumCPU)+1:N
    indx=INDXLISTS(i,:);
    temp1=bsxfun(@times,Alpha(i,:),newAtlasM2_1(:,indx));
    temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM2_2(:,indx));
    temp=temp1+temp2;
    temp=bsxfun(@times,temp,W(i,:));
    SubjM2(:,i)=sum(temp,2)/sum(W(i,:));
    %     temp=bsxfun(@times,temp,W(i,:)./S2(indx).^2);
    %     SubjM2(:,i)=sum(temp,2)/sum(W(i,:)./S2(indx).^2,2);
end

% Initialization of Sigma2
fprintf('Initializing Sigma2..\n');
% fprintf('Parallelizing processes...\n');
parfor j=1:Nt-mod(Nt,NumCPU)
    temp=zeros(L,length(R{j}));
    for i=1:length(R{j})
        temp(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM2_1(:,j)+...
            (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM2_2(:,j);
    end
    temp=SubjM2(:,R{j})-temp;
    S2(j)=sum(sum(temp.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
end
% fprintf('Back to single cpu...\n');
for j=Nt-mod(Nt,NumCPU)+1:N
    temp=zeros(L,length(R{j}));
    for i=1:length(R{j})
        temp(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM2_1(:,j)+...
            (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM2_2(:,j);
    end
    temp=SubjM2(:,R{j})-temp;
    S2(j)=sum(sum(temp.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
end
S2=sqrt(S2);
S2(S2==0)=min(S2(S2>0));
S2(isnan(S2)==1)=min(S2(S2>0));

% Initialization of Alpha, Alpha=0.5*ones(N,Ta^2);
% F=calc_alpha(x,atlasm1,atlasm2,subjm1,subjm2,sigma1,sigma2)
fprintf('Initializing alpha..\n');
% fprintf('Parallelizing processes...\n');
parfor i=1:N-mod(N,NumCPU)
    indx=INDXLISTS(i,:);
    subjm1=Subj(:,i);
    subjm2=SubjM2(:,i);
    %     tic
    for j=1:Ta
        atlsm1=[newAtlasM1_1(:,indx(j)) newAtlasM1_2(:,indx(j))];
        atlsm2=[newAtlasM2_1(:,indx(j)) newAtlasM2_2(:,indx(j))];
        sig1=S1(indx(j));
        sig2=S2(indx(j));
        Alpha(i,j)=calc_alpha(atlsm1,atlsm2,subjm1,subjm2,sig1,sig2);
        %         toc
    end
    %     toc
end
% fprintf('Back to single cpu...\n');
for i=N-mod(N,NumCPU)+1:N
    indx=INDXLISTS(i,:);
    subjm1=Subj(:,i);
    subjm2=SubjM2(:,i);
    %     tic
    for j=1:Ta
        atlsm1=[newAtlasM1_1(:,indx(j)) newAtlasM1_2(:,indx(j))];
        atlsm2=[newAtlasM2_1(:,indx(j)) newAtlasM2_2(:,indx(j))];
        sig1=S1(indx(j));
        sig2=S2(indx(j));
        Alpha(i,j)=calc_alpha(atlsm1,atlsm2,subjm1,subjm2,sig1,sig2);
    end
end
TTT=quantile(atlas_M{T+1}(:),0.99)/1000;
% END Alpha calculation
clear atlas_M;
err2=1E10;
iter=1;
maxiter=10;
% maxiter=3;
err=zeros(maxiter,2);

fprintf('Initial sigma and W found, starting iterations..\n');
fprintf('Convergence threshold set at %.4f.\n',TTT);
toc
tic
subj_M2=zeros(dim);
while err2>TTT && iter<=maxiter
    fprintf('Starting iteration %d : \n',iter);
    %     tic
    S1a=zeros(1,Nt);
    S2a=zeros(1,Nt);
    
    % iteration for W
    fprintf('   Finding W..\n');
    parfor i=1:N-mod(N,NumCPU)
        indx=INDXLISTS(i,:);
        temp1=bsxfun(@times,Alpha(i,:),newAtlasM1_1(:,indx));
        temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM1_2(:,indx));
        temp3=bsxfun(@times,Alpha(i,:),newAtlasM2_1(:,indx));
        temp4=bsxfun(@times,1-Alpha(i,:),newAtlasM2_2(:,indx));
        
        tempx=bsxfun(@minus,temp1+temp2,Subj(:,i));
        tempu=bsxfun(@minus,temp3+temp4,SubjM2(:,i));
        tempx=sum(tempx.^2,1);
        tempu=sum(tempu.^2,1);
        
        tempa=Alpha(i,:).*(1-Alpha(i,:));
        % May be multiply with tempa ?
        W(i,:)=tempa.*(1./(S1(indx))).*(1./(S2(indx))).*...
            exp(-tempx./(2*S1(indx).^2)).*exp(-tempu./(2*S2(indx).^2));
        W(i,:)=W(i,:)/sum(W(i,:));
    end
    for i=N-mod(N,NumCPU)+1:N
        indx=INDXLISTS(i,:);
        temp1=bsxfun(@times,Alpha(i,:),newAtlasM1_1(:,indx));
        temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM1_2(:,indx));
        temp3=bsxfun(@times,Alpha(i,:),newAtlasM2_1(:,indx));
        temp4=bsxfun(@times,1-Alpha(i,:),newAtlasM2_2(:,indx));
        
        tempx=bsxfun(@minus,temp1+temp2,Subj(:,i));
        tempu=bsxfun(@minus,temp3+temp4,SubjM2(:,i));
        tempx=sum(tempx.^2,1);
        tempu=sum(tempu.^2,1);
        
        tempa=Alpha(i,:).*(1-Alpha(i,:));
        % May be multiply with tempa ?
        W(i,:)=tempa.*(1./(S1(indx))).*(1./(S2(indx))).*...
            exp(-tempx./(2*S1(indx).^2)).*exp(-tempu./(2*S2(indx).^2));
        W(i,:)=W(i,:)/sum(W(i,:));
    end
    %     fprintf('NaNs in W= %d\n',sum(isnan(W(:))));
    W(isnan(W)==1)=0;
    
    % iteration for sigma1 and sigma2
    fprintf('   Finding sigma''s..\n');
    parfor j=1:Nt-mod(Nt,NumCPU)
        tempx=zeros(L*T,length(R{j}));
        tempu=zeros(L,length(R{j}));
        for i=1:length(R{j})
            tempx(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM1_1(:,j)+...
                (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM1_2(:,j);
            tempu(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM2_1(:,j)+...
                (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM2_2(:,j);
        end
        tempx=Subj(:,R{j})-tempx;
        tempu=SubjM2(:,R{j})-tempu;
        
        S1a(j)=sum(sum(tempx.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
        S2a(j)=sum(sum(tempu.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
        
    end
    for j=Nt-mod(Nt,NumCPU)+1:Nt
        tempx=zeros(L*T,length(R{j}));
        tempu=zeros(L,length(R{j}));
        for i=1:length(R{j})
            tempx(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM1_1(:,j)+...
                (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM1_2(:,j);
            tempu(:,i)=Alpha(R{j}(i),C{j}(i))*newAtlasM2_1(:,j)+...
                (1-Alpha(R{j}(i),C{j}(i)))*newAtlasM2_2(:,j);
        end
        tempx=Subj(:,R{j})-tempx;
        tempu=SubjM2(:,R{j})-tempu;
        
        S1a(j)=sum(sum(tempx.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
        S2a(j)=sum(sum(tempu.^2,1)*diag(W(R{j},C{j})))/sum(diag(W(R{j},C{j})));
        
    end
    S1a(S1a==0)=1E-4;S1a(isnan(S1a)==1)=1E-4;
    S2a(S2a==0)=1E-4;S2a(isnan(S2a)==1)=1E-4;
    %     S2a(S2a==0)=min(S2(S2>0));S2a(isnan(S2a)==1)=min(S2(S2>0));
    
    S1a=sqrt(S1a);
    S2a=sqrt(S2a);
    
    % Iterating over alpha in a parallel fashion
    fprintf('   Finding alpha..\n');
    
    %     matlabpool(NumCPU);
    parfor i=1:N-mod(N,NumCPU)
        indx=INDXLISTS(i,:);
        subjm1=Subj(:,i);
        subjm2=SubjM2(:,i);
        
        for j=1:Ta
            atlsm1=[newAtlasM1_1(:,indx(j)) newAtlasM1_2(:,indx(j))];
            atlsm2=[newAtlasM2_1(:,indx(j)) newAtlasM2_2(:,indx(j))];
            sig1=S1(indx(j));
            sig2=S2(indx(j));
            try
                Alpha(i,j)=calc_alpha(atlsm1,atlsm2,subjm1,subjm2,sig1,sig2);
            catch e
                fprintf('%s\n',e.message);
            end
        end
    end
    %     matlabpool close;
    for i=N-mod(N,NumCPU)+1:N
        indx=INDXLISTS(i,:);
        subjm1=Subj(:,i);
        subjm2=SubjM2(:,i);
        
        for j=1:Ta
            atlsm1=[newAtlasM1_1(:,indx(j)) newAtlasM1_2(:,indx(j))];
            atlsm2=[newAtlasM2_1(:,indx(j)) newAtlasM2_2(:,indx(j))];
            sig1=S1(indx(j));
            sig2=S2(indx(j));
            try
                Alpha(i,j)=calc_alpha(atlsm1,atlsm2,subjm1,subjm2,sig1,sig2);
            catch e
                fprintf('%s\n',e.message);
            end
        end
    end
    
    % find alternate contrast SubjM2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('   Finding subject M2..\n');
    parfor i=1:N-mod(N,NumCPU)
        indx=INDXLISTS(i,:);
        temp1=bsxfun(@times,Alpha(i,:),newAtlasM2_1(:,indx));
        temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM2_2(:,indx));
        temp=temp1+temp2;
        temp=bsxfun(@times,temp,W(i,:));
        tempSubjM2(:,i)=sum(temp,2)/sum(W(i,:));
        %         temp=bsxfun(@times,temp,W(i,:)./S2(indx).^2);
        %         tempSubjM2(:,i)=sum(temp,2)/sum(W(i,:)./S2(indx).^2,2);
        
    end
    for i=N-mod(N,NumCPU)+1:N
        indx=INDXLISTS(i,:);
        temp1=bsxfun(@times,Alpha(i,:),newAtlasM2_1(:,indx));
        temp2=bsxfun(@times,1-Alpha(i,:),newAtlasM2_2(:,indx));
        temp=temp1+temp2;
        temp=bsxfun(@times,temp,W(i,:));
        tempSubjM2(:,i)=sum(temp,2)/sum(W(i,:));
        %         temp=bsxfun(@times,temp,W(i,:)./S2(indx).^2);
        %         tempSubjM2(:,i)=sum(temp,2)/sum(W(i,:)./S2(indx).^2,2);
        
    end
    tempSubjM2(isnan(tempSubjM2)==1)=0;
%     err(iter,1)=mean(abs(S1-S1a));
%     err(iter,2)=mean(abs(S2-S2a));
    err(iter)=mean(abs(SubjM2(:)-tempSubjM2(:)));
    %     err(iter,2)=median(abs(S1-S1a));
    %     err(iter,3)=max(abs(S1-S1a));
    S1=S1a;
    S2=S2a;
    S1(isnan(S1)==1)=0;
    S2(isnan(S2)==1)=0;
    %     S1(S1<max(S1)/1000)=0;
    %     S2(S2<max(S2)/1000)=0;
    err2=err(iter);
    SubjM2=tempSubjM2;
    %     if sum(isnan(S2))>0
    %     end
    fprintf('Iter = %d, error = %.6f\n',iter,err(iter));
%     fprintf('Iter=%d,mean err 1=%.6f, mean err 3=%.6f\n',iter,err(iter,1),err(iter,3));
    iter=iter+1;
    %     toc
end
if iter<5
    fprintf('EM converged in %d iterations. Hmm. Threshold may be too high.\n',iter-1);
end
if iter==maxiter+1
    fprintf('EM converged in %d iterations. Hmm. Threshold may be too low.\n',iter-1);
end
clear temp*;

toc
fprintf('Iterations on Sigma and W finished. Generating synthesized volume..\n');
% errcount=0;
fprintf('Boundaries are forcefully made ZERO\n');
for t=1:T
subj_M{t}(1:dsize(1),:,:)=0;
subj_M{t}(:,1:dsize(2),:)=0;
subj_M{t}(:,:,1:dsize(3))=0;
subj_M{t}(dim(1)-dsize(1):dim(1),:,:)=0;
subj_M{t}(:,dim(2)-dsize(2):dim(2),:)=0;
subj_M{t}(:,:,dim(3)-dsize(3)+1:dim(3))=0;
end
tic
% for k=70:90
for k=dsize(3)+1:dim(3)-dsize(3)
    
    NumVox=sum(sum(subj_M{1}(dsize(1)+1:dim(1)-dsize(1),dsize(2)+1:dim(2)-dsize(2),k)>0));
    ImageIndex=zeros(1,dim(1)*dim(2));
    fprintf('Working on slice %d, number of voxels = %d...',k,NumVox);
    subpatch=zeros(3+L*T,NumVox);
    count=1;
    count2=1;
    for i=1:dim(1)
        for j=1:dim(2)
            if subj_M{1}(i,j,k)>0 && i>dsize(1) && i<=dim(1)-dsize(1) ...
                    && j>dsize(2) && j<=dim(2)-dsize(2)
                try
                    for t=1:T
                    subpatch((t-1)*L+1:t*L,count)=reshape(...
                        subj_M{t}(i-dsize(1):i+dsize(1),...
                        j-dsize(2):j+dsize(2),k-dsize(3):k+dsize(3)),[L 1]);
                    end
                catch e
                    fprintf('Probably boundary voxels. Ignore.\n');
                end
                subpatch(T*L+1:T*L+3,count)=[i j k]';
                count=count+1;
                ImageIndex(count2)=1;
            end
            count2=count2+1;
        end
    end
    
    tempIIB=cell(NumCPU,1);
    II=zeros(NumVox,Ta);
    Del=(NumVox-mod(NumVox,NumCPU))/NumCPU;
    parfor count=1:NumCPU
        tempIIB{count}=knnsearch(tree,subpatch(1:T*L,Del*(count-1)+1:Del*count)',...
            'dist','cityblock','K',Ta);
    end
    
    count=1;
    for i=1:NumCPU
        for j=1:size(tempIIB{i},1)
            II(count,:)=tempIIB{i}(j,:);
            count=count+1;
        end
    end
    for i=NumVox-mod(NumVox,NumCPU)+1:NumVox
        II(i,:)=knnsearch(tree,subpatch(1:T*L,i)','dist','cityblock','K',Ta);
    end
    clear tempIIB;
    
    
    
    subpatch2=zeros(1,NumVox);
    parfor i=1:NumVox-mod(NumVox,NumCPU)
        z=subpatch(1:T*L,i);
        
        indx=II(i,:);
        
        % For each subject patch, alpha and W are iterated, while
        % sigma's are kept from KDTree
        alpha=0.5*ones(1,Ta);
        iter=1;err=100*ones(1,15);
        err2=100;
        w=ones(1,Ta)/Ta;
        temp1=bsxfun(@times,alpha,newAtlasM2_1(:,indx));
        temp2=bsxfun(@times,1-alpha,newAtlasM2_2(:,indx));
        temp=temp1+temp2;
        temp=bsxfun(@times,temp,w./S2(indx).^2);
        z2=sum(temp,2)/sum(w./S2(indx).^2,2);
        while err2>0.01 && iter<6
            % iteration over alpha
            for t=1:Ta
                atlsm1=[newAtlasM1_1(:,indx(t)) newAtlasM1_2(:,indx(t))];
                atlsm2=[newAtlasM2_1(:,indx(t)) newAtlasM2_2(:,indx(t))];
                sig1=S1(indx(t));
                sig2=S2(indx(t));
                try
                    alpha(t)=calc_alpha(atlsm1,atlsm2,z,z2,sig1,sig2);
                catch e
                    fprintf('err\n');
                end
            end
            
            % iteration over w
            temp1=bsxfun(@times,alpha,newAtlasM1_1(:,indx));
            temp2=bsxfun(@times,1-alpha,newAtlasM1_2(:,indx));
            temp3=bsxfun(@times,alpha,newAtlasM2_1(:,indx));
            temp4=bsxfun(@times,1-alpha,newAtlasM2_2(:,indx));
            tempx=bsxfun(@minus,temp1+temp2,z);
            tempu=bsxfun(@minus,temp3+temp4,z2);
            tempx=sum(tempx.^2,1);
            tempu=sum(tempu.^2,1);
            tempa=alpha.*(1-alpha);
            w=tempa.*(1./(S1(indx))).*(1./(S2(indx))).*...
                exp(-tempx./(2*S1(indx).^2)).*exp(-tempu./(2*S2(indx).^2));
            
            w(isnan(w)==1)=0;
            w=w/sum(w);
            
            
            % iteration over z2
            temp1=bsxfun(@times,alpha,newAtlasM2_1(:,indx));
            temp2=bsxfun(@times,1-alpha,newAtlasM2_2(:,indx));
            temp=temp1+temp2;
            temp=bsxfun(@times,temp,w);
            z3=sum(temp,2)/sum(w,2);
            %             temp=bsxfun(@times,temp,w./S2(indx).^2);
            %             z3=sum(temp,2)/sum(w./S2(indx).^2,2);
            z3(isnan(z3)==1)=0;
            % convergence checking
            err(iter)=mean(abs(z2-z3));
            err2=err(iter);
            iter=iter+1;
            z2=z3;
        end
        subpatch2(i)=z2(P(1));
        
    end
    
    for i=NumVox-mod(NumVox,NumCPU)+1:NumVox
        z=subpatch(1:T*L,i);
        
        indx=II(i,:);
        
        % For each subject patch, alpha and W are iterated, while
        % sigma's are kept from KDTree
        alpha=0.5*ones(1,Ta);
        iter=1;err=1000;
        w=ones(1,Ta)/Ta;
        temp1=bsxfun(@times,alpha,newAtlasM2_1(:,indx));
        temp2=bsxfun(@times,1-alpha,newAtlasM2_2(:,indx));
        temp=temp1+temp2;
        temp=bsxfun(@times,temp,w./S2(indx).^2);
        z2=sum(temp,2)/sum(w./S2(indx).^2,2);
        while err>0.01 && iter<15
            % iteration over alpha
            for t=1:Ta
                atlsm1=[newAtlasM1_1(:,indx(t)) newAtlasM1_2(:,indx(t))];
                atlsm2=[newAtlasM2_1(:,indx(t)) newAtlasM2_2(:,indx(t))];
                sig1=S1(indx(t));
                sig2=S2(indx(t));
                try
                    alpha(t)=calc_alpha(atlsm1,atlsm2,z,z2,sig1,sig2);
                catch e
                    fprintf('err\n');
                end
            end
            
            % iteration over w
            temp1=bsxfun(@times,alpha,newAtlasM1_1(:,indx));
            temp2=bsxfun(@times,1-alpha,newAtlasM1_2(:,indx));
            temp3=bsxfun(@times,alpha,newAtlasM2_1(:,indx));
            temp4=bsxfun(@times,1-alpha,newAtlasM2_2(:,indx));
            tempx=bsxfun(@minus,temp1+temp2,z);
            tempu=bsxfun(@minus,temp3+temp4,z2);
            tempx=sum(tempx.^2,1);
            tempu=sum(tempu.^2,1);
            tempa=alpha.*(1-alpha);
            w=tempa.*(1./(S1(indx))).*(1./(S2(indx))).*...
                exp(-tempx./(2*S1(indx).^2)).*exp(-tempu./(2*S2(indx).^2));
            
            
            w=w/sum(w);
            w(isnan(w)==1)=0;
            
            % iteration over z2
            temp1=bsxfun(@times,alpha,newAtlasM2_1(:,indx));
            temp2=bsxfun(@times,1-alpha,newAtlasM2_2(:,indx));
            temp=temp1+temp2;
            temp=bsxfun(@times,temp,w);
            z3=sum(temp,2)/sum(w,2);
            %             temp=bsxfun(@times,temp,w./S2(indx).^2);
            %             z3=sum(temp,2)/sum(w./S2(indx).^2,2);
            z3(isnan(z3)==1)=z2(isnan(z3)==1);
            % convergence checking
            err=mean(abs(z2-z3));
            iter=iter+1;
            z2=z3;
        end
        
        subpatch2(i)=z2(P(1));
        
    end
    subpatch2=real(subpatch2);
    temp=zeros(1,dim(1)*dim(2));
    temp(ImageIndex==1)=subpatch2;
    subj_M2(:,:,k)=reshape(temp,[dim(2) dim(1)])';
    clear temp;
    
    fprintf('..Finished. \n');
    
end
toc
if flag
   delete(pl);
   rmdir(pooldirname,'s');
end
if param.isCT
    h1=fspecial3('gaussian',[3 3 3]);
    subj_M2=convn(subj_M2,h1,'same');
end
    
subj_M2=subj_M2.*(subj_M{1}>0);
subj_M2(subj_M2<imagethreshold)=imagethreshold;
if param.isCT
    subj_M2(subj_M{1}==0)=-1024;
end
    

fprintf('Writing output image volume %s\n');
try
    fprintf('Writing output volume\n');
    fprintf('%s\n',param.outname);
    if strcmpi(param.outname(end-2:end),'nii') 
        
        if strcmpi(subj.vol{1}(end-2:end),'nii')
            temp=load_untouch_nii(subj.vol{1});
            temp.img=subj_M2;
            temp.hdr.dime.datatype=16;
            temp.hdr.dime.bitpix=32;
            save_untouch_nii(temp,param.outname);
        else
        
            output2=make_nii(subj_M2,pixd,[],16);
            output2.fileprefix=extract_name(param.outname);
            save_nii(output2,param.outname);
            ret=1;
        end
    elseif strcmp(param.outname(end-2:end),'xml') 
                
           
        f.res=pixd;
        f.type='Float';
        f.orientation='Axial';        
        writeXml(subj_M2,f,param.outname);
        ret=1;
    end
catch e
    
    fprintf('Could not write NIFTI or RAW file, returning the result.\n');
    fprintf('%s\n',e.message);
    ret=subj_M2;
end




% to solve for alpha, atlas (Lx2 matrix) and subj (Lx2 matrix) are fed to
% this function.
function F=calc_alpha(atlsm1,atlsm2,subjm1,subjm2,sig1,sig2)
A=sum((atlsm1(:,1)-atlsm1(:,2)).^2)/sig1^2 + ...
    sum((atlsm2(:,1)-atlsm2(:,2)).^2)/sig2^2;
B=(atlsm1(:,2)-subjm1)'*(atlsm1(:,2)-atlsm1(:,1))/sig1^2 + ...
    (atlsm2(:,2)-subjm2)'*(atlsm2(:,2)-atlsm2(:,1))/sig2^2;
% F=(x^2)*(1-x)*A+x*(1-x)*B+2*x-1;
f=roots([-A A+B 2-B -1]);
F=f(f>=0 & f<=1);
if isempty(F)
    F=0.5;
else
    F=F(1);
end


