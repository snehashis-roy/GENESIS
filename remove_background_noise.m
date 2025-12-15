function ret=remove_background_noise(input,quant,output)
% function ret=remove_background_noise(INPUT,QUANT,OUTPUT)
% INPUT    An input volume (nifti nii/nii.gz or xml file or a 3D matrix)
% QUANT    Initial quantile for threshold. E.g. 40 means 40% of the image
%          histogram will be used as initial threshold. Enter as a character
%          array, e.g. '40'.
% OUTPUT   An output volume (nifti nii or xml), nifti preferable. If it is not
%          mentioned, the output is returned as a 3D matrix

if nargin==1
    output=[];
    quant=40;
    
elseif nargin==2
    if ~ischar(quant)
        fprintf('ERROR: Enter the Quantile value as string array. e.g. ''40''\n');
        ret=0;
        return;
    end
    quant=str2num(quant);    
    output=[];
elseif nargin==3  
    fprintf('Output will be written in \n');
    fprintf('%s\n',output);
    quant=str2num(quant);
end
fprintf('Quantile = %.2f %%\n',quant);

if ~isnumeric(input)
    if strcmpi(input(end-2:end),'nii') || strcmpi(input(end-5:end),'nii.gz')
        temp=load_untouch_nii(input);
        invol=double(temp.img);
        pixd=temp.hdr.dime.pixdim(2:4);
    elseif strcmpi(input(end-2:end),'xml')
        [invol,xmlparam]=ReadXml(input);
        pixd=xmlparam.res;
    else
        fprintf('Only XML and NII files are supported\n');
        ret=[];
        return;
    end
else
    invol=input;
    pixd=[1 1 ];
    fprintf('Input file size %d x %d x %d\n',size(invol,1),size(invol,2),size(invol,3));
end
invol=double(invol);
Q=quantile(invol(invol>0),quant/100);
fprintf('Initial threshold = %.2f\n',Q);
H=strel('ball',3,3);
start=invol;
for iter=1:4 
    temp=start;    
    initmask=uint8(temp>Q);
    mask=fillholes(initmask);
    mask=fillholes(mask);
    for t=1:3
        mask=imclose(mask,H);
    end
    mask=fillholes(mask);
    err=1-length(find(initmask~=mask & invol>0))/sum(invol(:)>0);
    fprintf('iter = %d, Q= %.2f, change = %.4f \n',iter,Q,err);
    start=invol.*double(mask);
    Q=1.10*Q; % 10% increment in each iteration. If very aggressive stripping, reduce to lower number
end

mask=remove_small_components(mask);
se=strel('ball',3,3);
mask=imclose(mask,se);

outvol=invol.*double(mask);
if isempty(output)
    ret=outvol;
    fprintf('Output is returned as %d x %d x %d matrix \n',size(ret,1),size(ret,2),size(ret,3));    
else
    if strcmpi(input(end-2:end),'nii') || strcmpi(input(end-5:end),'nii.gz')
%         output(end-2:end)='nii';
        temp=load_untouch_nii(input);
        temp.img=outvol;
        temp.hdr.dime.bitpix=32;
        temp.hdr.dime.datatype=16;
        save_untouch_nii(temp,output);
        ret=1;
    elseif strcmpi(input(end-2:end),'xml')
        if strcmpi(output(end-2:end),'xml')
            f.res=pixd;
            f.type='Float';
            f.orientation='Axial';
            writeXml(outvol,f,output);
            ret=1;
        elseif strcmpi(output(end-2:end),'nii')
            temp=make_nii(outvol,pixd,[],16);           
            save_nii(temp,output);
            ret=1;
        end
    end
end
    


function outmask=fillholes(inmask)
inmask=uint8(inmask);
U=length(unique(inmask(:)));
if U~=2
    fprintf('ERROR: Input mask must be binary.\n');
    outmask=inmask;
    return;
end
dim=size(inmask);
% fprintf('Z ');
for i=1:dim(3)
    inmask(:,:,i)=imfill(inmask(:,:,i),4,'holes');
end
% fprintf('X ');
for i=1:dim(1)
    inmask(i,:,:)=imfill(squeeze(inmask(i,:,:)),4,'holes');
end
% fprintf('Y ');
for i=1:dim(2)
    inmask(:,i,:)=imfill(squeeze(inmask(:,i,:)),4,'holes');
end
% fprintf('Z\n');
for i=1:dim(3)
    inmask(:,:,i)=imfill(inmask(:,:,i),4,'holes');
end
outmask=inmask;


function mask=remove_small_components(mask)

for k=1:size(mask,3)
    temp=squeeze(mask(:,:,k)); 
    cc=bwconncomp(temp,4);
    f=[];for t=1:cc.NumObjects f(t)=length(cc.PixelIdxList{t}); end;
    L=max(f);
    for t=1:cc.NumObjects
        if length(cc.PixelIdxList{t})<0.02*L
            temp(cc.PixelIdxList{t})=0;
        end
    end
    mask(:,:,k)=temp;
end
for k=1:size(mask,1)
    temp=squeeze(mask(k,:,:)); 
    cc=bwconncomp(temp,4);
    f=[];for t=1:cc.NumObjects f(t)=length(cc.PixelIdxList{t}); end;
    L=max(f);
    for t=1:cc.NumObjects
        if length(cc.PixelIdxList{t})<0.02*L
            temp(cc.PixelIdxList{t})=0;
        end
    end
    mask(k,:,:)=temp;
end
for k=1:size(mask,2)
    temp=squeeze(mask(:,k,:)); 
    cc=bwconncomp(temp,4);
    f=[];for t=1:cc.NumObjects f(t)=length(cc.PixelIdxList{t}); end;
    L=max(f);
    for t=1:cc.NumObjects
        if length(cc.PixelIdxList{t})<0.02*L
            temp(cc.PixelIdxList{t})=0;
        end
    end
    mask(:,k,:)=temp;
end

    


