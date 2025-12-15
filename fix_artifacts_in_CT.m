function fix_artifacts_in_CT(CTimage,brainmask,outputimage)

% fix_artifacts_in_CT(CTimage,brainmask,outputimage)

if nargin==2
    outputimage=CTimage;
    fprintf('Input image with be overwritten. Are you sure?\n');
    fprintf('I Will wait for 10 seconds. Press control+C to exit if you don''t want to overwrite.\n');
    pause(10);
end
ct=load_untouch_nii(CTimage);
mask=load_untouch_nii(brainmask);
H=strel('disk',3);
mask.img=imerode(mask.img,H);
mask.img=imerode(mask.img,H);
indx=find(mask.img>0);
x=ct.img(indx);
x(x<0)=10;
ct.img(indx)=x;

x=ct.img(indx);
x(x>120)=120;
ct.img(indx)=x;
fprintf('Writing %s\n',outputimage);
save_untouch_nii(ct,outputimage);


