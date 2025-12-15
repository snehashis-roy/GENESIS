function ret=HU2umap(input,output,kVp)
% function ret=HU2umap(inputvol,outputvol,kVp)
%	
% kVp: kVp used to generate CT image. (default = 120 for Siemens mCT)
%
% Based on Med. Phys. 33(4) April 2006
% "Method for transforming CT images for attenuation correction in PET/CT imaging"
%
% Before breaking point
% u = 9.6 x 10^-5 (HU+1000)
%
% After breaking point
% u = a*(HU+1000) + b
% 
% kVp = 120 is default (must be string, like '120')
% If output is empty [], the volume is returned


if ~isnumeric(input)
    if strcmpi(input(end-2:end),'nii')
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

inMat=double(invol);
if (nargin<2), kVp = 120; end

kVp=str2num(kVp);

switch kVp
case 80
	a=3.64*10^-5; b=6.26*10^-2; BP=1050;
	disp('kVp = 80');
case 100
	a=4.43*10^-5; b=5.44*10^-2; BP=1052;
	disp('kVp = 100');
case 110
	a=4.92*10^-5; b=4.88*10^-2; BP=1043;
	disp('kVp = 110');
case 120
	a=5.10*10^-5; b=4.71*10^-2; BP=1047;
	disp('kVp = 120');
case 130
	a=5.51*10^-5; b=4.27*10^-2; BP=1037;
	disp('kVp = 130');
case 140
	a=5.64*10^-5; b=4.08*10^-2; BP=1030;
	disp('kVp = 140');
otherwise
	disp('not supported');
	a=0.00*10^-5; b=0.01*10^-2; BP=1000;
end


inMat = inMat + 1000;
ind1 = inMat<=BP;
ind2 = inMat>BP;
outMat = inMat;
outMat(ind1) = 9.6*10^-5*outMat(ind1);
outMat(ind2) = a*outMat(ind2) + b;
outvol = outMat * 10^4;

if isempty(output)
    ret=outvol;
    fprintf('Output is returned as %d x %d x %d matrix \n',size(ret,1),size(ret,2),size(ret,3));    
else
    if strcmpi(input(end-2:end),'nii')
        output(end-2:end)='nii';
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
