
% [p,J]=find_peaks(filename,dim,type,thresh,option,delta)
% p is a vector having peaks, this code takes a volume and finds the maxima/peaks
% in the 1D histogram. It uses KDE1d and uses a default smoothing parameter
% of 1 for the KD estimator. If option=1, filename is a name of the file
% and dimension, type and thresh are needed. If option=2, filename is a vector
% and dim, type, could be empty vectors, however thresh could be a numeric.
% delta is the smoothing of KDE
% J is the height of distribution at the peaks
% 
% If delta==0, then it is found automatically
%
function [p,J]=find_peaks(filename,dim,type,thresh,option,delta)

if option==1
    if strcmp(filename(end-2:end),'nii') || strcmp(filename(end-5:end),'nii.gz')
        temp=load_nii(filename);
        x=double(temp.img);
        x=x(:);
    elseif strcmp(filename(end-2:end),'xml')
        x=ReadXml(filename);
        x=x(:);
    else
        fprintf('Only NIFTI and XML files are supported.\n');
        return;
    end      

elseif option==2
    x=filename;
else
    fprintf('option must be 1 or 2 : **ERROR**\n');
    p=[];
%     return ;
end
a= x>thresh;
x1=x(a);
x1=x1(:);

if delta==0
    [I, ~ ,~ ,~]=OutlierReduction(x1,0,0,[1 2]);
    
    x1(x1>I(2))=0;
    x1(x1<I(1))=0;

    delta=(I(2)-I(1))/80;
end
minn=min(x1);maxx=max(x1);

U=minn:delta:maxx; % default range of density estimator

ff=KDE1d(U,x1,delta);
% g=KDE1d(U,x1,0.5);figure;plot(g);grid minor;title(filename);
L=length(ff);

f1=gradient(ff);
f2=gradient(f1);
s=sign(f1);
s1=filter([1 1],1,s);    
count=1;
J=[];
for i=1:length(s1)
    if s1(i)==0
        if f2(i)<0
            p(count)=U(i);
            J(count)=ff(i);
            count=count+1;
        end
    end
end


    
    
