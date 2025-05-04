function J=imshift(I,u,inttype)
% USAGE    : J=imshift(I,u,inttype);
% FUNCTION : Shifts the image I by real(u) along the lines and by imag(u)
% along the columns. 'u' is an arbitrary complex vector field. If 'u' is a 
% scalar, all the pixels are shifted by the same amount.
%
% This function uses cubic OMOMS interpolation used in 'interp.m' with
% zero-padding image extension (can be changed).
%
% The optional argument inttype is a cell with two elements. inttype{1} can 
% take the values
%               * 'nearestneighbor'
%               * 'bilinear' 
%               * 'keys'
%               * 'cubicspline'
%               * 'cubicOMOMS' (default)
%               * 'shiftedlinear'
% and inttype{2} can take the values
%               * 'zpd'  (default)
%               * 'symh' (half-point symmetric extension)
%               * 'symw' (whole-point symmetric extension)
%               * 'ppd'  (periodization)
%
% DATE     : 23 November 2014
% AUTHOR   : Thierry Blu, mailto:thierry.blu@m4x.org

[M,N,P]=size(I);
if nargin==3
    if ischar(inttype)
        kernel=inttype;
        exttype='zpd';
    else
        if length(inttype)==1
            kernel=inttype{1};
            exttype='zpd';
        else
            kernel=inttype{1};
            exttype=inttype{2};
        end
    end
else
    kernel='cubicOMOMS';
    exttype='zpd';
end

if length(u)>1
    [M0,N0]=size(u);
    if M0==M&N0==N
        x=(1:M)'*ones(1,N);
        y=ones(M,1)*(1:N);
        J=zeros(M,N,P);
        for p=1:P
            J(:,:,p)=interp(x-real(u),y-imag(u),double(I(:,:,p)),{kernel exttype});
%             J(:,:,p)=interp(x-real(u),y-imag(u),double(I(:,:,p)),kernel);

        end
        J=cast(J,class(I));
    else
        error('Input image and flow dimensions do not match!')
    end
else
    integer_shift=(abs(round(u)-u)<=1e-6);

    if integer_shift
        u=round(u);
        xshift=abs(real(u));
        yshift=abs(imag(u));
        switch sign(real(u))
            case -1
                I1=[I;I(end-(1:xshift)+1,:)];
            case 0
                I1=I;
            case 1
                I1=[I((xshift:-1:1),:);I];
        end
        switch sign(imag(u))
            case -1
                I1=[I1 I1(:,end-(1:yshift)+1)];
            case 0
            case 1
                I1=[I1(:,(yshift:-1:1)) I1];
        end

        J=I1((1:M)+max(0,-real(u)),(1:N)+max(0,-imag(u)));
    else
        x=(1:M)'*ones(1,N);
        y=ones(M,1)*(1:N);
        J=zeros(M,N,P);
        for p=1:P
            J(:,:,p)=interp(x-real(u),y-imag(u),double(I(:,:,p)),{kernel exttype});
        end
        J=cast(J,class(I));
    end
end