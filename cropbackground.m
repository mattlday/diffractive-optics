function out=cropbackground(A,plot)
%% CROPBACKGROUND Beam profiles often have large background areas surrounding the beam, containing small amounts of power, this function crops the background so that the beam is dominant
% Inputs
% A - Array containing profile to be cropped
%
% Outputs
% out - Cropped array

s=size(A);
m=median(median(A)); % finds dominant value in the array, if the background is large enough, then the background terms should be the median.

if m>10
    disp('Warning: the median value of inputted array is higher than expected for an array dominated by background, check array is not already cropped')
end


[irow,icol]=find(A==max(max(A)));

irow=irow(1);
icol=icol(1);

%% North-east check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    jrow=jrow+1;
    jcol=jcol+1;
    
    if jrow==s(1) || jcol==s(2)
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(1)=jrow;
Jcol(1)=jcol;

%% East check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    jcol=jcol+1;
    
    if jrow==s(1) || jcol==s(2)
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(2)=jrow;
Jcol(2)=jcol;

%% South-east check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    jrow=jrow-1;
    jcol=jcol+1;
    
    if jrow==1 || jcol==s(2)
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(3)=jrow;
Jcol(3)=jcol;

%% South check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    jrow=jrow-1;
    %jcol=jcol+1;
    
    if jrow==1 || jcol==s(2)
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(4)=jrow;
Jcol(4)=jcol;

%% South-west check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    jrow=jrow-1;
    jcol=jcol-1;
    
    if jrow==1 || jcol==1
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(5)=jrow;
Jcol(5)=jcol;

%% West check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    %jrow=jrow-1;
    jcol=jcol-1;
    
    if jrow==1 || jcol==1
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(6)=jrow;
Jcol(6)=jcol;

%% North-west check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    jrow=jrow+1;
    jcol=jcol-1;
    
    if jrow==s(1) || jcol==1
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(7)=jrow;
Jcol(7)=jcol;

%% North check
probe=A(irow,icol);
jrow=irow;
jcol=icol;
while m<probe
    jrow=jrow+1;
    
    if jrow==s(1) || jcol==1
        break
    end
    
    probe=A(jrow,jcol);
    
end

Jrow(8)=jrow;
Jcol(8)=jcol;

out=A(min(Jrow):max(Jrow),min(Jcol):max(Jcol));

if plot==1
    figure
    imagesc(A)
    figure
    imagesc(out)
end
end