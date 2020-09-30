%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018 - Biomedical Signal Processing        %
% MSc Biomedical Engineering - upatras       %
% Elissaios Petrai                           %
% petrai AT ceid.upatras.gr                  %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [locs pks]=Rpeakfinder(x,minpeakdist,minpeakh)

% x is a vector input (generally a timecourse)
% minpeakdist is the minimum desired distance between peaks (optional, defaults to 1)
% minpeakh is the minimum height of a peak (optional)

    if size(x,2)==1, x=x'; end

    % Find all maxima and ties
    locs=find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>=x(3:end))+1;

    if nargin<2, minpeakdist=1; end % If no minpeakdist specified, default to 1.

    if nargin>2 % If there's a minpeakheight
        locs(x(locs)<=minpeakh)=[];
    end

    if minpeakdist>1
        while 1
            del=diff(locs)<minpeakdist;
            if ~any(del), break; end
            pks=x(locs);
            [garb mins]=min([pks(del) ; pks([false del])]); %#ok<ASGLU>
            deln=find(del);
            deln=[deln(mins==1) deln(mins==2)+1];
            locs(deln)=[];
        end
    end

    if nargout>1,
        pks=x(locs);
    end

end