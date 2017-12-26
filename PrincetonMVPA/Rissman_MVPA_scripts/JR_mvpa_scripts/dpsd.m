function [Ro,Rn,dp,crit,SSE,pts,flag]=dpsd(oldcum,newcum,varargin)
% 
% usage: [Ro,Rn,dp,crit,SSE,flag]=dpsd(old,new,'Rn'[opt],'trim'[opt],'plot'[opt])
%
% Estimates Recollection (Ro and Rn) and familiarity (d-prime), based on a
% dual-process signal detection model (Yonelinas, 2005), given response
% counts for old and new items.
%
% For N confidence levels, <old>, <new>, are each Nx1 arrays with each cell
% containing a response count for one confidence level (arranged from least
% confident, in cell 1, to most confident in cell N).
%
% (optional) Enter argument 'Rn' if you want to fit Rn. (Not sure how well
% this works; doesn't seem to plot properly.)
%
% (optional) Enter argument 'trim' if you want to trim points at the margin
% of model space (recommended).
%
% (optional) Enter argument 'print' to plot ROC.




% new1=new/sum(new); % compute probability distribution
% new2=flipud(new1); % arrange from most confident to least confident
% newcum=new2;%cumsum(new2); % compute cumulative probability distribution
% 
% old1=old/sum(old);
% old2=flipud(old1);
% oldcum= old2;%cumsum(old2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim off any points on an axis (conservative). Optional; off by default.
if strmatch ('trim',varargin)
    if or(newcum(1)==0,oldcum(1)==0)
        newcum=newcum(2:length(newcum));
        oldcum=oldcum(2:length(oldcum));
    end

    if or(newcum(length(newcum))==1,oldcum(length(oldcum))==1)
        newcum=newcum(1:length(newcum)-1);
        oldcum=oldcum(1:length(oldcum)-1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pts=length(oldcum); % How many data points did we end up using?

observed=[oldcum';newcum'];  %JR edit:  added transpose tics

starters=[.2 .2 .5 -2:4/(length(newcum)-1):2]; % Initial values for Ro, Rn, dp, N criteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pass starters to function solver fminsearch

options = optimset('Display','off', 'MaxFunEvals',100000,'TolX',1E-5,'TolFun',1E-6);
% set params for fminsearch, see help for fminsearch and optimset for
% details

[x, SSE, flag] = fminsearch(@predroc, starters', options);
% find the values for x that minimize fval for predroc fn (below)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unwrap final parameters (x)

Ro=abs(x(1));
if strmatch('Rn',varargin)
    Rn=abs(x(2));
else
    Rn=0;
end

dp=x(3);
crit=x(4:length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print ROC
if strmatch('print',varargin)
    chance=0:1;
    predfarP=0:.01:1;
    cP=norminv(predfarP,0,1);
    predhitP=min(1,predfarP+Ro+((1-Ro)*normcdf(cP,-dp,1))-(1-Rn)*normcdf(cP,0,1));

    plot(predfarP,predhitP,'b',newcum,oldcum,'bx',chance,chance,'k:')
    legend('predicted','observed','chance','Location','SouthEast')
    title(sprintf('Ro:%d\tRn:%d\td'':%d\tSSE%d',Ro,Rn,dp,SSE))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for fminsearch to optimize

    function SSE=predroc(x)
        Ro=abs(x(1));
        if strmatch('Rn',varargin)
            Rn=abs(x(2));
        else
            Rn=0;
        end
        dp=x(3);
        crit=x(4:length(x));

        m=1:pts;

        predhit(m,1)=Ro+(1-Ro)*normcdf(crit(m), -dp, 1); % Model predicted Hit distribution
        predfar(m,1)=(1-Rn)*normcdf(crit(m), 0, 1); % Model predicted FA distribution

        predicted=[predhit;predfar];


        SSE=sum((observed-predicted).^2);
    end
end

