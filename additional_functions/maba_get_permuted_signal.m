function [signal,permuted_signal]=maba_get_permuted_signal(signal,permutation_tps,force)
%
% [signal,permuted_signal]=maba_get_permuted_signal(signal,permutation_tps,force)
%
% Given a matrix "signal" (tps*dims) and a permutation 
% schema "permutation_tps" (tps*perms), function returns the 
% original and permuted signal. The function automatically 
% transposes input matrices to match dimensions (according to tps).
% In case of undefined solutions during transpose (e.g., tps==dims),
% the function return an error.
% To override the error, you need to set 'force' to 'yes', and the
% first dimension of both matrices will be considered as tps.
%
% authored by giacomo.handjaras@imtlucca.it
%

force_flag=false;
if ~exist('force','var')
	force_flag = false;
    force='no';
end
if (sum(ismember(lower(force),{'yes','y'}))>0); force_flag=true; end

%%% Empty matrix
permuted_signal=[];

%%% Let's check dimensionality of the data!
[a1]=size(signal);
[a2]=size(permutation_tps);
if (numel(a1)>2); error('signal should be 1D (tps*1) or 2D (tps*dims)!'); end
if (numel(a2)>2); error('permutation_tps MUST be 2D (tps*perms)!'); end
clear a1 a2

%%% Let's try to fix orientation!
if(force_flag==false)
    [a1,b1]=size(signal);
    [a2,b2]=size(permutation_tps);
    if(b1==a2); signal=signal'; end
    if(a1==b2); permutation_tps=permutation_tps'; end
    if(b1==b2); signal=signal'; permutation_tps=permutation_tps'; end
    if (a2==b2); error('Dimensions of permutation schema are identical: it is confusing!'); end
    if (a1==b1); error('Dimensions of signal are identical: it is confusing!'); end
    if (a1==b1 || a1==b2) && (a2==b1 || a2==b2); error('It is impossible to discriminate dims, tps and perms!'); end
    if (a1==a2 && b1==b2); error('It is impossible to discriminate dims, tps and perms!'); end
    clear a1 a2 b1 b2
end

if(force_flag==true)
    [a1,b1]=size(signal);
    [a2,b2]=size(permutation_tps);
    if (a1~=a2); error('The first dimensions of signal and permutation_tps MUST be identical!'); end
end    

%%% Let's permute the signal
tps=size(permutation_tps,1);
perms=size(permutation_tps,2);

%%% the idea is to apply the same perm schema to the columns of
%%% matrix signal without using a "for" (for each column...)
%%% So, we need to extract an offset that will be summed to the perm schema 
%%% and applied as vector indices to the signal matrix to re-sort tps

indx=[0:tps:((size(signal,2)*tps)-tps)];
indx=repmat(indx,tps,1);

permuted_signal=nan(tps,size(signal,2),perms);
for p=1:perms
    current_permutation_schema=repmat(permutation_tps(:,p),1,size(signal,2))+indx;
    permuted_signal(:,:,p)=signal(current_permutation_schema);
end

%%% In case of a one-column signal, let's remove the singleton dimension
permuted_signal=squeeze(permuted_signal);

end