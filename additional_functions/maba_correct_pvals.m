function [pvals_corr] = maba_correct_pvals(pvals,verbosity,method)
%
% [pvals_corr] = maba_correct_pvals(pvals,verbosity,method)
%
% 'pvals' is a vector of raw pvalues
% 'method' could be:    
%       a) 'fdr95' for fdr pdep (Bejnamini & Hochberg 1995)
%       b) 'fdr01' for fdr dep (Benjamini & Yekutieli 2001)
% Argument 'verbosity' handles debug msg on screen ('yes' | 'no')
%
% authored by giacomo.handjaras@imtlucca.it
%

verbose=false;
if ~exist('verbosity','var')
	verbosity = '';
end
if (strcmp(verbosity,'yes')); verbose=true; end

mustBeMember(method,{'fdr95','fdr01'}); 

pvals_corr=nan(size(pvals));
pvals_mask=find(~isnan(pvals));

if verbose; fprintf('%d results out of %d with a raw p <0.05\n',sum(pvals(pvals_mask)<0.05), numel(pvals_mask)); end

if (strcmp(method,'fdr95'))
[h, crit_p, adj_p]=fdr_bh(pvals(pvals_mask),0.05,'pdep','no');
pvals_corr(pvals_mask)=adj_p;
end

if (strcmp(method,'fdr01'))
[h, crit_p, adj_p]=fdr_bh(pvals(pvals_mask),0.05,'dep','no');
pvals_corr(pvals_mask)=adj_p;
end

if verbose; fprintf('Correction using [\b %s]\b: %d results out of %d with corrected p <0.05\n',method,sum(adj_p<0.05), numel(adj_p)); end

end