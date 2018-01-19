function[x] = slices(m,d)
%returns slices along the d^th (or last) dimension of m
if ~exist('d','var')
    d = ndims(m);
end
x = arrayfun(@slice_matrix,repmat({m},1,size(m,d)),d*ones(1,size(m,d)),1:size(m,d),'UniformOutput',false);

function[x] = slice_matrix(m,d,i)
%return the ith slice along the d^th dimension of m
if iscell(m) && (length(m) == 1)
    m = m{1};
end
x = eval(sprintf(['m(',indexer(d,ndims(m)),')'],i));

function[s] = indexer(i,n) %n: number of dimensions in m; i: index to change
n = max([i n]);
if max([n i]) < 1
    s = '';
    return;
end

%make a string of n ":"'s, separated by commas
s = join(',',repmat({':'},1,n));

%replace the ith ':' with '%d'
inds = strfind(s,':');
sep = inds(i);
s = [s(1:sep-1),'%d',s(sep+1:end)];
