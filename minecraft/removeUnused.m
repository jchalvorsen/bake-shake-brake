function [p tetr] = removeUnused(p, tetr),
% function [p tetr, tri] = removeUnused(p, tetr, tri),
% 
% description:
%    removes all nodal points not appearing in the element array 
%
% arguments:
%    - p      nodal points  (matrix of size Npts x 4 )
%    - tetr   element array (matrix of size (Nel x 4) or (Nel x 6) )
%    - tri    boundary elements (matrix of size (Nel x 3) or (Nel x 5) )
%   
% returns:
%    - p      the stripped down nodal array (smaller size)
%    - tetr   updates nodal indexes of the element array (updated values)
%    - tri    updates nodal indexes of the boundary elements (updated values)

% author: Kjetil A. Johannessen
% last edit: November 2013

tmp  = tetr(:,1:4);
N    = size(p,1);
used = zeros(N,1);
used(tmp(:)) = 1;
uu = sort(find(used == 1)); 

offset = zeros(1,N);
totOffset = 0;
k = 1;
for i=1:N,
	if(k <= numel(uu) && uu(k) ~= i),	
		totOffset = totOffset + 1;
	else
		k = k + 1;
	end
	offset(i) = totOffset;
end
p = p(uu,:);
Nel = size(tetr, 1);
for i=1:Nel,
	tetr(i,1:4) = tetr(i,1:4) - offset(tetr(i,1:4));
end

% Nel = size(tri, 1);
% for i=1:Nel,
% 	tri(i,1:3) = tri(i,1:3) - offset(tri(i,1:3));
% end
