function [pts el] = getLargestConnectedDomain(nodes, elements),
% function [pts el] = getLargestConnectedDomain(nodes, elements),
% 
% description:
%      generate an element description of the largest connected domain
%      described by an arbitrary set of points and elements.
%
% arguments:
%   - nodes     the nodal coodinates (size = nx3)
%   - elements  the hexahedral element description (size = mx10)
% returns:
%   - pts       nodal points of largest connected domain
%   - el        elements of largest connected domain

% author: Kjetil A. Johannessen
% last edit: November 2012

% load('elements.m')
% load('nodes.m')

elements(:,1:8) = elements(:,1:8) + 1;

N    = size(elements, 1);
Npt  = size(nodes, 1);
tetr = zeros(N*6, 4);

v = 8;

i = zeros(N*v^2, 1);
j = zeros(N*v^2, 1);
s = zeros(N*v^2, 1);

disp 'target range:'
N
for k = 1:N
    n = elements(k,1:8);
    l = (k-1)*v^2+1:k*v^2;
    [j(l) i(l)] = meshgrid(n, n);                   %Correct indices
    s(l) = 1;
	if mod(k, 5000) == 0,
		k
	end
end

C = sparse(i, j, s, Npt, Npt);                     %Construct matrix

disp 'sorting results'

start = zeros(Npt+1,1);
C = sortrows([i,j]);

disp 'creating start index array'
lastId = C(1,1);
start(lastId) = 1;
for k=1:numel(i),
	if(C(k,1) ~= lastId),
		lastId = C(k,1);
		start(lastId) = k;
	end
end

used = zeros(Npt,1);
stack = java.util.Stack();

allSize = zeros(Npt,1);

nPartitions = 0;
isDone = false;

while ~isDone,
	disp 'looking for parition '
	nPartitions+1
	didFindIt = false;
	for i=1:Npt,
		if used(i) == 0,
			disp 'found it'
			disp 'Total nodes covered:'
			sum(used ~= 0)
			stack.push(i);
			nPartitions = nPartitions + 1;
			used(i)     = nPartitions;
			thisSize    = 1;
			didFindIt   = true;
			break;
		end
	end
	isDone = ~didFindIt;
		
	while stack.size() > 0
		curNode = stack.pop();

		if(start(curNode) == 0),
			break;
		end

		s0 = start(curNode);
		s1 = start(curNode+1)-1;

		for k=s0:s1,
			if( used(C(k,2)) == 0  ),
				stack.push(C(k,2));
				used(C(k,2)) = nPartitions;
				thisSize = thisSize + 1;
			end
		end
		if mod(thisSize, 5000)==4999,
			thisSize
		end
		if mod(stack.size(), 1000)==999,
			disp 'huge stack'
			stack.size()
		end
	end

	allSize(nPartitions) = thisSize;
end
allSize = allSize(1:nPartitions)

m = max(allSize);
bestParition = find(allSize == m);

uu = sort(find(used == bestParition)); 

offset = zeros(Npt,1);
totOffset = 0;
k = 1;
for i=1:Npt,
	if(k <= numel(uu) && uu(k) ~= i),	
		totOffset = totOffset + 1;
	else
		k = k + 1;
	end
	offset(i) = totOffset;
end
	
properEl = zeros(N,1);
elCount  = 0;

for i=1:N,
	if(used(elements(i,2)) == bestParition),
		k = elements(i,1:8);
		off = offset(k);
		elements(i,1:8) = elements(i,1:8) - off';
		elCount = elCount + 1;
		properEl(elCount) = i;
	end
end

properEl = properEl(1:elCount);
el  = elements(properEl,:);
pts = nodes(uu,:);
save('el.m', 'el', '-ascii');
save('pts.m', 'pts', '-ascii');
