function tetr = hex2tetr(hex)
% function tetr = hex2tetr(hex)
% 
% description:
%      generate a tetrahedral element representation from a hexahedral one
%
% arguments:
%   - hex      hexahedral elements described by 8 corner points and two meta values (size = nx10)
% returns:
%   - tetr     tetrahedral elements described by 4 corner points and two meta values (size = (6n)x6 )

% author: Kjetil A. Johannessen
% last edit: October 2012

N = size(hex,1);

q2t = [1,2,6,7; % quad to tetrahedral
       1,5,6,7;
	   1,2,3,7;
	   1,3,4,7;
	   1,4,7,8;
	   1,5,7,8];


tetr = zeros(6*N, 6);
for i=1:N,
	for j=1:6,
		tetr(6*(i-1)+j,:) = hex(i,[q2t(j,:), 9,10]);
	end
end
