function writeVTF(p, tetr, u, filename),
% function writeVTF(p, tetr, u, filename),
% description:
%    writes 3d volumetric data to Ceetron GLview native format (.vtf)
%
% arguments:
%   - p        point cloud (nx3 matrix of n (x,y,z)-coordinates)
%   - tetr     tetrahedral elements. Index to the four corners of element i given in row i.
%   - u        primary solution. Vector of n components
%   - filename name of the file to store the results
% returns:
%   none

% author: Kjetil A. Johannessen
% last edit: September 2012


if nargin < 4,
	filename = 'out.vtf';
	disp 'writing result to "out.vtf"';
end

fid = fopen(filename, 'wt');
fprintf(fid, '*VTF-1.00\n');
fprintf(fid, '\n');
fprintf(fid, '*INTERNALSTRING 40001\n');
fprintf(fid, 'VTF Writer Version info:\n');
fprintf(fid, 'APP_INFO: GLview Express Writer: 1.1-12\n');
fprintf(fid, 'GLVIEW_API_VER: 2.1-22\n');
time = clock;
fprintf(fid, 'EXPORT_DATE: %d-%d-%d %02d:%02d:%02d\n', int32(time));
fprintf(fid, '\n');

fprintf(fid, '*NODES 1\n');
fprintf(fid, '%f %f %f\n', p');
fprintf(fid, '\n');


fprintf(fid, '*ELEMENTS 1\n');
fprintf(fid, '%%NODES #1\n');
fprintf(fid, '%%NAME "Patch 1"\n');
fprintf(fid, '%%NO_ID\n');
fprintf(fid, '%%MAP_NODE_INDICES\n');
fprintf(fid, '%%PART_ID \n');
fprintf(fid, '%%TETRAHEDRONS\n');
fprintf(fid, '%d %d %d %d\n', tetr');
fprintf(fid, '\n');

fprintf(fid, '*RESULTS 2 \n');
fprintf(fid, '%%NO_ID \n');
fprintf(fid, '%%DIMENSION 1 \n');
fprintf(fid, '%%PER_NODE #1 \n');
fprintf(fid, '%f \n', u);
fprintf(fid, '\n');

fprintf(fid, '*GLVIEWGEOMETRY 1\n');
fprintf(fid, '%%STEP 1\n');
fprintf(fid, '%%ELEMENTS\n');
fprintf(fid, '1 \n');
fprintf(fid, '\n');

fprintf(fid, '*GLVIEWSCALAR 1\n');
fprintf(fid, '%%NAME "Primary u"\n');
fprintf(fid, '%%STEP 1\n');
fprintf(fid, '2 \n');
fprintf(fid, '\n');

fclose(fid);


end