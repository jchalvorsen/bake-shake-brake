function writeVTF2(p, tetr, varargin),
% function writeVTF2(p, tetr, <name>, value, <name>, value, ...),
%
% description:
%    writes 3d volumetric data to Ceetron GLview native format (.vtf)
%
% arguments:
%   - p        point cloud (nx3 matrix of n (x,y,z)-coordinates)
%   - tetr     tetrahedral elements. Index to the four corners of element i given in row i.
%   - <name>   name of result set, a descriptive string
%   - value    values of that particular results
%
%  The interpretation of the result values are based on their size:
%     n x 1    - vectors are scalar results       (ex: temperature field)
%     n x 3    - matrices are vector results      (ex: displacement field)
%     n x m    - matrices are scalar time results (ex: time-dependent temperature)
%     n x 3 x m- 3D-matrices are multiple vector results (ex: eigenmodes)
%
%  Except for the following special <name> tags:
%    'FileName' - name of file to store results in
%    'Time'     - scalar values of the time iteration
%    
% returns:
%   none
%   
% EXAMPLES:
%   % store a scalar temperature field 
%   load tetmesh;
%   u = sqrt(sum(X.^2,2));
%   writeVTF2(X, tet, 'Temperature', u, 'FileName', 'myTemp.vtf');
%    
%   % store a time-varying temperature field
%   load tetmesh;
%   t  = linspace(0,2*pi, 50);
%   u  = sqrt(sum(X.^2,2));
%   Ut = u * t;
%   writeVTF2(X, tet, 'Temperature', Ut, 'Time', t, 'FileName', 'myTempTime.vtf');
%    
%   % store a vector displacement field
%   load tetmesh;
%   uz = 10*sin(2*pi*X(:,3)/max(X(:,3)));
%   uy = sqrt(sum(X.^2,2));
%   ux = zeros(size(uz));
%   U  = [ux, uy, uz];
%   writeVTF2(X, tet, 'Displacement', U, 'FileName', 'myDisplacement.vtf');
%    
%   % store a set of 3 eigenmodes
%   load tetmesh;
%   n  = size(X,1);
%   uz = 10*sin(2*pi*X(:,3)/max(X(:,3)));
%   uy = sqrt(sum(X.^2,2));
%   ux = zeros(size(uz));
%   allU = zeros(n,3,3);
%   allU(:,:,1) = [ux, uy, uz];
%   allU(:,:,2) = [uz, uy, ux];
%   allU(:,:,3) = [uz, uz, uz];
%   writeVTF2(X, tet, 'Eigenmodes', allU, 'FileName', 'myEigenmodes.vtf');
%    
%   % store a vector displacement field along with a scalar von mises stress
%   load tetmesh;
%   uz    = 10*sin(2*pi*X(:,3)/max(X(:,3)));
%   uy    = sqrt(sum(X.^2,2));
%   ux    = zeros(size(uz));
%   sigma = 10*sin(4*pi*X(:,1)/max(X(:,1)));
%   U     = [ux, uy, uz];
%   writeVTF2(X, tet, 'Displacement', U, 'Von Mises Stress', sigma, 'FileName', 'myDisplacement.vtf');
%

% author: Kjetil A. Johannessen
% modified by: Geir Bogfjellmo
% last edit: November 2014

filename = '';
timeSteps     = [];

for i=1:2:numel(varargin),
	if(strcmp(varargin{i}, 'FileName')),
		filename = varargin{i+1};
	elseif(strcmp(varargin{i}, 'Time')),
		timeSteps = varargin{i+1};
	end
end

if(numel(timeSteps) == 0),
	timeSteps = 1:size(p,1);
end

if strcmp(filename, ''),
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


res_id = 2;
for i=1:2:numel(varargin),
	if(strcmp(varargin{i}, 'FileName') || strcmp(varargin{i}, 'Time') )
		continue;
	end
	u = varargin{i+1};
	n = size(u);

	if(numel(n) > 2)  % time vector field
		for j=1:n(3),
			fprintf(fid, '*RESULTS %d \n', res_id);
			fprintf(fid, '%%NO_ID \n');
			fprintf(fid, '%%DIMENSION 3 \n');
			fprintf(fid, '%%PER_NODE #1 \n');
			fprintf(fid, '%f %f %f \n', u(:,:,j)');
			fprintf(fid, '\n');
			res_id = res_id + 1;
		end
	elseif(n(2) > 3) % time scalar field
		for j=1:n(2),
			fprintf(fid, '*RESULTS %d \n', res_id);
			fprintf(fid, '%%NO_ID \n');
			fprintf(fid, '%%DIMENSION 1 \n');
			fprintf(fid, '%%PER_NODE #1 \n');
			fprintf(fid, '%f \n', u(:,j));
			fprintf(fid, '\n');
			res_id = res_id + 1;
		end
	elseif(n(2) == 3) % vector field
		fprintf(fid, '*RESULTS %d \n', res_id);
		fprintf(fid, '%%NO_ID \n');
		fprintf(fid, '%%DIMENSION 3 \n');
		fprintf(fid, '%%PER_NODE #1 \n');
		fprintf(fid, '%f %f %f \n', u');
		res_id = res_id + 1;
	elseif(n(2) == 1) % scalar field
		fprintf(fid, '*RESULTS %d \n', res_id);
		fprintf(fid, '%%NO_ID \n');
		fprintf(fid, '%%DIMENSION 1 \n');
		fprintf(fid, '%%PER_NODE #1 \n');
		fprintf(fid, '%f \n', u);
		res_id = res_id + 1;
	else
		info = '';
		sprintf(info, 'Nonvalid dimensions (%d,%d) of solution field %s', n(1), n(2), varargin{i});
		disp(info);
		fclose(fid);
		return;
	end
	fprintf(fid, '\n');
end

fprintf(fid, '*GLVIEWGEOMETRY 1\n');
fprintf(fid, '%%STEP 1\n');
fprintf(fid, '%%ELEMENTS\n');
fprintf(fid, '1 \n');
fprintf(fid, '\n');

res_id = 2;
for i=1:2:numel(varargin),
	if(strcmp(varargin{i}, 'FileName') || strcmp(varargin{i}, 'Time') )
		continue;
	end
	u = varargin{i+1};
	n = size(u);

	if(numel(n) > 2)  % time vector field
		fprintf(fid, '*GLVIEWVECTOR %d\n', i);
		fprintf(fid, '%%NAME "%s"\n', varargin{i});
		for j=1:n(3),
			fprintf(fid, '%%STEP %d\n', j);
			fprintf(fid, '%%STEPNAME %d\n', timeSteps(j));
			fprintf(fid, '%d \n', res_id);
			res_id = res_id + 1;
		end
		fprintf(fid, '\n');
	elseif(n(2) > 3) % time scalar field
		fprintf(fid, '*GLVIEWVECTOR %d\n', i);
		fprintf(fid, '%%NAME "%s"\n', varargin{i});
		for j=1:n(2),
			fprintf(fid, '%%STEP %d\n', j);
			fprintf(fid, '%%STEPNAME %d\n', timeSteps(j));
			fprintf(fid, '%d \n', res_id);
			res_id = res_id + 1;
		end
		fprintf(fid, '\n');
	elseif(n(2) == 3) % vector field
		fprintf(fid, '*GLVIEWVECTOR %d\n', i);
		fprintf(fid, '%%NAME "%s"\n', varargin{i});
		fprintf(fid, '%%STEP 1\n');
		fprintf(fid, '%d \n', res_id);
		fprintf(fid, '\n');
		res_id = res_id + 1;
	elseif(n(2) == 1) % scalar field
		fprintf(fid, '*GLVIEWVECTOR %d\n', i);
		fprintf(fid, '%%NAME "%s"\n', varargin{i});
		fprintf(fid, '%%STEP 1\n');
		fprintf(fid, '%d \n', res_id);
		fprintf(fid, '\n');
		res_id = res_id + 1;
	end
end

fclose(fid);

end