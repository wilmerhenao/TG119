%assuming these were the gantry angles and couch angles
ga = 0:72:359
ca = ga*0;

Dij = 0;

%%
for i=1:length(ga)
  
  fname = ['Gantry' num2str(ga(i)) '_Couch' num2str(ca(i)) '_D.mat']; 
  load(fname)

  [b,j,d] = find (D);
%   if i<10
%       beamID = ['beam00' num2str(i)];
%   elseif i<100
%       beamID = ['beam0' num2str(i)];
%   else
%       beamID = ['beam' num2str(i)];
%   end
%   fid = fopen([beamID 'bixDs.bin'],'w');
%   fwrite(fid,[j';b';d'],'single');
%   fclose(fid);
  
  
  %num beamlets at angle i
  nBPB(i) = size(D,2);
  nDIJSPB(i) = length(j);
  numVoxels = size(D,1);

  
  
end

%% Set correct names

nBeams = length(nBPB);
beams = [0:nBeams-1];
dataFolder = '/media/troy/datadrive/Data/DataProject/TG119/';



%%
structnames = {'Core'; 'OuterTarget';'BODY'};
nStruct = length(structnames);
target = {0;1;0};
meanUB = {10;0;0};
doseLB = {0;20;0};
doseUB = {21;30;0};

for s=1:length(structnames)
    structures{s}{1} = structnames{s};
    structures{s}{2} = target{s};
    structures{s}{3} = meanUB{s};
    structures{s}{4} = doseLB{s};
    structures{s}{5} = doseUB{s};    
end


%% get reduced voxel list


structs = [1:nStruct];
bigZ = zeros(numVoxels,1);
total = 0;
for i=1:length(structs)
  s = structs(i);
  fname = [structnames{s} '_VOILIST.mat'];
  load(fname)
  
  V{s} = v;
  bigZ(v) = 1;
  length(v)
end
nVox = sum(bigZ)
voxelAssignment = zeros(nVox,1);
counter = 1;
for i=1:numVoxels
    if(bigZ(i)>0)
        voxelAssignment(counter) = i;
        counter = counter+1;
    end
end
counter

originalVoxels = zeros(numVoxels,1);
for i=1:nVox
    originalVoxels(voxelAssignment(i)) = i;
end

for s =1:nStruct
   structures{s}{6} = originalVoxels(V{s})
    
end



% singleStruct = zeros(nVox,1);
% 
% for i=1:length(structsPriority)
%     singleStruct(originalVoxels(V{structsPriority(i)})) = 2^(structsPriority(i)-1);    
% end
%% build new sparse matrices

for i=1:length(ga)
  
  fname = ['Gantry' num2str(ga(i)) '_Couch' num2str(ca(i)) '_D.mat']; 
  load(fname)

  [b,j,d] = find (D);
%   if i<10
%       beamID = ['beam00' num2str(i)];
%   elseif i<100
%       beamID = ['beam0' num2str(i)];
%   else
%       beamID = ['beam' num2str(i)];
%   end
%   fid = fopen([beamID 'bixDs.bin'],'w');
%   fwrite(fid,[j';b';d'],'single');
%   fclose(fid);
   newb = originalVoxels(b);
  
   max(newb)
   Dnew = sparse(newb,j,d,nVox,size(D,2));
   fnamenew = ['Gantry' num2str(ga(i)) '_Couch' num2str(ca(i)) '_Dnew.mat']; 
   save(fnamenew,'Dnew');
  
end


%%

save('tg119pythondata.mat','ga','ca','structures','nVox','nDIJSPB','nBPB','dataFolder','nBeams','nStruct','originalVoxels','voxelAssignment','-v7')



%% Here is where to read in data
%(make sure to have ga and ca set the same as above and CERR open and loaded)

% load in intensitiesOut.mat
load('intensitiesOut.mat')

craftViewOptDoseWithDijFilesName(ga,ca,x','testDoseName','/home/wilmer/Documents/TG119/');






