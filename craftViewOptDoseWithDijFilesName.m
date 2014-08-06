function craftViewOptDoseWithDijFilesName(ga,ca,xx,doseName,dataDir)
%function craftViewOptDoseWithDijFiles(s,ga,ca,xx)

global planC stateS
indexS = planC{end};

[nx,ny,nz] = size(planC{3}.scanArray)

d = 0;

bctr=1;

for i=1:length(ga)

  fname = [dataDir 'Gantry' num2str(ga(i)) '_Couch' num2str(ca(i)) '_D.mat']; 
  load(fname)

  [nv,nb] = size(D);
  
  myx = xx(bctr:bctr+nb-1);

  myd = D*myx;
  
  d = d + myd;

  bctr = bctr+nb;
end

dose3D = reshape(d,nx,ny,nz);

name=doseName;
size(dose3D)
showIMDose(dose3D,name)
