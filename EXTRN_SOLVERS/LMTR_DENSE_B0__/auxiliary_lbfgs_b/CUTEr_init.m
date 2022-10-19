function CUTEr_init()

if ispc
  disp('Running on Windows');
  disp('No installation available!');
end

if isunix,
  disp('Running on Unix');
  
  % adapt the next two lines for your CUTEr installation folder
  addpath('/sw/opt/CUTEr/130204-r150/cuter2/CUTEr.large.pc.lnx.gfo/bin');
  addpath('/sw/opt/CUTEr/130204-r150/cuter2/common/src/matlab');
end