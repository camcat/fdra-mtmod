function ok = testmex();
%function ok = testmex();

ok = makemex('private');

if (ok) disp('Good news, could compile mex files. Fast stress calculations are available.')
else disp('Could not compile mex files. Stress calculations will be slow.')
end
