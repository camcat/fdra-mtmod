path(path,'input');

%create input structure
f=fs_bp3();
%run simulation
%note: if you run this test more than once, the filename needs to be changed; alternatively, you can set f.overwrite=1 to allow overwriting
fdra(f);
%load simulation output:
q=Util('qload',f.saveFn);
%load reference simulation output:
q0=Util('qload','output/ref');

semilogy(q0.t, max(q0.v)); hold on
semilogy(q.t, max(q.v)); hold on
xlabel('Time (s)')
ylabel('V_{max} (m/s)')
legend('Reference output','Output')

