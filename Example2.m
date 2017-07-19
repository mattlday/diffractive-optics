function Example2()
%DOE_EXAMPLE2 initiates a desired DOE geometry, and does a quick simulation to
%demonstrate the Gaussian beam does not follow the designed geometry. An
%optimisation routine is then initiated that corrects for the finite size
%of the source beam, generating new DOE parameters. Another simulation is
%then performed to demonstrate the correct behaviour before the surface
%profile is exported for focused ion beam etching

r=60e-06;
N=10;
object=[0,0,1000].*1e-06;
image=[0,0,500].*1e-06;
direc='f';
material_sub='NBK7';
design_wavelength=405e-09;

% Initiate DOE, notice the geometry
fprintf('Initiating surface relief DOE object\n')
D1=surfacerelief(r,N,object,image,direc,material_sub,material_sub,design_wavelength);
% Let's generate the phase profile so we can simulate it
fprintf('Generating phase profile\n')
D1=D1.gen_p2p_phase_profile();

% This is just prepping the parallel pool to speed up simulation
p=gcp('nocreate');

if isempty(p)
    parpool;
end
% Just a quick simulation of the initisity as you move away from the DOE
% plane
fprintf('Quickly simulating DOE, should be less than a minute\n')
D1.quick_sim('design',1,3.75e-06,D1.design_wavelength,D1.object_coords(3),linspace(10,2000,100).*1e-06);
drawnow
pause(20)
% Uh oh, that doesn't even look to be focusing at all! Let's correct that
fprintf('That does not look right at all, let us correct that\n')
D1=D1.optimisation_1D(3.75e-06);
fprintf('Generating the optimised DOE profile\n')
D1=D1.gen_p2p_phase_profile();
D1=D1.gen_surface_profile();


fprintf('Quickly simulating the optimised DOE profile\n')
D1.quick_sim('design',1,D1.op_source_w0,D1.op_source_wavelength,D1.op_source_dist_from_DOE,linspace(10,1000,100).*1e-06);
drawnow

fprintf('Looks great, now saving it to file for fabrication\n')
D1.save_FIB_map('sdtDOE_example2_FIBmap');

end