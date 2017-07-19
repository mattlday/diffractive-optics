function Example3()
%DOE_EXAMPLE3 initiates a DOE geometry, displays the geometry, takes a
%brief pause to allow the geometry to be observed and then generates an
%electron beam lithography map to fabricate the DOE in PMMA

r=60e-06;
N=10;
object=[0,0,1000].*1e-06;
image=[0,0,300].*1e-06;
direc='f';
material_sub='NBK7';
material_DOE='PMMA';
design_wavelength=405e-09;

fprintf('Initiating sdtDOE object\n')

D1=surfacerelief(r,N,object,image,direc,material_sub,material_DOE,design_wavelength);

fprintf('Plotting design geometry to visualise location of source and focal point\n')

D1.plot_geometry;
drawnow

pause(20)

fprintf('Generating phase profile and corresponding surface relief in PMMA \n')
D1=D1.gen_p2p_phase_profile();
D1=D1.gen_surface_profile();

fprintf('Generating EBL dose map, including correcting for dose behaviour with feature size \n')

D1=D1.gen_EBL_map(1,1,400e-09);
D1.save_EBL_map('DOE_example3_EBLmap','full')

end