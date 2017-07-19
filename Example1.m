function Example1()
%sdtDOE_EXAMPLE1 initiates an sdtDOE object, generates its phase profile and simulates its focal plane and
%calculates its beam waist. For details of each function, type 'help
%sdtDOE.function.

r=60e-06; % radius of DOE
N=10; % number of points per micron
object=[0,0,1000].*1e-06;
image=[0,0,200].*1e-06;
direc='f';
material_sub='FS';
design_wavelength=405e-09;

fprintf('Initiating sdtDOE object\n')

D1=sdtDOE(r,N,object,image,direc,material_sub,design_wavelength); % initiate sdtDOE object

fprintf('Generating phase profile from point to point model\n')

D1=D1.gen_p2p_phase_profile();

fprintf('Plotting phase profile\n')

D1.plot_phase_profile('Phase profile')
drawnow

p=gcp('nocreate');

if isempty(p)
    fprintf('Starting parallel pool for simuations\n')
    parpool;
end

fprintf('Benchmarking your computer to estimate simulation time\n')
D1=D1.sim_benchmark();

fprintf('Preparing to simulate the focal plane, this is a good time to make a cup of tea \n')
D1=D1.simulation2D(25,1,3.75e-06,D1.design_wavelength,D1.object_coords(3),7e-06,500e-06);

fprintf('Fitting Gaussian to simulated beam profile\n')
fit=sdtDOE.gaussfit2D(D1.sim_intensity,(D1.sim_x(2) - D1.sim_x(1)),5,1);

fprintf('Do not worry, the matrix is only singular because the error on the fit is zero as we are fitting to an exact Gaussian\n')
fprintf('Waist at focal plane: %.3f um\n',2*fit.sigma_major.*1e06)

figure
imagesc(D1.sim_x.*1e06,D1.sim_y.*1e06,D1.sim_intensity(:,:,1))
axis square
xlabel('x [microns]')
ylabel('y [microns]')
set(gcf,'color','white')
end

