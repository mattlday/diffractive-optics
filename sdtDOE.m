classdef sdtDOE
    %sdtDOE is a class for generating, simulating and analysing diffractive
    %optical elements based on scalar diffraction theory. It is a
    %superclass to surfacerelief, which uses sdtDOE for design of its phase
    %profiles and analysis of measured data. All simulations are currently
    %based on Gaussian beams, future versions plan to include general
    %electric field intensity input.
    %
    % The units of the code are in SI, plotting is displayed in microns -
    % however the code applies to any scale of DOE (be careful about going
    % subwavelength though).
    %
    %  Methods:
    %
    %  Description of methods, for more information on any method type use
    %  "help sdtDOE.method"
    %
    %  CONSTRUCTOR
    %
    %  obj=sdtDOE(radius, profile_resolution, object_coords, image_coords,
    %  direction, material_substrate, material_DOE, design_wavelength)
    %
    %  INPUT:
    %  radius - the radius of the designed DOE
    %  profile_resolution - points per unit length for generating profile
    %  object_coords - [x1,y1,z1]
    %  image_coords - [x2,y2,z2]
    %  direction - 'f' (forward) or 'b' (backward), 'f' places object point
    %  inside substrate, 'b' places image point inside substrate
    %  material_substrate - string, corresponding to dictionary reference
    %  to material, to find out included materials, use the list_materials
    %  method. This is the material of the substrate the DOE is placed on.
    %  material_DOE - material the DOE is fabricated in
    %  design_wavelength - the wavelength of the light source being used,
    %  keep units consistent
    %  OUPUT:
    %  obj - a data structure for containing fields for all information about a DOE
    %  from design to measurement analysis
    %
    %  OPTIMISATION
    %
    %  obj=optimisation_1D(obj,op_source_w0) - uses the inputted design
    %  image and object point to design a DOE that focuses a Gaussian beam
    %  at the correct image point - compensating for the finite size of
    %  the source
    %
    %  PHASE PROFILE GENERATION
    %
    %  obj=gen_p2p_phase_profile(obj) - generates a point-to-point (p2p) phase profile that
    %  converts a point source at the object point to a point at the
    %  image point.
    %
    %  obj=gen_binarygrating_phase_profile(obj, pitch, contrast) -
    %  generates a binary grating phase profile
    %
    %  obj=gen_blazedgrating_phase_profile(obj, pitch, contrast) -
    %  generates a blazed grating phase profile
    %
    %  obj=collimate_output(obj,source_w0,image_z_coords_limits,number_scan_points,number_of_col_samples,collimation_range)
    %  - generates a phase profile that collimates a specified source beam
    %
    %  obj=gsip(obj,desired_image,hologram_radius,source_w0,iterations) -
    %  generates a phase profile using the Gerchberg-Saxton iterative
    %  procedure to generate a desired image in the focal plane of the DOE
    %
    %  SIMULATION
    %  [obj,t]=simulation1D(obj,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,Z)
    %  - simulates a single intensity point along the optical axis at each specified Z value
    %
    %  [obj,t]=simulation2D(obj,sim_axis,sim_no_points,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,sim_radius,Z,sim_centre,benchmark_bool)
    %  - Simulates beam profiles in one axis at different Z positions. If
    %  only using one z position, the simulation plane can be shifted to
    %  probe off-axis behaviour
    %
    %  [obj,t]=simulation3D(obj,sim_no_points,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,sim_radius,Z,sim_centre,benchmark_bool)
    %  - same as simulation2D but simulates 2D plane at each z value. For
    %  high res simulations, an HPC is required.
    %
    %  obj=sim_benchmark(obj) - determines the number of integrations per
    %  second the current computer can perform, useful for predicting HPC
    %  time or avoiding long simulations on a PC
    %
    %  time=sim_howlong(obj,sim_no_points,Z) - gives an estimate for how
    %  long a simulation of a sim_no_points by sim_no_points at length(Z)
    %  number of positions
    %
    %  I=quick_sim(obj, type, sim_source_amplitude, sim_source_w0,
    %  sim_source_wavelength, sim_source_dist_from_DOE,Z) - a simulation of a single intensity
    %  sample from each z position allows for quick DOE behaviour to be
    %  determined
    %
    %  IMPORTING MEASURED DATA
    %  sdtDOE is also a container for measured data, both of the field
    %  propgation and of the physical profiles
    %
    %  obj=import_and_process_data(obj,BG_file,keyphrase,pixel_size,magnification,
    %  measured_source_dist_from_DOE,
    %  measured_source_wavelength,measured_source_w0,start_z,z_increments)
    %  - RECOMMENDED, imports and processes measured data without storing
    %  all data in the object as this leads to unstable file sizes
    %
    %  obj=import_measured_beamprofiles(obj,keyphrase,pixel_size,magnifcation,measured_source_dist_from_DOE,measured_source_wavelenth,measured_source_w0,start_z,z_increments)
    %  - NOT RECOMMENDED, LARGE FILE SIZES AND UNSTABLE, data importer designed to import files ending in _xx where x are
    %  numerica digits. The importer sorts the numbers naturally from 01 to
    %  99, corresponding to the different z positions
    %
    %  obj=import_DE_beamprofiles(obj,filename_source,filename_DOE) -
    %  Diffraction efficiences are measured as the ratio of power in the
    %  source and output beams, this importer prepares the two beam profiles
    %  from file
    %
    % obj=import_BG(obj,background_filename) - imports a background file
    % for use with background_subtraction. Subtracting background off data
    % dramatically speeds up and improves quality of data processing.
    %
    % obj=transfer_data_from_object(obj,sdtDOE_obj,datatype) - for
    % transfering simulation data from one sdtDOE object to another.
    %
    %  ANALYSIS
    %
    %  obj=background_subtraction(obj,type) - subtracts background from
    %  measured data, highly recommended before calculating the diffraction
    %  efficiency.
    %
    %  obj=calc_diffraction_efficiency(obj) - calculates the diffraction
    %  efficency of the DOE provided the source and output beam profiles
    %  have been imported
    %
    %  obj=process_measured_beamprofiles(obj) - USE IMPORT_AND_PROCESS_DATA INSTEAD, fits a 2D Gaussian to each
    %  imported beam profile and returns vectors of Gaussian parameters for
    %  easy plotting (can be time consuming, includes progress bar)
    %
    %  obj=process_sim_beamprofiles(obj) - same as
    %  process_measured_beamprofiles but works on simulated beam profiles,
    %  in the future these functions should be COMBINED
    %
    %  I=extract_central_intensity(obj) - extracts the intensity along the
    %  optical axis of measured data for easy comparison to quick_sim
    %
    %  PLOTTING
    %  Some custom plotting functions have been written to speed up
    %  generating publication quality plots
    %
    %  plot_geometry(obj) - displays the DOE geometry for visualisation
    %
    %  plot_animated_propagation(obj,type,filename,titlestr) - ONLY WORKS
    %  WHEN FULL DATA IS IMPORTED, WILL BE FIXED IN LATER RELEASE
    %
    %  plot_propagation(obj,type,axis_select) - a universal function for
    %  visualising either design simulation ('sim'), fabricated simulation
    %  ('fab'), or measured ('measured') data as a 2D stich of 1D beam
    %  profile slices along either the 'x' or 'y' axis
    %
    %  plot_phase_profile(obj,titlestr) - plots the phae profile
    %
    %  plot_waist_compare(obj) - plots all the present beam radius data for
    %  comparison
    %
    %  STATIC
    %
    %  fit=gaussfit2D(A,ps,threshold,plot_bool) - fits a 2D gaussian to the
    %  profile A
    %
    %  I_pin=pinhole(I,ps,n) - places a pinhole of radius n times the beam
    %  radius of the beam profile
    %
    %  I_pin=pinhole_radius(I,ps,r) - same as pinhole but can place a
    %  pinhole of arbitrary radius
    %
    %  u=gauss(A,lambda,w0,x,y,z) - generates a Gaussian profile, used for
    %  simulation
    %
    %  I=hdr_profile(Iset,Rset,exposure_times,ps) - stitches together a set
    %  of exposures made at different exposure times
    %
    %  r=RMSD(data,ideal) - root-mean square deviation calculator
    %
    %  w=waist(z,w0,lambda,material) - the gaussian beam radius equation
    %
    %---------------------------------------------------------------
    %  AUTHOR: Matthew Day
    %  AFFILIATIONS: Quantum Engineering Centre for Doctoral, Training,
    %  University of Bristol, UK and The National Physical Laboratory, UK
    %  Date of current version release: 19/07/2017
    %---------------------------------------------------------------
    
    properties
        
        radius
        profile_resolution
        object_coords
        image_coords
        op_image_coords
        direction
        material_substrate
        design_wavelength
        phase_profile
        profile_x
        profile_y
        op_source_wavelength
        op_source_w0
        op_source_amplitude
        op_source_dist_from_DOE
        sim_source_dist_from_DOE
        sim_source_wavelength
        sim_source_w0
        sim_source_amplitude
        sim_x
        sim_y
        sim_z
        sim_intensity
        sim_approximate_beamradius
        sim_majorbeamwaist
        sim_minorbeamwaist
        sim_majorbeamwaist_95conf
        sim_minorbeamwaist_95conf
        sim_x0
        sim_y0
        sim_fullfitdata
        benchmark
        measured_source_dist_from_DOE
        measured_source_wavelength
        measured_source_w0
        measured_x
        measured_y
        measured_z
        CCD_pixelsize
        magnification
        measured_intensity
        measured_intensity_central
        measured_intensity_xslice
        measured_intensity_yslice
        measured_majorbeamwaist
        measured_minorbeamwaist
        measured_majorbeamwaist_95conf
        measured_minorbeamwaist_95conf
        measured_x0
        measured_y0
        measured_fullfitdata
        measured_beamprofile_source_w0
        measured_beamprofile_DOE_w0
        background
        diffraction_efficiency
    end
    
    methods
        
        %% Constructor
        
        function obj=sdtDOE(radius,profile_resolution,object_coords,image_coords,direction,material_substrate,design_wavelength)
            obj.radius=radius;
            obj.profile_resolution=profile_resolution;
            obj.object_coords=object_coords;
            obj.image_coords=image_coords;
            obj.op_image_coords=image_coords;
            obj.direction=direction;
            obj.material_substrate=material_substrate;
            obj.design_wavelength=design_wavelength;
            obj.profile_x=linspace(-radius,radius,int32(2*profile_resolution*radius*1e06));
            obj.profile_y= obj.profile_x;
        end
        
        %% Finite source 1D optimiser
        
        function obj=optimisation_1D(obj,op_source_w0)
            % OPTIMISATION_1D is an interactive method that guides the user
            % through a procedure to find the optimal point to point DOE image point
            % that correctly focusses a Gaussian beam to a desired Gaussian
            % image point.
            % INPUTS:
            % obj - the sdtDOE instance to optimise
            % op_source_w0 - the minimum waist radius of the Gaussian
            % source, where the minimum waist radius is placed at the
            % designed object point
            % OUTPUTS:
            % obj - containing optimised design parameters that are then
            % automatically used throughout the design procedure (i.e.
            % generating phase profiles)
            %
            % obj=optimisation_1D(obj,op_source_w0)
            
            obj.op_source_w0=op_source_w0;
            obj.op_source_amplitude=1;
            obj.op_source_dist_from_DOE=obj.object_coords(3);
            obj.op_source_wavelength=obj.design_wavelength;
            
            % ri=obj.refractive_index_substrate;
            d=obj.direction;
            lambda=obj.design_wavelength;
            x1=obj.object_coords(1); y1=obj.object_coords(2); z1=obj.object_coords(3);
            x2=obj.image_coords(1); y2=obj.image_coords(2); z2=obj.image_coords(3)+10e-06;
            w=op_source_w0;
            A=obj.op_source_amplitude;
            n=refractive_index(obj.material_substrate,lambda);
            
            function f = integral(d2)
                %% INTEGRAL First calculates the diffraction integral for an
                %%inputted d2, squares it to get the intensity and then outputs the
                %%inverse of this intensity. The code is based on
                %%p2p_DOE_xaxisanalysis.m, refer there for detailed comments.
                %%Differences will be highlighted below.
                %% Calculate the reduced wavelength inside the material
                lambda_s=lambda/n;
                
                %% Assign whether the substrate is before or after the DOE
                if obj.direction=='f'
                    Lambda1=lambda_s;
                    Lambda2=lambda;
                elseif obj.direction=='b'
                    Lambda1=lambda;
                    Lambda2=lambda_s;
                end
                
                if d=='f'
                    n1=n;
                    n2=1;
                    
                elseif d=='b'
                    n1=1;
                    n2=n;
                    
                end
                
                
                if d=='f' && z1<(n1/n2)*d2
                    k=0;
                elseif d=='f' && z1>=(n1/n2)*d2
                    k=1;
                elseif d=='b' && z1<(n1/n2)*d2
                    k=0;
                elseif d=='b' && z1>=(n1/n2)*d2
                    k=1;
                end
                
                %% DOE aperture coordinate space
                
                x=obj.profile_x;
                y=obj.profile_y;
                [X,Y]=meshgrid(x,y);
                
                %% Transmittance function
                
                Phi=(-1)^k * ((2*pi*n2)/lambda *sqrt((((d2*x2)/z2) - X).^2+((d2*y2)/z2-Y).^2+d2.^2) - (2*pi*n1)/lambda *sqrt((x1-X).^2+(y1-Y).^2+z1.^2));
                
                T=exp(1i.*Phi);
                
                
                %% Field before DOE
                
                E0=sdtDOE.gauss(A,Lambda1,w,X,Y,z1); % E0 stands for intial field
                
                
                %% Aperture shadow
                
                Ea=E0.*T;
                
                
                %% Solve the diffraction integral
                
                
                P=exp(1i.*(2*pi)./Lambda2.* sqrt(z2.^2+ (X-x2).^2 + (Y-y2).^2))./(z2.^2+ (X-x2).^2+ (Y-y2).^2);
                integrand=Ea.*P;
                E=z2/(1i.*Lambda2)*trapz(y,trapz(x,integrand));
                
                Intensity=abs(E).^2;
                
                %f=1/Intensity;
                f=-1*Intensity;
            end
            
            
            function d0=dinitial(z1,z2,w,lambda)
                %% Calculates an initial guess for d2 using theory of focussed Gaussian beams
                d0= 1/(1/z2 + ((pi*w^2)/(z1*lambda))^2/(z1*(((pi*w^2)/(z1*lambda))^2 + 1)));
            end
            
            bool2='y'; % Sets a choice the user will make later to 'yes' such that the while loop initiates
            i=1; % initialise a loop index
            d0=dinitial(z1,z2,w,lambda); %call the initial d2 value
            
            while bool2=='y' %while the user wishes to perform optimisations
                
                if i>1 %for the first run we use calculated d0, for later runs we provide the option to change d0
                    modify=input('Enter amount to change d0 by: ');
                    d0=d0+modify;
                end
                
                %[d2,fval]=fminsearch(@integral,d0,optimset('TolFun',1e-7,'TolX',1e-7)); % perform the minimisation routine using d0 as a starting point
                
                [d2,fval]=fminsearch(@integral,d0,optimset('TolFun',1e-5,'TolX',1e-5,'Display','iter','PlotFcns',@optimplotfval));
                
                xd2=(d2*x2)/z2; %define xd2 from found d2
                yd2=(d2*y2)/z2; % define yd2 from found d2
                
                % Report results to user
                result1=sprintf('Found (xd2,yd2,d2)=( %.3f, %.3f, %.3f) um with maximal intensity at desired focus of %.3f from an initial d0=%.3f um.',xd2.*1e06,yd2.*1e06,d2.*1e06,-fval,d0.*1e06);
                disp(result1)
                
                % Provide the option to view the optimisation landscape, that is to
                % plot the result of integral(d2) for a range of d2 and display the
                % found minimum.
                bool1=input('Would you like to display the optimisation landscape? [y/n]: ','s');
                
                if bool1=='y'
                    Z=linspace(d2*0.5,d2*2,100);
                    landscape=zeros(1,length(Z));
                    for i=1:length(Z)
                        landscape(i)=integral(Z(i));
                    end
                    figure
                    plot(Z,landscape);
                    xlabel('DOE design image point, d2 [um]');
                    ylabel('Negative intensity at desired z2 [a.u.]');
                    set(gcf,'Color','white');
                    hold on
                    h2=plot(d2,fval,'rx');
                    legend(h2,'Found minima');
                end
                
                bool2=input('Would you like to retry optimisation? [y/n]: ','s');
                
                i=i+1;
            end
            
            
            obj.op_image_coords=[xd2,yd2,d2];
            obj=obj.gen_p2p_phase_profile;
            
            bool3=input('Would you like to quickly simulate the output field of the optimised DOE?[y/n]: ','s');
            
            if bool3=='y'
                %Zsim=((3/2*obj.image_coords(3))/100):((3/2*obj.image_coords(3))/100):(3/2*obj.image_coords(3));
                %%% FIX THIS
                Zsim=(10:10:800).*1e-06; % THIS WORKS
                obj=obj.simulation1D(1,obj.op_source_w0,obj.design_wavelength,obj.op_source_dist_from_DOE,Zsim);
                %obj.quick_sim('design',1,3.75e-06,405e-09,obj.op_source_dist_from_DOE,(10:10:800).*1e-06)
                
                
            end
            
            
            
        end
        
        
        %% Phase profile generation
        
        function obj=gen_p2p_phase_profile(obj)
            %GEN_P2P_PHASE_PROFILE generates a phase profile based on a
            %point to point model such that the resulting DOE converts the
            %plane waves from a point source to another set of plane ways
            %converging to a different point
            %INPUTS:
            % obj - initiated sdtDOE object
            %OUTPUTS:
            % obj - sdtDOE object containing phase profile in
            % obj.phase_profile
            %
            % obj=gen_p2p_phase_profile(obj)
            
            %% Variable shorthand for object properties
            d=obj.direction;
            lambda=obj.design_wavelength;
            x1=obj.object_coords(1); y1=obj.object_coords(2); z1=obj.object_coords(3);
            x2=obj.op_image_coords(1); y2=obj.op_image_coords(2); z2=obj.op_image_coords(3);
            
            %% Calculate the refractive index of the substrate
            n=refractive_index(obj.material_substrate,obj.design_wavelength);
            
            %% If the direction is forwards, make the first refractive index that of the substrate and vice versa
            if d=='f'
                n1=n;
                n2=1;
                
            elseif d=='b'
                n1=1;
                n2=n;
            end
            
            %% Conditions on keeping the correct sign of the phase profile
            if d=='f' && z1<(n1/n2)*z2
                k=0;
            elseif d=='f' && z1>=(n1/n2)*z2
                k=1;
            elseif d=='b' && z1<(n1/n2)*z2
                k=0;
            elseif d=='b' && z1>=(n1/n2)*z2
                k=1;
            end
            
            %% The functional form of the point to point phase function
            function phi=phase_function(x,y)
                phi=(-1)^k * ((2*pi*n2)/lambda *sqrt((x-x2).^2+(y-y2).^2+z2.^2) - (2*pi*n1)/lambda *sqrt((x-x1).^2+(y-y1).^2+z1.^2));
            end
            
            %% Generate coordinates
            [X,Y]=meshgrid(obj.profile_x,obj.profile_y);
            
            %% Generate phase profile
            obj.phase_profile=phase_function(X,Y);
            
        end
        
        function obj=gen_binarygrating_phase_profile(obj, pitch, contrast)
            % GEN_BINARYGRATING_PHASE_PROFILE generates a phase profile of
            % a binary grating
            % INPUTS
            % obj - sdtDOE object
            % pitch - width of grating element
            % contrast - height of grating elements as fraction of 2*pi
            % OUTPUTS
            % obj - containing filled obj.phase_profile
            %
            % obj=gen_binarygrating_phase_profile(obj, pitch, contrast)
            n_element=int32((pitch*1e6)*obj.profile_resolution);
            
            m=floor(length(obj.profile_x)/(2*n_element));
            
            n1=1;
            n2=n_element;
            
            
            grating=zeros(1,length(obj.profile_x));
            
            for i=1:m
                
                grating(n1:n2)=ones(1,n_element);
                n3=n2+1;
                n4=n3+n_element-1;
                grating(n3:n4)=zeros(1,n_element);
                n1=n4+1;
                n2=n1+n_element-1;
                
            end
            
            %grating((end+1):(end+mod(length(obj.profile_x),n_element)-1))=ones(1,mod(length(obj.profile_x),n_element));
            
            for i=1:length(obj.profile_y)
                
                G(i,:)=grating;
                
            end
            
            
            
            G=(2*pi*contrast.*G)-2*pi;
            
            obj.phase_profile=G;
        end
        
        function obj=gen_blazedgrating_phase_profile(obj, pitch, contrast)
            % GEN_BINARYGRATING_PHASE_PROFILE generates a phase profile of
            % a blazed grating
            % INPUTS
            % obj - sdtDOE object
            % pitch - width of grating element
            % contrast - max height of grating elements as fraction of 2*pi
            % OUTPUTS
            % obj - containing filled obj.phase_profile
            %
            % obj=gen_binarygrating_phase_profile(obj, pitch, contrast)
            
            
            n_element=int32((pitch*1e6)*obj.profile_resolution);
            
            m=floor(length(obj.profile_x)/(n_element));
            
            
            element=linspace(contrast*2*pi,0,n_element);
            
            grating=repmat(element,1,m);
            
            mod(length(obj.profile_x),n_element);
            
            grating((end+1):(end+mod(length(obj.profile_x),n_element)))=element(1:mod(length(obj.profile_x),n_element));
            
           
            
            if length(grating)<length(obj.profile_x)
               
                grating((end+1):length(grating)+(length(grating)-length(obj.profile_x)))=zeros(1,(length(grating)-length(obj.profile_x)));
            elseif length(grating)>length(obj.profile_x)    
            
                grating=grating(1:length(obj.profile_x));
            end

            
            
            for i=1:length(obj.profile_x)
                
                G(i,:)=grating;
                
            end
            
            
            
            
            obj.phase_profile=G;
        end
        
         function obj=collimate_output(obj,source_w0,image_z_coords_limits,number_scan_points,number_of_col_samples,collimation_range)
            % COLLIMATE_OUTPUT generates a phase profile that takes a
            % diverging Gaussian beam and collimates it in the direction of
            % the design image points of the DOE. The code is
            % computationally intensive.
            %INPUTS:
            % obj - initiated sdtDOE object
            % source_w0 - beam waist of source Gaussian beam
            % image_z_coords_limits - the range of values the optimised
            % value of the z coordinate image point (trial and error,
            % smaller the range the better)
            % number_scan_points - number of points simulated to generate
            % FOM, higher the better but computational cost, can be as low
            % as three to narrow down z_coord_limits.
            % number_of_col_samples - number of z_coords to try in one scan
            % collimation range - two z coordinates over which collimation
            % is desired [start_z, end_z]
            %OUTPUTS:
            % obj - sdtDOE object containing phase profile in
            % obj.phase_profile
            %
            % obj=collimate_output(obj,source_w0,image_z_coords_limits,number_scan_points,number_of_col_samples,collimation_range)
            
            
            %% The possible design image coordinates that would collimate the beam
            z2=linspace(image_z_coords_limits(1),image_z_coords_limits(2),number_scan_points);
            
            
            Z=linspace(collimation_range(1),collimation_range(2),number_of_col_samples);
            
            disp('Performing coarse scan of design parameters')
            
            for j=1:number_scan_points
                obj.op_image_coords(1)=obj.image_coords(1)/obj.image_coords(3).*z2(j);
                obj.op_image_coords(2)=obj.image_coords(2)/obj.image_coords(3).*z2(j);
                obj.op_image_coords(3)=z2(j);
                obj=obj.gen_p2p_phase_profile;
                obj=obj.simulation1D(1,source_w0,obj.design_wavelength,obj.object_coords(3),Z);
                R(j)=sdtDOE.RMSD(obj.sim_intensity./obj.sim_intensity(1),ones(1,number_of_col_samples));
                
            end
            
            figure
            plot(z2.*1e6,R)
            
            disp('Initiating fine scan of design parameters')
            
            [~,idx]=min(R);
            
            if idx<2
                idx_lower=1;
            else
                idx_lower=idx-2;
            end
            
            if (idx+2)>length(R)
                idx_higher=length(R);
            else
                idx_higher=idx+2;
            end
            
            z2fine=linspace(z2(idx_lower),z2(idx_higher),number_scan_points);
            
            disp('Performing fine scan of design parameters')
            
            for j=1:number_scan_points
                obj.op_image_coords(1)=obj.image_coords(1)/obj.image_coords(3).*z2fine(j);
                obj.op_image_coords(2)=obj.image_coords(2)/obj.image_coords(3).*z2fine(j);
                obj.op_image_coords(3)=z2fine(j);
                obj=obj.gen_p2p_phase_profile;
                obj=obj.simulation1D(1,source_w0,obj.design_wavelength,obj.object_coords(3),Z);
                Rfine(j)=sdtDOE.RMSD(obj.sim_intensity./obj.sim_intensity(1),ones(1,number_of_col_samples));
                
            end
            
            [~,idx2]=min(Rfine);
            
            figure
            plot(z2fine.*1e6,Rfine)
            
            fprintf('The design image point is %.2f microns \n',z2fine(idx2).*1e6)
            
            obj.op_image_coords(1)=obj.image_coords(1)/obj.image_coords(3).*z2fine(idx2);
            obj.op_image_coords(2)=obj.image_coords(2)/obj.image_coords(3).*z2fine(idx2);
            obj.op_image_coords(3)=z2fine(idx2);
            
         end
        
         function obj=gsip(obj,desired_image,hologram_radius,source_w0,iterations)
             % GSIP performs the Gerchberg-Saxton iterative procedure to
             % iteratively design a DOE that generates a desired image in
             % the focal plane. As the method generates a design in the
             % focal plane and not the far field it cannot use fft, and so
             % it is computationally intensive. A seed phase profile must
             % already have been generated that creates a Gaussian beam of
             % approximately the right size in the plane of interest.
      
             rh=hologram_radius;
             lambda=obj.design_wavelength;
             dx=lambda/2;
             zB=obj.object_coords(3);
             L=obj.op_image_coords(3);
            
             %DOE plane
             x=obj.profile_x;
             y=obj.profile_y;
             [X,Y]=meshgrid(x,y);
             
             %Hologram plane
             u=-rh:dx:rh;
             v=fliplr(u);
             [U,V]=meshgrid(u,v);
             I0=desired_image;
             I0=imresize(I0,[length(u),length(u)]);
            
             %% Generate field before the DOE
             
             E0=sdtDOE.gauss(1,lambda/refractive_index(obj.material_substrate,lambda),source_w0,X,Y,zB);
             A0=abs(E0);
             
             figure
             imagesc(x,y,abs(E0).^2)
             title('Input beam')
             axis equal tight
             set(gcf,'color','white')
             drawnow
             %% Initialise phase hologram array
             
             Ph0=obj.phase_profile;
             Ph0=imresize(Ph0,[length(x),length(x)]);
             
             %% Initial transmittance function
             
             T=exp(1i*Ph0);
             E1=E0.*T;
             
             %% While loop starts
             
             for fidx=1:iterations
                 
                 %% Modify E1 to more resemble input beam 
                 
                 if fidx>1
                     E1=A0.*(E1./abs(E1));
                 end
                 
                 %% Propagate to hologram plane
                 tic
                 E2=zeros(length(v),length(u));
                 for i=1:length(u)
                     progress=(i/length(u))*100;
                     fprintf('Construct hologram plane: %.2f percent complete\n',progress)
                     parfor j=1:length(v)
                         Prop=exp(1i.*(2*pi)./lambda.* sqrt(L.^2+ (X-u(i)).^2 + (Y-v(j)).^2))./(L.^2+ (X-u(i)).^2+ (Y-v(j)).^2);
                         integrand=E1.*Prop;
                         E2(j,i)=L/(1i.*lambda)*trapz(y,trapz(x,integrand));
                     end
                 end
                 t=toc
                 
                 figure
                 imagesc(u,v,abs(E2).^2)
                 title('Output beam, focal plane')
                 axis equal tight
                 set(gcf,'color','white')
                 set(gca,'ydir','normal')
                 drawnow
                 
                 E2prop=E2;
                 
                 %% Error function goes here
                 difference=(abs(E2)-sqrt(I0)).^2;
                 epsilon(fidx)=trapz(u,trapz(v,difference));
                 
                 %% Modify E2 to more closely resemble I0
                 
                 E2=sqrt(I0).*(E2./abs(E2));
                 
                 % figure
                 % imagesc(u,v,abs(E2).^2)
                 % title('Modified beam for back propagation')
                 % axis equal
                 % set(gca,'ydir','normal')
                 % drawnow
                 
                 
                 %% Propagate back to DOE plane
                 
                 Prop=[];
                 
                 for i=1:length(x)
                     progress=(i/length(x))*100;
                     fprintf('Reconstruct DOE plane: %.2f percent complete\n',progress)
                     parfor j=1:length(y)
                         Prop=exp(-1i.*(2*pi)./lambda.* sqrt(L.^2 + (U-x(i)).^2 + (V-y(j)).^2))./(L.^2 + (U-x(i)).^2 + (V-y(j)).^2);
                         integrand=E2.*Prop;
                         E1(j,i)=L/(1i.*lambda)*trapz(u,trapz(v,integrand));
                     end
                 end
                 
                 % figure
                 % imagesc(x,y,abs(E1).^2)
                 % axis equal
                 % set(gca,'ydir','normal')
                 % title('DOE plane')
                 
                 Phi_mod=angle(E1./E0);
                 
                 % figure
                 % imagesc(x,y,angle(E1./E0))
                 % axis equal
                 % set(gca,'ydir','normal')
                 % title('Required phase profile')
                 
                
                 
             end
             
             obj.phase_profile=Phi_mod;
             
             figure
             scatter(1:length(epsilon),epsilon)
             
             
         end
        
        %% Simulation
        
        function [obj,t]=simulation1D(obj,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,Z)
            %SIMULATION1D simulates the 1D beam intensity along the optical axis as defined by the DOE design coordinates
            %INPUTS:
            % obj - a sdtDOE class object with a generated phase profile
            % sim_source_amplitude - amplitude of the source Gaussian beam
            % sim_source_w0 - the beam waist of the source beam
            % sim_source_wavelength - the wavelength of the
            % source beam in vacuum
            % sim_source_dist_from_DOE - the distance between the source
            % beam waist and the DOE, typically set at the design
            % object point
            % Z - a vector of z positions
            %OUTPUTS:
            % obj - the sdtDOE object containing simulation data in
            % obj.sim_intensity
            % t - the length of time the simulation took, for benchmarking
            % purposes
            %
            % [obj,t]=simulation2D(obj,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,Z)
            
            
            if isempty(obj.phase_profile)==1
                phase_bool=input('No phase profile has been generated, would you like to generate a basic point to point phase profile? [y/n]:','s');
                if strcmp(phase_bool,'y')
                    obj=obj.gen_p2p_phase_profile;
                    fprintf('Phase profile generated, continuing with simulation\n')
                elseif strcmp(phase_bool,'n')
                    error('No simulation can continue without a phase profile, please generate a suitable profile and retry simulation')
                end
            end
            
            %% Store inputs in the sdtDOE object for the record
            obj.sim_source_amplitude=sim_source_amplitude;
            obj.sim_source_w0=sim_source_w0;
            obj.sim_source_wavelength=sim_source_wavelength;
            obj.sim_source_dist_from_DOE=sim_source_dist_from_DOE;
            obj.sim_z=Z;
            
            %% Defining shorthand variable name
            lambda=sim_source_wavelength;
            
            %% Calculate the refractive index of the material
            n=refractive_index(obj.material_substrate,lambda);
            
            %% Calculate the reduced wavelength inside the material
            lambda_s=lambda/n;
            
            %% Assign whether the substrate is before or after the DOE
            if obj.direction=='f'
                Lambda1=lambda_s;
                Lambda2=lambda;
            elseif obj.direction=='b'
                Lambda1=lambda;
                Lambda2=lambda_s;
            end
            
            %% Transmittance function
            
            T=exp(1i.*obj.phase_profile);
            
            
            %% DOE aperture coordinate space
            
            x=obj.profile_x;
            y=obj.profile_y;
            [X,Y]=meshgrid(x,y);
            
            %% Output field coordinate space
            
            p=obj.image_coords(1)/obj.image_coords(3).*Z;
            q=obj.image_coords(2)/obj.image_coords(3).*Z;
            
            
            obj.sim_x=p;
            obj.sim_y=q;
            
            
            %% Field before DOE
            
            E0=sdtDOE.gauss(sim_source_amplitude,Lambda1,sim_source_w0,X,Y,sim_source_dist_from_DOE); % E0 stands for intial field
            %E0=ones(length(x))./((X.*1e6).^2 + (Y.*1e6).^2 +(sim_source_dist_from_DOE.*1e6).^2);
            
            %% Aperture shadow
            
            Ea=E0.*T;
            
            %% Output field
            
            tic
            %% Solve the diffraction integral for each p,q,z in the output field
            for k=1:length(Z)
                z=Z(k); %to avoid broadcast variable
                
                P=exp(1i.*(2*pi)./Lambda2.* sqrt(z.^2+ (X-p(k)).^2 + (Y-q(k)).^2))./(z.^2+ (X-p(k)).^2+ (Y-q(k)).^2);
                integrand=Ea.*P;
                E(k)=z/(1i.*Lambda2)*trapz(y,trapz(x,integrand));
            end
            
            I=abs(E).^2;
            t=toc;
            
            obj.sim_intensity=I;
            obj.sim_approximate_beamradius= sdtDOE.waist(sim_source_dist_from_DOE,sim_source_w0,sim_source_wavelength,obj.material_substrate) *sqrt(max(max(abs(E0).^2))./I);
            
        end
        
        function [obj,t]=simulation2D(obj,sim_axis,sim_no_points,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,sim_radius,Z,sim_centre,benchmark_bool)
            %SIMULATION2D simulates the 1D beam profile at a range of z
            %positions after the DOE given a Gaussian source as input
            %INPUTS:
            % obj - a sdtDOE class object with a generated phase profile
            % sim_no_points - the dimension of the output array (i.e.
            % sim_no_points=20 simulates a 20x20 array of points at each z
            % position)
            % sim_source_amplitude - amplitude of the source Gaussian beam,
            % simulation optimised for a value of 1
            % sim_source_w0 - the minimum beam waist of the source beam
            % sim_source_wavelength - the (monochromatic) wavelength of the
            % source beam
            % sim_source_dist_from_DOE - the distance between the source
            % minimum beam waist and the DOE, typically set at the design
            % object point
            % sim_radius - the physical radius of the simulation plane, set
            % as small as possible while still containing the output beam
            % at its greatest extent
            % Z - a vector of z positions
            % benchmark_bool - not used by the user, can safely ignore
            %OUTPUTS:
            % obj - the sdtDOE object containing simulation data in
            % obj.sim_intensity
            % t - the length of time the simulation took, for benchmarking
            % purposes
            %
            % [obj,t]=simulation3D(obj,sim_no_points,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,sim_radius,Z,benchmark_bool)
            
            %% if the user has ignored the benchmark_bool input, set it to the input that generates the correct behaviour
            if nargin<9
                benchmark_bool=0;
            end
            
            %% if the benchmark procedure has been run, then estimate simulation time, if not, suggest running it
            if isempty(obj.benchmark)==0
                fprintf('Estimated simulation time is %.4f seconds (%.3f minutes, or %.2f hours)\n',obj.benchmark*(sim_no_points^2)*length(Z),(obj.benchmark*(sim_no_points^2)*length(Z))/60, (obj.benchmark*(sim_no_points^2)*length(Z))/3600)
            elseif benchmark_bool==0
                %fprintf('Cannot estimate simulation time, to enable this feature please use the sim_benchmark method\n')
            else
                % do nothing (not a critical feature)
            end
            
            if isempty(obj.phase_profile)==1
                phase_bool=input('No phase profile has been generated, would you like to generate a basic point to point phase profile? [y/n]:','s');
                if strcmp(phase_bool,'y')
                    obj=obj.gen_p2p_phase_profile;
                    fprintf('Phase profile generated, continuing with simulation\n')
                elseif strcmp(phase_bool,'n')
                    error('No simulation can continue without a phase profile, please generate a suitable profile and retry simulation')
                end
            end
            
            %% Store inputs in the sdtDOE object for the record
            obj.sim_source_amplitude=sim_source_amplitude;
            obj.sim_source_w0=sim_source_w0;
            obj.sim_source_wavelength=sim_source_wavelength;
            obj.sim_source_dist_from_DOE=sim_source_dist_from_DOE;
            obj.sim_z=Z;
            
            %% Defining shorthand variable name
            lambda=sim_source_wavelength;
            
            %% Calculate the refractive index of the material
            n=refractive_index(obj.material_substrate,lambda);
            
            %% Calculate the reduced wavelength inside the material
            lambda_s=lambda/n;
            
            %% Assign whether the substrate is before or after the DOE
            if obj.direction=='f'
                Lambda1=lambda_s;
                Lambda2=lambda;
            elseif obj.direction=='b'
                Lambda1=lambda;
                Lambda2=lambda_s;
            end
            
            
            
            
            %% DOE aperture coordinate space
            
            res=sim_no_points;
            
            p=linspace(-sim_radius,sim_radius,res); obj.sim_x=p;
            q=fliplr(p); obj.sim_y=q;
            
            if strcmp(sim_axis,'x')==1
                
                x=obj.profile_x;
                y=0;
              
                
                int_dim=x;
                
                p=linspace(-sim_radius,sim_radius,res)+sim_centre(1); obj.sim_x=p;
                q=0+sim_centre(2); obj.sim_y=q;
                
                T=exp(1i.*middle_row(obj.phase_profile));
               
                
            elseif strcmp(sim_axis,'y')==1
                x=0;
                y=obj.profile_y;
                
                int_dim=y;
                
                p=sim_centre(1); obj.sim_x=p;
                q=fliplr(linspace(-sim_radius,sim_radius,res))+sim_centre(2); obj.sim_y=q;
                
                T=exp(1i.*middle_col(obj.phase_profile));
                
            else
                x=obj.profile_x;
                y=0;
                
                int_dim=x;
                
                p=linspace(-sim_radius,sim_radius,res); obj.sim_x=p;
                q=0; obj.sim_y=q;
                
                T=exp(1i.*middle_row(obj.phase_profile));
            end
            
            [X,Y]=meshgrid(x,y);
            
            %% Output field coordinate space
            
            
            
            
            %% Field before DOE
            
            E0=sdtDOE.gauss(sim_source_amplitude,Lambda1,sim_source_w0,X,Y,sim_source_dist_from_DOE); % E0 stands for intial field
            
            
            %% Aperture shadow
            
            Ea=E0.*T;
            
            %% Output field
            
            I=zeros(length(q),length(p));
            
            tic
            %% Solve the diffraction integral for each p,q,z in the output field
            for k=1:length(Z)
                z=Z(k); %to avoid broadcast variable
                E=zeros(length(q),length(p));
                for i=1:length(p)
                    ptemp=p(i); %to avoid broadcast variable
                    for j=1:length(q)
                        P=exp(1i.*(2*pi)./Lambda2.* sqrt(z.^2+ (X-ptemp).^2 + (Y-q(j)).^2))./(z.^2+ (X-ptemp).^2+ (Y-q(j)).^2);
                        integrand=Ea.*P;
                        
                        E(j,i)=z/(1i.*Lambda2)*trapz(int_dim,integrand);
                    end
                end
                
                if strcmp(sim_axis,'x')==1
                    I(k,:)=abs(E).^2;
                elseif strcmp(sim_axis,'y')==1
                    I(:,k)=abs(E).^2;
                end
            end
            t=toc;
            
            obj.sim_intensity=I;
            
            
        end
        
        function [obj,t]=simulation3D(obj,sim_no_points,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,sim_radius,Z,sim_centre,benchmark_bool)
            %SIMULATION3D simulates the 2D beam profile at a range of z
            %positions after the DOE given a Gaussian source as input
            %INPUTS:
            % obj - a sdtDOE class object with a generated phase profile
            % sim_no_points - the dimension of the output array (i.e.
            % sim_no_points=20 simulates a 20x20 array of points at each z
            % position)
            % sim_source_amplitude - amplitude of the source Gaussian beam,
            % simulation optimised for a value of 1
            % sim_source_w0 - the minimum beam waist of the source beam
            % sim_source_wavelength - the (monochromatic) wavelength of the
            % source beam
            % sim_source_dist_from_DOE - the distance between the source
            % minimum beam waist and the DOE, typically set at the design
            % object point
            % sim_radius - the physical radius of the simulation plane, set
            % as small as possible while still containing the output beam
            % at its greatest extent
            % Z - a vector of z positions
            % benchmark_bool - not used by the user, can safely ignore
            %OUTPUTS:
            % obj - the sdtDOE object containing simulation data in
            % obj.sim_intensity
            % t - the length of time the simulation took, for benchmarking
            % purposes
            %
            % [obj,t]=simulation3D(obj,sim_no_points,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,sim_radius,Z,benchmark_bool)
            
            %% if the user has ignored the benchmark_bool input, set it to the input that generates the correct behaviour
            if nargin<9
                benchmark_bool=0;
            end
            
            %% if the benchmark procedure has been run, then estimate simulation time, if not, suggest running it
            if isempty(obj.benchmark)==0
                fprintf('Estimated simulation time is %.4f seconds (%.3f minutes, or %.2f hours)\n',obj.benchmark*(sim_no_points^2)*length(Z),(obj.benchmark*(sim_no_points^2)*length(Z))/60, (obj.benchmark*(sim_no_points^2)*length(Z))/3600)
            elseif benchmark_bool==0
                %fprintf('Cannot estimate simulation time, to enable this feature please use the sim_benchmark method\n')
            else
                % do nothing (not a critical feature)
            end
            
            if isempty(obj.phase_profile)==1
                phase_bool=input('No phase profile has been generated, would you like to generate a basic point to point phase profile? [y/n]:','s');
                if strcmp(phase_bool,'y')
                    obj=obj.gen_p2p_phase_profile;
                    fprintf('Phase profile generated, continuing with simulation\n')
                elseif strcmp(phase_bool,'n')
                    error('No simulation can continue without a phase profile, please generate a suitable profile and retry simulation')
                end
            end
            
            %% Store inputs in the sdtDOE object for the record
            obj.sim_source_amplitude=sim_source_amplitude;
            obj.sim_source_w0=sim_source_w0;
            obj.sim_source_wavelength=sim_source_wavelength;
            obj.sim_source_dist_from_DOE=sim_source_dist_from_DOE;
            obj.sim_z=Z;
            
            %% Defining shorthand variable name
            lambda=sim_source_wavelength;
            
            %% Calculate the refractive index of the material
            n=refractive_index(obj.material_substrate,lambda);
            
            %% Calculate the reduced wavelength inside the material
            lambda_s=lambda/n;
            
            %% Assign whether the substrate is before or after the DOE
            if obj.direction=='f'
                Lambda1=lambda_s;
                Lambda2=lambda;
            elseif obj.direction=='b'
                Lambda1=lambda;
                Lambda2=lambda_s;
            end
            
            %% Transmittance function
            
            T=exp(1i.*obj.phase_profile);
            
            
            %% DOE aperture coordinate space
            
            x=obj.profile_x;
            y=obj.profile_y;
            [X,Y]=meshgrid(x,y);
            
            %% Output field coordinate space
            
            res=sim_no_points;
            p=linspace(-sim_radius,sim_radius,res)+sim_centre(1); obj.sim_x=p;
            q=fliplr(linspace(-sim_radius,sim_radius,res))+sim_centre(2); obj.sim_y=q;
            
            
            %% Field before DOE
            
            E0=sdtDOE.gauss(sim_source_amplitude,Lambda1,sim_source_w0,X,Y,sim_source_dist_from_DOE); % E0 stands for intial field
            
            
            %% Aperture shadow
            
            Ea=E0.*T;
            
            %% Output field
            
            I=zeros(length(q),length(p),length(Z));
            
            tic
            %% Solve the diffraction integral for each p,q,z in the output field
            for k=1:length(Z)
                z=Z(k); %to avoid broadcast variable
                E=zeros(length(q),length(p));
                for i=1:length(p)
                    ptemp=p(i); %to avoid broadcast variable
                    parfor j=1:length(q)
                        P=exp(1i.*(2*pi)./Lambda2.* sqrt(z.^2+ (X-ptemp).^2 + (Y-q(j)).^2))./(z.^2+ (X-ptemp).^2+ (Y-q(j)).^2);
                        integrand=Ea.*P;
                        E(j,i)=z/(1i.*Lambda2)*trapz(y,trapz(x,integrand));
                    end
                end
                I(:,:,k)=abs(E).^2;
            end
            t=toc;
            
            obj.sim_intensity=I;
            
            
        end
          
        function obj=sim_benchmark(obj)
            %SIM_BENCHMARK calculates the time per integration to simulate
            %the output field of the DOE, allowing simulation time
            %estimates to be made in simulation2D
            % INPUTS:
            % obj - a sdtDOE object
            % OUTPUTS:
            % obj - sdtDOE object with benchmark time recorded in
            % obj.benchmark
            %
            % obj=sim_benchmark(obj)
            
            N=20;
            [~,t]=simulation3D(obj,N,1,3.75,obj.design_wavelength,obj.object_coords(3),20,obj.op_image_coords(3)/2,1);
            
            obj.benchmark=t/(N^2);
            
        end
        
        function time=sim_howlong(obj,sim_no_points,Z)
            %SIM_HOWLONG calculates the length of time a particular
            %simulation will take, given the simulation parameters
            % INPUTS:
            % obj - a sdtDOE object with a stored obj.benchmark value
            % sim_no_points - the array dimension of the output field
            % Z - the z position vector used for simulation
            % OUTPUTS:
            % time - the length of time, in seconds the simulation will
            % take, result also printed to display
            %
            % time=sim_howlong(obj,sim_no_points,Z)
            
            time=obj.benchmark*(sim_no_points^2)*length(Z);
            fprintf('Estimated simulation time is %.4f seconds (%.3f minutes, or %.2f hours)\n',time,time/60, time/3600)
        end
        
        function I=quick_sim(obj,type,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,Z)
            % Redundant, use simulation1D instead
            % QUICK_SIM simulates one intensity point per z value, yielding
            %fast simulations of DOE behaviour
            % INPUTS
            % obj - a sdtDOE class object with a generated phase profile
            % type - 'design', 'fab' or 'both' depending on which profile is to be
            % simulated
            % sim_source_amplitude - amplitude of the source Gaussian beam,
            % simulation optimised for a value of 1
            % sim_source_w0 - the minimum beam waist of the source beam
            % sim_source_wavelength - the (monochromatic) wavelength of the
            % source beam
            % sim_source_dist_from_DOE - the distance between the source
            % minimum beam waist and the DOE, typically set at the design
            % object point
            % Z - a vector of z position values
            % OUPUTS
            % I - vector of intensity values at each z position
            %
            % I=quick_sim(obj,type,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,Z)
            
            if strcmp(type,'design')==1
                
                %% Use simulation2D to calculate the intensity propagation
                temp=simulation3D(obj,1,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,0,Z);
                
                I=zeros(1,length(Z));
                
                %% Store the intensity values in vector format
                for i=1:length(Z)
                    I(i)=temp.sim_intensity(1,1,i);
                end
                
                
                %% Plot for easy analaysis
                %                 figure
                %                 plot(Z.*1e06,I)
                %                 xlabel('Distance from DOE [\mum]')
                %                 ylabel('Intensity [a.u.]')
                %                 set(gcf,'color','white')
                %                 title('Design')
                %
                %% The rest follows the same way, but for the other options
            elseif strcmp(type,'fab')==1
                temp=fab_simulation3D(obj,1,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,0,Z);
                
                I=zeros(1,length(Z));
                
                for i=1:length(Z)
                    I(i)=temp.sim_fab_intensity(1,1,i);
                end
                
                figure
                plot(Z.*1e06,I)
                xlabel('Distance from DOE [\mum]')
                ylabel('Intensity [a.u.]')
                set(gcf,'color','white')
                title('Fabricated')
                
                
            elseif strcmp(type,'both')==1
                
                temp1=simulation3D(obj,1,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,0,Z);
                
                I1=zeros(1,length(Z));
                
                for i=1:length(Z)
                    I1(i)=temp1.sim_intensity(1,1,i);
                end
                
                temp2=fab_simulation3D(obj,1,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,0,Z);
                
                I2=zeros(1,length(Z));
                
                for i=1:length(Z)
                    I2(i)=temp2.sim_fab_intensity(1,1,i);
                end
                
                %                 figure
                %                 plot(Z,I1)
                %                 hold on
                %                 plot(Z,I2)
                %                 legend('Design','Fabricated')
                %                 xlabel('Distance from DOE [\mum]')
                %                 ylabel('Intensity [a.u.]')
                %                 set(gcf,'color','white')
                
            end
            
        end
        
        function [obj,I]=angular_spectrum_sim(obj,sim_source_amplitude,sim_source_w0,sim_source_wavelength,sim_source_dist_from_DOE,Z)
            %% NOT OPERATIONAL
            
            if isempty(obj.phase_profile)==1
                phase_bool=input('No phase profile has been generated, would you like to generate a basic point to point phase profile? [y/n]:','s');
                if strcmp(phase_bool,'y')
                    obj=obj.gen_p2p_phase_profile;
                    fprintf('Phase profile generated, continuing with simulation\n')
                elseif strcmp(phase_bool,'n')
                    error('No simulation can continue without a phase profile, please generate a suitable profile and retry simulation')
                end
            end
            
            fprintf('Meshing simulation\n')
            
            %% Store inputs in the sdtDOE object for the record
            obj.sim_source_amplitude=sim_source_amplitude;
            obj.sim_source_w0=sim_source_w0;
            obj.sim_source_wavelength=sim_source_wavelength;
            obj.sim_source_dist_from_DOE=sim_source_dist_from_DOE;
            obj.sim_z=Z;
            
            
            %% Defining shorthand variable name
            lambda=sim_source_wavelength;
            
            %% Calculate the refractive index of the material
            n=refractive_index(obj.material_substrate,lambda);
            
            %% Calculate the reduced wavelength inside the material
            lambda_s=lambda/n;
            
            %% Assign whether the substrate is before or after the DOE
            if obj.direction=='f'
                Lambda1=lambda_s;
                Lambda2=lambda;
            elseif obj.direction=='b'
                Lambda1=lambda;
                Lambda2=lambda_s;
            end
            
            
            % Materials and EM
            %             lambda0 = 422e-9; % wavelength [m]
            %             n1 = 1; % refractive index.
            %             lambda = lambda0/n1; % wavelength [m]
            
            k0 = 2*pi/Lambda1; %Wavevector [rad/m]
            k1 = 2*pi/Lambda2; %Wavevector [rad/m]
            % geometry
            xmx = obj.radius; % max. +ve value of transverse coordinate in source place
            % check: Shen and Wang, Appl. Opt. 45, 1102-1110 (2006)
            % discusses suitable choice of xmx
            dx = Lambda2/2; % spatial sampling
            N = 2^(1+nextpow2(xmx./dx));  % For FFT, grid is to be of size 2^m,
            % "nextpow2" pads the vector with zeros is
            % needed
            x =(-N/2:N/2-1).*dx;         % Generate x axis
            i_x0=N/2+1;                   % index of x = 0
            [X,Y] = meshgrid(x,x);
            % X-Y Array
            % mesh for reciprocal space
            kxmax = 2*pi./dx;
            dkx = kxmax./N;
            kx =(-N/2:N/2-1).*dkx;
            [KX,KY]=meshgrid(kx,kx);
            %Define Propagator
            %-----------------
            % NOTE: MUST BE POSITIVE EXPONENT TO ENSURE Kx^2+ky^2>k leads to
            % exponentially decaying evanescent waves!!
            % H -> propagator in reciprocal space
            % H(z) =  exp(-jk_z.z) = exp(-j*sqrt[k1^2-kx^2-ky^2]*z)
            % the fftshift function becomes useful at the IFFT2 stage.
            kz = sqrt(k1.^2-KX.^2-KY.^2);
            H=@(z) fftshift(exp(+1j*kz.*z)); %
            
            % % plot transverse wavenumbers (kx, ky) and their sum^2
            % cc = contour(kx/k1,kx/k1, (KX.^2+KY.^2)/k1^2); colorbar
            % title('(k^2_x+k^2_y)/k^2')
            %
            % hold on %
            % cc = contour(kx/k1,kx/k1, (KX.^2+KY.^2)/k1^2,1,'Linewidth',2)
            % hold off
            %Source or object field Field
            %---------------------
            % we assume a gaussian beam
            w0 = sim_source_w0; % waist of beam
            % E0 = zeros(size(X,Y)); %Initially a plane wave along z
            E0 = sdtDOE.gauss(1,Lambda1,w0,X,Y,0);
            % Angular spectra of object field
            A0 = fft2((E0));
            
            % % test - take compate IFT{FT{E0}} with E0
            % IFT_A0 = ifft2(A0);
            % plot(x,E0(:,i_x0),'-k' ,x,IFT_A0(:,i_x0),'or')
            % legend('E0','IFT(FT[E0])')
            % propagate field to a point z1.
            
            % Gaussian beam expression (paraxial)
            E1 = sdtDOE.gauss(1,Lambda1,w0,X,Y,sim_source_dist_from_DOE);
            
            % Angular spectra of object field
            A1 = A0.*H(sim_source_dist_from_DOE);
            iA1 = ifft2(A1);
            % test - compare IFT{FT{E1}} with E1
            figure;
            plot(1e6*x,real(E1(:,i_x0)),'-k' ,1e6*x,real(iA1(:,i_x0)),'.r'...
                ,1e6*x,imag(E1(:,i_x0)),'--k' ,1e6*x,imag(iA1(:,i_x0)),':r')
            legend('Re {E1}','Re {IFT(A1)}','Im {E1}','Im {IFT(A1)}')
            title(['Source to lens: z = ',num2str(1e6*sim_source_dist_from_DOE),' \mum']);
            xlabel('x [\mum]')
            
            
            figure
            plot(1e6*x,abs(E1(:,i_x0)).^2,'-k',1e6*x,abs(iA1(:,i_x0)).^2,'.r')
            legend('Intensity, E1','Intensity, IFT(A1)')
            title(['Source to lens: z = ',num2str(1e6*sim_source_dist_from_DOE),' \mum']);
            xlabel('x [\mum]')
            % ok up to here! %
            %% propagation through a phase lens:
            phi = resizem(obj.phase_profile,size(X)); % initialise phase of lens
            
            figure
            plot(1e6*x,mod(phi(:,i_x0),2*pi))
            
            
            t = exp(1j*(phi));
            % ****************************************************
            % optical field after the lens:
            E2 = iA1.*t;
            %             E2=E1.*t;
            A2 = fft2(E2);
            
            
            
            % A2 = fftshift(A2);
            
            
            %%  2D propagate to z = z2
            % mesh for propagation:
            
            fprintf('Propagating\n')
            
            %Initialise Output Intensity Array
            I=zeros(N,length(Z));
            %Propagation
            for pp = 1: length(Z)
                
                if mod(pp/length(Z)*100,5)==0
                    fprintf('%i percent complete\n',pp/length(Z)*100)
                end
                
                A3 = A2.*H(Z(pp));
                iA3 = ifft2(A3);
                I(:,pp)=iA3(i_x0,:).*conj(iA3(i_x0,:));
            end
            
            
        end

        
        %% IMPORT
        
        function obj=import_and_process_data(obj,BG_file,keyphrase,pixel_size,magnification, measured_source_dist_from_DOE, measured_source_wavelength,measured_source_w0,start_z,z_increments)
            %IMPORT_AND_PROCESS_DATA imports saved .tiff or .csv from the current
            %directory of measured beam profiles that contain a keyphrase (e.g. 'beamprofile') and
            %constructs an array of profiles based on the natural order of
            %the numeric digits at the end of each file name (e.g.
            %'beamprofile01.tiff' would be first, 'beamprofile02.tiff' second
            %etc.). This allows automatic file saving to be used in data
            %collection, while preserving the z position order of the data.
            % Data is background corrected and each frame is fitted with a
            % 2D Gaussian to extract its properties. The full 3D matrix of
            % beam profile is not saved so as to keep file sizes manageable
            % and stable.
            % INPUTS:
            % obj - a sdtDOE object that corresponds in geometry and design to
            % the measured profile
            % BG_file - the file that contains the background data (.csv or
            % .tiff)
            % keyphrase - a common identifier in all filenames to be
            % imported (e.g. if the file is 'DOE1_beamprofile_01.bmp' a
            % keyphrase could be 'DOE1_' as long as other data set files do
            % not contain this phrase)
            % pixel_size - the pixel dimension of the CCD used to capture
            % the beam profile
            % magnification - the calibrated magnification of the
            % microscope used to enlarge the beam profile
            % measured_source_dist_from_DOE - the measured position of the minimum
            % beam waist of the source beam before the DOE
            % measured_source_wavelength - the central/dominant wavelength of the source
            % beam
            % measured_source_w0 - the source minimum beam waist, measured
            % by the beam profiler
            % start_z - start position of data taking (i.e. which position '01' corresponds to) the DOE plane is z=0
            % z_increments - the step size between z positions
            % OUTPUTS
            % obj - sdtDOE object with all processed measured data and measurement parameters stored
            %
            % obj=import_measured_beamprofiles(obj,keyphrase,pixel_size,magnification, measured_source_dist_from_DOE, measured_source_wavelength,measured_source_w0,start_z,z_increments)
            
            
            [~,~,ext]=fileparts(BG_file); %obtain the extension of the background file (changes import procedure)
            
            %% Import background
            if strcmp('.csv',ext)
                obj.background=double(csvread(BG_file));
            elseif strcmp('.tiff',ext)
                obj.background=double(imread(BG_file));
            end
            s=size(obj.background);
            
            % Save measurement settings to file
            obj.CCD_pixelsize=pixel_size;
            obj.magnification=magnification;
            obj.measured_source_dist_from_DOE=measured_source_dist_from_DOE;
            obj.measured_source_wavelength=measured_source_wavelength;
            obj.measured_source_w0=measured_source_w0;
            
            obj.measured_x=linspace(-s(2)*pixel_size/(2*magnification),s(2)*pixel_size/(2*magnification),s(2));
            obj.measured_y=linspace(-s(1)*pixel_size/(2*magnification),s(1)*pixel_size/(2*magnification),s(1));
            
            %% Import each beam profile and process in order or z value
            files=dir; %find all files in current director
            filelist={files.name}'; % extract list of filenames
            filelist=sort_nat(filelist); % sort_nat places filename02 before filename10, not the case in sort - see help sort_nat for more details
            j=1;
            for i=1:length(filelist)
                d=strfind(filelist(i),keyphrase); % does the current file match the keyphrase?
                if isempty(d{1})==0 % if so place it in the beam profile array
                    
                    %% Import current frame and background subtract it
                    [~,~,ext]=fileparts(char(filelist(i)));
                    
                    if strcmp('.csv',ext)
                        Itemp=double(csvread(char(filelist(i))));
                    elseif strcmp('.tiff',ext)
                        Itemp=double(imread(char(filelist(i))));
                    end
                    
                    Itemp(Itemp<mean(mean(obj.background)))=0;
                    Itemp(Itemp>mean(mean(obj.background)))=Itemp(Itemp>mean(mean(obj.background)))-mean(mean(obj.background));
                    
                    %% Extract the central intensity for plotting
                    obj.measured_intensity_central(j)=max(max(Itemp));
                    
                    %% 2D Gaussian fit
                    
                    fit=sdtDOE.gaussfit2D(Itemp./max(max(Itemp)),pixel_size/magnification*1e06,5,0);
                    
                    major(j)=2*fit.sigma_major*1e-06;
                    minor(j)=2*fit.sigma_minor*1e-06;
                    major95(:,j)=2*fit.sigma_major_95conf.*1e-06;
                    minor95(:,j)=2*fit.sigma_minor_95conf.*1e-06;
                    x0(j)=fit.x0*1e-06;
                    y0(j)=fit.y0*1e-06;
                    
                    xindex=find_closest_member(obj.measured_y,y0(j)-max(obj.measured_y));
                    obj.measured_intensity_xslice(j,:)=Itemp(xindex,:);
                    
                    yindex=find_closest_member(obj.measured_x,x0(j)-max(obj.measured_x));
                    obj.measured_intensity_yslice(j,:)=Itemp(:,yindex);
                    
                    
                    j=j+1;
                    fprintf('Imported and processed file %i \n',j);
                end
                
            end
            
            %% Save profile parameters to file
            obj.measured_majorbeamwaist=major;
            obj.measured_minorbeamwaist=minor;
            obj.measured_majorbeamwaist_95conf=major95;
            obj.measured_minorbeamwaist_95conf=minor95;
            obj.measured_x0=x0;
            obj.measured_y0=y0;
            
            
            obj.measured_z=start_z:z_increments:((j-2)*z_increments);
            
            
        end
        
        function obj=import_measured_beamprofiles(obj,keyphrase,pixel_size,magnification, measured_source_dist_from_DOE, measured_source_wavelength,measured_source_w0,start_z,z_increments)
            %IMPORT_MEASURED_BEAMPROFILES imports saved bitmaps from the current
            %directory of measured beam profiles that contain a keyphrase (e.g. 'beamprofile') and
            %constructs an array of profiles based on the natural order of
            %the numeric digits at the end of each file name (e.g.
            %'beamprofile01.bmp' would be first, 'beamprofile02.bmp' second
            %etc.). This allows automatic file saving to be used in data
            %collection, while preserving the z position order of the data.
            %
            % IT IS ADVISED NOT TO USE IMPORT_MEASURED_BEAMPROFILES as it
            % generates unstable files due to large file sizes.
            % INPUTS:
            % obj - a sdtDOE object that corresponds in geometry and design to
            % the measured profile
            % keyphrase - a common identifier in all filenames to be
            % imported (e.g. if the file is 'DOE1_beamprofile_01.bmp' a
            % keyphrase could be 'DOE1_' as long as other data set files do
            % not contain this phrase)
            % pixel_size - the pixel dimension of the CCD used to capture
            % the beam profile
            % magnification - the calibrated magnification of the
            % microscope used to enlarge the beam profile
            % measured_source_dist_from_DOE - the measured position of the minimum
            % beam waist of the source beam before the DOE
            % measured_source_wavelength - the central/dominant wavelength of the source
            % beam
            % measured_source_w0 - the source minimum beam waist, measured
            % by the beam profiler
            % start_z - start position of data taking (i.e. which position '01' corresponds to) the DOE plane is z=0
            % z_increments - the step size between z positions
            % OUTPUTS
            % obj - sdtDOE object with all measured data and measurement parameters stored
            %
            % obj=import_measured_beamprofiles(obj,keyphrase,pixel_size,magnification, measured_source_dist_from_DOE, measured_source_wavelength,measured_source_w0,start_z,z_increments)
            
            disp('This function is not recommended as it generates large .mat files and leads to unstable behaviour, only use if importing small numbers of beam profiles. import_and_process_data is recommended.')
            
            files=dir; %find all files in current director
            filelist={files.name}'; % extract list of filenames
            filelist=sort_nat(filelist); % sort_nat places filename02 before filename10, not the case in sort - see help sort_nat for more details
            
            j=1;
            for i=1:length(filelist)
                d=strfind(filelist(i),keyphrase); % does the current file match the keyphrase?
                if isempty(d{1})==0 % if so place it in the beam profile array
                    obj.measured_intensity(:,:,j)=double(imread(char(filelist(i)))); % import bitmap as array of double values
                    j=j+1;
                end
                
            end
            
            s=size(obj.measured_intensity(:,:,1)); % for generating the CCD dimension arrays
            
            %% Record to object
            obj.CCD_pixelsize=pixel_size;
            obj.magnification=magnification;
            obj.measured_source_dist_from_DOE=measured_source_dist_from_DOE;
            obj.measured_source_wavelength=measured_source_wavelength;
            obj.measured_source_w0=measured_source_w0;
            obj.measured_z=start_z:z_increments:((j-2)*z_increments);
            obj.measured_x=linspace(-s(2)*pixel_size/(2*magnification),s(2)*pixel_size/(2*magnification),s(2));
            obj.measured_y=linspace(-s(1)*pixel_size/(2*magnification),s(1)*pixel_size/(2*magnification),s(1));
        end
        
        function obj=import_DE_beamprofiles(obj,filename_source,filename_DOE)
            %IMPORT_DE_BEAMPROFILEs imports the two beam profiles used to
            %calculate the diffraction efficiency of the DOE, the source
            %beam at its minimum beam waist and the minimum beam waist of
            %the DOE focal point
            % INPUTS:
            % obj - sdtDOE object corresponding to the fabricated DOE
            % filename_source - filename of the beam profile corresponding
            % to the source beam
            % filename_DOE - filename of the beam profile corresponding to
            % the minimum beam waist of the DOE output field
            % OUTPUTS:
            % obj - sdtDOE object with two beam profiles saved
            %
            % obj=import_DE_beamprofiles(obj,filename_source,filename_DOE)
            
            [~,~,ext]=fileparts(char(filename_source));
            
            if strcmp('.csv',ext)
                profile0=double(csvread(filename_source));
            elseif strcmp('.tiff',ext)
                profile0=double(imread(filename_source));
            end
            
            obj.measured_beamprofile_source_w0=profile0(:,:,1);
            
            [~,~,ext]=fileparts(char(filename_DOE));
            
            if strcmp('.csv',ext)
                profile1=double(csvread(filename_DOE));
            elseif strcmp('.tiff',ext)
                profile1=double(imread(filename_DOE));
            end
            obj.measured_beamprofile_DOE_w0=profile1(:,:,1);
            
        end
        
        function obj=import_BG(obj,background_filename)
            % IMPORT_BG allows the background file to be imported if it is
            % not already present. This is useful for background
            % subtracting the profiles used to calculate the diffraction
            % efficiency.
            % INPUTS
            % obj - sdtDOE object
            % background_filename - the filename of the background (.csv or
            % .tiff)
            % OUTPUT
            % obj - sdtDOE object with imported background
            %
            % obj=import_BG(obj,background_filename)
            
            [~,~,ext]=fileparts(char(background_filename));
            
            if strcmp('.csv',ext)
                obj.background=double(csvread(background_filename));
            elseif strcmp('.tiff',ext)
                obj.background=double(imread(background_filename));
            end
            
        end
        
        function obj=transfer_data_from_object(obj,sdtDOE_obj,datatype)
            % TRANSFER_DATA_FROM_OBJECT transfers the simulation data from
            % one object to another, this is useful if the simulation has
            % been performed on an HPC and needs to be merged with an
            % object that has already had its measurements imported and
            % processed.
            % INPUTS
            % obj - the sdtDOE object to import into
            % sdtDOE_obj - the sdtDOE object to import from
            % datatype - currently only supports 'sim' data transfer
            % OUTPUTS
            % obj - the sdtDOE object with imported sim data
            %
            % obj=transfer_data_from_object(obj,sdtDOE_obj,datatype)
            
            if strcmp(datatype,'sim')
                obj.sim_source_dist_from_DOE=sdtDOE_obj.sim_source_dist_from_DOE;
                obj.sim_source_wavelength=sdtDOE_obj.sim_source_wavelength;
                obj.sim_source_w0=sdtDOE_obj.sim_source_w0;
                obj.sim_source_amplitude=sdtDOE_obj.sim_source_amplitude;
                obj.sim_x=sdtDOE_obj.sim_x;
                obj.sim_y=sdtDOE_obj.sim_y;
                obj.sim_z=sdtDOE_obj.sim_z;
                obj.sim_intensity=sdtDOE_obj.sim_intensity;
                obj.sim_majorbeamwaist=sdtDOE_obj.sim_majorbeamwaist;
                obj.sim_minorbeamwaist=sdtDOE_obj.sim_minorbeamwaist;
                obj.sim_majorbeamwaist_95conf=sdtDOE_obj.sim_majorbeamwaist_95conf;
                obj.sim_minorbeamwaist_95conf=sdtDOE_obj.sim_minorbeamwaist_95conf;
                obj.sim_x0=sdtDOE_obj.sim_x0;
                obj.sim_y0=sdtDOE_obj.sim_y0;
                obj.sim_fullfitdata=sdtDOE_obj.sim_fullfitdata;
                obj.benchmark=sdtDOE_obj.benchmark;
            end
            
        end
        
        %% Analysis
        
        function obj=background_subtraction(obj,type)
            
            % BACKGROUND_SUBTRACTION subtracts the background either off
            % the diffraction efficiency profiles ('DE') or the measured
            % intensity profiles ('prop')
            
            % INPUTS
            % obj - sdtDOE object
            % type - 'DE' or 'prop'
            % OUTPUT
            % obj - sdtDOE object subtracted background data
            %
            % obj=background_subtraction(obj,type)
            
            
            if strcmp(type,'DE')==1
                
                obj.measured_beamprofile_source_w0(obj.measured_beamprofile_source_w0<mean(mean(obj.background)))=0;
                obj.measured_beamprofile_source_w0(obj.measured_beamprofile_source_w0>mean(mean(obj.background)))=obj.measured_beamprofile_source_w0(obj.measured_beamprofile_source_w0>mean(mean(obj.background)))-mean(mean(obj.background));
                
                obj.measured_beamprofile_DOE_w0(obj.measured_beamprofile_DOE_w0<mean(mean(obj.background)))=0;
                obj.measured_beamprofile_DOE_w0(obj.measured_beamprofile_DOE_w0>mean(mean(obj.background)))=obj.measured_beamprofile_DOE_w0(obj.measured_beamprofile_DOE_w0>mean(mean(obj.background)))-mean(mean(obj.background));
                
                
            elseif strcmp(type, 'prop')==1
                
                obj.measured_intensity(obj.measured_intensity<mean(mean(obj.background)))=0;
                obj.measured_intensity(obj.measured_intensity>mean(mean(obj.background)))=obj.measured_intensity(obj.measured_intensity>mean(mean(obj.background)))-mean(mean(obj.background));
                
            end
        end
        
        function obj=calc_diffraction_efficiency(obj,n)
            %CALC_DIFFRACTION_EFFICIENCY uses the imported source and
            %output field beam profiles to calculate the diffraction
            %efficiency of the DOE. The background should be subtracted
            %beforehand.
            % INPUTS
            % obj - sdtDOE object with beam profiles imported using
            % import_DE_beamprofiles
            % OUTPUTS
            % obj - sdtDOE object with diffraction efficiency recorded in
            % obj.diffractio_efficiency
            %
            % obj=calc_diffraction_efficiency(obj)
            
            %% Place a numerical pinhole of radius 2 times the minimum beam waist of the beam profiles over the profiles to reject background
            I0=sdtDOE.pinhole(obj.measured_beamprofile_source_w0,obj.CCD_pixelsize/obj.magnification.*1e06,n);
            %I0=obj.measured_beamprofile_source_w0;
            I1=sdtDOE.pinhole(obj.measured_beamprofile_DOE_w0,obj.CCD_pixelsize/obj.magnification.*1e06,n);
            
            %% Calculate the power of the profiles
            P0=powerofbeam(I0,obj.CCD_pixelsize/obj.magnification);
            
            P1=powerofbeam(I1,obj.CCD_pixelsize/obj.magnification);
            
            
            obj.diffraction_efficiency=P1/P0 *100; %efficiency is the ratio of these two values
            
            fprintf('Diffraction efficiency is %.2f percent\n',obj.diffraction_efficiency)
            
        end
        
        
        function obj=process_measured_beamprofiles(obj)
            %PROCESS_MEASURED_BEAMPROFILES processes the imported measured
            %beam profiles by fitting 2D Gaussians to each profile and
            %extracting the Gaussian parameters. These parameters are
            %sorted into vectors for easy plotting and analysis.
            % THIS IS NOT RECOMMENDED. USE IMPORT_AND_PROCESS_DATA INSTEAD.
            % INPUTS
            % obj - sdtDOE object with imported data included
            % OUTPUTS
            % obj - sdtDOE object with vectors of profile parameters including
            % their confidence intervals
            %
            % obj=process_measured_beamprofiles(obj)
            
            
            if isempty(obj.measured_intensity)==1
                error('No measured data found, please import and try again')
            end
            
            
            s=size(obj.measured_intensity);
            
            major=zeros(1,s(3));
            minor=zeros(1,s(3));
            major95=zeros(2,s(3));
            minor95=zeros(2,s(3));
            x0=zeros(1,s(3));
            y0=zeros(1,s(3));
            fulldata=zeros(1,s(3));
            
            for i=1:s(3)
                fit=gaussfit2D(obj.measured_intensity(:,:,i),obj.CCD_pixelsize/obj.magnification,5,0);
                
                major(i)=2*fit.sigma_major;
                minor(i)=2*fit.sigma_minor;
                major95(:,i)=2*fit.sigma_major_95conf;
                minor95(:,i)=2*fit.sigma_minor_95conf;
                x0(i)=fit.x0;
                y0(i)=fit.y0;
                
                %fulldata(i)=fit;
                
                p=i/s(3)*100;
                fprintf('%.3f percent complete\n',p)
                
            end
            
            obj.measured_majorbeamwaist=major;
            obj.measured_minorbeamwaist=minor;
            obj.measured_majorbeamwaist_95conf=major95;
            obj.measured_minorbeamwaist_95conf=minor95;
            obj.measured_x0=x0;
            obj.measured_y0=y0;
            obj.measured_fullfitdata=fulldata;
            
            obj.measured_intensity_central=obj.extract_central_intensity();
            
            s=size(obj.measured_intensity);
            obj.measured_intensity_xslice=zeros(s(3),s(2));
            for i=1:length(obj.measured_z)
                index=find_closest_member(obj.measured_y,obj.measured_y0(i)-max(obj.measured_y));
                obj.measured_intensity_xslice(i,:)=obj.measured_intensity(index,:,i);
            end
            
            obj.measured_intensity_yslice=zeros(s(3),s(1));
            for i=1:length(obj.measured_z)
                index=find_closest_member(obj.measured_x,obj.measured_x0(i)-max(obj.measured_x));
                obj.measured_intensity_yslice(i,:)=obj.measured_intensity(:,index,i);
            end
            
            del=input('Would you like to delete raw measured data? (y/n)\n','s');
            
            if strcmp('y',del)==1
                obj.measured_intensity=[];
            else
                fprintf('Data has not been deleted, if you intended to do so in the future input y as a string, for now manually set obj.measured_intensity to empty')
            end
            
            
        end
        
        
        function obj=process_sim_beamprofiles(obj)
            %PROCESS_SIM_BEAMPROFILES processes the simulated
            %beam profiles by fitting 2D Gaussians to each profile and
            %extracting the Gaussian parameters. These parameters are
            %sorted into vectors for easy plotting and analysis.
            % INPUTS
            % obj - sdtDOE object with simulation data
            % OUTPUTS
            % obj - sdtDOE object with vectors of profile parameters including
            % their confidence intervals
            %
            % obj=process_measured_beamprofiles(obj)
            
            if isempty(obj.sim_intensity)==1 && strcmp(type,'design')==1
                error('No design simulation data found, please run simulation and try again')
            end
            
            
            s=size(obj.sim_intensity);
            x=obj.sim_x;
            I=obj.sim_intensity;
            
            
            if length(s)>2
                major=zeros(1,s(3));
                minor=zeros(1,s(3));
                major95=zeros(2,s(3));
                minor95=zeros(2,s(3));
                x0=zeros(1,s(3));
                y0=zeros(1,s(3));
                max_index=s(3);
                fulldata=zeros(1,s(3));
            else
                max_index=1;
            end
            
            for i=1:max_index
                
                ps=(x(end) - x(1))/length(x);
                
                fit=sdtDOE.gaussfit2D(I(:,:,i),ps,5,0);
                
                major(i)=2*fit.sigma_major;
                minor(i)=2*fit.sigma_minor;
                major95(:,i)=2*fit.sigma_major_95conf;
                minor95(:,i)=2*fit.sigma_minor_95conf;
                x0(i)=fit.x0;
                y0(i)=fit.y0;
                
                %fulldata(i)=fit;
                
                
                p=i/max_index*100;
                fprintf('%.3f percent complete\n',p)
                
            end
            
            
            obj.sim_majorbeamwaist=major;
            obj.sim_minorbeamwaist=minor;
            obj.sim_majorbeamwaist_95conf=major95;
            obj.sim_minorbeamwaist_95conf=minor95;
            obj.sim_x0=x0;
            obj.sim_y0=y0;
            %obj.sim_fullfitdata=fulldata;
            
            
        end
        
        function I=extract_central_intensity(obj)
            
            if isempty(obj.measured_intensity)==1
                error('No measured data present, please import and try again')
            end
            
            s=size(obj.measured_intensity);
            
            for i=1:s(3)
                
                I(i)=max(max(obj.measured_intensity(:,:,i)));
                
            end
            
        end
        
        
        %% Plots
        
        function plot_geometry(obj)
            
            r=obj.radius;
            
            corners=[r,r;r,-r;-r,-r;-r,r];
            
            figure
            plot3(obj.image_coords(1),obj.image_coords(3),obj.image_coords(2),'rx')
            hold on
            plot3(obj.object_coords(1),-abs(obj.object_coords(3)),obj.object_coords(2),'bo')
            legend('Image point','Object point')
            for i=1:4
                hold on
                plot3([obj.image_coords(1),corners(i,1)],[obj.image_coords(3),0],[obj.image_coords(2),corners(i,2)],'k:')
                hold on
                plot3([corners(i,1),obj.object_coords(1)],[0,-abs(obj.object_coords(3))],[corners(i,2),obj.object_coords(2)],'k:')
            end
            hold on
            
            fill3(corners(:,1),[0,0,0,0]',corners(:,2),'b')
            grid on
            alpha(0.3)
            
            
            set(gcf,'color','white')
            hold on
            % DRAW AXIS LINEs
            plot3(get(gca,'XLim'),[0 0],[0 0],'k');
            plot3([0 0],[0 0],get(gca,'ZLim'),'k');
            plot3([0 0],get(gca,'YLim'),[0 0],'k');
            
            % GET TICKS
            X=get(gca,'Xtick');
            Y=get(gca,'Ytick');
            Z=get(gca,'Ztick');
            
            % REMOVE TICKS
            set(gca,'Xtick',[]);
            set(gca,'Ytick',[]);
            set(gca,'Ztick',[]);
            
            % GET OFFSETS
            Xoff=diff(get(gca,'XLim'))./30;
            %Yoff=diff(get(gca,'YLim'))./30;
            Zoff=diff(get(gca,'ZLim'))./30;
            
            % DRAW TICKS
            %%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
            for i=1:length(X)
                plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
            end;
            for i=1:length(Y)
                plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
            end;
            for i=1:length(Z)
                plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
            end;
            
            set(gca,'ytick',[],'Ycolor','w','box','off')
            set(gca,'xtick',[],'Xcolor','w','box','off')
            set(gca,'ztick',[],'Zcolor','w','box','off')
            
            view([103 26])
            axis equal tight
            
        end
        
        function plot_animated_propagation(obj,type,filename,titlestr)
            %PLOT_ANIMATED_PROPAGATION displays either 'design', 'fab' or
            %'measured' beam profiles as an animation and saves to file
            % CURRENTLY BROKEN
            % INPUTS
            % obj - sdtDOE object with either design or fab simulation data,
            % or measured data
            % type - 'design', 'fab', 'measured'
            % filename - string
            % titlestr - string to be displayed as title
            % OUTPUTS
            % saved gif of animation
            %
            % plot_animated_propagation(obj,type,filename,titlestr)
            
            close all
            
            if strcmp(type,'design')==1
                x=obj.sim_x;
                y=obj.sim_y;
                z=obj.sim_z;
                I=obj.sim_intensity;
            elseif strcmp(type,'fab')
                x=obj.sim_fab_x;
                y=obj.sim_fab_y;
                z=obj.sim_fab_z;
                I=obj.sim_fab_intensity;
                
            elseif strcmp(type,'measured')
                x=obj.measured_x;
                y=obj.measured_y;
                z=obj.measured_z;
                I=obj.measured_intensity;
            else
                fprintf('Please use a type: design, fab or measured')
            end
            
            
            figure
            for i=1:length(z)
                
                text_str=sprintf('%.3i microns',z(i)); %text to display on image of the current z position
                imagesc(x,y,I(:,:,i)) % plot the current beam profile
                text(max(x)-10,max(x)-5,text_str,'color','white') % place the z position text
                title(titlestr) % place title
                set(gcf,'color','white')
                xlabel('x [\mum]')
                ylabel('y [\mum]')
                axis square
                colormap(flipud(brewermap(100,'YlGnBu')))
                drawnow % draw the figure
                frame = getframe(1); % get the current fram from the figure
                im = frame2im(frame); % convert the frame to an image
                [A,map] = rgb2ind(im,256); % map to value value 1 to 256
                if i == 1;
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1); % write first frame to file
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1); % append remaining frames to file
                end
            end
        end
        
        function plot_waist_compare(obj,z_crop)
            %PLOT_WAIST_COMPARE plots the beam waist radii from the
            %simulated and measured data depending on whether they are
            %present or not
            % INPUTS
            % obj - sdtDOE object with at least one beam waist vector present
            % z_crop - the distance from the DOE the plot cuts off at
            % OUTPUTS
            % plot of beam waist radii
            
            
            present=[isempty(obj.sim_majorbeamwaist);
                isempty(obj.measured_majorbeamwaist);
                ];
            
            if sum(present)==3
                error('No waist data to compare')
            end
            
            figure
            if present(1)==0
                tmp=[];
                tmp=abs(obj.sim_z-z_crop);
                [~,idx]=min(tmp);
                
                plot(obj.sim_z(1:idx),obj.sim_majorbeamwaist(1:idx),'k--','DisplayName','Design, simulated')
                hold on
            end
            
            if present(2)==0
                tmp=[];
                tmp=abs(obj.measured_z-z_crop);
                [~,idx]=min(tmp);
                
                scatter(obj.measured_z(1:idx),obj.measured_majorbeamwaist(1:idx),'rx','DisplayName','Fabricated, measured')
                hold on
            end
            set(gcf,'color','white')
            box on
            xlabel('Distance from DOE plane [\mum]')
            ylabel('Beam waist radius [\mum]')
            legend(gca,'show')
            
        end
        
        function plot_propagation(obj,type,axis_select)
            %PLOT_PROPAGATION takes 1D slices of either measured ('fab') or
            %simulated ('sim') beam profiles in the 'x' or 'y' direction and stitches them together in
            %z position order to visualise a beam propagation
            % INPUTS
            % obj - a sdtDOE object with data to display
            % type - 'sim' or 'fab'
            % axis_select - 'x' (horizontal) or 'y' (vertical) slices
            % OUTPUTS
            % plot of propagation
            %
            % plot_propagation(obj,type,axis_select)
            
            if strcmp(type,'sim')==1
                
                if isempty(obj.sim_intensity)==1
                    error('No simulation data to display')
                end
                
                s=size(obj.sim_intensity);
                I=zeros(s(3),s(2));
                
                for i=1:length(obj.sim_z)
                    I(i,:)=middle_row(obj.sim_intensity(:,:,i));
                end
                
                z=obj.sim_z;
                r=obj.sim_y;
                
                
            elseif strcmp(type,'measured')
                
                
                if axis_select=='x'
                    
                    if isempty(obj.measured_intensity_xslice)==1
                        error('No processed measured data to display')
                    end
                    
                    I=obj.measured_intensity_xslice;
                    z=obj.measured_z;
                    r=obj.measured_x;
                    
                    
                elseif axis_select=='y'
                    
                    if isempty(obj.measured_intensity_yslice)==1
                        error('No processed measured data to display')
                    end
                    
                    s=size(obj.measured_intensity);
                    
                    
                    I=obj.measured_intensity_yslice;
                    
                    z=obj.measured_z;
                    r=obj.measured_y;
                    
                else
                    error('Please enter character x or y for axis choice')
                end
            end
            
            
            figure
            imagesc(z.*1e06,r.*1e06,rot90(I))
            set(gcf,'color','white')
            set(gca,'fontsize',14)
            xlabel('Distance from DOE plane [\mum]')
            ylabel('Radius [\mum]')
            %axis equal tight
            colormap(flipud(brewermap(100,'YlGnBu')))
            
        end
        
        function plot_phase_profile(obj,titlestr)
            %PLOT_PHASE_PROFILE displays the phase profile
            %as long as it are present
            % INPUTS
            % obj - sdtDOE object with desired profile to plot present
            % type - 'phase', 'surface', 'fab'
            % titlestr - a string that is to be displayed as the title
            % OUTPUTS
            % plot of the chosen profile
            %
            % plot_profile(obj,type,titlestr)
            
            
            
            if isempty(obj.phase_profile)==1
                error('No phase profile present, please generate a try again')
            end
            
            figure
            imagesc(obj.profile_x.*1e06,obj.profile_y.*1e06,mod(obj.phase_profile,2*pi))
            set(gcf,'color','white')
            set(gca,'fontsize',22)
            set(gcf,'color','white')
            xlabel('x [\mum]')
            ylabel('y [\mum]')
            title(titlestr)
            h=colorbar();
            xlabel(h,'Phase [rads]')
            axis square
            colormap(flipud(brewermap(100,'YlGnBu')))
            set(gca,'ydir','normal')
            
            
        end
        
    end
    
    methods(Static)
        function [fit,Intensity_x,Intensity_y]=gaussfit2D(A,ps,threshold,plot_bool)
            %% Fits a 2D Gaussian to an array of inputted data, finding the principle axes of the Gaussian, their rotation to the array's coordinate system and the Gaussian parameters (Amplitude, x0, y0, sigmax, sigmay, theta)
            % Input
            %
            %   A - An array of data for the Gaussian to be fit to
            %   ps - The pixel size of the data, that is, the spacing between array
            %   points. If using data from a beam profiler for example, this is the
            %   quoted pixel size of the CCD.
            %   threshold - the script crops out background to make the fitting more robust, the threshold is the percent cut-off of maximum value to crop beam profile (use
            %   5-10[%] depending on how noisy the background is - contemplate background removal)
            %   plot - Set to 1 for a confirmation plot, set to anything else to
            %   suppress plotting.
            %
            % Output
            %
            %   fit - A data structure containing the fields:
            %            amplitude: The maximum value of the Gaussian
            %            x0: The position of the centre of the Gaussian in the horizontal direction of the coordinate system
            %            generated in the script. Use to compare relative distances
            %            with other fits.
            %            y0: The position of the centre of the Gaussian in the vertical direction of
            %            the coordinate system, as above.
            %            sigmax: The variance of the Gaussian in its principle x-axis
            %            sigmay: The variance of the Gaussian in its principle y-axis
            %            theta: The angle of the principle y-axis from vertical,
            %            ranging from 0 to pi/2. It has been found that this range of
            %            angles gives the most consistent output results. Ranging from
            %            -pi/4 to pi/4 can allow the definition of the principle axes
            %            to vary between fits of similar but slightly rotated
            %            Gaussians.
            %            X_95conf: the 95 confidence interval of the value X
            %
            % This script used mainD2GaussFitRot.m (author: G. Nootz) as a base.
            % Features added:
            % 1. Functionalised to take any input array, 2. Background cropping feature
            % to make the fitting more robust, 3. Automatic, robust initial guesses for
            % fitting parameters, 4. Calculation of the 95% confidence interval of fit
            % parameters, 5. Uses a more reliable and robust form of a 2D Gaussian
            % function.
            %
            % Author: Matthew Day
            % Date:  February 2016
            % Affiliation: University of Bristol
            % Email: matt.day@bristol.ac.uk
            
            A=A./max(max(A));
            
            %% Crop the input array to be centered on the Gaussian
            sa=size(A);
            p=linspace(0,ps*sa(2),sa(2)); % the coordinates of the full array
            q=linspace(0,ps*sa(1),sa(1));
            
            A(isnan(A))=0;
            Aordered=sort(A(:),'descend'); % generate a vector of the values of A, from maximum to minimum
            
            
            for i=1:10 % for the top 10 values of A, generate a list of indices
                [x_temp,y_temp]=find(A==Aordered(i));
                xindex(i)=x_temp(1);
                yindex(i)=y_temp(1);
            end
            
            ix_peak=ceil(median(xindex)); % take the median value of the indices (i.e. a majority vote), this will bias against very localised but high intensity noise
            iy_peak=ceil(median(yindex));
            
            Amax=A(ix_peak,iy_peak);
            yprobe=A(ix_peak,iy_peak);
            Iy=iy_peak;
            while yprobe>(threshold/100 * Amax) & Iy<=sa(2) % while the current probe value is greater than the threshold value (i.e. a percentage of the peak maximum), keep increasing probe position in the y direction
                yprobe=A(ix_peak,Iy);
                Iy=Iy+1;
            end
            xprobe=A(ix_peak,iy_peak);
            Ix=ix_peak;
            while xprobe>(threshold/100 * Amax) & Ix<=sa(1) %same but in x direction
                xprobe=A(Ix,iy_peak);
                Ix=Ix+1;
            end
            
            ix=ix_peak;
            iy=iy_peak;
            
            jymin=(ix-abs(Ix-ix));
            jymax=(ix+abs(Ix-ix));
            jxmin=(iy-abs(Iy-iy));
            jxmax=(iy+abs(Iy-iy));
            
            if jymin<=0 % if the minimum index is out of bounds, just make the minimum possible value
                jymin=1;
            end
            
            if jymax>sa(1) % if minimum index out of bounds, make it the maximum possible value
                jymax=sa(1);
            end
            
            if jxmin<=0
                jxmin=1;
            end
            
            if jxmax>sa(2)
                jxmax=sa(2);
            end
            
            C=A(jymin:jymax,jxmin:jxmax); % construct a 'crop' array with a centred Gaussian to make fitting easier
            sc=size(C);
            
            if sc(1)<2 || sc(2)<2
                C=A;
                x=p;
                y=q;
                fprintf('Cropping failed, fitting to full profile - will take longer than usual')
            else
                x=p(jxmin:jxmax); %coordinates of the crop array
                y=q(jymin:jymax);
            end
            
            
            %% Define fit parameters
            
            
            [ix,iy]=find(max(max(C))); % now assume the noise has been cropped out and that the maximum value of the array is the peak of the Gaussian
            
            angle_origin=0; % determines the basis of the angle parameter, 0 is vertical, and all angles will be defined with respect to it.
            
            % param0=[Amp,x0,sigmax,y0,sigmay,angle]
            param0=[max(max(C)), x(ix), (max(x)-min(x))/4 ,y(iy), (max(y)-min(y))/4, angle_origin]; %initial guesses for the parameters
            
            lb=[0.99*max(max(C)), x(1), (max(x)-min(x))/16, y(1), (max(y)-min(y))/16, angle_origin-pi/4]; % lower bounds on fitting parameters
            ub=[1.01*max(max(C)), x(end), (max(x)-min(x))/2,y(end),(max(y)-min(y))/2, angle_origin+pi/4]; % upper bounds
            
            [X,Y]=meshgrid(x,y); % define coordinates for the fitting procedure
            coords(:,:,1)=X;
            coords(:,:,2)=Y;
            
            %% Fitting procedure
            options=optimoptions(@lsqcurvefit, 'Display','none','TolPCG',0.0001);
            [param,~,residual,exitflag,~,~,J] = lsqcurvefit(@gaussian2D,param0,coords,C,lb,ub,options); % see documentation of lsqcurvefit for more details on its operation
            
            if exitflag<=0
                fprintf('Warning: solution not found, check input array\n')
            end
            
            if (param(3) | param(5))<= 4/10*ps
                fprintf('Warning: fit interpolated from 10 pixels or less\n')
            end
            
            ci=nlparci(param,residual,'jacobian',J); % calculates confidence intervals from the residals and jacobian
            
            %% Output structure definition
            
            fit.amplitude=param(1);
            fit.amplitude_95conf=[ci(1,:)];
            fit.x0=param(2);
            fit.x0_95conf=[ci(2,:)];
            fit.sigmax=param(3);
            fit.sigmax_95conf=[ci(3,:)];
            fit.y0=param(4);
            fit.y0_95conf=[ci(4,:)];
            fit.sigmay=param(5);
            fit.sigmay_95conf=[ci(5,:)];
            fit.theta=param(6);
            fit.theta_95conf=[ci(6,:)];
            
            % define which axis is major and which is minor
            if fit.sigmax>fit.sigmay
                fit.sigma_major=fit.sigmax;
                fit.sigma_major_95conf=fit.sigmax_95conf;
                fit.sigma_minor=fit.sigmay;
                fit.sigma_minor_95conf=fit.sigmay_95conf;
            elseif fit.sigmay>fit.sigmax
                fit.sigma_major=fit.sigmay;
                fit.sigma_major_95conf=fit.sigmay_95conf;
                fit.sigma_minor=fit.sigmax;
                fit.sigma_minor_95conf=fit.sigmax_95conf;
            elseif fit.sigmax==fit.sigmay
                fit.sigma_major=fit.sigmax;
                fit.sigma_major_95conf=fit.sigmax_95conf;
                fit.sigma_minor=fit.sigmay;
                fit.sigma_minor_95conf=fit.sigmay_95conf;
            end
            
            %% Check validity of fit by finding the overlap with inputted array
            
            G=gaussian2D(param,coords);
            fit.overlap=(trapz(x,trapz(y,G.*C))).^2./(trapz(x,trapz(y,G.^2)).*trapz(x,trapz(y,C.^2)));
            
            if fit.overlap<0.90
                fprintf('Warning: overlap between inputted Gaussian and fitted Gaussian is sub-optimal, check plot for validity\n')
            end
            
            
            
            [P,Q]=meshgrid(p,q); % define coordinates for the fitting procedure
            coords_full(:,:,1)=P;
            coords_full(:,:,2)=Q;
            
            G_full=gaussian2D(param,coords_full);
            sg=size(G_full);
            
            fprintf('The major MFD is %.2f microns and the minor MFD is %.2f microns\n',4*fit.sigma_major.*1e06,4*fit.sigma_minor.*1e06)
            fprintf('The overlap is %.2f percent\n',fit.overlap*100)
            
            
            
            [~,ix]=min(abs(p-fit.x0));
            [~,iy]=min(abs(q-fit.y0));
            
            %span=floor(ix/2);
            %span=floor(1.5*2*fit.sigma_major/ps);
            span=floor((21e-6)/ps);
            
            Intensity_x=A(iy,(ix-span):(ix+span));
            Intensity_y=A((iy-span):(iy+span),ix);
            
            %% Plotting
            
            if plot_bool==1
                A(isinf(A))=NaN;
                
                
                x_profile=linspace(-ps*length(x)/2,ps*length(x)/2,length(x)).*1e06;
                y_profile=linspace(-ps*length(y)/2,ps*length(y)/2,length(y)).*1e06;
                Fzoom=gaussian2D(param,coords);
                orientation=cot(fit.theta).*(x_profile);%+(fit.y0.*1e06 - (fit.x0.*1e06)/tan(fit.theta));
                
                
                figure
                subplot(1,2,1)
                imagesc(x_profile,y_profile,C);
                xlabel('x [\mum]')
                ylabel('y [\mum]')
                %axis([min(x) max(x) min(y) max(y)])
                axis square
                title('Measured profile')
                set(gcf,'Color','white')
                set(gca,'fontsize',13)
                colormap(flipud(brewermap(100,'YlGnBu')))
                subplot(1,2,2)
                imagesc(x_profile,y_profile,Fzoom);
                hold on
                %plot(x_profile,orientation.*1e06)
                plot(x_profile(ceil(length(x_profile)/4):ceil((3*length(x_profile))/4)),orientation(ceil(length(x_profile)/4):ceil((3*length(x_profile))/4)),'r')
                %axis([min(x) max(x) min(y) max(y)])
                axis square
                xlabel('x [\mum]')
                ylabel('y [\mum]')
                title('Fitted Gaussian')
                set(gca,'fontsize',13)
                legend('Major axis','Location','Southeast')
                
                
               
                
                
                r=linspace(-span.*ps,span.*ps,(2*span+1)).*1e06;
                figure
                plot(r,G_full(iy,(ix-span):(ix+span)),'k:','linewidth',1.5)
                hold on
                scatter(r,A(iy,(ix-span):(ix+span)),'o')
                axis([min(r) max(r) 0 1])
                title('x-axis')
                set(gcf,'color','white')
                xlabel('Radius [\mum]')
                ylabel('Relative intensity to peak')
                box on
                set(gca,'fontsize',14)
                legend('Gaussian fit','Measured')
                
                
                figure
                scatter(r,A((iy-span):(iy+span),ix),'o')
                hold on
                plot(r,G_full((iy-span):(iy+span),ix),'k:','linewidth',1.5)
                axis tight
                axis([min(r) max(r) 0 1])
                title('y-axis')
                set(gcf,'color','white')
                xlabel('Radius [\mum]')
                ylabel('Relative intensity to peak')
                box on
                set(gca,'fontsize',14)
                legend('Measured','Gaussian fit')
                
                
                
                figure
                semilogy(r,A((iy-span):(iy+span),ix),'o')
                hold on
                semilogy(r,G_full((iy-span):(iy+span),ix),'k:','linewidth',1.5)
                axis tight
                axis([min(r) max(r) 0.00001 1])
                title('x-axis')
                set(gcf,'color','white')
                xlabel('Radius [\mum]')
                ylabel('Relative intensity to peak')
                set(gca,'fontsize',14)
                legend('Measured','Gaussian fit')
                
                
                
                figure
                semilogy(r,G_full(iy,(ix-span):(ix+span)),'k:','linewidth',1.5)
                hold on
                semilogy(r,A(iy,(ix-span):(ix+span)),'o')
                axis([min(r) max(r) 0.00001 1])
                title('y-axis')
                set(gcf,'color','white')
                xlabel('Radius [\mum]')
                ylabel('Relative intensity to peak')
                set(gca,'fontsize',14)
                legend('Gaussian fit','Measured')
                
%                 figure
%                 imagesc(r,r,log10(A((iy-span):(iy+span),(ix-span):(ix+span))))
%                 colormap(flipud(brewermap(100,'YlGnBu')))
%                 caxis([-5,0])
%                 axis square
%                 xlabel('x [\mum]')
%                 ylabel('y [\mum]')
%                 set(gcf,'color','white')
%                 set(gca,'fontsize',18)
                
            end
            
        end
        
        function I_pin=pinhole(I,ps,n)
            % PINHOLE places a pinhole on a Gaussian beam of n times that
            % Gaussian beam's radius, as calculated by sdtDOE.gaussfit2D
            % INPUT
            % I - intensity map
            % ps - pixel size
            % n - number of beam radii for pinhole radius to be
            % OUTPUT
            % I_pin - I with all values outside the pinhole radius set to 0
            %
            % I_pin=pinhole(I,ps,n)
            
            imageSize = size(I);
            fit=sdtDOE.gaussfit2D(double(I)/max(max(double(I))),ps,30,0);
            
            x=linspace(0,ps*imageSize(2),imageSize(2));
            y=linspace(0,ps*imageSize(1),imageSize(1));
            
            % find closest array index to known centre
            
            cmp_x=abs(x-fit.x0);
            [vx,ix]=min(cmp_x);
            cmp_y=abs(x-fit.y0);
            [vy,iy]=min(cmp_y);
            
            
            
            ci = [iy, ix, (2*fit.sigma_major*n)/ps];    % center and radius of circle ([c_row, c_col, r])
            [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
            %             mask = uint16(((xx.^2 + yy.^2)<ci(3)^2));
            %             croppedImage = uint16(zeros(size(I)));
            %             croppedImage=uint16(I).*mask;
            
            ci(3)*ps*1e6;
            
            mask = (((xx.^2 + yy.^2)<ci(3)^2));
            croppedImage = (zeros(size(I)));
            croppedImage=I.*mask;
            
            I_pin=double(croppedImage);
            
            %             figure
            %             imagesc(x,y,croppedImage);
            %             colormap jet
            
            
            
        end
        
        function I_pin=pinhole_radius(I,ps,r)
            % PINHOLE_RADIUS places a pinhole on a Gaussian beam of radius r
            % INPUT
            % I - intensity map
            % ps - pixel size
            % r - pinhole radius
            % OUTPUT
            % I_pin - I with all values outside the pinhole radius set to 0
            %
            % I_pin=pinhole(I,ps,n)
            
            imageSize = size(I);
            fit=sdtDOE.gaussfit2D(double(I)/max(max(double(I))),ps,30,0);
            
            x=linspace(0,ps*imageSize(2),imageSize(2));
            y=linspace(0,ps*imageSize(1),imageSize(1));
            
            % find closest array index to known centre
            
            cmp_x=abs(x-fit.x0);
            [vx,ix]=min(cmp_x);
            cmp_y=abs(x-fit.y0);
            [vy,iy]=min(cmp_y);
            
            
            
            ci = [iy, ix, (r)/ps];    % center and radius of circle ([c_row, c_col, r])
            [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
            %             mask = uint16(((xx.^2 + yy.^2)<ci(3)^2));
            %             croppedImage = uint16(zeros(size(I)));
            %             croppedImage=uint16(I).*mask;
            
            ci(3)*ps*1e6;
            
            mask = (((xx.^2 + yy.^2)<ci(3)^2));
            croppedImage = (zeros(size(I)));
            croppedImage=I.*mask;
            
            I_pin=double(croppedImage);
            
            %             figure
            %             imagesc(x,y,croppedImage);
            %             colormap jet
            
            
            
        end
        
        function u=gauss(A,lambda,w0,x,y,z)
            %Gaussian beam calculator
            % Calculates the (complex) electric field amplitude of a Gaussian beam in
            % the (x,y) plane, a distance z from its focal point. The other parameters
            % are:
            % A: electric field amplitude at (x,y,z)=(0,0,0)
            % w: beam waist radius at z=0
            % lambda: wavelength of light (remember to correct for the medium)
            % The output is a single, complex number.
            
            %[x,y]=meshgrid(x,y);
            
            zr=(pi * w0^2)./lambda; % the Rayliegh range
            
            k= (2*pi)./lambda;
            
            function [w]=waist(z)
                w=w0.*sqrt(1+(z./zr).^2);
            end
            
            function [R]=radius_curvature(z)
                R=z.*(1+(zr./z).^2);
            end
            
            
            function [xi]=gouy(z)
                xi=atan(z/zr);
            end
            
            
            function [rho]=radius(x,y)
                rho=sqrt(x.^2 + y.^2);
            end
            
            function curv=curvature(x,y,z)
                if z==0
                    curv=0;
                else
                    curv=(radius(x,y).^2)./(2*radius_curvature(z));
                end
            end
            
            u = ( A.*(w0./waist(z)).*exp(-(radius(x,y).^2)./(waist(z).^2)).*exp(-1i*k*z - 1i.*k.*curvature(x,y,z) + 1i.*gouy(z)));
            
            u=imag(u) + 1i*real(u); 
            
        end
        
        function I=hdr_profile(Iset,Rset,exposure_times,ps)
            % HDR_PROFILE stitches a set of beam profiles of different
            % exposure times to allow for high dynamic range profiles to be
            % constructed
            % INPUTS
            % Iset - a 3D matrix of 2D matrices of beam profiles, with the
            % third index specifying the beam profile number
            % Rset - a vector of crop radii to cut off each beam profile
            % where the index of the vector relates to the third index of
            % Iset
            % exposure_times - the exposure time in ms of each profile
            % measurement
            % ps - pixel size, including magnification
            % OUTPUT
            % I - stitched HDR profile
            %
            % I=hdr_profile(Iset,Rset,exposure_times,ps)
            
            s=size(Iset);
            x=0:ps:((s(2)-1)*ps);
            y=0:ps:((s(1)-1)*ps);
            
            fit=sdtDOE.gaussfit2D(Iset(:,:,1),ps,10,0);
            
            cmp_x=abs(x-fit.x0);
            [~,ix_centre]=min(cmp_x);
            cmp_y=abs(x-fit.y0);
            [~,iy_centre]=min(cmp_y);
            
            [xx,yy] = ndgrid((1:s(1))-iy_centre,(1:s(2))-ix_centre);
            I=(zeros(s(1),s(2)));
            for i=1:s(3)
                
                Iset(:,:,i)=sdtDOE.exposure_calibration(Iset(:,:,i),exposure_times(i),'WinCamD');
                
                if i==1
             
                    mask = (((xx.^2 + yy.^2)<((Rset(i))/ps)^2));
                
                    Itemp(:,:,i) = (zeros(s(1),s(2)));
                    Itemp(:,:,i)=Iset(:,:,i).*mask;
                elseif i>1
                   
                    mask = (((xx.^2 + yy.^2)<((Rset(i))/ps)^2) & ((Rset(i-1))/ps)^2<(xx.^2 + yy.^2) );
                    
                    Itemp(:,:,i) = (zeros(s(1),s(2)));
                    Itemp(:,:,i)=Iset(:,:,i).*mask;
                
                end
                
                
                I=I+Itemp(:,:,i);
            end
            
            I=I./max(max(I));
            
        end
        
        function Ical=exposure_calibration(I,exposure_time,camera)
            % EXPOSURE_CALIBRATION adjusts an intensity image according to
            % the exposure time used on the camera. This is to correct for
            % nonlinear intensity response with exposure time.
            % INPUT
            % I - intensity image
            % exposure_time - in milliseconds
            % camera - string, only 'WinCamD' is currently stored
            % OUTPUT
            % Ical - corrected intensity image
            
            if strcmp(camera,'WinCamD')==1
                
                c=0.1149.*exposure_time.^(-1.006)+0.9031;
                Ical=I./c;
            end
        end
        
        function r=RMSD(data,ideal)
            % RMSD is a lazy shortcut for a simple formula
            
            r=sqrt(mean(data-ideal).^2);
            
        end
        
        function w=waist(z,w0,lambda,material)
            % WAIST calculates gaussian beam radius at a distance z from
            % the beam waist
            % INPUT
            % z - distance from beam waist
            % w0 - beam waist size
            % lambda - wavelength
            % material - string, use list_materials.m to see available
            % materials to be used as a medium
            % OUPUT
            % w - beam radius at z
            
            n=refractive_index(material,lambda);
            
            w=w0.*sqrt(1+(z./((pi.*w0.^2)./(lambda/n))).^2);
            
        end
        
        
    end
    
end

