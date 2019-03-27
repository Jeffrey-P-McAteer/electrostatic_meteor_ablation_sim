; EPPIC parameter file (also readable by IDL)
	title = 'Meteor Ablation 3D Z path'
;Mesh size and Grid size
	ndim_space = 3
	nsubdomains = 32
	nx = 16
        dx = 0.005
	ny = 512
        dy = 0.005
	nz = 512
        dz = 0.005
; Time step, number of time steps to run, and time-steps between output:
	dt = 5e-9
        nt = 54000
	nout = 3000
	hdf_output_arrays = 0
	nout_avg = 1
; Fraction of particles to output (10240 originally)
	npout = 1e6
; Define the dielectric of free space (default epsilon=8.8542e-12)
	eps = 8.85418728e-12
; Bz in MKS (Tesla)
        Bx = 40000E-9
        Ey0_external = -1.6
; Damping width: (try fwidth=3)
	fwidth = 0
; Damping steepness: (try fsteep=3)
	fsteep = 0
; Time steps beteen dumps
        iwrite = 70000
; Start from t=0 (iread=0) or dump (iread=1)
	iread = 0
; Limit the frequency of outputing divj:
	divj_out_subcycle = -1
; Neutral parameters:
;      vth_neutral = 0.1
;      m_neutral = 1.0
;      v0_neutral = 10
; No EFIELD
       efield_algorithm = 10
       efield_zero_dc = 0
; particle boundary_type_xyz: 0=periodic, 1=inject, 2=open
       boundary_type_x = 0
       boundary_type_y = 0
       boundary_type_z = 1
; field_boundary type open=-1, periodic=0, dirichlet=1, neumann=2 for INJECT
       field_boundary_type_xl = 0;
       field_boundary_type_xh = 0;
       field_boundary_type_yl = 0;
       field_boundary_type_yh = 0;
       field_boundary_type_zl = 2;
       field_boundary_type_zh = 2;
; field boundary value
       field_boundary_value_xl = 0;
       field_boundary_value_xh = 0;
       field_boundary_value_yl = 0;
       field_boundary_value_yh = 0;
       field_boundary_value_zl = 0;
       field_boundary_value_zh = 0;
; Unscale the density
       unscale_density = 1
       inject_prop_dc = 0
; Number of distributions to initialize (just neutrals to start)
       ndist = 3

; Distribution 0: Ablating Neutrals
	dist = 0
; number of particles, density, mass, charge, thermal velocity and drift speed
	npd0 = 1
	part_pad0 = 6000000
	n0d0 = 2e3
; neutral mass: Everything is normalized to this
        md0 = 3.8E-26
	qd0 = 0
; Subcycle these particles
	subcycle0 = 1
; No injection
        n0rhsd1 = 0
        n0lhsd1 = 0
; The coll_rate
	coll_rate0 = 1
	coll_start_time0 = 0
	coll_type0 = 10
; beta_model = 0 is unity, beta_model = 1 is jones model, 2 is vondrak
	beta_model0 = 2
	background_neutral_dens0 = 1E19
	crosssec_m_model0 = 0
; if crosssec_m_model = 1, coll_cross_section should not matter, 1 is jones model 5.1067e-20 (20km/s)
        coll_cross_section0 = 2.933e-20
        coll_create_a_id0 = 1
	coll_create_b_id0 = 2
; thermal velocity 10eV=593994 m/s or 342368 for 9.1e-30 kg electron
        coll_create_vthb = 342368
; Atmospheric Neutrals N2 (4.7E-26 kg) but 100 x lighter(for e-)
	massd_neutral0 = 4.7e-26
        vthd_neutral0 = 272.8
        vx0d_neutral0 = 0
        vy0d_neutral0 = 0
        vz0d_neutral0 = 40000
; Thermal velocities
	vxthd0 = 951
	vythd0 = 951
	vzthd0 = 951
	vx0d0 = 0
	vy0d0 = 0
	vz0d0 = 0
        vz0rhsd0 = 0
        vz0lhsd0 = 0
; Output parameters
	pnvx0 = 64
	pnvy0 = 64
	pnvz0 = 64
	pvxmin0 = -80000
	pvxmax0 = 80000
	pvymin0 = -80000
	pvymax0 = 80000
	pvzmin0 = -80000
	pvzmax0 = 80000
        denft_out_subcycle0 = 1
	den_out_subcycle0 = 1
	vdist_out_subcycle0 = 1
	part_out_subcycle0 = 1
	flux_out_subcycle0 = 1
	nvsqr_out_subcycle0 = 16
; Ablation parameters
	init_dist0 = 11
        creation_rate0 =  2.5e10
	create_radius0 = 0.0001
        create_posx0 = 0.5
        create_posy0 = 0.5
        create_posz0 = 0.5
	create_vth0 = 951
	create_v00 = 0

; Second distribution: Post-Collision Ions
	dist = 1
; number of particles, density, mass, charge, thermal velocity and drift speed
        npd1 = 500000
	part_pad1 = 10
	n0d1 = 2e3
; neutral mass: Everything is normalized to this
        md1 = 3.8E-26
	qd1 = 1.6e-19
; Subcycle these particles
	subcycle1 = 1
; Injection parameters (ions should be injected)
        n0rhsd1 = 500000
        n0lhsd1 = 500000
; The coll_rate
	coll_rate1 = 1
	coll_start_time1 = 0
	coll_type1 = 8
	background_neutral_dens1 = 1E19
	crosssec_m_model1 = 0
; if crosssec_m_model = 1, coll_cross_section should not matter, 1 is jones model
	coll_cross_section1 = 5.1067e-20
; Atmospheric Neutrals N2 (4.7E-26 kg)
        massd_neutral1 = 4.7E-26
        vthd_neutral1 = 272.8
        vx0d_neutral1 = 0
        vy0d_neutral1 = 0
        vz0d_neutral1 = 40000
; Thermal velocities
	vxthd1 = 303.4
	vythd1 = 303.4
	vzthd1 = 303.4
	vx0d1 = 0
	vy0d1 = 0
	vz0d1 = 40000
        vz0rhsd1 = 40000
        vz0lhsd1 = 40000
; Output parameters
	pnvx1 = 64
	pnvy1 = 64
	pnvz1 = 64
	pvxmin1 = -80000
	pvxmax1 = 80000
	pvymin1 = -80000
	pvymax1 = 80000
	pvzmin1 = -80000
        pvzmax1 = 80000
        denft_out_subcycle1 = 1
	den_out_subcycle1 = 1
	vdist_out_subcycle1 = 1
	part_out_subcycle1 = 1
	flux_out_subcycle1 = 1
	nvsqr_out_subcycle1 = 16
; Ablation parameters (usually 11 but we have background ionosphere (0=rand_flat)
	init_dist1 = 0

; Third distribution: post collision electrons (1/40th mass of ions)
	dist = 2
; number of particles, density, mass, charge, thermal velocity and drift speed
	npd2 = 500000
	part_pad2 = 10
	n0d2 = 2e3
; electron mass: 10x e- mass (9.1E-31kg -> 9.1E-30kg)
        md2 = 9.1E-30
	qd2 = -1.6E-19
; Subcycle these particles
	subcycle2 = 1
; Injection parameters (electrons should be injected)
        n0rhsd2 = 500000
        n0lhsd2 = 500000
; The coll_rate
        coll_rate2 = 1
        coll_start_time2 = 0
        coll_type2 = 8
        background_neutral_dens2 = 1E19
        crosssec_m_model2 = 2
; if crosssec_m_model2 = 1, coll_cross_section should not matter, 1 is jones model, 2 is electron Engelhardt 1964 model
        coll_cross_section2 = 1.0E-20
; Atmospheric Neutrals N2 (4.7E-26 kg) but 100 x lighter(for e-)
        massd_neutral2 = 9.4e-28
        vthd_neutral2 = 272.8
        vx0d_neutral2 = 0
        vy0d_neutral2 = 0
        vz0d_neutral2 = 40000
; Thermal velocities
	vxthd2 = 19600
	vythd2 = 19600
	vzthd2 = 19600
	vx0d2 = 0
	vy0d2 = 0
	vz0d2 = 40000
        vz0rhsd2 = 40000
        vz0lhsd2 = 40000
; Output parameters
	pnvx2 = 64
	pnvy2 = 64
	pnvz2 = 64
	pvxmin2 = -80000
	pvxmax2 = 80000
	pvymin2 = -80000
	pvymax2 = 80000
	pvzmin2 = -80000
	pvzmax2 = 80000
        denft_out_subcycle2 = 1
	den_out_subcycle2 = 1
	vdist_out_subcycle2 = 1
	part_out_subcycle2 = 1
	flux_out_subcycle2 = 1
	nvsqr_out_subcycle2 = 16
; Ablation parameters (usually 11 but we have background ionosphere (0=rand_flat)
	init_dist2 = 0

; Simulation to real units
    vsim_to_kmsec = 0.001
    lsqsim_to_msq = 1
