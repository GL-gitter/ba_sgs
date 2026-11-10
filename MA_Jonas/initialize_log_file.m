function initialize_log_file(mat, fileID)

fprintf(fileID, '#####################################################\n');
fprintf(fileID, '############## POISSON OPTIONS ######################\n');
fprintf(fileID, 'iteration error: \t %1.3e\n', mat.poisson.opt.iter_err);
fprintf(fileID, 'maximum iteration num:\t %d\n', mat.poisson.opt.iter_max);
fprintf(fileID, 'update coefficient:\t %.2f\n', mat.poisson.opt.alpha);
fprintf(fileID, ['solver-type:\t', mat.poisson.opt.solve, '\n']);

fprintf(fileID, '#####################################################\n');
fprintf(fileID, '############## WIGNER OPTIONS #######################\n');
fprintf(fileID, 'number of kvals:\t %d\n', mat.wigner.params.Nk);
fprintf(fileID, 'maximum kvals:\t %.3E\n', mat.wigner.params.k_max);
fprintf(fileID, 'amplitude of the cap:\t %.2E\n', mat.wigner.params.cap_ampl);
fprintf(fileID, 'ratio of the cap:\t %.2E\n', mat.wigner.params.cap_rel);
fprintf(fileID, 'exponent of the cap:\t %d\n', mat.wigner.params.cap_n);

fprintf(fileID, '#####################################################\n');
fprintf(fileID, '############## SIGMA OPTIONS ########################\n');
fprintf(fileID, 'number of basis-funcs:\t %d\n', mat.sigma.params.N);
fprintf(fileID, 'discretiaztion points Ny:\t %.3E\n', mat.sigma.params.Ny);
fprintf(fileID, 'amplitude of the cap:\t %.2E\n', mat.sigma.params.cap_ampl);
fprintf(fileID, 'ratio of the cap:\t %.2E\n', mat.sigma.params.cap_rel);
fprintf(fileID, 'exponent of the cap:\t %d\n', mat.sigma.params.cap_n);

fprintf(fileID, '#####################################################\n');
fprintf(fileID, 'discretization width x:\t %.2E\n', mat.dx);
fprintf(fileID, 'discretization width x:\t %.2E\n', mat.dy);

