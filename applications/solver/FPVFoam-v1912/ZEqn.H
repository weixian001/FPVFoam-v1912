	fvScalarMatrix ZEqn
	(
		(
		  fvm::ddt(rho, Zmix)
		+ fvm::div(phi, Zmix)
        - fvm::laplacian(turbulence->DZEff(), Zmix)
        )
    );

	ZEqn.relax();
	fvOptions.constrain(ZEqn);
	ZEqn.solve();
    
    Info<< "----------> Zmix min/max   = " << min(Zmix).value() << ", "
        << max(Zmix).value() << endl;

    //combustion->correct();
    
