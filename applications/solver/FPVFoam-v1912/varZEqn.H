	fvScalarMatrix varZEqn
        (
		(
                    fvm::ddt(rho, varZ)
                    + fvm::div(phi, varZ)
                    - fvm::laplacian(turbulence->DZEff(), varZ)
                    - 2.0*turbulence->Production()*magSqr(fvc::grad(Zmix))
                    + turbulence->Destruction()*varZ
                )
	);

        varZEqn.relax();
	fvOptions.constrain(varZEqn);
	varZEqn.solve();
	//fvOptions.correct(varZ);
	bound(varZ, 0.0);

    Info<< "----------> varZ min/max   = " << min(varZ).value() << ", "
        << max(varZ).value() << endl;

	
	

    combustion->correct();
    
