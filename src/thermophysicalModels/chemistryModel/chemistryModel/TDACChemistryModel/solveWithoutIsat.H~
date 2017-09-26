 forAll(rho, celli)
{
    //DIClock_.timeIncrement();
    //label celli(cellIndexTmp[ci]);

    scalar rhoi = rho[celli];
    scalar Ti = this->thermo().T()[celli];
    //scalar hi = this->thermo().hs()[celli] + hc[celli];
    scalar pi = this->thermo().p()[celli];

	// Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);
    scalarField Rphiq(this->nEqns() + nAdditionalEqn);
    //scalarField phiq(this->nEqns());
    
    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);
    
    for(label i=0; i<this->nSpecie_; i++)
    {
     	c[i] = rhoi*this->Y_[i][celli]/this->specieThermo_[i].W();
        c0[i] = c[i];
        phiq[i] = this->Y()[i][celli];
    }
    phiq[this->nSpecie()]=Ti;
    phiq[this->nSpecie()+1]=pi;
    
    // Initialise time progress
    scalar timeLeft = deltaT[celli];

    //time step and chemical time step
    //scalar t = t0;
    //scalar tauC = this->deltaTChem_[celli];
    //scalar dt = min(deltaT, tauC);
    //scalar timeLeft = deltaT;

    //scalar cTot = 0.0;

	// Not sure if this is necessary
    Rphiq = Zero;
    
    clockTime_.timeIncrement();
    if (reduced)
    {
    	// Reduce mechanism change the number of species (only active)
    	mechRed_->reduceMechanism(c, Ti, pi);
        nActiveSpecies += mechRed_->NsSimp();
        nAvg++;
        scalar timeIncr = clockTime_.timeIncrement();
        reduceMechCpuTime_ += timeIncr;
        timeTmp += timeIncr;
    }
    // Calculate the chemical source terms
    while (timeLeft > SMALL)
    {
    	scalar dt = timeLeft;
        if (reduced)
        {
        	// completeC_ used in the overridden ODE methods
            // to update only the active species
            completeC_ = c;
            // Solve the reduced set of ODE
            this->solve
            (
            	simplifiedC_, Ti, pi, dt, this->deltaTChem_[celli]
            );
            for (label i=0; i<NsDAC_; i++)
            {
            	c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
            }
    	}
        else
        {
        	this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
        }
        timeLeft -= dt;
    }

    {
        scalar timeIncr = clockTime_.timeIncrement();
        solveChemistryCpuTime_ += timeIncr;
    	timeTmp += timeIncr;
    }
    
	// When operations are done and if mechanism reduction is active,
    // the number of species (which also affects nEqns) is set back
    // to the total number of species (stored in the mechRed object)
    if (reduced)
    {
     	this->nSpecie_ = mechRed_->nSpecie();
    }   
    deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            
    // Set the RR vector (used in the solver)
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->RR_[i][celli] =
        (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
    }
}//End of loop over all cells



