clockTime cellResClockTime= clockTime();
for(label i4=startat;i4<stopat;i4++)
{
	label celli(Rxf[i4][0]);//check
	cellTime_[celli]=0;
	cellResClockTime.timeIncrement();
	scalarField c(this->nSpecie());
	scalarField c0(this->nSpecie());
	for(label i=0; i<this->nSpecie(); i++)
	{
		c[i] = Rxf[i4][8+i];
	}
	c0=c;
	scalar Ti = Rxf[i4][4];
	scalar pi = Rxf[i4][5];
//	scalar hi = Rxf[i4][6];
scalar deltaTr=Rxf[i4][6];
	scalar rhoi = Rxf[i4][7];
	scalarField phiq(this->nEqns() + nAdditionalEqn);
    scalarField Rphiq(this->nEqns() + nAdditionalEqn);

	for(label i=0; i<this->nSpecie(); i++)
	{
		phiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
	}
	phiq[this->nSpecie()]=Ti;
	phiq[this->nSpecie()+1]=pi;
	if (tabulation_->variableTimeStep())
    {
     	phiq[this->nSpecie() + 2] = Rxf[i4][6];
    }
    scalar tauC = Rxf[i4][2];

    scalar timeLeft = deltaTr;

	// Not sure if this is necessary
    Rphiq = Zero;
    
    clockTime_.timeIncrement();
    
    //The tabulation algorithm try to retrieve the mapping
    if (tabulation_->retrieve(phiq, Rphiq))
    {
        for (label i=0; i<this->nSpecie(); i++)
        {
        	c[i] = rhoi*Rphiq[i]/this->specieThermo_[i].W();
        }
        searchISATCpuTime_ += clockTime_.timeIncrement();
    }
    //Retrieve has failed. The mapping is computed
    else
    {
        // Store total time waiting to attribute to add or grow
        scalar timeTmp = clockTime_.timeIncrement();
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
                tauC=this->deltaTChem_[celli];
            }
            else
            {
                this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
            tauC=this->deltaTChem_[celli];
            }
            timeLeft -= dt;
        }
        {
            scalar timeIncr = clockTime_.timeIncrement();
            solveChemistryCpuTime_ += timeIncr;
            timeTmp += timeIncr;
        }
        // If tabulation is used, we add the information computed here to
        // the stored points (either expand or add)
        if (tabulation_->active())
        {
        	forAll(c, i)
            {
                Rphiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
            }
            if (tabulation_->variableTimeStep())
            {
                Rphiq[Rphiq.size()-3] = Ti;
                Rphiq[Rphiq.size()-2] = pi;
                Rphiq[Rphiq.size()-1] = deltaT[celli];
            }
            else
            {
                Rphiq[Rphiq.size()-2] = Ti;
                Rphiq[Rphiq.size()-1] = pi;
            }
            label growOrAdd =
                tabulation_->add(phiq, Rphiq, rhoi, deltaT[celli]);
            if (growOrAdd)
            {
                this->setTabulationResultsAdd(celli);
            	addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
            }
            else
            {
                this->setTabulationResultsGrow(celli);
                growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
            }
        }

        // When operations are done and if mechanism reduction is active,
        // the number of species (which also affects nEqns) is set back
        // to the total number of species (stored in the mechRed object)
        if (reduced)
        {
            this->nSpecie_ = mechRed_->nSpecie();
        }
        deltaTMin = min(this->deltaTChem_[celli], deltaTMin);
        
    }

    //Compute the rate of reaction according to dc=c-c0
    //In the CFD solver the following equation is solved:
    //d(Yi*rho)/dt +convection+diffusion = RR*turbulentCoeff(=1 if not used)
    //Therefore, the unit of RR should be [kg/(m3.s)]
    scalarField dc = c - c0;
    for(label i=0; i<this->nSpecie(); i++)
    {
        this->RR()[i][celli] = dc[i]*this->specieThermo_[i].W()/deltaTr;
    }

    cellTime_[celli]+=cellResClockTime.timeIncrement();
}//End of loop over Resident cells



