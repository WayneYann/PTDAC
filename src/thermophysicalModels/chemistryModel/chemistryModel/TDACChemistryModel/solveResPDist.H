totNWhileIter=0;
clockTime cellResClockTime= clockTime();
for(label iii=totSend[Pstream::myProcNo()];iii<cellIDsTotSize;iii++)// cellIDsTot[Pstream::myProcNo()][i]--> cell number 
{
    label celli(cellIDsTot[Pstream::myProcNo()][iii]);
    cellTime_[celli]=0;
    cellResClockTime.timeIncrement();
    scalar rhoi = rho[celli];
    scalar Ti = this->thermo().T()[celli];
   // scalar hi = this->thermo().hs()[celli] + hc[celli];
    scalar pi = this->thermo().p()[celli];

	scalarField phiq(this->nEqns() + nAdditionalEqn);
    scalarField Rphiq(this->nEqns() + nAdditionalEqn);
    for(label i=0; i<this->nSpecie(); i++)
    {
        phiq[i] = this->Y()[i][celli];
    }
    phiq[this->nSpecie()]=Ti;
	phiq[this->nSpecie()+1]=pi;
	if (tabulation_->variableTimeStep())
    {
     	phiq[this->nSpecie() + 2] = deltaT[celli];
    }

    //Species are stored in mass fraction in the cells
    //c arrays indicate the molar concentration of the species
    // c = (Y * rho)/W [kmol/m3]
    //phiq array store the mass fraction, the temperature and pressure
    //(i.e. the composition) of the query point
    scalarField c(this->nSpecie());
    scalarField c0(this->nSpecie());
    for(label i=0; i<this->nSpecie(); i++)
    {
        c[i] = rhoi*phiq[i]/this->specieThermo_[i].W();;
    }
    //store the initial molar concentration to compute dc=c-c0
    c0 = c;

    //time step and chemical time step
    scalar tauC = this->deltaTChem_[celli];
    scalar timeLeft = deltaT[celli];

    //retrieveClock_.timeIncrement();
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
        this->RR()[i][celli] = dc[i]*this->specieThermo_[i].W()/deltaT[celli];
    }

    cellTime_[celli]+=cellResClockTime.timeIncrement();
}//End of loop over Resident cells
