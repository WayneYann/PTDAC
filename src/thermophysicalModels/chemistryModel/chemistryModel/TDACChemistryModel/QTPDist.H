//Uncommented by ali on bahman 93
label celli(cellIndexTmp[ci]);

scalar rhoi = rho[celli];
scalar Ti = this->thermo().T()[celli];
//scalar hi = this->thermo().hs()[celli] + hc[celli];
scalar pi = this->thermo().p()[celli];

// Composition vector (Yi, T, p)
scalarField phiq(this->nEqns() + nAdditionalEqn);
scalarField Rphiq(this->nEqns() + nAdditionalEqn);
    
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
if (tabulation_->variableTimeStep())
{
	phiq[this->nSpecie() + 2] = deltaT[celli];
}
    
// Initialise time progress
scalar timeLeft = deltaT[celli];

//time step and chemical time step
scalar tauC = this->deltaTChem_[celli];
//scalar timeLeft = deltaT;

// Not sure if this is necessary
Rphiq = Zero;
    
clockTime_.timeIncrement();

//The tabulation algorithm try to retrieve the mapping
if (tabulation_->retrieve(phiq, Rphiq))
{
    o++;
    //nFound_ ++;
    for (label i=0; i<this->nSpecie(); i++)
    {
       	c[i] = rhoi*Rphiq[i]/this->specieThermo_[i].W();
    }
    searchISATCpuTime_ += clockTime_.timeIncrement();
		
	// Set the RR vector (used in the solver)
	for (label i=0; i<this->nSpecie_; i++)
	{
	    this->RR_[i][celli] =
	        (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
	}

}
else if(timeSteps_==tabulation_->numberOfInitDI()+1)
{
    oo++;
    
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
      
    // Set the RR vector (used in the solver)
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->RR_[i][celli] =
            (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
    }
}//END if iterator==1
else
{
    this->setTabulationResultsDistributed(celli);
    ss++;
    cellIDsTot[Pstream::myProcNo()].setSize(ss);
    cellIDsTot[Pstream::myProcNo()][ss-1]=celli;
    totcellTime[Pstream::myProcNo()].setSize(ss);
    totcellTime[Pstream::myProcNo()][ss-1]=cellTime_[celli];
}
