label totNWhileIter=0; 
label pack=0;
label bp=infoRecv[pack][1];
scalar deltaTMin1=0;
if(cTOTRecv.size()!=0)
{
    //bug fix
	//deltaTMin1=cTOTRecv[pack][this->nSpecie()+5];
	deltaTMin1=cTOTRecv[bp-1][this->nSpecie()+5];
}
clockTime cellClockTime= clockTime();
for(label i4=0;i4<cTOTRecv.size();i4++)
{
    label celli(i4);
    cellTimeRecv[celli]=0;
    cellClockTime.timeIncrement();
    scalarField c(this->nSpecie());
    scalarField c0(this->nSpecie());
    for(label i=0; i<this->nSpecie(); i++)
    {
        c[i] = cTOTRecv[celli][i];
    }
    c0=c;
    scalar Ti = cTOTRecv[celli][this->nSpecie()];
    scalar pi = cTOTRecv[celli][this->nSpecie()+1];
    //scalar hi = cTOTRecv[celli][this->nSpecie()+2];
    scalar deltaTnr=cTOTRecv[celli][this->nSpecie()+2];
    scalar rhoi = cTOTRecv[celli][this->nSpecie()+3];
    scalarField phiq(this->nEqns() + nAdditionalEqn);
    scalarField Rphiq(this->nEqns() + nAdditionalEqn);
    for(label i=0; i<this->nSpecie(); i++)
    {
        phiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
    }
    phiq[this->nSpecie()]=Ti;
    phiq[this->nSpecie()+1]=pi;
    //store the initial molar concentration to compute dc=c-c0
    c0 = c;
    //time step and chemical time step
    scalar tauC = cTOTRecv[celli][this->nSpecie()+4];
    DeltatchemRecv[celli] = tauC;
    scalar timeLeft = deltaTnr;
    //retrieveClock_.timeIncrement();
    
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
                simplifiedC_, Ti, pi, dt, tauC
            );
            for (label i=0; i<NsDAC_; i++)
            {
              	c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
            }
            //tauC=this->deltaTChem_[celli];
       }
       else
       {
           this->solve(c, Ti, pi, dt, tauC);
           //tauC=this->deltaTChem_[celli];
       }
            timeLeft -= dt;
            //Rxf_nr_buf[celli][2]=tauC;
    }
    {
        scalar timeIncr = clockTime_.timeIncrement();
        solveChemistryCpuTime_ += timeIncr;
        timeTmp += timeIncr;
    }
    if (reduced)
    {
        this->nSpecie_ = mechRed_->nSpecie();
    }
    deltaTMin1 = min(tauC, deltaTMin1);
    //Rphiq array store the mapping of the query point
    //scalarField Rphiq(this->nSpecie());
    //Transform c array containing the mapping in molar concentration [mol/m3]
    //to Rphiq array in mass fraction
    for(label i=0; i<this->nSpecie(); i++)
    {
       Rphiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
    }
    if (tabulation_->variableTimeStep())
    {
        Rphiq[Rphiq.size()-3] = Ti;
        Rphiq[Rphiq.size()-2] = pi;
        Rphiq[Rphiq.size()-1] = deltaTnr;
    }
    else
    {
     	Rphiq[Rphiq.size()-2] = Ti;
       	Rphiq[Rphiq.size()-1] = pi;
    }
    
   	for(label i=0; i<this->nSpecie(); i++)
	{
        RQhiRecv[celli][i] = Rphiq[i];
	}
	RQhiRecv[celli][this->nSpecie()]=Ti;
	RQhiRecv[celli][this->nSpecie()+1]=pi;
	RQhiRecv[celli][this->nSpecie()+2]=0;//flag to tell source proc that DI performed in remote proc so u can add the chempoint

	scalarField dc = c - c0;
    for(label i=0; i<this->nSpecie(); i++)
    {
        RRRecv[celli][i] = dc[i]*this->specieThermo_[i].W()/deltaTnr;
    }
    deltaTminRecv[pack]=deltaTMin1;
    if(celli+1==bp && pack+1<infoRecv.size())//infoRecv[pack][1]+infoRecv[pack][1])
    {
        deltaTMin1=cTOTRecv[celli+1][this->nSpecie()+5];
        pack++;
        bp+=infoRecv[pack][1];
    }
    cellTimeRecv[celli]+=cellClockTime.timeIncrement();
}//END for(label i=0;i<NonResidentSize;i++)
