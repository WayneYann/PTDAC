clockTime cellClockTime= clockTime();
scalar deltaTMin1=GREAT;
for(label i4=startat;i4<stopat;i4++)
{
	label celli(i4);//check
	cellTime_[celli]=0;
	cellClockTime.timeIncrement();
	scalarField c(this->nSpecie());
	scalarField c0(this->nSpecie());
	for(label i=0; i<this->nSpecie(); i++)
	{
		c[i] = Rxf_nr_buf[i4][8+i];
	}
	c0=c;
	scalar Ti = Rxf_nr_buf[i4][4];
	scalar pi = Rxf_nr_buf[i4][5];
	//scalar hi = Rxf_nr_buf[i4][6];
	scalar deltaTnr=Rxf_nr_buf[i4][6];
	scalar rhoi = Rxf_nr_buf[i4][7];
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
     	phiq[this->nSpecie() + 2] = Rxf_nr_buf[i4][6];
    }
    scalar tauC = Rxf_nr_buf[i4][2];
    scalar timeLeft = deltaTnr;
   
		Rxf_nr_buf[celli][8+2*this->nSpecie()]=0;//flag to tell source proc that DI performed in remote proc so u can add the chempoint
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
            Rxf_nr_buf[celli][2]=tauC;
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
 	    for(label i=0; i<this->nSpecie() ; i++)
		{
            //RQhiRecv[celli][i] = Rphiq[i];
			Rxf_nr_buf[celli][i+8+this->nSpecie()]=Rphiq[i];
		}
		Rxf_nr_buf[celli][ldxf-2] = Ti;
        Rxf_nr_buf[celli][ldxf-1] = pi;
		
		Rxf_nr_buf[celli][4]=phiq[this->nSpecie()];//Rxf_nr_buf[celli][4]=Ti;
		Rxf_nr_buf[celli][5]=phiq[this->nSpecie()+1];//Rxf_nr_buf[celli][5]=pi;
		Rxf_nr_buf[celli][6]=deltaTnr;//hi;
		Rxf_nr_buf[celli][7]=rhoi;
		
    //Compute the rate of reaction according to dc=c-c0
    //In the CFD solver the following equation is solved:
    //d(Yi*rho)/dt +convection+diffusion = RR*turbulentCoeff(=1 if not used)
    //Therefore, the unit of RR should be [kg/(m3.s)]
    scalarField dc = c - c0;
    for(label i=0; i<this->nSpecie(); i++)
    {
        Rxf_nr_buf[celli][i+8]=dc[i]*this->specieThermo_[i].W()/deltaTnr;
    }

    if(i4+1==stopat)//infoRecv[pack][1]+infoRecv[pack][1])
		{
			for(label i6=startat;i6<stopat;i6++) Rxf_nr_buf[i6][3]=deltaTMin1;
			//Rxf_nr_buf[celli][this->nSpecie()+5]=deltaTMin1;
			//deltaTMin1=Rxf_nr_buf[celli+1][3];			
		}
		Rxf_nr_buf[celli][1]+=cellClockTime.timeIncrement();	
}//End of loop over Resident cells
