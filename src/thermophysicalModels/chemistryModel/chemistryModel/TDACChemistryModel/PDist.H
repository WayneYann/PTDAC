//Info<<"MPI Variables Declaration"<<endl;
	//iniTime.timeIncrement();
	label ss=0;
	//matrix with rows=Processor ID and columns=Cell IDs of that processor which is not retrieved
	//the matrix first declared just like a matrix with 0 columns and then the size of the matrix is changed during the simulation in the retrieve2.H 
	//these cells were not retrieved from the local tables (in retrieve2.H i.e. QT part) and shpuld be distributed 
	//the ID i.e. number of the cell in mesh should be stored to be used after recevieng
	//this matrix is accessible just by 1 processor itself. 
	List<List<scalar> > cellIDsTot(Pstream::nProcs(), List<scalar>(0) ); 
	
	//matrix with rows=Processor ID and columns=time spent to resolve the cell in the previous step (i.e. cellTime_)
	//the matrix first declared just like a matrix with 0 columns and then the size of the matrix is changed during the simulation in the retrieve2.H 
	//this matrix is accessible just by 1 processor itself. 
	List<List<scalar> > totcellTime(Pstream::nProcs(), List<scalar>(0) );
	
	//matrix with number of rows=number of Processors and number of columns=4*number of packs to be sent by sender 
	//during the simulation number of rows may be reduced 
	//the final size of rows is se: number of sender processors
	//this matrix is used to store 
	//column0:Processor ID of sender
    //column1:Processor ID of receiver
	//column2:Number of cells that will be exchanged
	//column3:Pack ID, i.e. each sender may send more than 1 pack. there maybe 1 sender and 10 receiver. so sender sends 10 packs
	//so sendRecepie has 40 columns. 
	//e.g. row0,coulmn 4 is again processor ID of sender (which is equal to row0,column 0
	//column 5:Processor ID of receiver (which is NOT equal to c.1)
	//column 6:NO. cells and column 7:Pack ID which is now equal to 2.
	List<List<label> > sendRecepie(Pstream::nProcs(), List<label>(4) );
	
	//this is to store size of sendRecepie list i.e. number of rows
	//as stated before, at first it is equal to number of processors
	//during the simulation it is reduced to number of processors that are senders 
	label sendRecepieSize=Pstream::nProcs();
	
	//nFSC=number of flame surface cells for each processor
	//in this implementation it store total estimated time needed to resolve the cells that are not retrieved in each processor
	//i.e. sum(totCellTime)
	//i.e. nFSC[0]=sum(totCellTime[0][]) = sum of or total estimated time needed to resolve the cells that are not retrieved in proc. 0
	SortableList<scalar> nFSC(Pstream::nProcs());
	
	//A scalar to store mean of nFSC 
	//i.e. average over nFSC[i] or average of chemistry workload or average of time needed to resolve all cells remained after QT in all processors
	//MeanNFSC=sum(nFSC)/nFSC.size()
	scalar MeanNFSC=0;
	
	//to store Processors IDs after sorting
	//when nFSC is sorted nFSC[0] is not total estimated time needed to resolve the cells that are not retrieved in proc. 0
	//so we store processor IDs in Indices to know after sorting nFSC, nFSC[0] is related to which processor
	//so after sorting e.g. nFSC[0]=sum(totcellTime[Indices[0][]) not sum(totcellTime[0][])
	labelList Indices(Pstream::nProcs()); 
	
	//Total number of cells which is in flame surface in each processor
	//in this implementation: Total number of cells which is not retrieved in the QT part in each processor
	//i.e. cellIDsTot[Pstream::myProcNo()].size() or number of columns of each row
	label cellIDsTotSize=0;
	
	//Total time which is estimated that total number of cells of each proc. need to be resolved
	//i.e. sigma i (totcellTime[Pstream::myProcNo()][i])
	scalar totcellTimeSize=0;
	
	//labelList cellIDsResidentSize(Pstream::nProcs(),0);
	//label NonResidentSize=0;
	
	//A list whith the size of nProc to inform each processor that "are you sender or receiver?" in the PDist algorithm
	//0: processor that does nothing 
	//-1: if a processor is reciever 
	//1: if a process is sender
	//this is determined by average load or MeanNFSC; 
	//if nFSC[i]>=1.1*MeanNFSC then the processor is sender
	//if nFSC[i]<=0.9*MeanNFSC then the processor is receiver
	//else --> the processor does nothing
	List<scalar> amISenderOrReciever(Pstream::nProcs(),0.0);
	
	
	//these two lists are used to store total number of cells that each processor sends or recieves
	List<label> totSend(Pstream::nProcs(),0.0);
	List<label> totRecv(Pstream::nProcs(),0.0);
	
	label TotalCellsDI0=0;
	scalar deltaTminSendOrRecv=0;	
	/*
		cTOTSendOrRecv to store composition of cells that need to be distributed
		at first the size is (1,number of species +6)
		then e.g. for proc number 0 the rows will be increased to number of cells in each pack i.e. sendRecepie[0][2] or sendRecepie[0][6] and so on
		
		jj is cell counter
		for(label ii=0; ii<this->nSpecie(); ii++)
		{
			cTOTSendOrRecv[jj][ii] = c[ii]; c[ii] is concentration of specie ii in the cell jj
		}
		cTOTSendOrRecv[jj][this->nSpecie()]=this->thermo().T()[cellIDsTot[Pstream::myProcNo()][jjj]];
		cTOTSendOrRecv[jj][this->nSpecie()+1]=this->thermo().p()[cellIDsTot[Pstream::myProcNo()][jjj]];
		cTOTSendOrRecv[jj][this->nSpecie()+2]=this->thermo().hs()[cellIDsTot[Pstream::myProcNo()][jjj]] + hc[cellIDsTot[Pstream::myProcNo()][jjj]];
		cTOTSendOrRecv[jj][this->nSpecie()+3]=rho[cellIDsTot[Pstream::myProcNo()][jjj]];
		cTOTSendOrRecv[jj][this->nSpecie()+4]=this->deltaTChem_[cellIDsTot[Pstream::myProcNo()][jjj]];
		cTOTSendOrRecv[jj][this->nSpecie()+5]=deltaTMin;
	*/
	List<List<scalar> > cTOTSendOrRecv(1, List<scalar>(this->nSpecie()+6,0.0));
	
	//when cTOTSendOrRecv is distributed and received by a receiver proc. the composition is stored in cTOTRecv in each processor
	//so cTOTRecv has number of rows that is the size of totRecv or total number of cells that each proc receives
	//and the number of columns are equal to number of columns of cTOTSendOrRecv
	//see solveNonResPDist.H
	List<List<scalar> > cTOTRecv(totRecv[Pstream::myProcNo()], List<scalar>(this->nSpecie()+6,0.0));
	
	//when the reaction mapping is calculated on each processr for non resident cells (i.e. cells that are received), RQhiRecv is used to store these
	//N.Rows is equal to total number of cells that are recieved, i.e. totRecv
	//N.Columns is equal to dimension of composition space +1 -->to store concentrations, T, P and ...
	//see solveNonResPDist.H
	List<List<scalar> > RQhiRecv(totRecv[Pstream::myProcNo()], List<scalar>(this->nEqns()+1,0.0));  
	
	
	//List<scalar> normalisedGradTSendOrRecv(1,0.0);
	//List<List<scalar> > RRSendOrRecv(1, List<scalar>(this->nSpecie(),0.0));   
	//List<List<scalar> > RQhiSendOrRecv(1, List<scalar>(this->nEqns()+1,0.0));   
	
	
	
	
	
	
	
	//List<scalar> normalisedGradTRecv(totRecv[Pstream::myProcNo()],0.0);
	
	//this is to store RR or reaction rate of specie i 
	//when we have cTOTRecv and RQhiRecv so we can calculate RR on host processor
	/*
		scalarField dc = c - c0;
		for(label i=0; i<this->nSpecie(); i++)
		{
			RRRecv[celli][i] = dc[i]*Wi[i]/deltaT;
		}
	*/
	//see solveNonResPDist.H
	List<List<scalar> > RRRecv(totRecv[Pstream::myProcNo()], List<scalar>(this->nSpecie(),0.0));  
	
	//to store deltaTmin of received cells on host processr for non resident queries
	//see solveNonResPDist.H
	//at first the size is 1 and then will be increased to number of packs that each proc. recieves i.e. numRecv
	//packs means that each receiver proc. may receive cells from many proc. each sender packs the cells and sends the pack
	//so e.g. proc. 100 may receives from proc. 0 and proc.1
	List<scalar> deltaTminRecv(1,0.0);
	
	//to store deltatchem of received cells on host processr for non resident queries
	//DeltatchemRecv[celli] = tauC;
	//see solveNonResPDist.H
	List<scalar> DeltatchemRecv(totRecv[Pstream::myProcNo()],0.0);
	
	//to store time to resolve each cell on host processor
	//see solveNonResPDist.H
	List<scalar> cellTimeRecv(totRecv[Pstream::myProcNo()],0.0);
	
	//to store some info. of packs on host proc. after receiving cTOTSendOrRecv
	//this is a matrix with n.rows equal to numRecv or number of packs that each proc. receivves (at first is 1 and then will be increased)
	//infoRecv[numRecv][0]=sendRecepie[i][4*j]; i.e. the ID of the processor that has sended this pack
	//infoRecv[numRecv][1]=cTOTSendOrRecv.size(); i.e. total number of cells
	List<List<label> > infoRecv(1, List<label>(2,10000));  
	
	//it is not used
	List<scalar> Deltatchem(totSend[Pstream::myProcNo()],0.0);
	//List<scalar> cellTimeSendOrRecv(totSend[Pstream::myProcNo()],0.0);
	
	/*
	when host proc. resolve nonresident cells, the reaction mapping and all needed info. should be send back to the original proc (i.e. sender)
	so receivers pack the info in infoPack and send them to the original processors one by one
	the size of the matrix is as below
	n. columns= number of cells that each receiver has received in pack i
	n.rows
	=
	n. species (to store RR or reaction rate of specie i)
	+
	n. species + 1 + 1 i.e. n. eqns (to store RQhi or reaction mapping i.e. new concentration of n. species + T + P)
	+
	4 to store:
	infoPack[0][j]=cellTimeRecv[j]
	infoPack[1][j]=DeltatchemRecv[j];
	infoPack[2][j]=deltaTminRecv[j];
		
	*/
	List<List<scalar> > infoPack(this->nSpecie()+this->nEqns()+4, List<scalar>(1,0.0));   

	label sendSize=0;
	label numRecv=0;// number of packs that each processor recieves
	label lastIt=0;
	label lastIt3=0;
	label sendSize3=0;
	
	//some wierd variables!
	label o=0;
	label oo=0;
	
	//end of declarations 
	
	
	//qTClock.timeIncrement();
    //commented on bahman 93
	//forAll(rho, celli) //if you want to use nonRandom ISAT
	forAll(rho, ci)//if you want to use random ISAT
    {
    	clockTime_cellTime.timeIncrement();
        #include "QTPDist.H"
        cellTime_[celli]=clockTime_cellTime.timeIncrement();
    }//END for all cells
    //QTTime_+=qTClock.timeIncrement();

	cellIDsTotSize=cellIDsTot[Pstream::myProcNo()].size();
	Pout<<"n cells failed to retrieve ="<<cellIDsTotSize<<endl;
	
    //totcellTimeSize=sum(totcellTime[Pstream::myProcNo()]);
    for(label ti=0;ti<cellIDsTotSize;ti++)
    {
       totcellTimeSize+= totcellTime[Pstream::myProcNo()][ti];
    }


        //commClock_.timeIncrement();
        Pstream::gatherList(totcellTime);
        //communicateTime_+=commClock_.timeIncrement();
        //Pout<<"first communication lasts:"<<communicateTime_<<" s"<<endl;

        clockTime clockT2=clockTime();
        clockTime clockT6=clockTime();
        clockT2.timeIncrement();
        clockT6.timeIncrement();
        clockTime temp=clockTime();
        
        //the main PDist algorithm which is done by master node
        //1. identify sender or recievers procs
        //2. distribute the load based on cell solution time
        if (Pstream::myProcNo() == Pstream::masterNo())
        {
            //phase 1: identification of senders or recievers
            forAll(nFSC, c2i)
            {
                nFSC[c2i]=sum(totcellTime[c2i]);
            }
            //sort in increasing order
            nFSC.sort();
            Indices=nFSC.indices();
            MeanNFSC=sum(nFSC)/nFSC.size();
            label se=0; //number of sender procs
            label re=0; //number of reciever procs
            	
            forAll(nFSC, ci) // Processors whith nFSC>1.1 are senders and those with nFSC<0.9 are recievers
            {
                if(nFSC[ci]>=1.1*MeanNFSC)
                {
                    nFSC[ci] -=MeanNFSC;
                    //nFSC[ci] -=MeanNFSC;
                    //nFSC[ci] =nFSC[ci]*Ratio[Indices[ci]];
                    amISenderOrReciever[Indices[ci]]=1;
                    se++;
                }
                else if(nFSC[ci]<=0.9*MeanNFSC)
                {
                    nFSC[ci] -=MeanNFSC;
                    //nFSC[ci] -=MeanNFSC;
                    //nFSC[ci] =nFSC[ci]/Ratio[Indices[ci]];
                    amISenderOrReciever[Indices[ci]]=-1;
                    re++;
                }
                else
                {
                    nFSC[ci]=0;
                    amISenderOrReciever[Indices[ci]]=0;
                }
            }
            if (MeanNFSC==0)
            {
                amISenderOrReciever=0;
                se=0;
            }
            sendRecepie.setSize(se, List<label>(4) );
            sendRecepieSize=se;
            label ci=nFSC.size()-1;
            label cii=0;
            label sr=0;
            scalar res=0;
            
            //phase 2:
            //here is the programmed distribution algorithm based on previous time step solution time of each cell
            //sendRecepie will show each cell will go to which processor
            while(amISenderOrReciever[Indices[ci]]==1)
            {
                //if(amISenderOrReciever[Indices[ci]]!=1)	break;
                label noc=0;
                label e=0;
                sr=0;
                while(cii<nFSC.size())
                {
                    if(amISenderOrReciever[Indices[cii]]!=-1)
                    {
                        ci--;
                        break;
                    }
                    res=nFSC[ci]+nFSC[cii];
                    if(res>=0)
                    {
                        /*if(res<=0.1*MeanNFSC)
		    				{
		    					sendRecepie[nFSC.size()-1-ci][4*sr]=Indices[ci];
		    					sendRecepie[nFSC.size()-1-ci][4*sr+1]=Indices[cii];
		    					label tnoc=0;
		    					label ni=0;
		    					while(ni<=fabs(nFSC[ci]))
		    					{
		    						ni+=totcellTime[Indices[ci]][noc];
		    						noc++;
		    						tnoc++;
		    					}
		    					sendRecepie[nFSC.size()-1-ci][4*sr+2]=tnoc;
		    					//sendRecepie[nFSC.size()-1-ci][4*sr+2]=nFSC[ci];
		    					sendRecepie[nFSC.size()-1-ci][4*sr+3]=sr+1;
							//totSend[Indices[ci]]+=tnoc;
							//totRecv[Indices[cii]]+=tnoc;
		    					ci--;
		    					cii++;
		    					break;
		    				}
		    				else
		    				{*/
                        sendRecepie[nFSC.size()-1-ci][4*sr]=Indices[ci];
                        sendRecepie[nFSC.size()-1-ci][4*sr+1]=Indices[cii];
                        label tnoc=0;
                        scalar ni=0;
                        while(ni<=fabs(nFSC[cii]) && noc<totcellTime[Indices[ci]].size())
                        {
                            ni+=totcellTime[Indices[ci]][noc];
                            nFSC[ci]-=totcellTime[Indices[ci]][noc];
                            noc=noc+1;
                            tnoc=tnoc+1;
                        }
                        sendRecepie[nFSC.size()-1-ci][4*sr+2]=tnoc;
                        sendRecepie[nFSC.size()-1-ci][4*sr+3]=sr+1;

                        if(cii+1<re && noc<totcellTime[Indices[ci]].size()) //if there is othere reciever(s)
                        {
                            //nFSC[ci]=res;
                            //totSend[Indices[ci]]+=tnoc;
                            //totRecv[Indices[cii]]+=tnoc;
                            cii++;
                            sr++;
                            sendRecepie[nFSC.size()-1-ci].setSize(4*(sr+1));
                        }

                        else if(cii+1>=re || noc>=totcellTime[Indices[ci]].size())//if there is no other reciever(s)
                        {
                            //sendRecepie[nFSC.size()-1-ci][4*sr+2]=nFSC[ci];
                            //nFSC[ci]=res;
                            //totSend[Indices[ci]]+=tnoc;
                            //totRecv[Indices[cii]]+=tnoc;
                            ci--;
                            break;
                        }
                        //}
                    }
                    else if(res<0)
                    {
                        /*if(fabs(res)<=0.1*MeanNFSC)
		    				{
		    					sendRecepie[nFSC.size()-1-ci][4*sr]=Indices[ci];
		    					sendRecepie[nFSC.size()-1-ci][4*sr+1]=Indices[cii];
		    					//sendRecepie[nFSC.size()-1-ci][4*sr+2]=nFSC[ci];
		    					label tnoc=0;
		    					label ni=0;
		    					while(ni<=fabs(nFSC[ci]))
		    					{
		    						ni+=totcellTime[Indices[ci]][noc];
		    						noc++;
		    						tnoc++;
		    					}
		    					sendRecepie[nFSC.size()-1-ci][4*sr+2]=tnoc;
		    					sendRecepie[nFSC.size()-1-ci][4*sr+3]=sr+1;
							//totSend[Indices[ci]]+=tnoc;
							//totRecv[Indices[cii]]+=tnoc;
		    					ci--;
		    					cii++;
		    					break;
		    				}
		    				else
		    				{*/
                        sendRecepie[nFSC.size()-1-ci][4*sr]=Indices[ci];
                        sendRecepie[nFSC.size()-1-ci][4*sr+1]=Indices[cii];
                        //sendRecepie[nFSC.size()-1-ci][4*sr+2]=nFSC[ci];
                        label tnoc=0;
                        scalar ni=0;
                        while(ni<=fabs(nFSC[ci]) && noc<totcellTime[Indices[ci]].size())
                        {

                            ni+=totcellTime[Indices[ci]][noc];
                            nFSC[cii]+=totcellTime[Indices[ci]][noc];
                            noc=noc+1;
                            tnoc=tnoc+1;

                        }
                        sendRecepie[nFSC.size()-1-ci][4*sr+2]=tnoc;
                        sendRecepie[nFSC.size()-1-ci][4*sr+3]=sr+1;
                        //nFSC[cii]=res;
                        //totSend[Indices[ci]]+=tnoc;
                        //totRecv[Indices[cii]]+=tnoc;
                        ci--;
                        break;
                        //}
                    }
                }
            }// END while
            //commClock_.timeIncrement();
            
            //How many number of queries each processor will recieve?
            //this is broadcasted by master node  
            for(int slave=Pstream::firstSlave();slave<=Pstream::lastSlave();slave++)
            {
                //OPstream toSlave(Pstream::scheduled, slave);
                //use this for OFDev
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                //Pstream::commsTypes::nonBlocking
                toSlave << sendRecepieSize;
            }
            //communicateTime_+=commClock_.timeIncrement();
        }
        else
        {
        	//slaves see how many queries they will recieve...
        
            //commClock_.timeIncrement();
           // temp.timeIncrement();
            
            //IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
            //use this for OFDev
            IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
            fromMaster >> sendRecepieSize;
            //communicateTime_+=commClock_.timeIncrement();
        }
        if (Pstream::myProcNo() == Pstream::masterNo())
        {
        	//master sends the list of queries for each proc
        	
        	
            //commClock_.timeIncrement();
            //temp.timeIncrement();
            for(int slave=Pstream::firstSlave();slave<=Pstream::lastSlave();slave++)
            {
            	//OPstream toSlave(Pstream::scheduled, slave);
            	//use this for OFDev
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << sendRecepie;
            }
           // communicateTime_+=commClock_.timeIncrement();
        }
        else
        {
        	//slaves recieve the list of queries they shoudl handle 
        	
            //commClock_.timeIncrement();
            //temp.timeIncrement();
            sendRecepie.setSize(sendRecepieSize);
            
            //IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
            //use this for OFDev
            IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
            fromMaster >> sendRecepie;
           // communicateTime_+=commClock_.timeIncrement();
        }
        for(label i4=0;i4<sendRecepie.size();i4++)
        {
            for(label j4=0;j4<sendRecepie[i4][sendRecepie[i4].size()-1];j4++)
            {
                totSend[sendRecepie[i4][4*j4]]+=sendRecepie[i4][4*j4+2];
                totRecv[sendRecepie[i4][4*j4+1]]+=sendRecepie[i4][4*j4+2];
            }
        }
        //START COMMUNICATION BETWEEN PROCESSORS TO EVALUATE THOSE CELLS THAT HAVE NORMALIZED GRADIENT TEMPRETURE BIGGER THAN 0.3
        //if the normalized gradient in cell was bigger  than 0.3 but there was no need to commiunication, direct integration set place
        //i.e. DI set to one
        cTOTRecv.setSize(totRecv[Pstream::myProcNo()],List<scalar>(this->nSpecie()+6,0.0));
        RRRecv.setSize(totRecv[Pstream::myProcNo()],List<scalar>(this->nSpecie(),0.0));
        RQhiRecv.setSize(totRecv[Pstream::myProcNo()],List<scalar>(this->nEqns()+1,0.0));
        DeltatchemRecv.setSize(totRecv[Pstream::myProcNo()],0.0);
        cellTimeRecv.setSize(totRecv[Pstream::myProcNo()],0.0);
        Deltatchem.setSize(totSend[Pstream::myProcNo()],0.0);
        //commClock_.timeIncrement();
        for(label i=0;i<sendRecepieSize;i++)//i counts number of processors which sends their cells
        {
            label j=0;// j counts number of packs that each processor send to others
            sendSize=0;
            while(j<sendRecepie[i][sendRecepie[i].size()-1])
            {
                cTOTSendOrRecv.setSize(sendRecepie[i][4*j+2], List<scalar>(this->nSpecie()+6,0.0));
                if(Pstream::myProcNo()==sendRecepie[i][4*j])
                {
                    for(label jj=0;jj<sendRecepie[i][4*j+2];jj++) //jj -->cell Ids which should solve in other processor
                    {
                        //TotalCellsDI0++;
                        label jjj=jj+sendSize;
                        scalarField phiq(this->nEqns());
                        scalar rhoi = rho[cellIDsTot[Pstream::myProcNo()][jjj]];
                        for(label ii=0; ii<this->nSpecie(); ii++)
                        {
                            phiq[ii] = this->Y()[ii][cellIDsTot[Pstream::myProcNo()][jjj]];
                        }
                        scalarField c(this->nSpecie());
                        for(label ii=0; ii<this->nSpecie(); ii++)
                        {
                            c[ii] = rhoi*phiq[ii]/this->specieThermo_[ii].W();
                        }
                        for(label ii=0; ii<this->nSpecie(); ii++)
                        {
                            cTOTSendOrRecv[jj][ii] = c[ii];
                        }
                        cTOTSendOrRecv[jj][this->nSpecie()]=this->thermo().T()[cellIDsTot[Pstream::myProcNo()][jjj]];
                        cTOTSendOrRecv[jj][this->nSpecie()+1]=this->thermo().p()[cellIDsTot[Pstream::myProcNo()][jjj]];
                        //cTOTSendOrRecv[jj][this->nSpecie()+2]=this->thermo().hs()[cellIDsTot[Pstream::myProcNo()][jjj]] + hc[cellIDsTot[Pstream::myProcNo()][jjj]];
                        cTOTSendOrRecv[jj][this->nSpecie()+2]=deltaT[cellIDsTot[Pstream::myProcNo()][jjj]];
                        cTOTSendOrRecv[jj][this->nSpecie()+3]=rho[cellIDsTot[Pstream::myProcNo()][jjj]];
                        cTOTSendOrRecv[jj][this->nSpecie()+4]=this->deltaTChem_[cellIDsTot[Pstream::myProcNo()][jjj]];
                        cTOTSendOrRecv[jj][this->nSpecie()+5]=deltaTMin;
                    }
                    sendSize+=sendRecepie[i][4*j+2];
                    OPstream toReciever0(Pstream::commsTypes::scheduled,sendRecepie[i][4*j+1]);
                    toReciever0 << cTOTSendOrRecv;
                }
                else if(Pstream::myProcNo()==sendRecepie[i][4*j+1])
                {

                    IPstream fromSender0(Pstream::commsTypes::scheduled, sendRecepie[i][4*j]);
                    fromSender0 >> cTOTSendOrRecv;
                    //numRecv is set to 0 at first to the size is numRecv+1
					infoRecv.setSize(numRecv+1,List<label>(2));
                    deltaTminRecv.setSize(numRecv+1);
                    infoRecv[numRecv][0]=sendRecepie[i][4*j];
                    infoRecv[numRecv][1]=cTOTSendOrRecv.size();
                    for(label it=0;it<cTOTSendOrRecv.size();it++)
                    {
                        label it3=it+lastIt;
                        for(label it2=0;it2<this->nSpecie()+6;it2++)
                        {
                            cTOTRecv[it3][it2]=cTOTSendOrRecv[it][it2];
                        }
                    }
                    lastIt+=cTOTSendOrRecv.size();
                    numRecv++;
                }
                j++;
            }
        }
        //communicateTime_+=commClock_.timeIncrement();

        //clockTime temp2=clockTime();
        //temp2.timeIncrement();

        //Solution of nonresident cells
        #include "solveNonResPDist.H"

        //Solution of resident Cells  with Its Own ISAT Table
        #include "solveResPDist.H"
        //sendBack
        // commClock_.timeIncrement();
        //temp.timeIncrement();
        for(label i=0;i<sendRecepieSize;i++)
        {
            label j=0;
            label lastp1=0;
           // clockTime sendtime=clockTime();
           // sendtime.timeIncrement();
            for(label i1=0;i1<infoRecv.size();i1++)
            {
                if(infoRecv[i1][0]==sendRecepie[i][4*j])
                {
                    //Making infoPack
                    forAll(infoPack, fieldI)
                    {
                        infoPack[fieldI].setSize(infoRecv[i1][1]);
                    }

                    infoPack[2][0]=deltaTminRecv[i1];
                    for(label ii=lastp1;ii<lastp1+infoRecv[i1][1];ii++)
                    {
                        infoPack[0][ii-lastp1]=cellTimeRecv[ii];

                        infoPack[1][ii-lastp1]=DeltatchemRecv[ii];

                        for(label jj=0; jj<this->nSpecie(); jj++)//jj--> species number
                        {
                            infoPack[jj+3][ii-lastp1] = RRRecv[ii][jj];
                        }
                        for(label jj=0; jj<this->nEqns()+1; jj++)//jj--> species number
                        {
                            infoPack[jj+3+this->nSpecie()][ii-lastp1] = RQhiRecv[ii][jj];
                        }
                    }

                    OPstream toReciever2(Pstream::commsTypes::scheduled, sendRecepie[i][4*j]);
                    toReciever2 << infoPack;
                }
                lastp1+=infoRecv[i1][1];
            }

            label lastp=0;
            while(j<sendRecepie[i][sendRecepie[i].size()-1])
            {
                if(Pstream::myProcNo()==sendRecepie[i][4*j])
                {
                    clockTime recvtime=clockTime();
                    recvtime.timeIncrement();
                    forAll(infoPack, fieldI)
                    {
                        infoPack[fieldI].setSize(sendRecepie[i][4*j+2]);
                    }

                    IPstream fromSender2(Pstream::commsTypes::scheduled, sendRecepie[i][4*j+1]);
                    fromSender2 >> infoPack;

                    deltaTMin=min(deltaTMin,infoPack[2][0]);

                    for(label ii=lastp;ii<lastp+sendRecepie[i][4*j+2];ii++)
                    {
                        //this->cellTime_[cellIDsTot[Pstream::myProcNo()][ii]] = infoPack[0][ii-lastp];
                        this->cellTime_[cellIDsTot[Pstream::myProcNo()][ii]] += infoPack[0][ii-lastp];
                        this->deltaTChem_[cellIDsTot[Pstream::myProcNo()][ii]] = infoPack[1][ii-lastp];
                        for(label jj=0; jj<this->nSpecie(); jj++)//jj--> species number
                        {
                            this->RR()[jj][cellIDsTot[Pstream::myProcNo()][ii]] = infoPack[jj+3][ii-lastp];
                        }
                        //Add this leaf
                        if(infoPack[3+this->nSpecie()+this->nEqns()][ii-lastp]!=0)
                        {
                            scalar rhoi = rho[cellIDsTot[Pstream::myProcNo()][ii]];
                            // no bug
							scalar Ti = this->thermo().T()[cellIDsTot[Pstream::myProcNo()][ii]];
                            scalar pi = this->thermo().p()[cellIDsTot[Pstream::myProcNo()][ii]];
                            //scalar Ti = infoPack[3+this->nSpecie()+this->nEqns()-2][ii-lastp];
							//scalar pi = infoPack[3+this->nSpecie()+this->nEqns()-1][ii-lastp];
							scalarField phiq(this->nEqns() + nAdditionalEqn);
							scalarField Rphiq(this->nEqns() + nAdditionalEqn);
							Rphiq=Zero;
                            for(label in=0; in<this->nSpecie(); in++)
                            {
                                phiq[in] = this->Y()[in][cellIDsTot[Pstream::myProcNo()][ii]];
                            }
                            phiq[this->nSpecie()]=Ti;
                            phiq[this->nSpecie()+1]=pi;
                            if (tabulation_->variableTimeStep())
							{
							 	phiq[this->nSpecie() + 2] = deltaT[cellIDsTot[Pstream::myProcNo()][ii]];
							}
                            for(label in=0; in<this->nSpecie(); in++)
                            {
                                Rphiq[in]=infoPack[in+3+this->nSpecie()][ii-lastp];
                            }
                            if (tabulation_->variableTimeStep())
							{
								Rphiq[Rphiq.size()-3] = infoPack[3+this->nSpecie()+this->nEqns()-2][ii-lastp];
								Rphiq[Rphiq.size()-2] = infoPack[3+this->nSpecie()+this->nEqns()-1][ii-lastp];
								Rphiq[Rphiq.size()-1] = deltaT[cellIDsTot[Pstream::myProcNo()][ii]];
							}
							else
							{
								Rphiq[Rphiq.size()-2] = infoPack[3+this->nSpecie()+this->nEqns()-2][ii-lastp];
								Rphiq[Rphiq.size()-1] = infoPack[3+this->nSpecie()+this->nEqns()-1][ii-lastp];
							}
							label growOrAdd =
								tabulation_->add(phiq, Rphiq, rhoi, deltaT[cellIDsTot[Pstream::myProcNo()][ii]]);
							if (growOrAdd)
							{
								this->setTabulationResultsAdd(cellIDsTot[Pstream::myProcNo()][ii]);
								//addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
							}
							else
							{
								this->setTabulationResultsGrow(cellIDsTot[Pstream::myProcNo()][ii]);
								//growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
							}
                        }
                    }
                    lastp+=sendRecepie[i][4*j+2];
                }
                j++;
            }
        }
        //communicateTime_+=commClock_.timeIncrement();

        //Display information about ISAT
        //scalar foundRatio =  (static_cast<scalar> (nFound_))/meshSize;
        //Pout << "Tabulation found " << foundRatio*100 << "% of the cells in the binary tree" << endl;
        //Pout << "Chemistry library size = " ;
        //Pout << tabPtr_->size() << endl;

        //Pout<<"Number of DI for Resident Cells in each step is "<<nDIResidentCells<<endl;
        //Pout<<"Number of DI for NonResident Cells in each step is"<<nDINonResidentCells<<endl;
        //Pout<<"Number of DI Total in each step is "<<nDIResidentCells+nDINonResidentCells<<endl;
		
		//Pout<<"communicationTime="<<communicateTime_<<endl;
