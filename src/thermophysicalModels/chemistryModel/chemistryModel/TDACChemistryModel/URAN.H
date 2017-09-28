//notRet is number of cells that are not retrieved at QT stage 
//and need to be distributed
label notRet=0;
scalar deltaTminSendOrRecv=0;	

//needed somewhere I am not sure
label o=0;
label oo=0;

//nx is the dimension of vector which contains information of each query
label nx=2*this->nSpecie()+9+1+1;
//xf is the matrix of information
List<List<scalar> > xf(nx, List<scalar>(notRet,0.0));


//commented on bahman 93

//first let's do quick try!
//forAll(rho, celli)
forAll(rho, ci)
{


	clockTime_cellTime.timeIncrement();
	#include "QTURAN.H"
	cellTime_[celli]=clockTime_cellTime.timeIncrement();


}
//to test
Pout<<"notRet = "<<notRet<<endl;

//again many variables declarations for URAN method!
label notRet_un=notRet;
label nf=nx;//should be corrected
label ldxf=nx;//should be corrected
label k_pos=0;//should be corrected
label myrank=Pstream::myProcNo();
label nproc=Pstream::nProcs();
label ntogo=0; //should be corrected
label msgsize=0;
label nnr=0;
label nlp=0;
label nout=0;
label start_marker=0;
label start_marker1=0;
label count_nr=0;
label count_d=0;
label startat=0;
label stopat=0;
label index=0;

List<label> sourceof(notRet_un,0);
List<label> sourceof1(notRet_un,0);
List<label> p1(notRet_un,0);
List<label> n_nonres (nproc,0);
List<label> n_nonres_recv (nproc,0);
List<label> n_distrib (nproc,0);
List<label> n_distrib_recv (nproc,0);
List<label> start_nr_buf (nproc,0);
List<label> start_nr_buf1 (nproc,0);
List<label> start_d_buf  (nproc,0);
List<label> start_d_buf1  (nproc,0);
MPI_Request req_d[nproc];//List<MPI_Request> req_d (nproc,0);
MPI_Request req_nr1[nproc];//List<MPI_Request> req_nr1 (nproc,0);
MPI_Request req_nr2[nproc];//List<MPI_Request> req_nr2 (nproc,0);

double Rxf[notRet_un][nx];

for(label i = 0;i< notRet_un;i++)
{
	sourceof[i]  = i;
	sourceof1[i] = i;
}
//Directly from line 284 of x2f_mpi of Pope
forAll(p1,proci) p1[proci]=-1;
ntogo      = notRet_un;
x2f_uran_p( nproc,  myrank, ntogo, notRet_un, p1, n_nonres, n_distrib );

forAll(n_nonres_recv,proci) n_nonres_recv[proci]=n_nonres[proci];
n_nonres_recv[myrank]=0;
forAll(n_distrib_recv,proci) n_distrib_recv[proci]=n_distrib[proci];
n_distrib_recv[myrank]=0;


nnr = sum( n_nonres ) - n_nonres[myrank];

double Rxf_nr_buf[nnr][ldxf];

for(label i9=0;i9<nnr;i9++)
{
	for(label j9=0;j9<ldxf;j9++)
	{
		Rxf_nr_buf[i9][j9]=0;
	}
}

nout = sum( n_distrib ) - n_distrib[myrank];

nlp  = notRet_un - nout;

forAll(start_nr_buf,proci) start_nr_buf[proci]=0;
forAll(start_nr_buf1,proci) start_nr_buf1[proci]=0;
forAll(start_d_buf,proci) start_d_buf[proci]=0;
forAll(start_d_buf1,proci) start_d_buf1[proci]=0;

for(label i8=0;i8<nproc;i8++) req_d[i8]=MPI_REQUEST_NULL; 
for(label i8=0;i8<nproc;i8++) req_nr1[i8]=MPI_REQUEST_NULL; 
for(label i8=0;i8<nproc;i8++) req_nr2[i8]=MPI_REQUEST_NULL; 

// start communication: post nonblocking receives for all expected messages
start_marker  = 0;//check
start_marker1 = 0;
count_nr      = 0;
for(label j=0;j<nproc;j++)
{
	if ( n_nonres_recv[j]!=0 )// correct this
	{
		MPI_Irecv( &Rxf_nr_buf[start_marker][0], ldxf*n_nonres[j],	MPI_DOUBLE, j, 201+j, MPI_COMM_WORLD, &req_nr1[j]);

		start_nr_buf[j] = start_marker;
		start_marker = start_marker + n_nonres[j];
		count_nr = count_nr + 1;
	}
}

start_marker  = 0;
start_marker1 = 0;
count_d       = 0;
for(label j=0;j<nproc;j++)
{
	if ( n_distrib[j]!=0 )
	{
		start_d_buf[j] = start_marker;
		start_marker = start_marker + n_distrib[j];
		if ( j != myrank ) count_d = count_d + 1;
	}
}

// rearrange data to make an output buffer in xf, using a bucket sort algorithm
// the number of buckets equals nproc plus 1 for those already done by quick try
// sourceof is an array to link sorted locations in xf back to source locations
// for nout < i <= notRet, sourceof(i) identifies a source for local processing
// note, each vector to be fed as input to the x2f function will be of length nx

bucket_sort( notRet_un, ldxf, nx, xf, p1, nproc, n_distrib, sourceof, sourceof1 ) ;

for(label ii=0;ii<nx;ii++)
{
	for(label jj=0;jj<notRet_un;jj++)
	{
		Rxf[jj][ii]=xf[ii][jj];
	}
}

// buffer is set: post nonblocking sends AND receives for all outgoing messages
for(label j=0;j<nproc;j++)
{
	if (  n_distrib[j]!=0  &&  j != myrank  )	
	{
		MPI_Isend( &Rxf[start_d_buf[j]][0], ldxf*n_distrib[j],MPI_DOUBLE, j, 201+myrank, MPI_COMM_WORLD, &req_d[j] );

		MPI_Request_free( &req_d[j] );
		MPI_Irecv( &Rxf[start_d_buf[j]][0], ldxf*n_distrib[j], MPI_DOUBLE, j, 301+j, MPI_COMM_WORLD, &req_d[j] );
	}
}
Pout<<"@@@@ solve RESIDENT cells @@@@"<<endl;
//=========== phase 3: local processing of resident vectors ========================
// identify chunk of sorted xf array holding vectors for local processing
startat = start_d_buf[myrank];
stopat  = start_d_buf[myrank] + n_distrib[myrank];


//let's first solve for resident queries
if( startat!= stopat)
{
	#include "solveResURAN.H"
}

//=========== phase 4: processing of non-resident vectors ==========================
Pout<<"@@@@ solve NON resident cells @@@@"<<endl;
index = 0;//check if index output starts from 1 or 0
for(label ict=0;ict<count_nr;ict++)
{
	MPI_Waitany( nproc, req_nr1, &index, MPI_STATUS_IGNORE );

	startat = start_nr_buf[index];
	stopat  = start_nr_buf[index] + n_nonres[index];
	#include "solveNonResURAN2.H"
	MPI_Isend( &Rxf_nr_buf[startat][0], ldxf*n_nonres[index],  MPI_DOUBLE, index, 301+myrank, MPI_COMM_WORLD, &req_nr2[index] );


}
//=========== phase 5: receive f vectors from other processors =====================
Pout<<"@@@@ recive back Cells @@@@"<<endl;
for(label ict=0;ict<count_d;ict++)
{
	MPI_Waitany( nproc, req_d, &index, MPI_STATUS_IGNORE );

	startat = start_d_buf[index];
	stopat  = start_d_buf[index] + n_distrib[index];
	//Pout<<"Phase 5 index = "<<index<<endl;
	//Pout<<"Phase 5 :startat="<<startat<<"	stopat = "<<stopat<<endl;

	for(label i5=startat;i5<stopat;i5++)
	{

		this->cellTime_[(Rxf[i5][0])] += Rxf[i5][1];				
		this->deltaTChem_[(Rxf[i5][0])]=Rxf[i5][2];	
		deltaTMin=min(deltaTMin,Rxf[i5][3]);//check
		scalar Ti =Rxf[i5][4];
		scalar pi =Rxf[i5][5];
		//scalar hi =Rxf[i5][6];
		scalar rhoi = Rxf[i5][7];
		for(label jj=0; jj<this->nSpecie(); jj++)
		{
			this->RR()[jj][(Rxf[i5][0])]=Rxf[i5][jj+8];
		}
		if(Rxf[i5][8+2*this->nSpecie()]==1)//check the flag then add the chempoint
		{
			//clockTime_.timeIncrement();
			//compute A and Add this leaf
			
			scalarField phiq(this->nEqns() + nAdditionalEqn);
			scalarField Rphiq(this->nEqns() + nAdditionalEqn);
			Rphiq=Zero;

			for(label in=0; in<this->nSpecie(); in++)
			{
				phiq[in] = this->Y()[in][(Rxf[i5][0])];
			}
			phiq[this->nSpecie()]=Ti;
			phiq[this->nSpecie()+1]=pi;
			if (tabulation_->variableTimeStep())
			{
			 	phiq[this->nSpecie() + 2] = Rxf[i5][6];
			}
			for(label in=0; in<this->nSpecie(); in++)
			{
				Rphiq[in]=Rxf[i5][in+8+this->nSpecie()];
			}
			if (tabulation_->variableTimeStep())
            {
                Rphiq[Rphiq.size()-3] = Rxf[i5][8+1+2*this->nSpecie()];
                Rphiq[Rphiq.size()-2] = Rxf[i5][8+2+2*this->nSpecie()];
                Rphiq[Rphiq.size()-1] =  Rxf[i5][6];
            }
            else
            {
                Rphiq[Rphiq.size()-2] = Rxf[i5][8+1+2*this->nSpecie()];
                Rphiq[Rphiq.size()-1] =Rxf[i5][8+2+2*this->nSpecie()];
            }
		    label growOrAdd =
                tabulation_->add(phiq, Rphiq, rhoi, Rxf[i5][6]);
            if (growOrAdd)
            {
                this->setTabulationResultsAdd(Rxf[i5][0]);
            	//addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
            }
            else
            {
                this->setTabulationResultsGrow(Rxf[i5][0]);
                //growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
            }
		}

	}

}
//   ntogo = 0
//for(label i=0;i<notRet_un;i++)
//{
//	if ( xf[k_pos][i]>0 ) 
//	{
//		p1[i]=0;
//	}
//	else 
//	{
//	    xf[k_pos][i] = -xf[k_pos][i]
//       p1[i] = 1
//        ntogo  = ntogo + 1
//	}
//}
//bucket_total[0] = ntogo
//call bucket_sort( notRet_un, ldxf, ldxf, xf, p1, 1, bucket_total, sourceof, sourceof1 )
//notRet_un = ntogo + notRet_un - sum(n_distrib);

// mop up the nonblocking calls that sent back computed, nonresident f-vectors

for(label ict=0;ict<count_nr;ict++)
{
	MPI_Waitany( nproc, req_nr2, &index, MPI_STATUS_IGNORE );
}
//call pigeonhole_sort( notRet, ldxf, nf, xf, sourceof, sortedto)
//is it needed???
//END URAN
