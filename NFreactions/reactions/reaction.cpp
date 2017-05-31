#include "reaction.hh"
#include "../../NFutil/setting.hh" //razi added for debugging purpose, last update 2017-3-29

using namespace std;
using namespace NFcore;




FunctionalRxnClass::FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet, System *s) :
	BasicRxnClass(name,1,"",transformationSet,s)
{
	if ((RAZI_DEBUG) & (SHOW_FIRE | CREATE_REACTION)) {cout<<"\n\tFunctional RXN:"<<name<<" with global function:"<< gf->getNiceName() <<" is created.\n";  mypause(-1);}
	this->cf=0;
	this->gf=gf;
	for(int vr=0; vr<gf->getNumOfVarRefs(); vr++) {
		if(gf->getVarRefType(vr)=="Observable") {
			Observable *obs = s->getObservableByName(gf->getVarRefName(vr));
			obs->addDependentRxn(this);
		} else {
			cerr<<"When creating a FunctionalRxnClass of name: "+name+" you provided a function that\n";
			cerr<<"depends on an observable type that I can't yet handle! (which is "+gf->getVarRefType(vr)+"\n";
			cerr<<"try using type: 'MoleculeObservable' for now.\n";
			cerr<<"quiting..."<<endl; exit(1);
		}
	}
}

FunctionalRxnClass::FunctionalRxnClass(string name, CompositeFunction *cf, TransformationSet *transformationSet, System *s) :
	BasicRxnClass(name,1, "", transformationSet,s)
{
	if ((RAZI_DEBUG) & (SHOW_FIRE | CREATE_REACTION)) {cout<<"\n\tFunctional RXN:"<<name<<" with composite function:"<< cf->getName() <<" is created.\n";  mypause(-1);}
	this->gf=0;
	this->cf=cf;
	this->cf->setGlobalObservableDependency(this,s);
}


FunctionalRxnClass::~FunctionalRxnClass() {};

double FunctionalRxnClass::update_a() {
	//cout<<"udpating a"<<endl;
	if(this->onTheFlyObservables==false) {
		cerr<<"Warning!!  You have on the fly observables turned off, but you are using functional\n";
		cerr<<"reactions which depend on observables.  Therefore, you cannot turn off onTheFlyObservables!\n";
		cerr<<"exiting now."<<endl;
		exit(1);
	}



	//	cout<<"here"<<endl;
	if(gf!=0) {
	//	cout<<"in here"<<endl;
		a=FuncFactory::Eval(gf->p);
	} else if(cf!=0) {
		int * reactantCounts = new int[this->n_reactants];
		for(unsigned int r=0; r<n_reactants; r++) {
			reactantCounts[r] = (int)getReactantCount(r);
		}
		a=cf->evaluateOn(0,0, reactantCounts, n_reactants);
		delete [] reactantCounts;
	//	cout<<"and here"<<endl;
	} else {
		cout<<"Error!  Functional rxn is not properly initialized, but is being used!"<<endl;
		exit(1);
	}

	if(a<0) {
		cout<<"Warning!!  The function you provided for functional rxn: '"<<name<<"' evaluates\n";
		cout<<"to a value less than zero!  You cannot have a negative propensity!";
		cout<<"here is the offending function: \n";
		gf->printDetails();
		cout<<"\nhere is the offending reaction: \n";
		this->printDetails();
		cout<<"\n\nquitting."<<endl;
		exit(1);
	}


	// check here for the total rate flag - if this is set to true, then
	// use the rate exactly as given by the function, but if it is false,
	// then we have to multiply here by the reactant counts
	if(!this->totalRateFlag) {
		for(unsigned int i=0; i<n_reactants; i++)
			a*=(double)getCorrectedReactantCount(i);
	}
	else
	{
		// Check that we have at least one set of reactants!
		for(unsigned int i=0; i<n_reactants; i++) {
			if(getCorrectedReactantCount(i)==0) {
				a=0.0;
				break;
				//cout<<"Warning!  Function evaluates to positive rate for a reaction, but"<<endl;
				//cout<<"one of the reactant lists is empty!"<<endl;
				//this->printDetails();
				//cf->printDetails(reactantTemplates[0]->getMoleculeType()->getSystem());
				//exit(1);
			}
		}
	}
	
	return a;
}

void FunctionalRxnClass::printDetails() const {

	string trate = "off";
	if(this->totalRateFlag) trate = "on";

	if(gf!=0)
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<gf->getNiceName()<<"="<<FuncFactory::Eval(gf->p)<<",  a="<<a<<", fired="<<fireCounter<<" times, TotalRate="<<trate<<" )"<<endl;
	else if(cf!=0) {
		int * reactantCounts = new int[this->n_reactants];
		for(unsigned int r=0; r<n_reactants; r++) {
			reactantCounts[r]=getReactantCount(r);
		}
		double value=cf->evaluateOn(0,0, reactantCounts, n_reactants);
		delete [] reactantCounts;
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<cf->getName()<<"="<<value<<",  a="<<a<<", fired="<<fireCounter<<" times, TotalRate="<<trate<<" )"<<endl;

	}

	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}









MMRxnClass::MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet,System *s) :
	BasicRxnClass(name,1,"",transformationSet,s)
{
	if ((RAZI_DEBUG) & (SHOW_FIRE | CREATE_REACTION)) {cout<<"\n\tMM-RXN:"<<name<<" is created.\n";  mypause(-1);}

	this->Km = Km;
	this->kcat = kcat;
	this->sFree=0;
	if(n_reactants!=2) {
		cerr<<"You have tried to create a reaction with a Michaelis-Menten rate law (named: '"+name+"'\n')";
		cerr<<"but you don't have the correct number of reactants!  Michaelis-Menten reactions require\n";
		cerr<<"exactly 2 reactants.  A substrate (always given first) and an enzyme (always given second)\n";
		cerr<<"Read your tutorial next time... now I will quit."<<endl;
		exit(1);
	}
}
MMRxnClass::~MMRxnClass() {};

double MMRxnClass::update_a()
{
	double S = (double)getCorrectedReactantCount(0);
	double E = (double)getCorrectedReactantCount(1);
	sFree=0.5*( (S-Km-E) + pow((pow( (S-Km-E),2.0) + 4.0*Km*S),  0.5) );
	a=kcat*sFree*E/(Km+sFree);
	return a;
}

void MMRxnClass::printDetails() const {
	cout<<"ReactionClass: " << name <<"  ( Km="<<Km<<", kcat="<<kcat<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}









BasicRxnClass::BasicRxnClass(string name, double baseRate, string baseRateName, TransformationSet *transformationSet, System *s) :
	ReactionClass(name,baseRate,baseRateName,transformationSet,s)
{
	if ((RAZI_DEBUG) & (CREATE_REACTION)) {cout<<"\n\tBasic RXN:"<<name<<" with base rate:"<< baseRate <<" is created.\n";  mypause(-1);}
	this->reactionType = BASIC_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	reactantLists = new ReactantList *[n_reactants];
	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists[r]=(new ReactantList(r,transformationSet,25));
}


BasicRxnClass::~BasicRxnClass()
{
    //cout<<"  -------------------------------\n  ----------------------------\n";
	//cout<<"Reaction: "<<name<<endl;
	//this->reactantLists[0]->printDetails();

	//this->reactantLists[0]->removeMappingSet(0);
	//this->reactantLists[0]->removeMappingSet(3);
	//this->reactantLists[0]->removeMappingSet(6);
	//this->reactantLists[0]->removeMappingSet(9);

	//cout<<endl<<endl<<endl;
	//this->reactantLists[0]->printDetails();




	if(DEBUG) cout<<"Destroying rxn: "<<name<<endl;

	for(unsigned int r=0; r<n_reactants; r++)
	{
		//delete reactantTemplates[r]; DO NOT DELETE HERE (MoleculeType has responsibility of
		//deleting all template molecules of its type now.
		delete reactantLists[r];
	}
	delete [] reactantLists;
}

void BasicRxnClass::init()
{
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}


void BasicRxnClass::prepareForSimulation()
{

}


int BasicRxnClass::checkForEquality(Molecule *m, MappingSet* ms, int rxnIndex, ReactantList* reactantList){
	/*
	Check if mapping set clashes with any of the mapping sets already in reactantList
	*/
	set<int> tempSet = m->getRxnListMappingSet(rxnIndex);
	for(set<int>::iterator it= tempSet.begin();it!= tempSet.end(); ++it){
		MappingSet* ms2 = reactantList->getMappingSet(*it);
		if(MappingSet::checkForEquality(ms,ms2)){
			return *it;
		}
	}
	return -1;


}

bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos)
{
	//First a bit of error checking, that you should skip unless we are debugging...
	//	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	//	{
	//		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
	//		exit(1);
	//	}

	//Get the specified reactantList
	rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	//cout<<" got mappingSetId: " << m->getRxnListMappingId(rxnIndex)<<" size: " <<rl->size()<<endl;
	//cout<< " testing whether to add molecule ";
	//m->printDetails();
	//cout<<" ... as a normal reaction "<<this->name<<endl;


	//If this reaction has multiple instances, we always remove them all!
	// then we remap because other mappings may have changed.  Yes, this may
	// be more ineffecient, but it is the fast implementation
	if(rl->getHasClonedMappings()) {
		while(m->getRxnListMappingId(rxnIndex)>=0) {
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->deleteRxnListMappingId(rxnIndex,m->getRxnListMappingId(rxnIndex));
			//m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	}

	/*
	//Here we get the standard update...
	if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
	{
		if(!reactantTemplates[reactantPos]->compare(m)) {
			//	cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}

	} else {
		//Try to map it!
		ms = rl->pushNextAvailableMappingSet();
		if(!reactantTemplates[reactantPos]->compare(m,rl,ms)) {
			//we must remove, if we did not match.  This will also remove
			//everything that was cloned off of the mapping set
			rl->removeMappingSet(ms->getId());
		} else {
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
	}
	*/

	//Here we get the standard update...
	while(m->getRxnListMappingId(rxnIndex)>=0) {
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->deleteRxnListMappingId(rxnIndex,m->getRxnListMappingId(rxnIndex));
		//m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}


	//Try to map it!
	ms = rl->pushNextAvailableMappingSet();
	symmetricMappingSet.clear();
	comparisonResult = reactantTemplates[reactantPos]->compare(m,rl,ms,false,&symmetricMappingSet);
	if(!comparisonResult) {
		//cout << "no mapping in normal reaction, remove"<<endl;
		//we must remove, if we did not match.  This will also remove
		//everything that was cloned off of the mapping set
		rl->removeMappingSet(ms->getId());
		//JJT: removes any symmetric mapping sets that might have been added since we are not using them
		for(vector<MappingSet *>::iterator it=symmetricMappingSet.begin();it!=symmetricMappingSet.end();++it){
			rl->removeMappingSet((*it)->getId());
		}
	} else {
		//cout << "should be in normal reaction, confirm push"<<endl;
		//ms->printDetails();
		
		//TODO: it is necessary to remove elements that are not used anymore from the rl as well as from the m
		//for that
		//m->setRxnListMappingId(rxnIndex,-1);

		if (symmetricMappingSet.size() > 0){
            rl->removeMappingSet(ms->getId());
			for(vector<MappingSet *>::iterator it=symmetricMappingSet.begin();it!=symmetricMappingSet.end();++it){
					//XXX: JJT this is a band-aid, symmetricMappingSet should not have repeated elements in the first place
					int mapIndex = checkForEquality(m,*it,rxnIndex,rl);
					if(mapIndex >= 0){
						rl->removeMappingSet((*it)->getId());
					}
					else{
						m->setRxnListMappingId(rxnIndex,(*it)->getId());
					}
            }
		}
		else{
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
		
	}

	return true;
}



void BasicRxnClass::remove(Molecule *m, unsigned int reactantPos)
{

	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}


	//Get the specified reactantList
	ReactantList *rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);


	if(isInRxn)
	{
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}
}




void BasicRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex)
{
	cerr<<"You are trying to notify a Basic Reaction of a rate Factor Change!!! You should only use this"<<endl;
	cerr<<"function for DORrxnClass rules!  For this offense, I must abort now."<<endl;
	exit(1);
}

double BasicRxnClass::update_a()
{
	// Use the total rate law convention (macroscopic rate)


#ifdef FIX_A   //Razi: Propensity of a reaction should be perhaps the minimum of reactant counts and not multiplication of it, later check
	double c;
	a = 1.0;
	a = getCorrectedReactantCount(0);
	if (n_reactants >= 1){
		for(unsigned int i=1; i<n_reactants; i++) {
			c=getCorrectedReactantCount(i);	a=(a<c)?a:c;    //min(a,getCorrectedReactantCount(i)
			mypause(0);
			cout<<"\nupdate a for reaction:"<< this->getName()<<"  i:"<<i<<" count:"<<c<< "  a:min(a,c):"<<a<<endl;
			mypause(0);
		}
	}
	a*=baseRate;
	cout<<"\tupdating finished: a for reaction after applying baserate is:"<< a<<endl;

#else

	if(this->totalRateFlag) {
		a=baseRate;
		for(unsigned int i=0; i<n_reactants; i++)
			if(getCorrectedReactantCount(i)==0) a=0.0;

	// Use the standard microscopic rate
	} else {
		a = 1.0;
		for(unsigned int i=0; i<n_reactants; i++) {
			a*=getCorrectedReactantCount(i);
		}
		a*=baseRate;
	}

#endif


	if ((RAZI_DEBUG & (SHOW_FIRE | CREATE_REACTION) && (this->getName().compare("XbPlusYbBind")==0))) {// && (this->getName().compare("XbPlusYbBind")==0) && 0){ //don't show anymore
		int i,j,k;
		cout <<"\tupdate_a() is called for Basic reaction:"<<this->name << " totalRateFlag:"<<totalRateFlag<<"   Base rate:"<<baseRate<<endl;
		for(i=0; i<n_reactants; i++){
			j = getCorrectedReactantCount(i);
			cout<<"\tReactantlist "<<i<<"  Population Type: "<< isPopulationType[i] <<" Corrected count: "<< j<<endl;
		}
		cout <<"\tFinal a is :"<<a<<endl<<endl;
		mypause(6000);
	}

	return a;
}


int BasicRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	return isPopulationType[reactantIndex] ?
			   reactantLists[reactantIndex]->getPopulation()
			 : reactantLists[reactantIndex]->size();
}


int BasicRxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	/*cerr << "  getCorrectedReactantCount rindex: " << reactantIndex << "  isPop? " << isPopulationType[reactantIndex] << endl;
	if ( isPopulationType[reactantIndex] )
	{
		cerr << "  corr:  " << identicalPopCountCorrection[reactantIndex];
		cerr << "  pop:   " << reactantLists[reactantIndex]->getPopulation() << endl;
		cerr << "  final: " << std::max( reactantLists[reactantIndex]->getPopulation()
	             - identicalPopCountCorrection[reactantIndex], 0 ) << endl;
	}
	else
	{
		cerr << "  count: " << reactantLists[reactantIndex]->size() << endl;
	}
	*/
	return isPopulationType[reactantIndex] ?
			   std::max( reactantLists[reactantIndex]->getPopulation()
			             - identicalPopCountCorrection[reactantIndex], 0 )
			 : reactantLists[reactantIndex]->size();
}

void BasicRxnClass::printFullDetails() const
{
	if ((RAZI_DEBUG) & (SHOW_FIRE | CREATE_REACTION)){
		unsigned int i=0;
//#ifdef RHS_FUNC
//			cout<<"BasicRxnClass: "<<name<<"  has "<< n_reactants<<" reactants and "<< n_mappingsets<< " mappingSets and"<< n_productTemplates <<" Output products." << endl;  //razi added
//#else
			cout<<"BasicRxnClass: "<<name<<"  has "<< n_reactants<<" reactants and "<< n_mappingsets<< " mappingSets." << endl;  //razi added
//#endif
			int n_printItems= (n_mappingsets<=3? n_mappingsets: 3); //Razi changed to show only a few
			for (i =0 ; i < n_printItems; i++){
				if (mappingSet[i]){
					try{
						cout<<" Check why there is an error in the following line for i:"<<i<<endl;
						//cout<<"A-mappsingSet Id: "<< mappingSet[i]->getId() << " # of mappings:" << mappingSet[i]->getNumOfMappings() << " complex Id:" << mappingSet[i]->getComplexID()<<endl;
					}catch(std::string err){ //int e string err
						cout<<"BasicRxnClass::printFullDetails: Error occurred for i:"<<i<<"  Err:"<< err<< endl; mypause(5000);//" Err:"<<err<<" !!!"<<endl;
					}
				}else
					cout<<"mappsingSet ["<<i<<"] for Reaction Class is not set yet." <<endl;
			}
			cout<<"------------------------------------------------------------------\n";
	}
	else cout<<"BasicRxnClass: "<<name<<endl; //razi commented
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists[i]->printDetails();
}


void BasicRxnClass::pickMappingSets(double random_A_number) const
{
	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire

	//Select a reactant from each list
	for(unsigned int i=0; i<n_reactants; i++)
	{
		if ( isPopulationType[i] ) {
			reactantLists[i]->pickRandomFromPopulation(mappingSet[i]);
		} else {
			reactantLists[i]->pickRandom(mappingSet[i]);
		}
	}
}


