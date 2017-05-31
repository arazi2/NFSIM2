#include "reaction.hh"
#include "../../NFutil/setting.hh" //razi added for debugging purpose, last update 2017-3-29

#define DEBUG_MESSAGE 0//was 0  Razi changed to show more debugging messages


using namespace std;
using namespace NFcore;



//should also accept list of local functions and list of PointerNames for each of the functions...
RHSRxnClass::RHSRxnClass(   //Razi: this version supports RHS functions
		string name,
		double baseRate,
		string baseRateName,
		TransformationSet *transformationSet,
		CompositeFunction *function,
		vector <string> &lfArgumentPointerNameList, System *s) :
	ReactionClass(name,baseRate,baseRateName,transformationSet,s)
{
	bool verbose=false;
	if (RAZI_DEBUG & CREATE_REACTION)
		verbose = system->getverbose();

	if (verbose) {cout<<"\n\tRHS RXN:"<<name<<" with RHS composite functions:"<< function->getName() << "and "<< n_reactants <<" reactants is created.\n"; mypause(-1);}

	this->reactionType = RHS_RXN;
	reactantLists = new ReactantList *[n_reactants];
	check_mappingSet = new MappingSet *[n_reactants];

	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists[r]=(new ReactantList(r,transformationSet,25));
	

	if (verbose) cout<<"RHSreaction: I developed everything up to here. More may be needed."; //mypause(-1); //exit(0);

	//Set the actual function
	this->cfo = function;

	//Initialize a to zero
	this->a=0; //later check

}




RHSRxnClass::~RHSRxnClass() {

	for(unsigned int r=0; r<n_reactants; r++) {
		delete reactantLists[r];
	}

	delete [] reactantLists;

/*
	delete [] argIndexIntoMappingSet;
	delete [] argMappedMolecule;
	delete [] argScope;
*/

}

void RHSRxnClass::init() {

	//Here we have to tell the molecules that they are part of this function
	//and for single molecule functions, we have to tell them also that they are in
	//this function, so they need to update thier value should they be transformed
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}




void RHSRxnClass::remove(Molecule *m, unsigned int reactantPos)
{

	cerr<<"irrelevant function remove for RHSRxnClass\n";
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





int RHSRxnClass::checkForEquality(Molecule *m, MappingSet* ms, int rxnIndex, ReactantList* reactantList){
	//Razi: Comes from BASIC RXN, check if any modifications are required
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




bool RHSRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos) {
	//based on Basic RXN
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





int RHSRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	//Razi: Check if this is needed
	return isPopulationType[reactantIndex] ?
		       reactantLists[reactantIndex]->getPopulation()
	         : reactantLists[reactantIndex]->size();
}


int RHSRxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	//Razi: Check if this is needed
	return isPopulationType[reactantIndex] ?
			   std::max( reactantLists[reactantIndex]->getPopulation()
			             - identicalPopCountCorrection[reactantIndex], 0 )
			 : reactantLists[reactantIndex]->size();
}


//This function takes a given mappingset and looks up the value of its local
//functions based on the local functions that were defined
double RHSRxnClass::evaluateLocalFunctions(MappingSet *ms)
{
	cerr<< "RHSRxnClass::evaluateLocalFunctions not developed yet.\n";
	return -1;
	//Razi: This may need substantial changes: Since we will need to apply the function on the reaction products
	// after firing a reaction, [and not on the input reactants]


	//Go through each function, and set the value of the function
	//this->argMappedMolecule
	//cout<<"\t\t\t\tRHSRxnClass::evaluateLocalFunctions()"<<endl;
	//cout<<"dor is reevaluating its function."<<endl;

	//Grab the molecules needed for the local function to evaluate
	for(int i=0; i<this->n_argMolecules; i++) {
		//cout<<"here."<<endl;
		//cout<<"\t\t\t\t\t"<<i<<": argMappedMolecule="<<argMappedMolecule[i]<<" argIndexIntoMappingSet="<<argIndexIntoMappingSet[i]<<endl;
		this->argMappedMolecule[i] = ms->get(this->argIndexIntoMappingSet[i])->getMolecule();
		//cout<<"\t\t\t\t\t"<<"argMappedMoleculeType="<<argMappedMolecule[i]->getMoleculeTypeName()<<endl;
		//cout<<"\t\t\t\t\t"<<"argMappedMoleculeScope="<<argScope[i]<<endl;
	}

	//cout<<"done setting molecules, so know calling the composite function evaluate method."<<endl;
	int * reactantCounts = new int[this->n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		reactantCounts[r]=reactantLists[r]->size();
	}

	double value = this->cfo->evaluateOn(argMappedMolecule,argScope, reactantCounts, n_reactants);
	delete [] reactantCounts;
	//cout<<"\t\t\t\t\t"<<"composite function value="<<value<<endl;

	return value;

	/*Molecule

	for(int i=0; i<(signed)lfList.size(); i++) {
		Molecule *molObject = ms->get(this->indexIntoMappingSet.at(i))->getMolecule();
		int index = lfList.at(i)->getIndexOfTypeIFunctionValue(molObject);
		this->localFunctionValue.at(i)=molObject->getLocalFunctionValue(index);
		//cout<<"found that local function: "<<getName()<<" evaluates to: " <<localFunctionValue.at(i)<<endl;
	}
	return this->localFunctionValue.at(0);
	*/
}


double RHSRxnClass::update_a() {
//Razi :Source is BasicRxnClass, check for modifications
bool verbose = false;
if (RAZI_DEBUG & SHOW_SIM)
	verbose = system->getverbose();

#ifdef FIX_A   //Razi: Propensity of a reaction should be perhaps the minimum of reactant counts and not multiplication of it, later check
	double c;
	a = 1.0;
	a = getCorrectedReactantCount(0);
	if (n_reactants >= 1){
		for(unsigned int i=1; i<n_reactants; i++) {
			c=getCorrectedReactantCount(i);	a=(a<c)?a:c;    //min(a,getCorrectedReactantCount(i)
			mypause(0);
			cout<<"\t update a for RHS reaction:"<< this->getName()<<"  i:"<<i<<" count:"<<c<< "  a:min(a,c):"<<a<<endl;
			mypause(0);
		}
	}
	a*=baseRate;
	cout<<"\t updating finished: a for reaction after applying base rate is:"<< a<<endl;

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


	if (verbose){// && (this->getName().compare("XbPlusYbBind")==0) && 0){ //don't show anymore
		int i,j,k;
		cout <<"\tupdate_a() is called for RHS reaction:"<<this->name << " totalRateFlag:"<<totalRateFlag<<"   Base rate:"<<baseRate<<endl;
		for(i=0; i<n_reactants; i++){
			j = getCorrectedReactantCount(i);
			cout<<"\tReactantlist "<<i<<"  Population Type: "<< isPopulationType[i] <<" Corrected count: "<< j<<endl;
		}
		cout <<"\tFinal a is :"<<a<<endl<<endl;		mypause(500);
	}
	return a;
}

void RHSRxnClass::pickMappingSets(double randNumber) const
{
	//Razi: Source is BasicRxnClass::pickMappingSets check for modifications

	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire
	for(unsigned int i=0; i<n_reactants; i++)
	{
		if ( isPopulationType[i] ) {
			reactantLists[i]->pickRandomFromPopulation(mappingSet[i]);
		} else {
			reactantLists[i]->pickRandom(mappingSet[i]);
		}
	}
}

void RHSRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) {
	cerr<<"notifyRateFactorChange is not developed or is an invalid function for RHS reaction!!!\n."; exit(1);
}


bool ReactionClass::checkReaction()   //clone reactants
{
	cout<<"ReactionClass:: checkReaction() not developed yet.\n"; exit(0);
	return true;
}





// old version, change list to vector


bool RHSRxnClass::checkReaction()   //razi: clone reactants //static const
{
	list <Molecule *>  origMs(100); //at most 100 connected molecules for each reactant
	list <Molecule *>  copyMs(100);
	vector<list <Molecule *> >  AllcopyMs(10, list <Molecule *>(100)); //at most 10 reactants per reaction rule
	vector<list <Molecule *> >  AllorigMs(10, list <Molecule *>(100));


	//return true;  //Razi: bypass check only for test, later uncomment
	bool verbose = false;
	if (RAZI_DEBUG & (SHOW_FIRE | SHOW_RHS))
		verbose = this->system->getverbose();
	//verbose = false;


	if (verbose) cout<<"RHSRxnClass:: checkReaction() is called.\n";

	Molecule * copyM;
	int maxDepth=10;
	double result=1;

//	try{

	//clone molecules and complexes
	//int start_id = Molecule::uniqueIdCount;
	int start_id = mappingSet[0]->get(0)->getMolecule()->getUniqueIdCount();
	start_id=start_id+1000; //make sure that test molecules do not overlap with actual molecules

	//AllcopyMs.clear(); AllorigMs.clear();


	if(verbose){mappingSet[0]->get(0)->getMolecule()->printDetails();mappingSet[1]->get(0)->getMolecule()->printDetails();}

	for(unsigned int k=0; k<n_reactants; k++) {

		AllcopyMs[k].clear();
		AllorigMs[k].clear();

		if (verbose){
			cout<<"Copying molecules: ["<<k+1<<"/"<<n_reactants<<"] with total number of "<< mappingSet[k]->getNumOfMappings() <<" mappings"<<endl; //mypause(-1);
		}
		//for(unsigned int p=0; p<mappingSet[k]->getNumOfMappings();p++) {
		Molecule *origM = mappingSet[k]->get(0)->getMolecule();
//cout<<"orig molecule before copying"; origM->printDetails(); mypause(-1);

		//find all other connected molecules
		origMs.clear(); copyMs.clear();
		origM->CopybreadthFirstSearch(origM, copyM, origMs, copyMs, maxDepth, start_id, verbose);
		start_id=start_id+copyMs.size();


		//productMptr.push_back(copyM);
		if (verbose) {
//origM->printDetails();
			cout<<"one molecule is cloned, name:"<< origM->getMoleculeTypeName() <<"   mol id: "<<  origM->getUniqueID() <<"-->"<< copyM->getUniqueID()<<"   complex id: "<<  origM->getComplexID() <<"-->"<< copyM->getComplexID()<<endl;
//cout<<"copy molecule after copying"; copyM->printDetails(); mypause(-1);
		}

		AllorigMs[k] = origMs; //keep pointer to all original molecules
		AllcopyMs[k] = copyMs; //keep pointer to all copy molecules

		ms = new MappingSet(mappingSet[k], copyMs);
		//mappingSet[k]->printDetails();ms->printDetails();

		check_mappingSet[k] = ms;

	}
	if (transformationSet->getNmappingSets() != n_reactants){
		cerr<<"Inconsistency between the number of reactants"<<n_reactants <<" and transformations:"<< transformationSet->getNmappingSets()<<"!!!";
		result = 0; //exit(0);
	}


	if (result){
	//cout<<"press a key to apply transformation...."; mypause(-1);

		//if(verbose){mappingSet[0]->get(0)->getMolecule()->printDetails();check_mappingSet[0]->get(0)->getMolecule()->printDetails();mappingSet[1]->get(0)->getMolecule()->printDetails();			check_mappingSet[1]->get(0)->getMolecule()->printDetails();} //mypause(-1);		}


		this->transformationSet->transform(check_mappingSet, true);  //apply the transformation on test molecules
	//cout<<"press a key to evaluate RHS func(1)...."; mypause(-1);

		DORproductIndex = this->transformationSet->RHSreactantIndex;
	//cout<<"DORproductIndex:"<<  DORproductIndex <<endl;

		//Razi: apply the composite function cfo on molecule defined by DORproductIndex
		Molecule * mol = check_mappingSet[DORproductIndex]->get(0)->getMolecule();

		if (!mol){ cerr<<"Error when running RHS reaction. No valid molecule to apply RHS function..."; exit(0);}
		int scope = LocalFunction::SPECIES; 		//razi: check if the scope is correct !!!!LocalFunction::MOLECULE   LocalFunction::RHS_SPECIES      LocalFunction::RHS_MOLECULE
	//cout<<"press a key to evaluate RHS func(2)...."; mypause(-1);



		// just for test, later delete
		result = cfo->evaluateOnProduct(mol, scope, CompositeFunction::EvalConditionalPart, false); //verbose=false
		if(result!=0) result=1; //razi: the result is either pass or fail
		if (verbose) {cout<<"RHS evaluation result: "<<result<<endl; }//mypause(-1);


		//razi: extra test, later delete
		/*
		 /if (result == 1){
			mappingSet[0]->get(0)->getMolecule()->printDetails();check_mappingSet[0]->get(0)->getMolecule()->printDetails();mappingSet[1]->get(0)->getMolecule()->printDetails();			check_mappingSet[1]->get(0)->getMolecule()->printDetails();
			mypause(-1);
		}*/


		//result = cfo->evaluateOnProduct(mol, scope, CompositeFunction::EvalConditionalPart, false);

	}

	//delete test molecules and mappings
	for(unsigned int k=0; k<n_reactants; k++){
		if (check_mappingSet[k]){
			if (verbose) cout<<"k:"<<k<< "  pointer:"<< check_mappingSet[k]<< " mol id:" <<check_mappingSet[k]->get(0)->getMolecule()->getUniqueID()<<endl;
			delete check_mappingSet[k];
		}else break;
	}
	for(unsigned int k=0; k<n_reactants; k++){
		if(AllorigMs[k].size()>0){
			for (list <Molecule *>::iterator ii=AllorigMs[k].begin(); ii!=AllorigMs[k].end();ii++)
				if (*ii) {(*ii)->copyMptr=0; (*ii)->setCopy(NULL); if (verbose) cout<<"copy cleared for molecule id:"<<(*ii)->getUniqueID()<<endl;} //forget the copies made when firing this reactions
		}
		//else{if (verbose) cout<<"Can not find AllorigMs["<<k<<"]\n"; break;}
	}
	for(unsigned int k=0; k<n_reactants; k++){
		if (verbose) cout<<"deleting all copy molecules for reactant:"<<k<<endl; mypause(10);
		if(AllcopyMs[k].size()>0){
			for (list <Molecule *>::iterator jj=AllcopyMs[k].begin(); jj!=AllcopyMs[k].end();jj++){
				if (*jj){
					//if (verbose) cout<<"deleting one copy molecule:" << (*jj)->getUniqueID()<<" for reactant:"<<k<<endl;
					try{delete (*jj);}	catch (...){ cout<<"Unknown Error when Firing RHS reaction!";}
				}
				//else {if (verbose)cout<<"  null pointer, no molecule is deleted"<<endl; mypause(100);}
			}
		}
		//else {if (verbose) cout<<"no molecule in AllcopyMs["<<k<<"]\n"; break;}
	}

	//if (verbose) {cout<<"RHS evaluation on test molecule is finished. result: "<<result<<endl;mypause(20);}

	if (result == 0)
		return false;
	else if (result==1)
		return true;
	else{
		cerr<<"Invalid RHS function for reaction:"<< this->getName() <<", evaluation result:" <<result<<endl; exit(0);
	}


//	}
//	catch (...)
//	{
//		cout<<"Unknown Error";
//	}
//	catch (std::exception& e) {
//	    	cout<<"Error in RHS reaction check: "<<e.what();
//	}
//	catch (const runtime_error& error)
//	{
//		cout<<"Runtime Error in RHS reaction check: "<<error.what();
//	}
}



void RHSRxnClass::printDetails() const
{
	cout<<"RHSRxnClass: " << name <<"  ( baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -|"<< this->getReactantCount(r)<<" mappings|\t";
		cout<<this->reactantTemplates[r]->getPatternString()<<"\n";
	}

/*
	if (RAZI_DEBUG & SHOW_FIRE) //Razi: print details of ractants: for some reasons it does not work, later check (e.g.null pointers)
		for(unsigned int i=0; i<n_reactants; i++){
			try{
				reactantLists[i]->printDetails();
			}catch(...){
				cout<<"Error occurred while printing details of the DOR reaction.\n";
			}
	}
	*/

	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}



