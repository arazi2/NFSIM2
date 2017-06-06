#include <iostream>
#include "NFcore.hh"
#include "../NFutil/setting.hh" //Razi added for debugging purpose, last update 2017-3-29
#include <queue>

using namespace std;
using namespace NFcore;

int Molecule::uniqueIdCount = 0;





// Molecule Constructor
//
//



Molecule::Molecule(MoleculeType * parentMoleculeType, int listId)
{
	if(DEBUG) cout<<"-creating molecule instance of type " << parentMoleculeType->getName() << endl;
	this->parentMoleculeType = parentMoleculeType;

	// set population type (1 for particle type, 0 for population type)
	this->population_count = ( parentMoleculeType->isPopulationType()  ?  0  :  1 );

	//First initialize the component states and bonds
	this->numOfComponents = parentMoleculeType->getNumOfComponents();
	this->component = new int [parentMoleculeType->getNumOfComponents()];
	for(int c=0; c<numOfComponents; c++)
		component[c] = parentMoleculeType->getDefaultComponentState(c);

	// initialize bond sites
	this->bond = new Molecule * [numOfComponents];
	this->indexOfBond = new int [numOfComponents];
	this->hasVisitedBond = new bool [numOfComponents];
	for(int b=0; b<numOfComponents; b++) {
		bond[b]=0; indexOfBond[b]=NOBOND;
		hasVisitedBond[b] = false;
	}


	hasVisitedMolecule = false;
	hasEvaluatedMolecule = false;
	isMatchedTo=0;
	rxnListMappingId2 = 0;
	nReactions = 0;
	useComplex = parentMoleculeType->getSystem()->isUsingComplex();
	isPrepared = false;
	isObservable = 0;
	localFunctionValues=0;
	//isDead = true;

	//register this molecule with moleculeType and get some ID values

	ID_complex = this->parentMoleculeType->createComplex(this);
	ID_type = this->parentMoleculeType->getTypeID();
	ID_unique = Molecule::uniqueIdCount++;
	this->listId = listId;
	isAliveInSim = false;


#ifdef RHS_FUNC //Razi:  Added to support RHS functions
	originalMptr = 0; //Razi: default is original molecule, not a cloned molecule for test
	copyMptr = 0;   //Razi: no copy exists yet
#endif

	if((RAZI_DEBUG & CREATE_MOLECULE)&& (listId<=3))	{
		std::cout<<"\tI just created a Molecule of Type:"<< this->parentMoleculeType->getName() <<" and a complex with ID:"<<ID_complex <<"  Type ID:" <<ID_type << "  Unique ID:" << ID_unique <<" List ID:" << listId<<" UsingComlexFlag:"<<useComplex<<"....\n";
	}
}







// Molecule Deconstructor
//
//
Molecule::~Molecule()
{
	if(DEBUG) cout <<"   -destroying molecule instance of type " << parentMoleculeType->getName() << endl;

//razi added to debug memory crash when deleting test molecule. the problem solved. cout<<"XX1-";
	delete [] bond;

//cout<<"XX2-";
	parentMoleculeType = 0;


	delete [] isObservable;
//cout<<"XX3-";
	delete [] component;
//cout<<"XX4-";
	delete [] indexOfBond;
//cout<<"XX5-";
	delete [] rxnListMappingId2;
//cout<<"XX6-";
	delete [] hasVisitedBond;
//cout<<"XX7-";

	if(localFunctionValues!=0)
		delete [] localFunctionValues;
//cout<<"XX8-";
}


void Molecule::prepareForSimulation()
{
	if(isPrepared) return;
	nReactions = parentMoleculeType->getReactionCount();
	this->rxnListMappingId2 = new set<int>[nReactions];

	isPrepared = true;

	//We do not belong to any observable... yet.
	isObservable=new int [parentMoleculeType->getNumOfMolObs()];
	for(int o=0;o<parentMoleculeType->getNumOfMolObs(); o++) {
		isObservable[o]=0;
	}


}






void Molecule::setUpLocalFunctionList()
{
	if (parentMoleculeType->getNumOfTypeIFunctions() > 0)
	{
		localFunctionValues=new double[parentMoleculeType->getNumOfTypeIFunctions()];
		for(int lf=0; lf<parentMoleculeType->getNumOfTypeIFunctions(); lf++) {
			localFunctionValues[lf]=0;
		}
	}
}



//Used so that this molecule can remember what its local function was
//evaluated to.  Only TypeI local functions are set up in this way
void Molecule::setLocalFunctionValue(double newValue,int localFunctionIndex) {
	if(localFunctionIndex<0 || localFunctionIndex>=parentMoleculeType->getNumOfTypeIFunctions()) {
		cout<<"Error in Molecule: trying to set the value of a local function, but the\n";
		cout<<"index provided was out of bounds!  I shall quit now."<<endl;
		std::exit(1);
	}
	//cout<<"here, mol: "<<this->getMoleculeTypeName()<<"_"<<this->getUniqueID()<<endl;
	//cout<<"localfunctionIndex given: "<< localFunctionIndex<<endl;
	//cout<<"n_localfunctionValue given: "<< localFunctionIndex<<endl;

	//cout<<"localFunctionValues="<<endl;
	//for(int k=0; k<parentMoleculeType->getNumOfTypeIFunctions(); k++) {
	//	cout<<"  "<<localFunctionValues[k]<<endl;
	//}
	//cout<<"done."<<endl;

	localFunctionValues[localFunctionIndex] = newValue;
}


double Molecule::getLocalFunctionValue(int localFunctionIndex) {
	if(localFunctionIndex<0 || localFunctionIndex>=parentMoleculeType->getNumOfTypeIFunctions()) {
		cout<<"Error in Molecule: trying to set the value of a local function, but the\n";
		cout<<"index provided was out of bounds!  I shall quit now."<<endl;
		std::exit(1);
	}
	return localFunctionValues[localFunctionIndex];
}


LocalFunction * Molecule::getLocalFunction(int localFunctionIndex) {
	if(localFunctionIndex<0 || localFunctionIndex>=parentMoleculeType->getNumOfTypeIFunctions()) {
			cout<<"Error in Molecule: trying to set the value of a local function, but the\n";
			cout<<"index provided was out of bounds!  I shall quit now."<<endl;
			std::exit(1);
		}
	return parentMoleculeType->getTypeILocalFunction(localFunctionIndex);
}



void Molecule::updateRxnMembership()
{
	parentMoleculeType->updateRxnMembership(this);
}

void Molecule::updateTypeIIFunctions()
{
	for (int i=0; i<parentMoleculeType->getNumOfTypeIIFunctions(); i++) {
		parentMoleculeType->getTypeIILocalFunction(i)->evaluateOn(this, LocalFunction::SPECIES);
	}
}

void Molecule::updateTypeIIFunctions( vector <Complex *> & productComplexes )
{
	for (int i=0; i<parentMoleculeType->getNumOfTypeIIFunctions(); i++) {
		vector <Complex *>::iterator complexIter;
		for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter ) {
			parentMoleculeType->getTypeIILocalFunction(i)->evaluateOn( (*complexIter)->getFirstMolecule(), LocalFunction::SPECIES);
		}
	}
}

void Molecule::updateDORRxnValues()
{
	ReactionClass *rxn=0; int rxnIndex=-1, rxnPos=-1;
	//cout<<"Looping over :"<<parentMoleculeType->getNumOfDORrxns()<<" dor rxns "<<endl;
	for(int i=0; i<parentMoleculeType->getNumOfDORrxns(); i++) {
		rxn=parentMoleculeType->getDORrxn(i);
		rxnIndex=parentMoleculeType->getDORrxnIndex(i);
		rxnPos=parentMoleculeType->getDORrxnPosition(i);
		//cout<<"\t\ti="<<i<<" rxn="<<rxn->getName()<<" rxnIndex="<<rxnIndex<<" rxnPos="<<rxnPos<<endl;

		if(isPrepared) {
			//If we are in this reaction, then we have to update our value...
			if(getRxnListMappingId(rxnIndex)>=0) {
				//iterate over all mappings
				set<int> tempSet = getRxnListMappingSet(rxnIndex);
				//iterate over all agent-mappings  for the same reaction
				for(set<int>::iterator it= tempSet.begin();it!= tempSet.end(); ++it){

				//Careful here!  remember to update the propensity of this
				//reaction in the system after we notify of the rate factor change!
					double oldA = rxn->get_a();
					rxn->notifyRateFactorChange(this,rxnPos,*it);
					parentMoleculeType->getSystem()->update_A_tot(rxn,oldA,rxn->update_a());
				}
			}
		}
	}
}

///////////////
//  MOLECULE_DEPENDENT_UPDATE_ADDITION
//void Molecule::addDependentUpdateMolecule(Molecule *m) {
//	for(molIter=dependentUpdateMolecules.begin();molIter!=dependentUpdateMolecules.end();molIter++)
//		if((*molIter)->getUniqueID()==m->getUniqueID())
//			return;
//	dependentUpdateMolecules.push_back(m);
//}
//void Molecule::removeDependentUpdateMolecule(Molecule *m) {
//	for(molIter=dependentUpdateMolecules.begin();molIter!=dependentUpdateMolecules.end();molIter++)
//		if((*molIter)->getUniqueID()==m->getUniqueID()) {
//			dependentUpdateMolecules.erase(molIter);
//		}
//}
////////////////



//void Molecule::updateDORs()
//{
//
//	for(int r=0; r<parentMoleculeType->getDORrxnCount(); r++)
//	{
//
//		ReactionClass * DORrxn = parentMoleculeType->getDORrxn(r);
//		int dorRxnIndex = parentMoleculeType->getDORreactantIndex(r);
//		int dorRxnPos = parentMoleculeType->getDORreactantPosition(r);
//
//	//	cout<<" identified DOR RXN index: "<<dorRxnIndex<<endl;
//	//	cout<<" identified DOR RXN pos: "<<dorRxnPos<<endl;
//		DORrxn->notifyRateFactorChange(this, dorRxnPos, rxnListMappingId[dorRxnIndex]);
//	}
//
//}

//double Molecule::getDORvalueFromGroup(string groupName, int valueIndex)
//{
////	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
////	{
////		if(groupName==(*listenerIter)->getGroupName())
////			return (*listenerIter)->getValue(valueIndex);
////	}
//
//	cerr<<"Error!! trying to get DOR value for a group, but no name match!"<<endl;
//	cerr<<"    Looking for group: "<<groupName<<" from molecule ";
//	cerr<<this->getMoleculeTypeName()<<"_"<<this->getUniqueID()<<endl;
//	exit(1);
//}



void Molecule::removeFromObservables()
{
	parentMoleculeType->removeFromObservables(this);
}
void Molecule::addToObservables()
{
	parentMoleculeType->addToObservables(this);
}


// set population
bool Molecule::setPopulation( int count )
{
	if ( isPopulationType()  &&  (count >= 0) )
	{
		population_count = count;
		return true;
	}
	else return false;
}

// get popualtion
int Molecule::getPopulation() const
{
	return population_count;
}

// increase population by one
bool Molecule::incrementPopulation()
{
	if ( isPopulationType() )
	{
		++population_count;
		return true;
	}
	else return false;
}

// decrease population by one
bool Molecule::decrementPopulation()
{
	if ( isPopulationType()  &&  (population_count > 0) )
	{
		--population_count;
		return true;
	}
	else return false;

}


void Molecule::setComponentState(int cIndex, int newValue)
{
	this->component[cIndex]=newValue;
	if (useComplex)
		// Need to manually unset canonical flag since we're not calling a Complex method
		getComplex()->unsetCanonical();

	//if(listeners.size()>0) cout<<"Molecule State has changed..."<<endl;
	//Let all the listeners know that the state of a molecule has changed...
	//for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
	//	(*listenerIter)->notify(this,stateIndex);
}
void Molecule::setComponentState(string cName, int newValue) {
	this->component[this->parentMoleculeType->getCompIndexFromName(cName)]=newValue;

	if (useComplex)
		// Need to manually unset canonical flag since we're not calling a Complex method
		getComplex()->unsetCanonical();

}


void Molecule::printDetails() {
	this->printDetails(cout);
}
void Molecule::printDetails(ostream &o)
{
	int degree = 0;
	o<<"++ Molecule instance of type: " << parentMoleculeType->getName();
	o<< " (uId="<<ID_unique << ", listId=" << listId << ", tId=" << ID_type << ", cId=" << ID_complex<<", degree="<<degree<<")"<<endl;
	o<<"      components: ";
	for(int c=0; c<numOfComponents; c++)
	{
		if(c!=0)o<<"                  ";
		o<< parentMoleculeType->getComponentName(c) <<"=";
		o<<parentMoleculeType->getComponentStateName(c,component[c]);
		o<<"\tbond=";
		if(bond[c]==NOBOND) o<<"empty";
		else {
			o<<bond[c]->getMoleculeTypeName()<<"_"<<bond[c]->getUniqueID();
			o<<"("<<bond[c]->getMoleculeType()->getComponentName(this->indexOfBond[c])<<")";
		}

		o<<endl;
	}

	o.flush();
	if(parentMoleculeType->getNumOfTypeIFunctions()>0) {
		o<<"      loc funcs:";
		for(int lf=0; lf<parentMoleculeType->getNumOfTypeIFunctions(); lf++) {
			if(lf!=0) o<<"                  ";
			o<<parentMoleculeType->getTypeILocalFunction(lf)->getNiceName();
			o<<"="<<localFunctionValues[lf]<<"\n";
		}
	}
}

//Get the number of molecules this molecule is bonded to
int Molecule::getDegree()
{
	int degree = 0;
	for(int c=0; c<numOfComponents; c++)
		if(bond[c]!=NOBOND) degree++;
	return degree;
}

// Get a label for this molecule or one of it's components (labels are not unique)
//   cIndex==-1  =>  get label for molecule, "m:typename"
//   cIndex>=0   =>  get label for component cIndex, "c:name~state"
string Molecule::getLabel ( int cIndex ) const
{
    string label("");
    if ( cIndex < 0 )
    {	// molecule label
        label += "m:" + getMoleculeTypeName();
    }
    else
    {	// component label
        label += "c:" + (   getMoleculeType()->isEquivalentComponent(cIndex)
        		          ? getMoleculeType()->getEquivalenceClassComponentNameFromComponentIndex(cIndex)
        		          : getMoleculeType()->getComponentName(cIndex) )
        		      + "~" + getMoleculeType()->getComponentStateName(cIndex, getComponentState(cIndex));
    }
    return label;
}


bool Molecule::isBindingSiteOpen(int cIndex) const
{
	if(bond[cIndex]==NOBOND) return true;
	return false;
}

bool Molecule::isBindingSiteBonded(int cIndex) const
{
	if(bond[cIndex]==NOBOND) return false;
	return true;
}

Molecule * Molecule::getBondedMolecule(int cIndex) const
{
	return bond[cIndex];
}

// given the component index, look up what we are bonded to.  Then
// in the molecule we are bonded to, look at what site we are bonded to
//
//  for instance
//
//    this(a!1).other(b!),   and we call this->getBondedMoleculeBindingSiteIndex(0)
//
//    where index 0 = this site a, then this function would return the
//    component index of b in molecule other.
//
int Molecule::getBondedMoleculeBindingSiteIndex(int cIndex) const
{
	return indexOfBond[cIndex];
}



void Molecule::bind(Molecule *m1, int cIndex1, Molecule *m2, int cIndex2)
{
	if(m1->bond[cIndex1]!=NOBOND || m2->bond[cIndex2]!=NOBOND) {
		cerr<<endl<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cerr<<"Your universal traversal limit was probably set too low, so some updates were not correct!\n\n";

		cerr<<"Trying to bond "<< m1->getMoleculeTypeName() << "_"<<m1->getUniqueID()<<"(";
		cerr<<m1->getMoleculeType()->getComponentName(cIndex1)<<") & ";
		cerr<< m2->getMoleculeTypeName()<<"_"<<m2->getUniqueID()<<"(";
		cerr<<m2->getMoleculeType()->getComponentName(cIndex2)<<")\n";
		cerr<<" to sites that are already occupied!  Check rxn rules!!\n";
		cerr<<"\n";
		m1->printDetails(cerr);
		m2->printDetails(cerr);
		std::exit(1);		//Razi commented due to an error when using RHS functions: later check==> std::exit(1);
	}

	m1->bond[cIndex1] = m2;
	m2->bond[cIndex2] = m1;

	m1->indexOfBond[cIndex1] = cIndex2;
	m2->indexOfBond[cIndex2] = cIndex1;

	//Handle Complexes
	if(m1->useComplex)
	{
		if(m1->getComplex()!=m2->getComplex())
		{
			// NOTE: mergeWithList will handle canonical flags
			m1->getComplex()->mergeWithList(m2->getComplex());
		}
		else
			// Need to manually unset canonical flag since we're not calling a Complex method
			m1->getComplex()->unsetCanonical();
	}
}

void Molecule::bind(Molecule *m1, string compName1, Molecule *m2, string compName2)
{
	int cIndex1 = m1->getMoleculeType()->getCompIndexFromName(compName1);
	int cIndex2 = m2->getMoleculeType()->getCompIndexFromName(compName2);
	Molecule::bind(m1, cIndex1, m2, cIndex2);
}


void Molecule::unbind(Molecule *m1, int cIndex)
{
	//get the other molecule bound to this site
	//cout<<"I am here. "<<bSiteIndex<<endl;
	Molecule *m2 = m1->bond[cIndex];
	if(m2==NULL)
	{
		cerr<<endl<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cerr<<"Your universal traversal limit was probably set too low, so some updates were not correct!"<<endl;
		cerr<<"Trying to unbind a binding site that is not bound!!  Check rxn rules, and traversal limits! Quitting."<<endl;
		cerr<<endl<<endl<<"The molecule is:"<<endl;
		m1->printDetails(cerr);
		cerr<<endl<<"The site trying to be unbound was: ";
		cerr<<m1->getMoleculeType()->getComponentName(cIndex)<<endl;
		//Razi commented due to an error when using RHS functions: later check==> std::exit(3);
	}


	int cIndex2 = m1->indexOfBond[cIndex];

	//break the bond (older compilers don't let you assign NOBOND to type molecule)
	m1->bond[cIndex] = 0; //NOBOND;
	m2->bond[cIndex2] = 0; //NOBOND;

	m1->indexOfBond[cIndex] = NOINDEX;
	m2->indexOfBond[cIndex2] = NOINDEX;

	//Handle Complexes
	if(m1->useComplex)
	{
		// NOTE: mergeWithList will handle canonical flags
		m1->getComplex()->updateComplexMembership(m1);
	}

	//cout<<" UnBinding!  mol1 complex: ";
	//m1->getComplex()->printDetails();
}

void Molecule::unbind(Molecule *m1, char * compName)
{
	int cIndex = m1->getMoleculeType()->getCompIndexFromName(compName);
	Molecule::unbind(m1,cIndex);
}








queue <Molecule *> Molecule::q;
queue <int> Molecule::d;
list <Molecule *>::iterator Molecule::molIter;

void Molecule::breadthFirstSearch(list <Molecule *> &members, Molecule *m, int depth)
{
	if(m==0) {
		cerr<<"Error in Molecule::breadthFirstSearch, m is null.\n";
		cerr<<"Likely an internal error where a MappingSet is on a list and\n";
		cerr<<"is not actually mapped to any molecule!";
		std::exit(3);
	}

	//Create the queues (for effeciency, now queues are a static attribute of Molecule...)
	//queue <Molecule *> q;
	//queue <int> d;
	int currentDepth = 0;

	//cout<<"traversing on:"<<endl;
	//m->printDetails();

	//First add this molecule
	q.push(m);
	members.push_back(m);
	d.push(currentDepth+1);
	m->hasVisitedMolecule=true;

	//Look at children until the queue is empty
	while(!q.empty())
	{
		//Get the next parent to look at (currentMolecule)
		Molecule *cM = q.front();
		currentDepth = d.front();
		q.pop();
		d.pop();

		//Make sure the depth does not exceed the limit we want to search
		if((depth!=ReactionClass::NO_LIMIT) && (currentDepth>=depth)) continue;

		//Loop through the bonds
		int cMax = cM->numOfComponents;
		for(int c=0; c<cMax; c++)
		{
			//cM->getComp
			if(cM->isBindingSiteBonded(c))
			{
				Molecule *neighbor = cM->getBondedMolecule(c);
				//cout<<"looking at neighbor: "<<endl;
				//neighbor->printDetails();
				if(!neighbor->hasVisitedMolecule)
				{
					neighbor->hasVisitedMolecule=true;
					members.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
					//cout<<"adding... to traversal list."<<endl;
				}
			}
		}
	}


	//clear the has visitedMolecule values
	for( molIter = members.begin(); molIter != members.end(); molIter++ )
  		(*molIter)->hasVisitedMolecule=false;
}





//
void Molecule::traverseBondedNeighborhood(list <Molecule *> &members, int traversalLimit)
{
	//always call breadth first search, it is a bit faster
	//if(traversalLimit>=0)
		Molecule::breadthFirstSearch(members, this, traversalLimit);
	//else
	//	this->depthFirstSearch(members);
}


//Isn't ever called really, but is availabe.  Note that it cannot use traversal limits
//because it is depth first
void Molecule::depthFirstSearch(list <Molecule *> &members)
{
	if(this->hasVisitedMolecule==true) {
		return;
	}

	this->hasVisitedMolecule=true;
	members.push_back(this);

	int cMax = this->numOfComponents;
	for(int c=0; c<cMax; c++)
	{
		if(hasVisitedBond[c]==true) continue;
		if(this->isBindingSiteBonded(c))
		{
			Molecule *neighbor = this->getBondedMolecule(c);
			neighbor->hasVisitedBond[indexOfBond[c]]=true;
			hasVisitedBond[c]=true;
			neighbor->depthFirstSearch(members);
		}
	}

	//clear things out
	hasVisitedMolecule = false;
	for(int c=0; c<numOfComponents; c++)
		hasVisitedBond[c] = false;
}



void Molecule::printMoleculeList(list <Molecule *> &members)
{
	cout<<"List of molecules contains: "<<endl;
	list <Molecule *>::iterator molIter;
	for( molIter = members.begin(); molIter != members.end(); molIter++ ) {
		cout<<"   -"<<(*molIter)->getMoleculeTypeName();
		cout<<"_u"<<(*molIter)->getUniqueID()<<endl;
	}
}





#ifdef RHS_FUNC //Razi:  Added to support RHS functions
void Molecule::setUniqueID(int ID){
	//this function sets unique ID, it should be used only for copy and temporary molecules, for regular molecules it is automatically set in constructor
	if (ID < Molecule::uniqueIdCount){
		cerr<<" Setting invalid(duplicate) molecule ID:" << ID << "< limit: "<< Molecule::uniqueIdCount; exit(0);
	}
	ID_unique = ID;
}


Molecule::Molecule(Molecule &obj){ //const Molecule &obj

	localFunctionValues;

	if (obj.copyMptr != 0){
		cerr<<"This molecule is already cloned:"<< obj.copyMptr<<"...\n"; exit(0);
	}

	//useComplex = obj.useComplex;
	useComplex=0; //razi: complex is not developed yet for the copy molecules
	originalMptr =(Molecule *) &obj;   //Razi: lets point to the original pointer

	//cout<<"cloning a test molecule: start";

	hasVisitedMolecule=obj.hasVisitedMolecule;
	hasEvaluatedMolecule=obj.hasEvaluatedMolecule;
	static const int NOSTATE = -1;
	static const int NOBOND = 0;
	static const int NOINDEX = -1;

	isPrepared=obj.isPrepared;
	isAliveInSim = obj.isAliveInSim;
	ID_type = obj.ID_type;

	//ID_complex=this->parentMoleculeType->createComplex(this); no need to add to the allcomplexes array
	//ID_complex=-1;
	ID_complex = obj.ID_complex;
//cout<<"Molecule::uniqueIdCount before copying:"<<Molecule::uniqueIdCount<<endl;

	//ID_unique = Molecule::uniqueIdCount++; //Razi: we do not need to increase uniqueIdCount for test molecule created here

	//cout<<"Molecule::uniqueIdCount after copying:"<<Molecule::uniqueIdCount<<endl; mypause(-1);


	parentMoleculeType = obj.parentMoleculeType;
	population_count = obj.population_count;

	numOfComponents=obj.numOfComponents;


	this->component = new int [numOfComponents];
	this->bond = new Molecule * [numOfComponents];
	this->indexOfBond = new int [numOfComponents];
	this->hasVisitedBond = new bool [numOfComponents];


	for (int i=0; i<numOfComponents; i++){
		component[i]=obj.component[i];   //later check
		bond[i]=0; //obj.bond[i];  //bond to the new relevant test molecule
		indexOfBond[i]=NOBOND;
		hasVisitedBond[i] = false;
	}

	isMatchedTo=0; //obj.isMatchedTo; //razi: this causes error when killing the test object
	rxnListMappingId2 = 0; //obj.rxnListMappingId2;  //razi: this causes error when killing the test object
	nReactions = obj.nReactions;
	isPrepared = obj.isPrepared;
	isObservable = 0; //obj.isObservable; //razi: this causes error when killing the test object
	localFunctionValues=new double[parentMoleculeType->getNumOfTypeIFunctions()];
	//localFunctionValues=0;
	//localFunctionValues=obj.localFunctionValues; //razi: this causes error when killing the test object

	//isDead = obj.isDead;

	//register this molecule with moleculeType and get some ID values

	listId = obj.listId;
	isAliveInSim = obj.isAliveInSim;

	//cout<<": end.\n"; mypause(-1);
}




//This function is developed to copy a molecule and its connected subnetwork
void Molecule::CopybreadthFirstSearch(Molecule *origM, Molecule * &copyM, list <Molecule *> &origMs, list <Molecule *> &copyMs, int maxDepth, int start_id, bool verbose)
//void Molecule::CopybreadthFirstSearch(Molecule *origM, Molecule * &copyM, vector <Molecule *> &origMs, vector <Molecule *> &copyMs, int maxDepth, int start_id)
{
	if(origM==0) {
		cerr<<"Error in Molecule::CopybreadthFirstSearch, m is null.\n";	std::exit(3);
	}
	int currentDepth = 0;


//cout<<"AAA-traversing on start"<<endl;mypause(-1);
	//m->printDetails();

	//First add this molecule
	q.push(origM);
	origMs.push_back(origM);

	//make a fresh copy or get an existing copy
	Molecule * copyM1;
	if (origM->getCopy()){
		copyM1 = origM->getCopy();
if (verbose) cout<<"AAA-get previous copy for original molecule.\n";
	}else{
		copyM1 = new Molecule(*origM);  //Razi: check
		copyM1->setUniqueID(start_id++);
		origM->setCopy(copyM1);  //Razi: let the source molecule remember its copies
if (verbose) cout<<"AAA-fresh copy for original molecule. ID:"<<start_id-1<<"="<<copyM1->getUniqueID()<<endl; //mypause(-1);
	}
	copyM=copyM1;

	origM->setCopy(copyM);   //Razi: let the source molecule remember its copies
	copyMs.push_back(copyM);

	//copyM=copyM1;
	d.push(currentDepth+1);
	origM->hasVisitedMolecule=true;


//cout<<"AAA-before while"; mypause(-1);
	Molecule * neighbor_copy;
	//Look at children until the queue is empty
	while(!q.empty())
	{

//cout<<"AAA-start while"; mypause(-1);

		//Get the next parent to look at (currentMolecule)
		Molecule *cM = q.front();
		Molecule * cM_copy = cM->getCopy();
		if (!cM_copy){
			cerr<<"Molecule::CopybreadthFirstSearch: No copy for a molecule alreaidy in the queue. I am quitting..."; exit(0);
		}

		currentDepth = d.front();
		q.pop();
		d.pop();

		//Make sure the depth does not exceed the limit we want to search
		if((maxDepth!=ReactionClass::NO_LIMIT) && (currentDepth>=maxDepth)) continue;

		//Loop through the bonds
		int cMax = cM->numOfComponents; int cIndex1,cIndex2;
		for(int c=0; c<cMax; c++)
		{
//cout<<"AAA-start while-for"; mypause(-1);

			//cM->getComp
			if(cM->isBindingSiteBonded(c))
			{
				Molecule *neighbor = cM->getBondedMolecule(c);
				//cout<<"looking at neighbor: "<<endl;
				//neighbor->printDetails();
				if(!neighbor->hasVisitedMolecule)
				{
					neighbor->hasVisitedMolecule=true;

					//make a fresh copy or get the pointer to the existing copy
					origMs.push_back(neighbor);
					if (neighbor->getCopy()){
						neighbor_copy = neighbor->getCopy();
//if (verbose) cout<<"BBB-get previous copy for a neighbor molecule.\n";
					}else{
						neighbor_copy = new Molecule(*neighbor);  //Razi: check
						neighbor_copy->setUniqueID(start_id++);
						neighbor->setCopy(neighbor_copy);   //Razi: let the source molecule remember its copies
//if (verbose) {cout<<"BBB-get a fresh copy for a neighbor molecule. ID:"<<start_id-1<<"="<<neighbor_copy->getUniqueID()<<endl; mypause(-1);}
					}
//cout<<"BBB-m id after copy: "<< m->getUniqueID()<<"  " <<endl; mypause(-1);


					cIndex1 = c;
					cIndex2 = cM->getBondedMoleculeBindingSiteIndex(cIndex1);
//if (verbose) {cout<<"Binding copy molecule when traversing m1 uID:" << cM_copy->getUniqueID() <<" index1:"<<c<<"  m2 uID:" << neighbor_copy->getUniqueID() <<"  index2:"<< cIndex2 <<endl; mypause(-1);}
					cM->bind(cM_copy, cIndex1, neighbor_copy, cIndex2); //Razi: make bonds for copy molecules
//if (verbose) {cout<<"Binding copy molecule finished successfully."<<endl; mypause(-1);}
					copyMs.push_back(neighbor_copy);


					q.push(neighbor);
					d.push(currentDepth+1);
					//cout<<"adding... to traversal list."<<endl;
				}
			}
		}
	}
//cout<<"AAA-start while end"; mypause(-1);


	//clear the has visitedMolecule values
	for( molIter = origMs.begin(); molIter != origMs.end(); molIter++ )
  		(*molIter)->hasVisitedMolecule=false;

//cout<<"AAA-middle"; mypause(-1);
//cout<<"copyMs.size():"<<copyMs.size()<<endl; mypause(-1);

	for( molIter = copyMs.begin(); molIter != copyMs.end(); molIter++ ){
		if (*molIter){
			try{
//			cout<<"good pointer 1: "; mypause(-1);
			(*molIter)->hasVisitedMolecule=false;
//			cout<<"good pointer 2"; mypause(-1);
			}catch (string err){
				cout<<"Error:"<<err;
			}
		}
	}
//cout<<"XXX-traversing on end"<<endl;mypause(-1);

}

#endif


