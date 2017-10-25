#include "NFcore.hh"
#include "../NFutil/setting.hh" //Razi added for debugging purpose, last update 2017-3-29

using namespace std;
using namespace NFcore;



ReactionClass::ReactionClass(string name, double baseRate, string baseRateParameterName, TransformationSet *transformationSet, System *s)
{
	bool verbose=false;
	if (RAZI_DEBUG & CREATE_REACTION)
		verbose = system->getverbose();


	//cout<<"\n\ncreating reaction "<<name<<endl;
	this->system=s;
	this-> tagged = false;

	totalRateFlag=false;
	isDimerStyle=false;
	//Setup the basic properties of this reactionClass
	this->name = name;
	this->baseRate = baseRate;
	this->baseRateParameterName=baseRateParameterName;
	this->fireCounter = 0;
	this->a = 0;
	this->traversalLimit = ReactionClass::NO_LIMIT;
	this->transformationSet = transformationSet;


	//Set up the template molecules from the transformationSet
	this->n_reactants   = transformationSet->getNreactants();
	this->n_mappingsets = transformationSet->getNmappingSets();
	this->reactantTemplates = new TemplateMolecule *[n_reactants];
	vector <TemplateMolecule*> tmList;
	vector <int> hasMapGenerator;


	for(unsigned int r=0; r<n_reactants; r++)
	{
		//The main reactant should be the one that is getting modified...
		//In other words, we select the reactant that has at least one map generator, and
		//to minimize mistakes, with the least sym sites...
		TemplateMolecule *curTemplate = transformationSet->getTemplateMolecule(r);
		TemplateMolecule::traverse(curTemplate,tmList,TemplateMolecule::FIND_ALL);

		//First, single out all the templates that have at least one map generator
		for(unsigned int i=0; i<tmList.size(); i++) {
			//cout<<"looking at:"<<endl;
			//tmList.at(i)->printDetails();

			if(tmList.at(i)->getN_mapGenerators()>0) {
				hasMapGenerator.push_back(i);
			}
		}

		//Find the one with the least sym comp bonds...


#ifdef FIX_BUG1   //Razi: This bypasses the condition of having less symmetric sites, because it causes some problem in evaluating local functions.
		//The local function tries to be evalaued on main reactant
		//For now we bypassing the condition of having less symmetric sites and chose the reactant that appears earlier on the transformationSet object and has at least one mapGenerator
		if (hasMapGenerator.size()<1) {
			cout<<"There is no molecule with at least one transformation for this reactant, Perhaps something is wrong, I am quitting!!!\n"; exit(1);
		}else{
				curTemplate = tmList.at(hasMapGenerator.at(0));
		}
#else
		int minSymSites = 999999;
		for(unsigned int k=0; k<hasMapGenerator.size(); k++) {
			if(tmList.at(hasMapGenerator.at(k))->getN_symCompBonds()<minSymSites) {
				curTemplate = tmList.at(hasMapGenerator.at(k));
				minSymSites = curTemplate->getN_symCompBonds();
			}
		}
#endif
		reactantTemplates[r] = curTemplate;
		tmList.clear(); hasMapGenerator.clear();
	}
	mappingSet = new MappingSet *[n_mappingsets];
#ifdef RHS_FUNC //Razi added to support RHS functions
	check_mappingSet = new MappingSet *[n_mappingsets];
#endif


	/* create blank mappingSets for the added molecules. These will be used
	 * to hold mappings to added molecules, which is useful for rules that create
	 *  molecules and then perform other transformations.  --Justin, 1Mar2011
	 */
	for ( unsigned int r = n_reactants; r < n_mappingsets; ++r )
	{
		mappingSet[r] = transformationSet->generateBlankMappingSet(r,0);
	}



	//Here, if we identify that there are disjoint sets in this pattern, from
	//the connected-to syntax, then we have to flag the ones that we actually
	//have to traverse down...
	for(unsigned int r=0; r<n_reactants; r++)
	{
		tmList.clear();

		// Get the connected set of molecules
		TemplateMolecule *curTemplate = reactantTemplates[r];
		TemplateMolecule::traverse(curTemplate,tmList,TemplateMolecule::FIND_ALL);

		//Label the unique sets, and only continue if we have more than one set
		vector <vector <TemplateMolecule *> > sets;
		vector <int> uniqueSetId;
		int setCount = TemplateMolecule::getNumDisjointSets(tmList,sets,uniqueSetId);
		if(setCount<=1) continue;


		//count the number of map generators (rxn centers) in each set
		vector <int> numMapGenerators;
		for(int s=0; s<setCount; s++) { numMapGenerators.push_back(0); }

		int curTemplateSetId = -1;
		for(unsigned int t=0;t<tmList.size();t++) {
			if(tmList.at(t)==curTemplate) {
				curTemplateSetId = uniqueSetId.at(t);
			}
			int n_maps = numMapGenerators.at(uniqueSetId.at(t));
			numMapGenerators.at(uniqueSetId.at(t)) = n_maps+tmList.at(t)->getN_mapGenerators();
		}

		// Debug output
		//cout<<"found "<<setCount<<" unique sets."<<endl;
		//cout<<"found that reactant molecule head is in set: "<<curTemplateSetId<<endl;
		//for(int s=0; s<setCount; s++) {
		//	cout<<"set: "<<s<<" has "<<numMapGenerators.at(s)<<" map generators."<<endl;
		//}
		//for(unsigned int i=0; i<tmList.size(); i++) {
		//	cout<<"looking at:"<<endl;
		//	tmList.at(i)->printDetails();
		//}


		//Lets rearrange the connected-to elements so that the one head is listed as
		//connected to all other molecules.  This will better suit our needs.

		// first, clear out the old connections
		for(unsigned int i=0; i<tmList.size(); i++) {
			tmList.at(i)->clearConnectedTo();
		}

		// add back the connections, but always through the head template
		int rxnCenterSets = 1;
		int curSet=0;
		for(unsigned int i=0; i<uniqueSetId.size(); i++) {
			if(uniqueSetId.at(i)==curTemplateSetId) {
				if(curSet==curTemplateSetId) curSet++;
				continue;
			}
			if(uniqueSetId.at(i)==curSet) {
				bool otherHasRxnCenter = false;
				if(numMapGenerators.at(curSet)>0) {
					otherHasRxnCenter=true;
					rxnCenterSets++;
				}
				TemplateMolecule *otherTemplate = tmList.at(i);
				int ctIndex1=curTemplate->getN_connectedTo();
				int ctIndex2=otherTemplate->getN_connectedTo();
				curTemplate->addConnectedTo(otherTemplate,ctIndex2,otherHasRxnCenter);
				otherTemplate->addConnectedTo(curTemplate,ctIndex1);
				curSet++;
			}
		}

		if(rxnCenterSets>2) {
			cout.flush();
			cerr<<"\n\n   Error in Reaction Rule: "<<name<<endl;
			cerr<<"   You created a reaction with a pattern that includes the connected-to\n";
			cerr<<"   syntax (ie: A().B()).  You included 3 or more disjoint sets of molecules\n";
			cerr<<"   where there are more than 2 sets with rxn centers.  This may work ok, \n";
			cerr<<"   but you really shouldn't ever do something this crazy, so I'm just going\n";
			cerr<<"   to stop you now.  Goodbye.\n"<<endl;
			exit(1);
		}


		//cout<<"++++++++++++++++"<<endl;
		//for(unsigned int i=0; i<tmList.size(); i++) {
		//	tmList.at(i)->printDetails();
		//}


		//Finally, clear out the data structures.
		for(unsigned int i=0; i<sets.size(); i++) sets.at(i).clear();
		sets.clear(); uniqueSetId.clear();
		numMapGenerators.clear();
	}


	//Check here to see if we have molecule types that are the same across different reactants
	//Because if so, we will give a warning
	if(n_reactants>2) cerr<<"Warning!! You created a reaction ("<< name <<") that has more than 2 reactants.  This has not been extensively tested!"<<endl;

	if(n_reactants==2)
	{	//If the reactants are of the same type, then we have to make a few special considerations
		if(reactantTemplates[0]->getMoleculeType()==reactantTemplates[1]->getMoleculeType())
		{
			cout<<endl;
			cout<<"Warning! You have a binding rxn (" << name << ") that allows a moleculeType to bind another of the same type."<<endl;
			cout<<"Make sure that is correct, because this can potentially make long polymers or large aggregates."<<endl;
			cout<<endl;
		}
	}
	/*
	if (verbose){// cls();
		for(unsigned int r=0; r<n_reactants; r++)
			try{cout<<"\n\n\nAAA3-DORreactantIndex:"<<r<<"   MolType:"<< reactantTemplates[r]->getMoleculeType()->getName()<<endl;}catch(...){cerr<<"Error in AAA3\n";}
	}
	*/

	if ( this->transformationSet->usingSymmetryFactor() )
	{	// new general method for handling reaction center symmetry
		baseRate *= this->transformationSet->getSymmetryFactor();
	}
	else
	{	// old method for handling symmetric binding and unbinding
		if(n_reactants==2)
		{
			//If the binding is symmetric
			if(transformationSet->hasSymBindingTransform()) {
				cout<<endl;
				cout<<"Warning! You have an binding rxn (" << name << ") that is symmetric."<<endl;
				cout<<"Make sure that is correct."<<endl;

				cout<<endl;
				baseRate = baseRate*0.5;  //We have to correct the rate to get the proper factor
				isDimerStyle=true;
			}
		}
		if(n_reactants==1)
		{
			if(transformationSet->hasSymUnbindingTransform())
			{
				cout<<endl;
				cout<<"Warning! You have an unbinding rxn (" << name << ") that is symmetric."<<endl;
				cout<<"Make sure that is correct."<<endl;
				cout<<endl;
				baseRate = baseRate*0.5;  //We have to correct the rate to get the proper factor
				isDimerStyle=true;
			}
		}
	}


	onTheFlyObservables=true;


	// check for population type reactants
	isPopulationType = new bool[n_reactants];
	for( unsigned int i=0; i < n_reactants; ++i )
	{
		isPopulationType[i] = reactantTemplates[i]->getMoleculeType()->isPopulationType();
	}


	// calculate discrete count corrections for symmetric population reactants
	//  e.g. number of reactant pairs = A*(A-1)/2.  Note that the factor of two
	//  is part of the symmetry factor above.
	identicalPopCountCorrection = new int[n_reactants];
	for ( int i=0; i < (int)n_reactants; ++i )
	{
		identicalPopCountCorrection[i] = 0;
		if ( isPopulationType[i] )
		{
			for ( int j=i-1; j >= 0; --j )
			{
				if ( reactantTemplates[i]->getMoleculeType() == reactantTemplates[j]->getMoleculeType() )
				{
					identicalPopCountCorrection[i] = identicalPopCountCorrection[j] + 1;
					break;
				}
			}
		}
	}
}



ReactionClass::~ReactionClass()
{
	delete [] reactantTemplates;
	delete transformationSet;
	for ( unsigned int r = n_reactants; r < n_mappingsets; ++r )
	{
		delete mappingSet[r];
	}

	delete [] mappingSet;
	delete [] isPopulationType;
	delete [] identicalPopCountCorrection;
}


void ReactionClass::setBaseRate(double newBaseRate,string newBaseRateName) {
	if ( this->transformationSet->usingSymmetryFactor() )
	{	this->baseRate = this->transformationSet->getSymmetryFactor() * newBaseRate;   }
	else if (isDimerStyle)
	{	this->baseRate = 0.5 * newBaseRate;   }
	else
	{	this->baseRate = newBaseRate;   }

	this->baseRateParameterName = newBaseRateName;
	update_a();
};


void ReactionClass::resetBaseRateFromSystemParamter() {

	if(!this->baseRateParameterName.empty()) {
		if ( transformationSet->usingSymmetryFactor() ) {
			this->baseRate = transformationSet->getSymmetryFactor() * system->getParameter(this->baseRateParameterName);
		}
		else {
			// TODO: do we need to handle DimerStyle here?? --Justin
			this->baseRate=system->getParameter(this->baseRateParameterName);
		}
		this->update_a();
	}

}



void ReactionClass::printDetails() const {
	cout<< name <<"  (id="<<this->rxnId<<", baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -|"<< this->getReactantCount(r)<<" mappings|\t";
		cout<<this->reactantTemplates[r]->getPatternString()<<"\n";
		//cout<<"head: "<<endl; this->reactantTemplates[r]->printDetails(cout);
		//reactantTemplates[r]->printDetails();
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
	cout<<"\n";
}


void ReactionClass::fire(double random_A_number) {

	bool verbose = false;
	if (RAZI_DEBUG & (SHOW_FIRE | SHOW_RHS))
		verbose = this->system->getverbose();


	//cout<<endl<<">FIRE "<<getName()<<endl;
	fireCounter++;

	if(verbose){
		cout<<"\n=======================================================================================\n";
		if (RAZI_DEBUG & SHOW_FIRE){ // & ((name.compare(TEST_REAC)==0) | (TEST_REAC == "ALL")) ){
			this->printFullDetails();  //mypause(-1);
		}else{
			this->printDetails();
		}
	}


	//cout<<"Now, you printed the details. Lets continue on picking a random molecule and fire the reaction.\n";


	// First randomly pick the reactants to fire by selecting the MappingSets
	this->pickMappingSets(random_A_number);



	if(verbose){ //Razi commented & ((getName()==TEST_REAC) | (TEST_REAC == "ALL")) ){
		cout<<"These molecules are finally chosen:\n";  //show complexes that are chosen
		for (unsigned int ii=0; ii< n_reactants; ii++){
			if (mappingSet[ii])
				cout<<"\tMapping id:"<< mappingSet[ii]->getId()<< "   complex id:"<< mappingSet[ii]->getComplexID()<<" molecule list id:"<< mappingSet[ii]->get(0)->getMolecule()->getMolListId() <<" molecule id:"<< mappingSet[ii]->get(0)->getMolecule()->getUniqueID() <<  " MT:"<< mappingSet[ii]->get(0)->getMolecule()->getMoleculeTypeName()<<endl;
			else
				cout<<"mappsingSet ["<<ii<<"] for Reaction Class is not set yet." <<endl;
		}

	}



	// Check reactants for correct molecularity:
	if ( ! transformationSet->checkMolecularity(mappingSet) ) {
		// wrong molecularity!  this is a NULL event
		++(System::NULL_EVENT_COUNTER);
		return;
	}


	// output something if the reaction was tagged
	if(tagged && verbose) {
		cout<<"#RT "<<this->rxnId<<" "<<this->system->getCurrentTime();
		for(unsigned int k=0; k<n_reactants; k++) {
			cout<<" [";
			for(unsigned int p=0; p<mappingSet[k]->getNumOfMappings();p++) {
				Molecule *mForTag = mappingSet[k]->get(p)->getMolecule();
				cout<<" "<<mForTag->getMoleculeTypeName()<<mForTag->getUniqueID();
			}
			cout<<" ]";
		}
		cout<<endl;
	}


#ifdef RHS_FUNC
	// Razi: the first check is to make sure that the products are not connected
	// Check #1: Check connections between the reaction products

	// Check #2: if RHS function is included, we need to keep firing reaction, until the output complexes satisfy the function
	// in order to do so, we keep picking reactants, make copies, check the output and repeat again
	// If retry until a predefined number of time defined as RHS_maxRETRY, setting this number to 1 results in true rates !
	// Steps:
	//    1-make a copy of selected complexes
	//    2-apply the transformation on the cloned objects
	//    3-check the RHS function


	if (!this->checkReaction()){
		if (verbose || 1) cout<<"\t\tReaction:" << name <<" did not pass the output product test !!!\n";
		return;
	}else{
		if (verbose || 1) cout<<"output check passed ....\n";
	}



#endif



	// Generate the set of possible products that we need to update
	// (excluding new molecules, we'll get those later --Justin)
	this->transformationSet->getListOfProducts(mappingSet,products,traversalLimit);
	if(verbose){
		cout<<"Product molecules are: ";
		for( molIter = products.begin(); molIter != products.end(); molIter++ ) {
			if (*molIter) cout<<"M:"<<(*molIter)->getMoleculeTypeName()<<"_Lid"<< (*molIter)->getMolListId()<<",    "; else cout<< " invalid, ";
		} cout<<endl;
	}


	// display product molecules for debugging..
	//for( molIter = products.begin(); molIter != products.end(); molIter++ ) {
	//	cout<<">>molecule: "<<(*molIter)->getMoleculeTypeName()<<endl;
	//	(*molIter)->printDetails();
	//	cout<<"<<"<<endl;
	//}


	// Loop through the products (excluding added molecules) and remove from observables
	if (this->onTheFlyObservables) {

		// molecule observables..
		for ( molIter = products.begin(); molIter != products.end(); molIter++ )
			(*molIter)->removeFromObservables();

		// species observables..
		if(system->getNumOfSpeciesObs()>0) {
			// we can find reactant complexes by following mappingSets to target molecules
			int matches = 0;
			Complex * c;
			for ( unsigned int k=0; k<transformationSet->getNreactants(); k++) {
				// get complexID and check if we've already updated that complex
				int complexId = mappingSet[k]->get(0)->getMolecule()->getComplexID();
				if ( std::find( updatedComplexes.begin(), updatedComplexes.end(), complexId ) == updatedComplexes.end() ) {
					// complex has not been updated, so do it now.
					updatedComplexes.push_back(complexId);
					c = mappingSet[k]->get(0)->getMolecule()->getComplex();
					for(int i=0; i<system->getNumOfSpeciesObs(); i++) {
						matches = system->getSpeciesObs(i)->isObservable(c);
						system->getSpeciesObs(i)->straightSubtract(matches);
					}
				}
			}

			// grab added molecules that are represented as populations and remove from observables
			for ( int k=0; k<transformationSet->getNumOfAddMoleculeTransforms(); k++)
			{
				Molecule * addmol = transformationSet->getPopulationPointer((unsigned int)k);
				if ( addmol == NULL ) continue;

				// get complexID and check if we've already updated that complex
				int complexId = addmol->getComplexID();
				if ( std::find( updatedComplexes.begin(), updatedComplexes.end(), complexId ) == updatedComplexes.end() ) {
					// complex has not been updated, so do it now.
					updatedComplexes.push_back(complexId);
					c = addmol->getComplex();
					for (int i=0; i < system->getNumOfSpeciesObs(); i++) {
						matches = system->getSpeciesObs(i)->isObservable(c);
						system->getSpeciesObs(i)->straightSubtract(matches);
					}
				}
			}
			updatedComplexes.clear();
		}
	}


	// Through the MappingSet, transform all the molecules as necessary
	//  This will also create new molecules, as required.  As a side effect,
	//  deleted molecules will be removed from observables.
	this->transformationSet->transform(this->mappingSet);


	// Add newly created molecules to the list of products
	this->transformationSet->getListOfAddedMolecules(mappingSet,products,traversalLimit);


	// if complex bookkeeping is on, find all product complexes
	// (this is useful for updating Species Observables and TypeII functions, so keep the info handy).
	// NOTE: this is a brute force approach: check complex of each molecule. there may be a more
	//  elegant way to do this, but it's tricky to get it right.
	if (system->isUsingComplex()) {
		Complex * complex;
		for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
			// skip dead molecules
			if ( ! (*molIter)->isAlive() ) continue;
			// get complexID and check if we've already updated that complex
			complex = (*molIter)->getComplex();
			if ( std::find( productComplexes.begin(), productComplexes.end(), complex ) == productComplexes.end() )
				productComplexes.push_back(complex);
		}
	}


	// If we're handling observables on the fly, tell each molecule to add itself to observables.
	if (onTheFlyObservables) {

		// molecule observables..
		for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
			// skip dead molecules
			if ( ! (*molIter)->isAlive() ) continue;
			(*molIter)->addToObservables();
		}

		// species observables..
		if (system->getNumOfSpeciesObs()>0) {
			Complex * c;
			int matches;
			// we can assume that complex bookkeeping is enabled..
			for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter ) {
				// update all species observables for this complex
				c = *complexIter;
				matches = 0;
				for ( int i=0; i < system->getNumOfSpeciesObs(); i++ ) {
					matches = system->getSpeciesObs(i)->isObservable(c);
					system->getSpeciesObs(i)->straightAdd(matches);
				}
			}

			// NOTE: we don't need to handle added population types separately since they are
			//  among the product molecules
		}
	}


	// Now update reaction membership, functions, and update any DOR Groups
	//  also, gather a list of typeII dependencies that will require updating
	typeII_products.clear();
	for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
		Molecule * mol = *molIter;
		MoleculeType * mt = mol->getMoleculeType();

		// If this moleculeType has typeII dependencies, add it to the list
		// (do this for alive and dead molecules, since molecule deletion may influence
		//    the value of a local function)
		if ( mt->getNumOfTypeIIFunctions() > 0 ) {
			if ( std::find( typeII_products.begin(), typeII_products.end(), mt ) == typeII_products.end() )
				typeII_products.push_back( mt );
		}

		//Update this molcule's reaction membership
		//  NOTE: as a side-effect, DORreactions that depend on molecule-scoped local functions
		//   (typeI relationship) will be updated as long as UTL is set appropriately.
		if ( mol->isAlive() )
			mol->updateRxnMembership();
	}


	// update complex-scoped local functions for typeII dependencies
	// NOTE: as a side-effect, dependent DOR reactions (via typeI molecule dependencies) will be updated
	if (system->getEvaluateComplexScopedLocalFunctions()) {
		// for each typeII product molecule, update all dependent local functions
		if (system->isUsingComplex()) {
			// this is the easy way: update all typeI molecules on each complex
			for ( typeII_iter = typeII_products.begin(); typeII_iter != typeII_products.end(); ++typeII_iter ) {
				MoleculeType * mt = *typeII_iter;
				for (int i=0; i < mt->getNumOfTypeIIFunctions(); i++) {
					for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter )
						mt->getTypeIILocalFunction(i)->evaluateOn( *complexIter );
				}
			}
		}
		else {
			// this is the hard way: find a representative molecule from each connected set
			//  and evaluate TypeII functions on that representative.
			list <Molecule *> allMols;
			Molecule * mol;
			for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
				mol = *molIter;
				if ( std::find( allMols.begin(), allMols.end(), mol ) == allMols.end() ) {
					// remember everything connected to this molecule
					//  (so we don't evaluate this connected set multiple times)
					mol->traverseBondedNeighborhood( allMols, ReactionClass::NO_LIMIT );
					// evaluate typeII local functions on this connected set
					for ( typeII_iter = typeII_products.begin(); typeII_iter != typeII_products.end(); ++typeII_iter ) {
						MoleculeType * mt = *typeII_iter;
						for (int i=0; i<mt->getNumOfTypeIIFunctions(); i++)
							mt->getTypeIILocalFunction(i)->evaluateOn( mol, LocalFunction::SPECIES );
					}
				}
			}
		}
	} // done updating complex-scoped local functions


	//display final product molecules for debugging..

#if RAZI_DEBUG == 0x0FFF
	if(system->getverbose() && 0 && (RAZI_DEBUG & SHOW_RHS))
	for( molIter = products.begin(); molIter != products.end(); molIter++ ) {
		cout<<">>molecule: "<<(*molIter)->getMoleculeTypeName()<<endl;
	 	(*molIter)->printDetails();
	  	cout<<"<<"<<endl;
	} mypause(1000);
#endif


	//Tidy up
	products.clear();
	productComplexes.clear();
}




#ifdef RHS_FUNC

bool ReactionClass::checkReaction()   //clone reactants
{

	//return true;  //Razi: bypass check only for test, later uncomment
	bool verbose = false;
	if (RAZI_DEBUG & (SHOW_FIRE | SHOW_RHS))
		verbose = this->system->getverbose();
	//verbose = false;
	if (verbose) cout<<"RxnClass:: checkReaction() is called.\n";



	bool check_products =system->get_check_products();
	if ((this->reactionType != RHS_RXN) && (!check_products)){
		//cout<<"skip product checking. neither RHS function exist, nor ring check is active.\n";
		return true;
	}
	//proceed with checking output reactants


	list <Molecule *>  origMs(100); //at most 100 connected molecules for each reactant
	list <Molecule *>  copyMs(100);
	vector<list <Molecule *> >  AllcopyMs(10, list <Molecule *>(100)); //at most 10 reactants per reaction rule
	vector<list <Molecule *> >  AllorigMs(10, list <Molecule *>(100));
	Molecule * copyM;
	int maxDepth=10;
	bool result=true;


//	try{

	//clone molecules and complexes
	//int start_id = Molecule::uniqueIdCount+1000;
	int start_id = mappingSet[0]->get(0)->getMolecule()->getUniqueIdCount();
	start_id=start_id+1000; //make sure that test molecules do not overlap with actual molecules



	if(verbose){mappingSet[0]->get(0)->getMolecule()->printDetails();mappingSet[1]->get(0)->getMolecule()->printDetails();}

	//razi: lets first copy all reactants and their connected molecules
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



	//problem here
	if (transformationSet->getNmappingSets() != n_reactants){
		cerr<<"Inconsistency between the number of reactants"<<n_reactants <<" and transformations:"<< transformationSet->getNmappingSets()<<"!!!";
		result = false; //exit(0);
	}


	//proceed with checking the reaction
	if (result){
	//cout<<"press a key to apply transformation...."; mypause(-1);

		//if(verbose){mappingSet[0]->get(0)->getMolecule()->printDetails();check_mappingSet[0]->get(0)->getMolecule()->printDetails();mappingSet[1]->get(0)->getMolecule()->printDetails();			check_mappingSet[1]->get(0)->getMolecule()->printDetails();} //mypause(-1);		}






		//in this block, check if the products are disjoint molecules
		if (check_products){
			result = this->transformationSet->transform(check_mappingSet, true, false);  //apply the transformation on test molecules
			if ((verbose |1) && (!result)){
				cout<<"The reaction does not pass the output molecularity check. Perhaps an alternative connection exist after unbinfding molecules...\n";
			}
			if ((verbose |1) && (result)){
				cout<<"The reaction passed the output molecularity check.\n";
			}

		}else{
			//apply the reaction to the test molecules
			this->transformationSet->transform(check_mappingSet, true, false);  //apply the transformation on test molecules
		}




		//razi: now it is time to check the RHS function
		if (result && (this->reactionType == RHS_RXN)){
			//cout<<"press a key to evaluate RHS func(1)...."; mypause(-1);

			DORproductIndex = this->transformationSet->RHSreactantIndex;
		//cout<<"DORproductIndex:"<<  DORproductIndex <<endl;
			//Razi: apply the composite function cfo on molecule defined by DORproductIndex
			Molecule * mol = check_mappingSet[DORproductIndex]->get(0)->getMolecule();

			if (!mol){ cerr<<"Error when running RHS reaction. No valid molecule to apply RHS function..."; exit(0);}
			int scope = LocalFunction::SPECIES; 		//razi: check if the scope is correct !!!!LocalFunction::MOLECULE   LocalFunction::RHS_SPECIES      LocalFunction::RHS_MOLECULE
		//cout<<"press a key to evaluate RHS func(2)...."; mypause(-1);


			// just for test, later delete
			double res = cfo->evaluateOnProduct(mol, scope, CompositeFunction::EvalConditionalPart, false); //verbose=false
			if(res==0) {
				if (verbose|1) cout<<"\t\tRHS Reaction:" << name <<" DID NOT pass the RHS test. Continue firing ....\n";
				result=false;
			}else {

				if (verbose|1) cout<<"\t\tRHS Reaction:" << name <<" passed the RHS test. Continue firing ....\n";
				result=true; //razi: the result is either pass or fail
			}

		}
	}









	//razi: cleanup
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

	return result;
}
#endif





