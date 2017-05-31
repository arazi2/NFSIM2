
#include "transformationSet.hh"
#include "../../NFutil/setting.hh" //razi added for debugging purpose, last update 2017-3-29
using namespace NFcore;




list <Molecule *> TransformationSet::deleteList;
list <Molecule *> TransformationSet::updateAfterDeleteList;
list <Molecule *>::iterator TransformationSet::it;

TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates)
{
	this->hasSymUnbinding=false;
	this->hasSymBinding = false;

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->n_addmol  = 0;

#ifdef RHS_FUNC  //razi added to support RHS functions
	this->n_productTemplates = 0;
	this->includeRHSFunc=false;
	this->RHSreactantIndex=-1;
#endif

	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	this->addmol = new TemplateMolecule *[n_addmol];

	// complex bookkeeping is off by default
	this->complex_bookkeeping = false;

	// for now, symmetry factor is off by default
	this->useSymmetryFactor = false;
	this->symmetryFactor = 1.0;

	// check collisions is off by default
	this->check_collisions = false;

	//Set up our transformation vectors


	if(RAZI_DEBUG & RUN_REACTIONS){
		if (n_reactants>1){
			cout<<"\tTransformationSet object created for Templates["<<n_reactants<<"]="<<reactantTemplates[0]->getPatternString()<< ", "<<reactantTemplates[1]->getPatternString(); cin.get();
		}
		else{
			if (n_reactants==1){
				cout<<"\tTransformationSet object created for Template["<<n_reactants<<"]="<<reactantTemplates[0]->getPatternString(); cin.get();
			}
		}
		//string TM=" ";		for (i=0; i<n_reactants; i++) {TM= std::strcat(TM, reactantTemplates[i]->getMoleculeTypeName()); TM=strcat(TM, ", ");}
	}

	this->transformations = new vector <Transformation *> [n_reactants];
	finalized = false;
}

#ifdef RHS_FUNC  //razi added to support RHS functions
TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates,
        vector <TemplateMolecule *> addMoleculeTemplates,
        vector <TemplateMolecule *> OutputTemplates)
{

	this->hasSymUnbinding = false;
	this->hasSymBinding   = false;

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	this->includeRHSFunc=false; //by default, functions are LHS

	//Razi: Lets add output product Templates
	this->n_productTemplates = OutputTemplates.size();
	this->productTemplates = new TemplateMolecule *[n_productTemplates];
	for(unsigned int r=0; r<n_productTemplates; r++)
		this->productTemplates[r] = OutputTemplates.at(r);

	if (addMoleculeTemplates.size()>0){		//Remember our add molecules
		this->n_addmol = addMoleculeTemplates.size();
		this->addmol = new TemplateMolecule *[n_addmol];
		for(unsigned int r=0; r<n_addmol; r++)
			this->addmol[r] = addMoleculeTemplates.at(r);

		//Set up our transformation vectors
		this->transformations = new vector <Transformation *> [ this->getNmappingSets() ];
		cout<<"\t"<<n_addmol<< "\tNew molecules are added, when creating transformation set [new format].\n";
	}
	else{
		this->n_addmol  = 0;
		this->addmol = new TemplateMolecule *[n_addmol];
		this->transformations = new vector <Transformation *> [n_reactants];

		//cout<<"\tNo new molecule is added, when creating transformation set [new format].\n";
	}
	// complex bookkeeping is off by default
	this->complex_bookkeeping = false;

	// for now, symmetry factor is off by default
	this->useSymmetryFactor = false;
	this->symmetryFactor = 1.0;

	// check collisions is off by default
	this->check_collisions = false;


	finalized = false;

}
#endif


TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates,
		                             vector <TemplateMolecule *> addMoleculeTemplates )
{
	this->hasSymUnbinding = false;
	this->hasSymBinding   = false;

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	//Remember our add molecules
	this->n_addmol = addMoleculeTemplates.size();
	this->addmol = new TemplateMolecule *[n_addmol];
	for(unsigned int r=0; r<n_addmol; r++)
		this->addmol[r] = addMoleculeTemplates.at(r);

	// complex bookkeeping is off by default
	this->complex_bookkeeping = false;

	// for now, symmetry factor is off by default
	this->useSymmetryFactor = false;
	this->symmetryFactor = 1.0;

	// check collisions is off by default
	this->check_collisions = false;

	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [ this->getNmappingSets() ];
	finalized = false;
}


TransformationSet::~TransformationSet()
{
	for(unsigned int r=0; r<getNmappingSets(); r++)  {
		Transformation *t;
		while(transformations[r].size()>0)
		{
			t = transformations[r].back();
			transformations[r].pop_back();
			delete t;
		}
	}

	Transformation *t;
	while(addMoleculeTransformations.size()>0)
	{
		t = addMoleculeTransformations.back();
		addMoleculeTransformations.pop_back();
		delete t;
	}

	while(addSpeciesTransformations.size()>0)
	{
		t = addSpeciesTransformations.back();
		addSpeciesTransformations.pop_back();
		delete t;
	}

	delete [] transformations;
	delete [] reactants;
	delete [] addmol;

	this->n_reactants = 0;
	this->n_addmol = 0;
}

#ifdef RHS_FUNC //razi added
void TransformationSet::printDetails(){
	cout<<"\n--------------------------------------------------------------------------\n";
	cout<<"\tPrinting details of transformation set...."<<endl;
	cout<<"--------------------------------------------------------------------------\n";
	if(includeRHSFunc) cout<<"\tRHS func exists.\n"; else   cout<<"\tRHS func does not exist.\n";

	cout<<"\tReactants: ";
	for (int i=0; i< n_reactants; i++) cout<<reactants[i]->getPatternString()<<",  ";  cout<<endl;

	cout<<"\tProducts: ";
	for (int i=0; i< n_productTemplates; i++) cout<<productTemplates[i]->getPatternString()<<",  ";  cout<<endl;
	cout<<"--------------------------------------------------------------------------\n";
}

TemplateMolecule * TransformationSet::getTemplateMolecule(unsigned int reactantIndex)
{
	return getTemplateMolecule(reactantIndex, false);
}

TemplateMolecule * TransformationSet::getTemplateMolecule(unsigned int reactantIndex, bool RHSfunc) const
{
	if ((reactantIndex < n_reactants ) && (!RHSfunc))
	{
		return reactants[reactantIndex];
	}
	if ( reactantIndex < getNmappingSets() && (!RHSfunc))
	{
		return addmol[reactantIndex-n_reactants];
	}
	if ((reactantIndex < n_productTemplates ) && (RHSfunc))
	{
		return productTemplates[reactantIndex];
	}
	cerr<<"requesting invalid molecule template from a transformationSet object.\n";
	return NULL;
}



//Razi added: return the number of transformations associated with input and output templates to support RHS functions
int TransformationSet:: getNumOfTransformations(int reactantIndex){
	return getNumOfTransformations(reactantIndex, false);
}



//Razi added: return the number of transformations associated with input and output templates to support RHS functions
int TransformationSet:: getNumOfTransformations(int reactantIndex, bool RHSfunc) const{
//	if (RHSfunc)
//		return product_transformations[reactantIndex].size();
//	else
		return transformations[reactantIndex].size();
};



Transformation * TransformationSet::getTransformation(int reactantIndex, int index) {
	return getTransformation(reactantIndex, index, false);
}

//Razi added: return transformations associated with input and output templates to support RHS functions
Transformation * TransformationSet::getTransformation(int reactantIndex, int index, bool RHSfunc) const{
//	if (RHSfunc)
//	return product_transformations[reactantIndex].at(index);
//	else
		return transformations[reactantIndex].at(index);
};


#else
TemplateMolecule *
TransformationSet::getTemplateMolecule( unsigned int reactantIndex ) const
{
	if ( reactantIndex < n_reactants )
	{
		return reactants[reactantIndex];
	}
	else if ( reactantIndex < getNmappingSets() )
	{
		return addmol[reactantIndex-n_reactants];
	}
}
#endif


bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string cName, int finalStateValue)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addStateChangeTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	//cout<<"Adding state change transform to value: "<<finalStateValue<<endl;
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genStateChangeTransform(cIndex, finalStateValue);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);


	if (RAZI_DEBUG & RUN_REACTIONS){
		cout<<"\tStateChange Transformation for ReactId:"<< reactantIndex<<"  TM:"<< t->getPatternString()<<",  Comp:" << cName <<", Val:" <<finalStateValue << ", Total Trans for This Reactant:"<<transformations[reactantIndex].size();
		cin.get();
	}

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}

#ifdef RHS_FUNC //Razi: added to support RHS functions
bool TransformationSet::addLocalFunctionReference(TemplateMolecule *t, string PointerName, int scope)
{
	return addLocalFunctionReference(t, PointerName, scope, false); //Razi: Default is LHS function
}

bool TransformationSet::addLocalFunctionReference(TemplateMolecule *t, string PointerName, int scope, bool RHSfunc)
{

	//Razi: scope here means either LocalFunction::SPECIES or LocalFunction::MOLECULE (not local or global)
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	//Razi: For functions, where the argument is one of the output product templates, ...
	// Find the reactant whose first molecule matches the first molecule of the corresponding product
	//For isntance: X.Y + Z  -> u:X + Y + Z       [the reactant X.Y is relevant]
	//ReactantIndex includes the relevant reactant index [see find()]
	int ReactantIndex = find(t, RHSfunc);

	if(ReactantIndex==-1) {
		cerr<<"Couldn't find the Output template you gave me!  In transformation set - addStateChangeTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	if (RHSfunc){
		if (RAZI_DEBUG & SHOW_SIM) {
			cout<<"\n\tRHS Reaction: found Local Func for ReactId:"<< ReactantIndex<<"  TM:"<< t->getPatternString()<<",  PointerName:" << PointerName <<", Scope:" <<scope<<endl ;}

		this->RHSreactantIndex=ReactantIndex; //Razi: contains the reactant index with RHS function
		includeRHSFunc=true;
		return true;
/*
		cout <<"transformationSet: I don't know how to add local function reference for RHS functions.\n"; mypause(-1);
		//Transformation *transformation = TransformationFactory::genLocalFunctionReference(PointerName,scope,t);
		//transformations[reactantIndex].push_back(transformation);
		Transformation *transformation = TransformationFactory::genLocalFunctionReference(PointerName,scope,t);
		product_transformations[ReactantIndex].push_back(transformation);

		if (RAZI_DEBUG  & RUN_REACTIONS){
			cout<<"\n\tAdd RHS Local Func for ProductId:"<< ReactantIndex<<"  TM:"<< t->getPatternString()<<",  PointerName:" << PointerName <<", Scope:" <<scope << ", Total Trans for This product :"<<product_transformations[ReactantIndex].size();
			cin.get();
		}


		MapGenerator *mg = new MapGenerator(product_transformations[ReactantIndex].size()-1);
		t->addMapGenerator(mg, 1);
		return true;
*/

	}else{
		if(includeRHSFunc) cerr<<"The function type is LHS, but RHS flag is already set. Perhaps two functions are included for a reaction which is not allowed.\n"; mypause(10000);
		includeRHSFunc=false; //make sure it is set to false
		Transformation *transformation = TransformationFactory::genLocalFunctionReference(PointerName,scope,t);
		transformations[ReactantIndex].push_back(transformation);

		if (RAZI_DEBUG  & (CREATE_FUNC | RUN_REACTIONS)){
			cout<<"\n\tAdd Local Func for ReactId:"<< ReactantIndex<<"  TM:"<< t->getPatternString()<<",  PointerName:" << PointerName <<", Scope:" <<scope << ", Total Trans for This Reactant:"<<transformations[ReactantIndex].size();
		}


		MapGenerator *mg = new MapGenerator(transformations[ReactantIndex].size()-1);
		t->addMapGenerator(mg);
		return true;
	}
}
#else
bool TransformationSet::addLocalFunctionReference(TemplateMolecule *t, string PointerName, int scope)
{
	//Razi: scope here means either LocalFunction::SPECIES or LocalFunction::MOLECULE (not local or global)
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addStateChangeTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	Transformation *transformation = TransformationFactory::genLocalFunctionReference(PointerName,scope,t);
	transformations[reactantIndex].push_back(transformation);

	if (RAZI_DEBUG  & RUN_REACTIONS){
		cout<<"\n\tAdd Local Func for ReactId:"<< reactantIndex<<"  TM:"<< t->getPatternString()<<",  PointerName:" << PointerName <<", Scope:" <<scope << ", Total Trans for This Reactant:"<<transformations[reactantIndex].size();
		cin.get();
	}


	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}
#endif

bool TransformationSet::addIncrementStateTransform(TemplateMolecule *t, string cName)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addIncrementStateTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genIncrementStateTransform(cIndex);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}
bool TransformationSet::addDecrementStateTransform(TemplateMolecule *t, string cName)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDecrementStateTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genDecrementStateTransform(cIndex);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}



bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string cName, string finalStateValue)
{
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	int fStateValue = t->getMoleculeType()->getStateValueFromName(cIndex,finalStateValue);
	return TransformationSet::addStateChangeTransform(t,cName, fStateValue);
}


bool TransformationSet::addBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	//Again, first find the reactants that the binding pertains to
	int reactantIndex1 = find(t1);
	int reactantIndex2 = find(t2);
	if(reactantIndex2==-1 || reactantIndex2==-1) {
		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	//Find the index of the respective binding sites
	unsigned int cIndex1 = t1->getMoleculeType()->getCompIndexFromName(bSiteName1);
	unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);


	//Check for symmetric binding
	bool isSymmetric = TemplateMolecule::checkSymmetry(t1,t2,bSiteName1,bSiteName2);
	if( isSymmetric )
		hasSymBinding = true;


	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
	//equal to the size.
	Transformation *transformation1;
	if(reactantIndex1==reactantIndex2)
		transformation1 = TransformationFactory::genBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1);
	else
		transformation1 = TransformationFactory::genBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size());

	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2);


	if (RAZI_DEBUG & RUN_REACTIONS){
		cout<<"\n\tAdd Binding Transform for ReactIds:"<< reactantIndex1 <<", "<<reactantIndex2<<"  TMs:"<< t1->getPatternString()<<", "<< t2->getPatternString()<<",  Sites:" << bSiteName1<<", "<<bSiteName2 <<", SiteIndexes: "<<cIndex1 <<", "<< cIndex2<< ", Total Trans for These Reactant:"<<1+transformations[reactantIndex1].size()<<", "<<1+transformations[reactantIndex2].size();
		cin.get();
	}

	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);

	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex2].size()-1);
	t2->addMapGenerator(mg2);

	return true;
}



bool TransformationSet::addNewMoleculeBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	//Again, first find the reactants that the binding pertains to
	int reactantIndex1 = find(t1);
	int reactantIndex2 = find(t2);
	if(reactantIndex2==-1 || reactantIndex2==-1) {
		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	//Find the index of the respective binding sites
	unsigned int cIndex1 = t1->getMoleculeType()->getCompIndexFromName(bSiteName1);
	unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);


	//Check for symmetric binding
	bool isSymmetric = TemplateMolecule::checkSymmetry(t1,t2,bSiteName1,bSiteName2);
	if( isSymmetric )
		hasSymBinding = true;


	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
	//equal to the size.
	Transformation *transformation1;
	if(reactantIndex1==reactantIndex2)
		transformation1 = TransformationFactory::genNewMoleculeBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1);
	else
		transformation1 = TransformationFactory::genNewMoleculeBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size());

	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2);

	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);

	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex2].size()-1);
	t2->addMapGenerator(mg2);

	return true;
}




// deprecated
//
//bool TransformationSet::addBindingSeparateComplexTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
//{
//	cout<<"adding separate complex binding"<<endl;
//	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
//	//Again, first find the reactants that the binding pertains to
//	int reactantIndex1 = find(t1);
//	int reactantIndex2 = find(t2);
//	if(reactantIndex2==-1 || reactantIndex2==-1) {
//		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!\n";
//		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
//		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
//		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
//		return false;
//	}
//
//	//Find the index of the respective binding sites
//	unsigned int cIndex1 = t1->getMoleculeType()->getCompIndexFromName(bSiteName1);
//	unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);
//
//
//	//Check for symmetric binding
//	bool isSymmetric = TemplateMolecule::checkSymmetry(t1,t2,bSiteName1,bSiteName2);
//	if( isSymmetric )
//		hasSymBinding = true;
//
//
//
//	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
//	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
//	//equal to the size.
//	Transformation *transformation1;
//	if(reactantIndex1==reactantIndex2)
//		transformation1 = TransformationFactory::genBindingSeparateComplexTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1);
//	else
//		transformation1 = TransformationFactory::genBindingSeparateComplexTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size());
//
//	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2);
//
//	transformations[reactantIndex1].push_back(transformation1);
//	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
//	t1->addMapGenerator(mg1);
//
//	transformations[reactantIndex2].push_back(transformation2);
//	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex1].size()-1);
//	t2->addMapGenerator(mg2);
//
//	return true;
//}



bool TransformationSet::addUnbindingTransform(TemplateMolecule *t, string bSiteName, TemplateMolecule *t2, string bSiteName2)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	TemplateMolecule *tToTransform = 0;
	if(t==0 && t2==0) {
		cerr<<"Error in transformation set! when creating unbinding transform!"<<endl;
		cerr<<"Both molecules you gave me are null!\n";
		return false;
	} else if(t2==0) {
		tToTransform = t;
	} else if(t==0) {
		tToTransform = t2;
	} else {
		// they are both real, so randomly pick t1
		tToTransform=t;

		//Check for symmetric unbinding
		bool isSymmetric = TemplateMolecule::checkSymmetryAroundBond(t,t2,bSiteName,bSiteName2);
		if( isSymmetric ) {
			hasSymUnbinding = true;
		}

	}

	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(tToTransform);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you might get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	unsigned int cIndex = tToTransform->getMoleculeType()->getCompIndexFromName(bSiteName);
	Transformation *transformation = TransformationFactory::genUnbindingTransform(cIndex);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	if (RAZI_DEBUG & RUN_REACTIONS){
		cout<<"\n\tAdd UnBinding Transform for ReactId:"<< reactantIndex <<",  TM:"<<tToTransform->getPatternString()<< ",  out of TMs:"<< t->getPatternString()<<", "<< t2->getPatternString()<<",  Sites:" << bSiteName<<", "<<bSiteName2 <<", Selected SiteIndex: "<<cIndex <<", Total Trans for Thess Reactant:"<<transformations[reactantIndex].size();
		cin.get();
	}

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	tToTransform->addMapGenerator(mg);

	return true;
}





/*!
	Adds a delete rule to the given TemplateMolecule.
	@author Michael Sneddon
*/
bool TransformationSet::addDeleteMolecule(TemplateMolecule *t, int deletionType) {
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDeleteMolecule!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}
	Transformation *transformation = TransformationFactory::genRemoveMoleculeTransform(deletionType);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}


/*!
	Adds a decrement population rule to the given TemplateMolecule.
	@author Justin Hogg
*/
bool TransformationSet::addDecrementPopulation(TemplateMolecule *t)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDecrementPopulation!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}
	Transformation *transformation = TransformationFactory::genDecrementPopulationTransform();

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}

bool TransformationSet::addAddSpecies( SpeciesCreator *sc )
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	// We don't need a polymorphic transform because AddTransforms are handled separately!
	//  But we do need to call some methods specific to AddMoleculeTransform.
	//  So we're modified TransformationFactory to return the specific object type  --JUstin
	AddSpeciesTransform * transformation = TransformationFactory::genAddSpeciesTransform( sc );

	// 3) Add the transformation object to the TransformationSet
	addSpeciesTransformations.push_back( transformation );

	// 3) No map generators needed for an add species!
	return true;
}


bool TransformationSet::addAddMolecule( MoleculeCreator *mc )
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	// We don't need a polymorphic transform because AddTransforms are handled separately!
	//  But we do need to call some methods specific to AddMoleculeTransform.
	//  So we're modified TransformationFactory to return the specific object type  --JUstin
	AddMoleculeTransform * transformation = TransformationFactory::genAddMoleculeTransform( mc );

	// 3) Add the transformation object to the TransformationSet
	addMoleculeTransformations.push_back( transformation );

	// 3) No map generators needed for an add molecule!
	return true;
}



#ifdef RHS_FUNC  //Razi added to support RHS functions
int TransformationSet::find(TemplateMolecule *t){
	return find(t, 0);
}
int TransformationSet::find(TemplateMolecule *t, bool RHSfunc)
{
//cout<<"TransformationSet::find called for "<< t->getPatternString()<< "flag: "<<RHSfunc<<endl;

	if(finalized) { cerr<<"TransformationSet cannot search for a templateMolecule once it has been finalized!"<<endl; exit(1); }
	int findindex2 = -1; int findIndex = -1;

	//Razi: first make sure that the output patterns exists (is valid)
	if (RHSfunc){ //the component belongs to the product templates
		for(unsigned int r=0; r<n_productTemplates; r++)  {
			//cout<<"product:  "<<this->productTemplates[r]->getPatternString() <<endl;
			if(this->productTemplates[r]->contains(t)) {
				if(findindex2==-1) {
					findindex2 = r;
				}
				else {
					cerr<<"Found duplicate output template molecule in two reaction lists!!  (in transformationSet)."<<endl;
					exit(1);
				}
			}
		}

		if (findindex2 != -1){
			//Razi: Find the relevant reactant
			//The assumption is the same order of molecules in the input and out of the reaction

			int molid = this->productTemplates[findindex2]->getMoleculeType()->getTypeID(); //this is the molecule we are looking for
			string moltext = this->productTemplates[findindex2]->getMoleculeType()->getName(); //this is the molecule we are looking for
			//cout<<"TransformationSet::find looking to find relevant reactant for: "<<moltext<<endl;


			//Razi: find the order of appearance in output products
			int order2=1;
			for(int i=0; (i<=this->n_productTemplates) && (i < findindex2) ; i++){
				string pattern= this->productTemplates[i]->getPatternString();
				int pos=0;
				while ((pos = pattern.find(moltext, pos)) < pattern.length())
				{
					order2++;
				    pos += moltext.length();
				}
			}

			//Razi: scan the input reactants until finding the molecules [at the same order]
			int cnt=0;
			for(int i=0; (i<=this->n_reactants)  ; i++){
				string pattern= this->reactants[i]->getPatternString();
				int pos=0;
				while ((pos = pattern.find(moltext, pos)) < pattern.length())
				{
					cnt++;
				    pos += moltext.length();
				    if (cnt==order2){ ///found the relevant reactant
				    	if (RAZI_DEBUG & (CREATE_REACTION|SHOW_FIRE)) cout<<"transformationset:find successfully, molecule:"<<moltext <<" reatant["<<i <<"]:"<<pattern<<"   product["<< findindex2 <<"]:" << this->productTemplates[findindex2]->getPatternString() <<endl;
				    	findIndex = i; return findIndex;
				    }
				}
			}
		}
	}
	else{  //the component belongs to the reactant templates
		for(unsigned int r=0; r<n_reactants; r++)  {
			if(this->reactants[r]->contains(t)) {
				if(findIndex==-1) {
					findIndex = r;
				}
				else {
					cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
					exit(1);
				}
			}
		}
		// also check add molecule templates
		for(unsigned int r=0; r<n_addmol; r++)  {
			if(this->addmol[r]->contains(t)) {
				if(findIndex==-1) {
					findIndex = r + n_reactants;
				}
				else {
					cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
					exit(1);
				}
			}
		}
	}
	return findIndex;
}
#else
int TransformationSet::find(TemplateMolecule *t)
{
	if(finalized) { cerr<<"TransformationSet cannot search for a templateMolecule once it has been finalized!"<<endl; exit(1); }
	int findIndex = -1;
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(this->reactants[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r;
			}
			else {
				cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
				exit(1);
			}
		}
	}
	// also check add molecule templates
	for(unsigned int r=0; r<n_addmol; r++)  {
		if(this->addmol[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r + n_reactants;
			}
			else {
				cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
				exit(1);
			}
		}
	}
	return findIndex;
}
#endif


#ifdef RHS_FUNC //Razi added to support RHS functions
bool TransformationSet::transform(MappingSet **mappingSets, bool testmode)
{
	if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	// addMolecule transforms are applied before other transforms so the molecules exist
	//  for potential modification by other transforms.

//cout<<"AAA1 ";mypause(-1);
//cout<<"addMoleculeTransformations.size():"<<addMoleculeTransformations.size()<<"   addSpeciesTransformations.size():"<<addSpeciesTransformations.size()<<endl;


	if(!testmode){
		int size = addMoleculeTransformations.size();
		if(size>0) {
			for(int i=0; i<size; i++) {
				addMoleculeTransformations.at(i)->apply_and_map( mappingSets[n_reactants+i]  );
			}
		}

		// apply addSpecies transforms, so we have all the molecules out there
		size = addSpeciesTransformations.size();
		if(size>0) {
			for(int i=0; i<size; i++) {
				addSpeciesTransformations.at(i)->apply(NULL,NULL);
			}
		}
	}

//cout<<"AAA2 nMappingSet:"<<getNmappingSets()<<"  ";mypause(-1);
	// loop over reactants and added molecules, apply transforms to each
	for(unsigned int r=0; r<getNmappingSets(); r++)
	{
		MappingSet *ms = mappingSets[r];
		for ( unsigned int t=0;  t<transformations[r].size();  t++ )
		{
//cout<<"BBB2 r:"<<r<<"  t:"<<t<<"  ";mypause(-1);
			if(transformations[r].at(t)->getType()==(int)TransformationFactory::REMOVE )
			{	// handle deletions
//cout<<"BBB2-1 r:"<<r<<"  t:"<<t<<"  "; mypause(-1);
				if (!testmode) {
					Molecule * mol = ms->get(t)->getMolecule();
					if ( transformations[r].at(t)->getRemovalType()==(int)TransformationFactory::COMPLETE_SPECIES_REMOVAL )
					{	// complex deletion: flag connected molecules for deletion
						mol->traverseBondedNeighborhood(deleteList,ReactionClass::NO_LIMIT);
					}
					else
					{	// molecule deletion: flag this molecule for deletion
						deleteList.push_back( mol );
					}
				}
			}
			else
			{	// handle other transforms
				//if (RAZI_DEBUG & SHOW_FIRE) {cout<<"Transform:  reactant:"<<r<<"  transformation:"<<t<<"  transformation type:"<< transformations[r].at(t)->getType()<<".\n";//ms->get(t)->printDetails(); mypause(-1);}

				transformations[r].at(t)->apply(ms->get(t), mappingSets);
			}
		}
	}


	if (!testmode) {
		//Each molecule that is on the delete list must be dealt with
		Molecule * mol;
		for( it = deleteList.begin(); it!=deleteList.end(); it++)
		{
			mol = *it;
			mol->getMoleculeType()->removeMoleculeFromRunningSystem(mol);
		}
		deleteList.clear();
	}
	return true;
}

bool TransformationSet::transform(MappingSet **mappingSets){
	return TransformationSet::transform(mappingSets, false);
}

#else
bool TransformationSet::transform(MappingSet **mappingSets)
{
	if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	/*
	 * NOTE: Check for "null conditions" was moved to ReactionClass::fire. This allows rejection of a reaction
	 * prior to removing molecules from observables. A general method TransformationSet::checkMolecularity has
	 * been implemented to check for incorrect molecularity or reaction center conflicts. --Justin
	 */


	// addMolecule transforms are applied before other transforms so the molecules exist
	//  for potential modification by other transforms.

//cout<<"AAA1 ";mypause(-1);
	int size = addMoleculeTransformations.size();
	if(size>0) {
		for(int i=0; i<size; i++) {
			addMoleculeTransformations.at(i)->apply_and_map( mappingSets[n_reactants+i]  );
		}
	}

	// apply addSpecies transforms, so we have all the molecules out there
	size = addSpeciesTransformations.size();
	if(size>0) {
		for(int i=0; i<size; i++) {
			addSpeciesTransformations.at(i)->apply(NULL,NULL);
		}
	}

//cout<<"AAA2 nMappingSet:"<<getNmappingSets()<<"  ";mypause(-1);
	// loop over reactants and added molecules, apply transforms to each
	for(unsigned int r=0; r<getNmappingSets(); r++)
	{
		MappingSet *ms = mappingSets[r];
		for ( unsigned int t=0;  t<transformations[r].size();  t++ )
		{
//cout<<"BBB2 r:"<<r<<"  t:"<<t<<"  ";mypause(-1);
			if( transformations[r].at(t)->getType()==(int)TransformationFactory::REMOVE )
			{	// handle deletions
//cout<<"BBB2-1 r:"<<r<<"  t:"<<t<<"  "; mypause(-1);
				Molecule * mol = ms->get(t)->getMolecule();
				if ( transformations[r].at(t)->getRemovalType()==(int)TransformationFactory::COMPLETE_SPECIES_REMOVAL )
				{	// complex deletion: flag connected molecules for deletion
					mol->traverseBondedNeighborhood(deleteList,ReactionClass::NO_LIMIT);
				}
				else
				{	// molecule deletion: flag this molecule for deletion
					deleteList.push_back( mol );
				}
			}
			else
			{	// handle other transforms
				//if (RAZI_DEBUG & SHOW_FIRE) {cout<<"Transform:  reactant:"<<r<<"  transformation:"<<t<<"  transformation type:"<< transformations[r].at(t)->getType()<<".\n";//ms->get(t)->printDetails(); mypause(-1);}

				transformations[r].at(t)->apply(ms->get(t), mappingSets);
			}
		}
	}


	//Each molecule that is on the delete list must be dealt with
	Molecule * mol;
	for( it = deleteList.begin(); it!=deleteList.end(); it++)
	{
		mol = *it;
		mol->getMoleculeType()->removeMoleculeFromRunningSystem(mol);
	}
	deleteList.clear();

	return true;
}

#endif






bool TransformationSet::checkMolecularity( MappingSet ** mappingSets )
{
	if ( n_reactants < 2 )
	{	// unimolecular, so there's nothing to check
		return true;
	}
	else if ( complex_bookkeeping )
	{	// verify that each reactant pattern points to a unique complex
		complex_ids.clear();
		for ( unsigned int ir = 0;  ir < n_reactants;  ++ir )
		{
			// skip populations
			if ( reactants[ir]->getMoleculeType()->isPopulationType() ) continue;

			complex_id = mappingSets[ir]->getComplexID();
			complex_id_iter = std::find( complex_ids.begin(), complex_ids.end(), complex_id );
			if ( complex_id_iter == complex_ids.end() )
			{
				complex_ids.push_back( complex_id );
			}
			else
			{   // two reactant patterns matched the same complex!
				return false;
			}
		}
		return true;
	}
	else if ( check_collisions )
	{	// we won't do a proper check for molecularity, but we should ensure that mappingSets
		//  point to non-overlapping reaction centers.
		for ( collision_pair_iter = collision_pairs.begin(); collision_pair_iter != collision_pairs.end(); ++collision_pair_iter )
		{
			if ( MappingSet::checkForCollisions( mappingSets[ (*collision_pair_iter).first  ],
					                             mappingSets[ (*collision_pair_iter).second ]  ) )
			{	// reaction centers overlap!
				return false;
			}
		}
		return true;
	}
	else
	{   // do nothing
		return true;
	}
}


bool TransformationSet::getListOfProducts(MappingSet **mappingSets, list <Molecule *> &products, int traversalLimit)
{
	//if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	list <Molecule *>::iterator molIter;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		// if we are deleting the entire complex, we don't have to track molecules in this complex
		if (mappingSets[r]->hasSpeciesDeletionTransform()) continue;

		//cout<<"Traversing:"<<endl;
		//mappingSets[r]->get(0)->getMolecule()->printDetails();
		//mappingSets[r]->get(0)->getMolecule()->traverseBondedNeighborhood(products,traversalLimit);

		/*
		 * I thought that making sure we don't go over the same molecule multiple
		 * times would make the code faster - but this is rarely used for most rxn
		 * systems, so it is commented out for now.  But actually, we do have to check
		 * because if have the same molecule in here twice, then it can mess up our
		 * observable lists...  --michael */
		else
		{
			// For each of the molecules that we possibly affect, traverse the neighborhood
			// Q: Is it sufficient to just look at the first mapping?
			// A: It should be if the traversal limit is set high enough, at least for
			// all standard reactions.  I'm wondering now, though, if it is enough in
			// all cases where you would use the connected-to syntax.  I think so, but
			// someone should test it.  --michael 9Mar2011
			Molecule * molecule = mappingSets[r]->get(0)->getMolecule();

			// is this molecule already on the product list?
			if ( std::find( products.begin(), products.end(), molecule ) == products.end() )
			{	// Traverse neighbor and add molecules to list
				molecule->traverseBondedNeighborhood(products,traversalLimit);
				//molecule->traverseBondedNeighborhoodForUpdate(products,traversalLimit);
			}
		}
	}

	// Next, find added molecules that are treated as populations.
	//  Populations molecules have to be removed from observables, then incremented,
	//  and then added back to the observables (Add molecules treated as particles are handled later)
	vector <AddMoleculeTransform *>::iterator addmol_iter;
	for ( addmol_iter = addMoleculeTransformations.begin();
			addmol_iter != addMoleculeTransformations.end();  ++addmol_iter )
	{
		// get molecule creator
		AddMoleculeTransform * addmol = *addmol_iter;
		if ( !(addmol->isPopulationType()) ) continue;

		// Get the population molecule pointer
		Molecule * molecule = addmol->get_population_pointer();

		// is this molecule already on the product list?
		if ( std::find( products.begin(), products.end(), molecule ) == products.end() )
		{	// Add molecule to list
			products.push_back( molecule );
		}
	}

	//cout<<"All together, we have: "<<products.size()<<endl;
	return true;
}


Molecule * TransformationSet::getPopulationPointer( unsigned int r ) const
{
	return addMoleculeTransformations.at(r)->isPopulationType()
			? addMoleculeTransformations.at(r)->get_population_pointer()
		    : NULL;
}


bool TransformationSet::getListOfAddedMolecules(MappingSet **mappingSets, list <Molecule *> &products, int traversalLimit)
{
	//if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	// Add new molecules (particle type) to the list of products
	list <Molecule *>::iterator molIter;
	for (unsigned int r=n_reactants; r<getNmappingSets(); r++)
	{
		//For each of the molecules that we possibly affect, traverse the neighborhood
		// NOTE: in this instance, it's okay to only look at the first mapping
		Molecule * molecule = mappingSets[r]->get(0)->getMolecule();

		// Skip populations
		if ( molecule->isPopulationType() ) continue;

		// Is the molecule already in the products list?  If not, add to list.
		if ( std::find( products.begin(), products.end(), molecule ) == products.end() )
		{	// Add molecule to list.
			products.push_back( molecule );
			// NOTE: we don't need to traverse neighbors. All new molecules will be put in this
			//  list separately and old molecules that bind to new molecules will be traversed elsewhere
		}
	}

	return true;
}


MappingSet *TransformationSet::generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId)
{
	if(!finalized) { cerr<<"TransformationSet cannot generate blank mapping if it is not finalized!"<<endl; exit(1); }
	if( reactantIndex>=getNmappingSets() ) {
		cerr<<"Gave me (a transformation Set) a reactant index that was too high!"<<endl;
		exit(1);
	}
	return new MappingSet(mappingSetId, transformations[reactantIndex]);
}

void TransformationSet::finalize()
{
	//Be sure to add at least a blank transformation to every reactant if there is no transformation
	//specified so that we count the reactants even if we don't do anything to it.
	for(unsigned int r=0; r<getNmappingSets(); r++)  {
		if(transformations[r].size()==0) {
			transformations[r].push_back(TransformationFactory::genEmptyTransform());
			MapGenerator *mg = new MapGenerator(transformations[r].size()-1);
#ifdef RHS_FUNC
			getTemplateMolecule(r,false)->addMapGenerator(mg, false);   //Razi: check if similar process is needed for product reactants
#else
			getTemplateMolecule(r)->addMapGenerator(mg);
#endif
			}
	}

	// Determine if we need to do any reactant center overlap checks.
	// Currently we check a necessary (but not sufficient) condition for the
	//   possibility of reactant center overlap: are there a common molecule types in
	//   a pair of reactant templates. In the future, we could check a necessary and sufficent condition
	//   (e.g. pattern overlap) to avoid extra work.

	//Razi note: collision_pairs include mappings between reactants (non population types) that share common molecule types
	if ( (n_reactants>1)  &&  !complex_bookkeeping )
	{
		vector <TemplateMolecule *> tmList1;
		vector <TemplateMolecule *> tmList2;
		vector <TemplateMolecule *>::iterator tm_iter;

		vector <MoleculeType *> moltypes;
		vector <MoleculeType *>::iterator  found_iter;

		for ( unsigned int ir1 = 0;  ir1 < n_reactants;  ++ir1 )
		{
			// skip populations
			if ( reactants[ir1]->getMoleculeType()->isPopulationType() ) continue;

			tmList1.clear();
			moltypes.clear();

			// get all the template molecules in reactant pattern ir1
			TemplateMolecule::traverse(reactants[ir1], tmList1, false);

			// collect the molecule types included in the pattern
			for ( tm_iter = tmList1.begin(); tm_iter != tmList1.end();  ++tm_iter )
			{
				moltypes.push_back( (*tm_iter)->getMoleculeType() );
			}

			for ( unsigned  int ir2 = ir1 + 1;  ir2 < n_reactants;  ++ir2 )
			{
				// skip populations
				if ( reactants[ir2]->getMoleculeType()->isPopulationType() ) continue;

				tmList2.clear();

				// get all the template molecules in reactant pattern ir2
				TemplateMolecule::traverse(reactants[ir2], tmList2, false);

				// check if any moleculeTypes collide
				for ( tm_iter = tmList2.begin(); tm_iter != tmList2.end();  ++tm_iter )
				{
					found_iter = std::find( moltypes.begin(), moltypes.end(), (*tm_iter)->getMoleculeType() );
					if ( found_iter != moltypes.end() )
					{
						// we found a molecule type that is common to patterns ir1 and ir2!
						check_collisions = true;
						collision_pairs.push_back( pair<int,int>(ir1,ir2) );
						break;
					}
				}
			}
		}

	}

	finalized = true;
}
