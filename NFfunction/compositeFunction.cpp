/*
 * compositeFunction.cpp
 *
 *  Created on: Dec 4, 2008
 *      Author: msneddon
 */

#include "NFfunction.hh"
#include "../NFutil/setting.hh" //Razi added for debugging purpose, last update 2017-3-29



using namespace std;
using namespace NFcore;
using namespace mu;


CompositeFunction::CompositeFunction(System *s,
					string name,
					string expression,
					vector <string> &functions,
					vector <string> &argNames,
					vector <string> &paramNames)
{

	//cout<<"creating composite function"<<endl;
	this->name = name;
	this->originalExpression=expression;
	this->parsedExpression="";
#ifdef RHS_FUNC  //Razi: This block is added to support RHS functions
	this->RHS = false; //false;
#endif


	this->n_params=paramNames.size();
	this->paramNames = new string[n_params];
	for(unsigned int i=0; i<n_params; i++) {
		this->paramNames[i]=paramNames.at(i);
	}

	this->n_args=argNames.size();
	this->argNames=new string[n_args];
	for(unsigned int i=0; i<n_args; i++) {
		this->argNames[i]=argNames.at(i);
	}


	this->n_allFuncs=functions.size();
	allFuncNames = new string[n_allFuncs];
	for(unsigned int i=0; i<n_allFuncs; i++) {
		this->allFuncNames[i] = functions.at(i);
	}

	p=0;

	if (RAZI_DEBUG & CREATE_FUNC){
		cout<<"\tCreating Composite Function for:\t"<<name<<" expression:" << expression<<endl;
		cout<<"\t\tParams: \t";		for(unsigned int an = 0; an < this->n_params; an++) cout<<this->paramNames[an]<< ", "; cout<<endl;
		cout<<"\t\tArguments: \t";		for(unsigned int an = 0; an < this->n_args; an++) cout<<this->argNames[an]<< ", "; cout<<endl;
		cout<<"\t\tFunctions: \t";		for(unsigned int an = 0; an < this->n_allFuncs; an++) cout<<this->allFuncNames[an]<< ", "; cout<<endl;
	}
}


CompositeFunction::~CompositeFunction()
{
	delete [] allFuncNames;

	delete [] argNames;
	delete [] paramNames;

	delete [] gfNames;
	delete [] gfs;
	delete [] gfValues;

	delete [] lfNames;
	delete [] lfs;

	delete [] refLfInds;
	delete [] refLfRefNames;
	delete [] refLfScopes;
	delete [] refLfValues;

	delete [] reactantCount;

	if(p!=NULL) delete p;

}


void CompositeFunction::setGlobalObservableDependency(ReactionClass *r, System *s) {

	for(int i=0; i<this->n_gfs; i++) {
		GlobalFunction *gf=gfs[i];

		for(int vr=0; vr<gf->getNumOfVarRefs(); vr++) {
			if(gf->getVarRefType(vr)=="Observable") {
				Observable *obs = s->getObservableByName(gf->getVarRefName(vr));
				obs->addDependentRxn(r);
			} else {
				cerr<<"When creating a FunctionalRxnClass of name: "+r->getName()+" you provided a function that\n";
				cerr<<"depends on an observable type that I can't yet handle! (which is "+gf->getVarRefType(vr)+"\n";
				cerr<<"try using type: 'MoleculeObservable' for now.\n";
				cerr<<"quiting..."<<endl; exit(1);
			}
		}
	}
}

//call this immediately after you have read in all the functions, but before preparing for the simulation
//or reading in reactions
void CompositeFunction::finalizeInitialization(System *s)
{

	//first, find the global functions by name
	vector <GlobalFunction *> gf_tempVector;
	for(unsigned int i=0; i<n_allFuncs; i++) {
		GlobalFunction *gf=s->getGlobalFunctionByName(allFuncNames[i]);
		if(gf!=0) gf_tempVector.push_back(gf);
	}

	this->n_gfs = gf_tempVector.size();
	this->gfNames = new string[n_gfs];
	this->gfValues = new double[n_gfs];
	this->gfs = new GlobalFunction * [n_gfs];
	for(int i=0; i<n_gfs; i++) {
		gfNames[i] = gf_tempVector.at(i)->getName();
		gfValues[i] = 0;
		gfs[i] = gf_tempVector.at(i);
	}


	//now the local functions...
	vector <LocalFunction *> lf_tempVector;
	for(unsigned int i=0; i<n_allFuncs; i++) {
		LocalFunction *lf=s->getLocalFunctionByName(allFuncNames[i]);
		if(lf!=0) lf_tempVector.push_back(lf);
	}

	this->n_lfs = lf_tempVector.size();
	this->lfNames = new string[n_lfs];
	this->lfs = new LocalFunction * [n_lfs];
	for(int i=0; i<n_lfs; i++) {
		lfNames[i] = lf_tempVector.at(i)->getName();
		lfs[i] = lf_tempVector.at(i);
	}






	//Now we have to do the dirty work of parsing out the function expression
	//so that we can properly get at the variables and functions we need
	this->parsedExpression=originalExpression;

	//First
	for(int f=0; f<n_gfs; f++) {
		string::size_type sPos=parsedExpression.find(gfNames[f]);
		for( ; sPos!=string::npos; sPos=parsedExpression.find(gfNames[f],sPos+1)) {

			string::size_type openPar = parsedExpression.find_first_of('(',sPos);
			string::size_type closePar = parsedExpression.find_first_of(')',sPos);
			if(openPar!=string::npos && closePar!=string::npos) {
				if(closePar>openPar) { //if we got here, we found a valid parenthesis to look at
					string inBetween = parsedExpression.substr(openPar+1,closePar-openPar-1);
					NFutil::trim(inBetween);

					if(inBetween.size()==0) {
						parsedExpression.replace(openPar,closePar-openPar+1,"");
					}
				}
			}
		}
	}

	if (s->getverbose() && (RAZI_DEBUG & CREATE_FUNC)){
		cout<<"Finalizing Composite Functions. The parsed expression is: "<<parsedExpression<<endl;
	}

	///////// do the same for local functions here (can be a bit tricky, because different
	///////// composite functions will have different numbers of arguments
	vector <int> lfIndexValues;
	vector <string> lfReferenceName;
	vector <int> lfScope;


	for(int f=0; f<n_lfs; f++) {
		string::size_type sPos=parsedExpression.find(lfs[f]->getName());
		for( ; sPos!=string::npos; sPos=parsedExpression.find(lfNames[f],sPos+1)) {

			string::size_type openPar = parsedExpression.find_first_of('(',sPos);
			string::size_type closePar = parsedExpression.find_first_of(')',sPos);
			if(openPar!=string::npos && closePar!=string::npos) {
				if(closePar>openPar) { //if we got here, we found a valid parenthesis to look at
					string possibleArg = parsedExpression.substr(openPar+1,closePar-openPar-1);
					NFutil::trim(possibleArg);

					for(unsigned int a=0; a<n_args; a++) {

						if(possibleArg==argNames[a]) {

#ifdef RHS_FUNC
							if (s->getverbose() && (RAZI_DEBUG & CREATE_FUNC)) cout<<" Initializing a function possibleArg:" << possibleArg<< "   argNames[a]:"<< argNames[a]<<endl;
#endif
							string identifier = "_"+argNames[a];
							parsedExpression.replace(openPar,closePar-openPar+1,identifier);

							bool found = false;
							for(unsigned int x=0; x<lfIndexValues.size(); x++) {
								if(lfReferenceName.at(x)==(lfNames[f]+identifier)){
									found=true;
								}
							}
							if(!found) {
								lfIndexValues.push_back(f);
								lfReferenceName.push_back(lfNames[f]+identifier);
#ifdef RHS_FUNC  //Razi: This block is added to support RHS functions
								lfScope.push_back(a);    //change the scope to RHSLocal
#else
								lfScope.push_back(a);
#endif

							}

							break; //break cause we're done with this scope...
						}
					}
				}
			}
		}
	}



	this->n_refLfs=lfIndexValues.size();
	this->refLfInds = new int[n_refLfs];
	this->refLfRefNames = new string[n_refLfs];
	this->refLfScopes = new int[n_refLfs];
	this->refLfValues = new double[n_refLfs];
	for(unsigned int i=0; i<lfIndexValues.size(); i++) {
		this->refLfInds[i]=lfIndexValues.at(i);
		this->refLfRefNames[i]=lfReferenceName.at(i);
		this->refLfScopes[i] = lfScope.at(i);
		this->refLfValues[i]=0;
	}



	// Parse out the ability to get reactant counts in composite reactions

#ifdef RHS_FUNC  //Razi: This block is added to support RHS functions
	if (RHS){return;}   //Razi: This is for regular LHS composite functions
#endif
	//Razi: This is for regular LHS composite functions
	int maxReactantIndex = 0;
	string::size_type sPos=parsedExpression.find("reactant_");
	for( ; sPos!=string::npos; sPos=parsedExpression.find("reactant_",sPos+1)) {
		string::size_type openPar = parsedExpression.find_first_of('(',sPos);
		string::size_type closePar = parsedExpression.find_first_of(')',sPos);
		if(openPar!=string::npos && closePar!=string::npos) {
			if(closePar>openPar) { //if we got here, we found a valid parenthesis to look at
				string inBetween = parsedExpression.substr(openPar+1,closePar-openPar-1);
				NFutil::trim(inBetween);
				if(inBetween.size()==0) {
					parsedExpression.replace(openPar,closePar-openPar+1,"");
				}
			}
		}

		string numOneDigit = parsedExpression.substr(sPos+9,1);
		string numTwoDigits = parsedExpression.substr(sPos+9,2);
		NFutil::trim(numTwoDigits);



		int iOneDigit;
		try {
			iOneDigit = NFutil::convertToInt(numOneDigit);
		} catch (std::runtime_error e) {
			cerr<<"When referencing a reactant, you must include the reactant number"<<endl;
			cerr<<e.what()<<endl;
			exit(1);
		}

		bool isTwoDigitNumber = true;
		if(numTwoDigits.size()<2) {
			isTwoDigitNumber = false;
		} else {
			//cout<<endl<<numTwoDigits<<endl;
			try {
				NFutil::convertToInt(numTwoDigits);
			} catch (std::runtime_error e) {
				isTwoDigitNumber = false;
			}
		}

		if(isTwoDigitNumber) {
			cerr<<"When referencing a reactant, you can only reference reactant numbers up to 9."<<endl;
			exit(1);
		}

		if(iOneDigit>maxReactantIndex) maxReactantIndex = iOneDigit;
	}

//razi: later del
	//cout<<"now the expression is finally: "<<parsedExpression<<endl;	cout<<" max reactant index: "<<maxReactantIndex<<endl;

	this->n_reactantCounts = maxReactantIndex;
	reactantCount = new double[n_reactantCounts];
	for(int r=0; r<n_reactantCounts; r++) {
		reactantCount[r]=0;
	}

	//cout<<"prepared: "<<this->name<<endl;
	//cout<<"expression= "<<this->parsedExpression<<endl;
}

int CompositeFunction::getNumOfArgs() const {
	return this->n_args;
}
string CompositeFunction::getArgName(int aIndex) const {
	return this->argNames[aIndex];
}



void CompositeFunction::updateParameters(System *s)
{
	//cout<<"Updating parameters for function: "<<name<<endl;
	for(unsigned int i=0; i<n_params; i++) {
		p->DefineConst(paramNames[i],s->getParameter(paramNames[i]));
	}
}

void CompositeFunction::prepareForSimulation(System *s)
{
	bool verbose = false;
	if (RAZI_DEBUG & CREATE_FUNC) verbose = s->getverbose();

	try {
		if (verbose){  cout<<"Preparing composite function for simulation. Func:"<< this->name<<endl; mypause(-1);}

		p=FuncFactory::create();
		for(int f=0; f<n_gfs; f++) {
			p->DefineVar(gfNames[f],&gfValues[f]);
			if (verbose)  cout<<"\tSet Global Func:"<< gfNames[f] <<" Value:" << gfValues[f]<<endl;
		}

		//Define local function variables here...
#ifdef RHS_FUNC  //Razi: This block is added to support RHS functions
		if(RHS){
			for(int f=0; f<this->n_refLfs; f++) {
				p->DefineVar(refLfRefNames[f],0);
				if (verbose)  cout<<"\tSet Local Func:"<<refLfRefNames[f] <<" Value:" << 0<<endl;
			}
		}else
#endif
		//Razi: Applies to LHS functions
		for(int f=0; f<this->n_refLfs; f++) {
			p->DefineVar(refLfRefNames[f],&refLfValues[f]);
			if (verbose)  cout<<"\tSet Local Func:"<<refLfRefNames[f] <<" Value:" << refLfValues[f]<<endl;
		}

		for(unsigned int i=0; i<n_params; i++) {
			p->DefineConst(paramNames[i],s->getParameter(paramNames[i]));
			if (verbose)  cout<<"\tParam:"<<paramNames[i] <<" Value:" << s->getParameter(paramNames[i]) <<endl;
		}

#ifdef RHS_FUNC  //Razi: This block is added to support RHS functions
		if(!RHS)  //run only for LHS functions
#endif
		for(int r=0; r<n_reactantCounts; r++) {
			string reactantStr = "reactant_"+NFutil::toString((r+1));
			p->DefineVar(reactantStr,&reactantCount[r]);
			if (verbose)  cout<<"\tReactant:"<<reactantStr <<" Count:" << reactantCount[r] <<endl;
		}

		p->SetExpr(this->parsedExpression);
		if(verbose){ cout<<"\tExpression:"<<this->parsedExpression<<endl;}
	}
	catch (mu::Parser::exception_type &e)
	{
		cout<<"Error preparing function "<<name<<" in class CompositeFunction!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		exit(1);
	}

	//cout<<"preparing composite function.."<<this->name<<endl;
//	exit(0);
}


void CompositeFunction::printDetails(System *s) {

#ifdef RHS_FUNC
	if(RHS) cout<<" RHS Composite : '"<< this->name << "()'"<<endl; else  cout<<" LHS Composite: '"<< this->name << "()'"<<endl;
#else
	cout<<"Composite Function: '"<< this->name << "()'"<<endl;
#endif
	cout<<" = "<<this->originalExpression<<endl;
	cout<<" parsed expression = "<<this->parsedExpression<<endl;


	cout<<"   -Function References:"<<endl;
	for(int f=0; f<n_gfs; f++) {
		gfValues[f]=FuncFactory::Eval(gfs[f]->p);
		cout<<"         global function: "<<gfNames[f]<<" = "<<gfValues[f]<<endl;

		gfs[f]->printDetails(s);
	}
	for(int f=0; f<n_lfs; f++) {
		cout<<"         local function: "<<lfs[f]->getNiceName()<<endl;
	}


	if(n_args>0) {
		cout<<"   -Arguments:"<<endl;
		for(unsigned int a=0; a<n_args; a++)
			cout<<"         "<<argNames[a]<<endl;
	}

	if(n_params>0) {
		cout<<"   -Constant Parameters:"<<endl;
		for(unsigned int i=0; i<n_params; i++) {
			cout<<"         "<<paramNames[i]<<" = " << s->getParameter(paramNames[i])<<endl;
		}
	}

	if(p!=0)
		cout<<"   Function last evaluated to: "<<FuncFactory::Eval(p)<<endl;

#ifdef RHS_FUNC //Razi added to support RHS function
	else
		cout<<"   The Function has evaluated yet"<<endl;
#endif


	if (!(RAZI_DEBUG & (CREATE_FUNC | SHOW_FIRE))) return;

	if(n_lfs>0) {
		for(int i=0; i<n_refLfs; i++) {
			//cout<<" embedded local fun:"<<lfs[refLfInds[i]]->getNiceName()<<" with scope: "<<refLfScopes[i]<< endl;
			cout<<" Composite Func: embedded local fun:"<<lfs[refLfInds[i]]->getNiceName()<<" with scope: "<<refLfScopes[refLfInds[i]] <<endl;
		}
	}
	cout<< "Number of reactants is: "<< this->n_reactantCounts <<endl;
//	cout<<"trying something new..."<<endl;
//	Molecule ** molList = new Molecule *[2];
//	molList[0] = s->getMoleculeTypeByName("Receptor")->getMolecule(0);
//	molList[1] = s->getMoleculeTypeByName("Receptor")->getMolecule(1);
//	int *scope = new int[1];
//	scope[0]=0;
//	scope[1]=1;
//
//	double x = this->evaluateOn(molList,scope);
//	cout<<"got final value: "<<x<<endl;
//
//	exit(1);
}



void CompositeFunction::addTypeIMoleculeDependency(MoleculeType *mt) {
#ifdef RHS_FUNC //Razi added to support RHS function
	if (RHS) return;  //irrelevant for RHS function
#endif
	for(int i=0; i<n_lfs; i++) {
		// add typeI dependency, which means this local function influences
		//  the propensity of some DOR reaction for which mt is the head template molecule.
		lfs[i]->addTypeIMoleculeDependency(mt);
		if ( refLfScopes[i]==LocalFunction::SPECIES ) {
			// enable complex-scoped evaluation for this local fcn!
			lfs[i]->setEvaluateComplexScope( true );
		}
	}
}


double CompositeFunction::evaluateOn(Molecule **molList, int *scope, int *curReactantCounts, int n_reactants) {
bool verbose=false;
#ifdef RHS_FUNC //Razi added to support RHS function
	if (RHS){ cerr<<" Error,  calling wrong function evaluateOn for RHS function..."<<endl; return -1;}  //irrelevant for RHS function
	//if (RHS && (RAZI_DEBUG & SHOW_FIRE)) verbose=true;
#endif
	if (verbose) cout << "CompositeFunction::evaluateOn() for "<<name<<endl;

	//1 evaluate all global functions
	if (verbose) cout << "n_gfs=" << n_gfs << "     ";
	for(int f=0; f<n_gfs; f++) {
		gfValues[f]=FuncFactory::Eval(gfs[f]->p);
	}

	//2 evaluate all local functions
	if (verbose) cout << "n_lfs=" << n_lfs << "     ";
	if(n_lfs>0) {
		if (verbose) cout << "scope[0]=" << scope[0] << "     ";
		if (verbose) cout << "molList[0]" << molList[0]->getMoleculeTypeName() << "     ";
		if(molList!=0 && scope!=0) {

			if (verbose) cout << "n_refLfs=" << n_refLfs << endl;
			for(int i=0; i<n_refLfs; i++) {
				if (verbose) cout<<"--- evaluating: "<<lfs[refLfInds[i]]->getNiceName()<<" with scope: "<<scope[refLfScopes[i]]<< "     ";
				this->refLfValues[i] = this->lfs[refLfInds[i]]->getValue(molList[refLfScopes[i]],scope[refLfScopes[i]]);
				if (verbose) cout<<"answer: "<<this->refLfValues[i]<<"     ";
			}

			//for (n_refLfs)  set the value by calling the correct local function to evaluate on the specified scope
			//which we reference through the given molList.
			//this->refLfValues[i] = this->lfs[refLfInds[i]]->evaluateOn(molList[refLfScopes[i]],scope[refLfScopes[i]]);

		} else {

			cout<<"Error evaluating composite function: "<<name<<endl;
			if (molList==0) cout<<"This function depends on local functions, but you gave no molecules when calling this function.  Time to quit."<<endl;
			if (scope==0) cout<<"This function depends on local functions, but you gave no scope when calling this function.  Time to quit."<<endl;
			exit(1);
		}

	}
	if (verbose) cout<<endl;

	//3 update reactant counts as need be
	if(n_reactants<this->n_reactantCounts) {
		cerr<<"Not given enough reactants for this composite function!"<<this->name<<endl;
	}
	for(int r=0; r<n_reactantCounts; r++) {
		reactantCount[r]=curReactantCounts[r];
	}
	//evaluate this function
	return FuncFactory::Eval(p);

}



#ifdef RHS_FUNC //Razi added to support RHS function

void CompositeFunction::setRHSFlag(bool val) {
	RHS=val;

	cout<<"\n\nClean up the unnecessary links (e.g. type I, type II molecule type) ....\n\n\n"; mypause(-1);

}; //set to true for RHS, default: false


//Razi: directly apply composite function on a test molecule [test product of a reaction]
double CompositeFunction::evaluateOnProduct(Molecule *mol, int scope, int evaluationType, bool verbose) {

	if (verbose) {cout << "RHS CompositeFunction::evaluateOn() for "<<name<<endl;}// mypause(-1);};


	if ((evaluationType!=EvalGlobalPart) && (evaluationType!=EvalConditionalPart) && (evaluationType!=EvalBothParts)){
		cout<<"Invalid evaluation type when evaluating RHS function.\n"; exit(1);
	}


	//1 evaluate all global functions
	if (verbose) cout << "n_gfs=" << n_gfs << endl;
	for(int f=0; f<n_gfs; f++) {
		if (evaluationType==EvalConditionalPart)
			gfValues[f]=1; //only evaluate the conditional part, so all global parts is set to 1
		else
			gfValues[f]=FuncFactory::Eval(gfs[f]->p);
	}

	//2 evaluate all local functions
	if (verbose) cout << "n_lfs=" << n_lfs << endl;
	if(n_lfs>0) {
		if (verbose) cout << "scope[0]=" << scope << endl;
		//if (verbose) cout << "molList[0]" << molList[0]->getMoleculeTypeName() << endl;
		if (verbose) cout<<"evaluating composite function with local dependencies."<<endl;

		if (evaluationType==EvalGlobalPart){
			for(int i=0; i<n_refLfs; i++)
				this->refLfValues[i] = 1;   //just evaluate global parts and set all local functions to 1

		}else{

			if(mol!=0 && ((scope == LocalFunction::SPECIES)|| (mol!=0 && scope == LocalFunction::MOLECULE)))
			{
				if (verbose) cout << "n_refLfs=" << n_refLfs << endl;
				for(int i=0; i<n_refLfs; i++) {
					//if ((RAZI_DEBUG) & SHOW_FIRE) cout<<"--- evaluating: "<<lfs[refLfInds[i]]->getNiceName()<<" with scope: "<<scope[refLfScopes[i]]<<endl;
					if (verbose) cout<<"Local func i:"<< i << "  refLfInds[i]:"<<refLfInds[i]<<"   "<< lfs[refLfInds[i]]->getName();


					//this->refLfValues[i] = this->lfs[refLfInds[i]]->getValue(mol,scope);
					this->refLfValues[i] = this->lfs[refLfInds[i]]->evaluateOnTestMol(mol,scope, verbose);
					if (verbose) cout<<"  answer: "<<this->refLfValues[i]<<endl;
				}
				//for (n_refLfs)  set the value by calling the correct local function to evaluate on the specified scope
				//which we reference through the given molList.
				//this->refLfValues[i] = this->lfs[refLfInds[i]]->evaluateOn(molList[refLfScopes[i]],scope[refLfScopes[i]]);

			}else {
				if (mol==0) cout<<"Error evaluating composite function: "<<name<<", No molecule is specified."<<endl;
				else cout<<"This function depends on local functions, but you gave no scope when calling this function.  Time to quit."<<endl;
				exit(1);
			}
		}
	}

	//evaluate this function
	double result = FuncFactory::Eval(p);
	if (verbose) cout<<"Final result when evaluating a RHS function is: "<< result<<endl;
	return result;
}

#endif
