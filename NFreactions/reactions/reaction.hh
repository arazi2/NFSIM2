#ifndef BASICREACTIONS_HH_
#define BASICREACTIONS_HH_


#include "../NFreactions.hh"




using namespace std;

namespace NFcore
{


#ifdef RHS_FUNC //Razi added to support RHS functions


	//Razi: This class includes reactions with local functions applied to the right hand side
    //The function is assumed to be a composite fucntion, where the global part defines the rate and the local function typically
   // is a condition on the output products templates.
	// Later, we may need to develop a class that supports functions applied to both input reactants and output products
	class RHSRxnClass : public ReactionClass {
		public:
		RHSRxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					CompositeFunction *function,
					vector <string> &lfArgumentPointerNameList,
					System *s);
			virtual ~RHSRxnClass();


			virtual void init();
			virtual void prepareForSimulation() {};



			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();

			//Razi: Check if needed: virtual int getDORreactantPosition() const { return DORreactantIndex; };

			//JJT: checks if there's an existing mapping set in *m equal to *ms that maps to this reaction
			//virtual int checkForCollision(Molecule *m, MappingSet* ms,int rxnIndex);
			//JJT: checks if there's an existing mapping set in *m equal to *ms that maps to this reaction
			//Razi: Check if needed [comes from basic reaction class]:
			virtual int checkForEquality(Molecule *m, MappingSet* ms,int rxnIndex, ReactantList*);



			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;

			virtual void printDetails() const;
			virtual void printFullDetails() const {};

			//void directAddForDebugging(Molecule *m);
			//void printTreeForDebugging();

			//static void test1(System *s);

			/*moved to reactionClass
			//Razi added to handle RHS functions
			unsigned int n_productTemplates;
			virtual bool checkReaction();   //clone reactants
			MappingSet ** check_mappingSet; //include pointers to test molecules
			*/

		protected:

			virtual double evaluateLocalFunctions(MappingSet *ms);

			virtual void pickMappingSets(double randNumber) const;

			ReactantList **reactantLists;
			ReactantList **productLists;

//razi: moved to reactionClass			MappingSet *ms;
//razi: moved to reactionClass			CompositeFunction *cfo;
			int n_argMolecules;
			int * argIndexIntoMappingSet;
			Molecule ** argMappedMolecule;
			int * argScope;


//			bool includeRHSFunc;

			//Parameters to keep track of local functions
//razi moved to reactionClass			int DORproductIndex;
			ReactantList *rl;  //Razi: From Basic Rxn


			//vector <int> argIndexIntoMappingSet;

			//vector <LocalFunction *> lfList;
			//vector <int> indexIntoMappingSet;
			//vector <double> localFunctionValue;

	};

#endif //RHS_FUNC


	class BasicRxnClass : public ReactionClass {
		public:
			BasicRxnClass(string name, double baseRate, string baseRateName, TransformationSet *transformationSet, System *s);
			virtual ~BasicRxnClass();

			virtual void init();
			virtual void prepareForSimulation();

			//JJT: checks if there's an existing mapping set in *m equal to *ms that maps to this reaction
			virtual int checkForEquality(Molecule *m, MappingSet* ms,int rxnIndex, ReactantList*);

			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();
			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;

			virtual void printFullDetails() const;

		protected:
			virtual void pickMappingSets(double randNumber) const;

			ReactantList **reactantLists;

//#ifdef RHS_FUNC //Razi added to support RHS functions, 			//check if needed
//			ReactantList ** productLists;
//#endif
			ReactantList *rl;
			MappingSet *ms;
	};


	class FunctionalRxnClass : public BasicRxnClass {

		public:
			FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet, System *s);
			FunctionalRxnClass(string name, CompositeFunction *cf, TransformationSet *transformationSet, System *s);

			virtual ~FunctionalRxnClass();

			virtual double update_a();
			virtual void printDetails() const;

		protected:
			GlobalFunction *gf;
			CompositeFunction *cf;
	};

	class MMRxnClass : public BasicRxnClass {

		public:
			MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet, System *s);
			virtual ~MMRxnClass();

			virtual double update_a();
			virtual void printDetails() const;

		protected:
			double Km;
			double kcat;
			double sFree;
	};


	class DORRxnClass : public ReactionClass {
		public:
			DORRxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					CompositeFunction *function,
					vector <string> &lfArgumentPointerNameList,
					System *s);
			virtual ~DORRxnClass();

			virtual void init();
			virtual void prepareForSimulation() {};
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();

			virtual int getDORreactantPosition() const { return DORreactantIndex; };

			//JJT: checks if there's an existing mapping set in *m equal to *ms that maps to this reaction
			virtual int checkForCollision(Molecule *m, MappingSet* ms,int rxnIndex);

			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;

			virtual void printDetails() const;
			virtual void printFullDetails() const {};

			void directAddForDebugging(Molecule *m);
			void printTreeForDebugging();

			static void test1(System *s);

		protected:

			virtual double evaluateLocalFunctions(MappingSet *ms);

			virtual void pickMappingSets(double randNumber) const;

			ReactantList **reactantLists;
			ReactantTree *reactantTree;

			MappingSet *ms;


			CompositeFunction *cf;

/*
#ifdef RHS_FUNC //Razi added to support RHS functions
			//check if needed
			ReactantList **productLists;
			CompositeFunction *cfo;
			bool includeRHSFunc;
#endif
*/

			//Parameters to keep track of local functions
			int DORreactantIndex;

			int n_argMolecules;
			int * argIndexIntoMappingSet;
			Molecule ** argMappedMolecule;
			int * argScope;


			//vector <int> argIndexIntoMappingSet;



			//vector <LocalFunction *> lfList;
			//vector <int> indexIntoMappingSet;
			//vector <double> localFunctionValue;

	};

	/* A reaction class with DOR calculations on two reactants.
	 * The rate function must be factored as  f(x,y) = g(x)*h(y)
	 * */
	class DOR2RxnClass : public ReactionClass {
		public:
			DOR2RxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					CompositeFunction *function1,
					CompositeFunction *function2,
					vector <string> &lfArgumentPointerNameList1,
					vector <string> &lfArgumentPointerNameList2,
					System *s);
			virtual ~DOR2RxnClass();

			virtual void init();
			virtual void prepareForSimulation() {};
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();

			virtual int getDORreactantPosition()  const { return DORreactantIndex1; };
			virtual int getDORreactantPosition2() const { return DORreactantIndex2; };

			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;

			virtual void printDetails() const;
			virtual void printFullDetails() const {};

			void directAddForDebugging(Molecule *m);
			void printTreeForDebugging();

			static void test1(System *s);

		protected:

			virtual double evaluateLocalFunctions1(MappingSet *ms);
			virtual double evaluateLocalFunctions2(MappingSet *ms);

			virtual void pickMappingSets(double randNumber) const;

			ReactantList **reactantLists;
			ReactantTree *reactantTree1;
			ReactantTree *reactantTree2;

			// TODO: figure out if we need one of two of these
			MappingSet *ms;

			CompositeFunction *cf1;
			CompositeFunction *cf2;

/*
#ifdef RHS_FUNC //Razi added to support RHS functions
			//check if needed
			ReactantList **productLists;
			CompositeFunction *cfo1;
			CompositeFunction *cfo2;
#endif
*/

			//Parameters to keep track of local functions
			int DORreactantIndex1;
			int DORreactantIndex2;

			int n_argMolecules1;
			int n_argMolecules2;
			int * argIndexIntoMappingSet1;
			int * argIndexIntoMappingSet2;
			Molecule ** argMappedMolecule1;
			Molecule ** argMappedMolecule2;
			int * argScope1;
			int * argScope2;

	};

}




#endif /*BASICREACTIONS_HH_*/
