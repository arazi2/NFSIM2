/*
 * razi.hh
 *
 *  Created on: Jun 19, 2016
 *      Author: Abolfazl Razi
 */
#include <iostream>
#include <ctime>


#ifndef SETTING_HH_
#define SETTING_HH_




#ifdef _WIN32
#include <windows.h>
#elif __linux__
#elif __unix__
#elif __APPLE__
#endif


//namespace RZ{


//selecting extra debugging messages to show
#define CREATE_MOLECULE		0x0001
#define CREATE_SPECIES		0x0002
#define CREATE_FUNC			0x0004
#define CREATE_TEMPLATES	0x0008
#define CREATE_REACTION		0x0010
#define CREATE_OBS			0x0020
#define RUN_REACTIONS		0x0040
#define READ_FILE 			0x0080
#define SHOW_ARGS 			0x0100
#define SHOW_SIM			0x0200
#define SHOW_FIRE			0x0400
#define CHECK_LBG 			0x0800  //activate debugging message for local function bug
#define SHOW_PARAMS			0x1000
#define SHOW_RHS			0x2000

//#define RAZI_DEBUG    		0x0FFF  //activate everything
//#define RAZI_DEBUG 0//(READ_FILE|CHECK_LBG) //SHOW_FIRE//(READ_FILE|SHOW_SIM)//0xFFFF //(READ_FILE |CREATE_REACTION | CREATE_FUNC | CREATE_OBS)  //(READ_FILE | CREATE_FUNC | CREATE_OBS)
//#define RAZI_DEBUG CREATE_TEMPLATES | CREATE_MOLECULE | CREATE_TEMPLATES | CREATE_REACTION //SHOW_FIRE//(READ_FILE|SHOW_SIM)//0xFFFF //(READ_FILE |CREATE_REACTION | CREATE_FUNC | CREATE_OBS)  //(READ_FILE | CREATE_FUNC | CREATE_OBS)
#define RAZI_DEBUG 			(SHOW_FIRE | SHOW_RHS | SHOW_SIM) //| CREATE_REACTION| CREATE_FUNC


#define TEST_REAC   "XWBindCond1"//"XYBind"//"Bug1"//"XaYbBind" //Razi:Show more messages for reactions: Use "ALL" to show all reactions
#define RAZI_VER   1 //Razi: activates some debugging messages
#define FIX_BUG1 1   //Razi: this applies some changes to fix the local function bug
//#define FIX_OBS_COUNT //Razi: Activate later after checking the results. This macro corrects the problem of counting a molecule twice against an observable: EX: X(P~Phos) for OBS: X(p~Phos,w!1).W(a!1,b!2).X(p~Phos,w!2)

#define RHS_FUNC 1  //Razi: Activates function evaluation in the output
#define RHS_maxRETRY  1  //define number of retries for reactions with RHS functions


//#define FIX_A     1 //Razi: Change propensity update for reactions, and calculating reaction fire time [dt]
/* Razi: Later check, but it seems that there is a problem in calculating propensities. In the original version the propensity of a reaction is calculated
 * as base rate * n1*n2*..nk, where ni is the number of  the ith  reactant
 * This is problematic and perhaps a better option would be baserate*min(n1,n2,....n_k)
 * or baserate*(n1,n2,....n_k)^(1/k)
 * consider two reactions X(p~0)->X(p~1) 0.1
 * and X(y)+Y(x)->X(y!1).Y(x!1)   0.1
 * if the number of X and Y species is 100, the propensities of two reactions with the current calculations are
 * p1=0.2*100=50 and p2=0.2*100*100=5000>100, which is more than the available reactants X(y) and Y(x).
 * Obviously it is not correct for the second reaction [one may expect a lower rate for the second reaction or a rate < [100*.2]]
 *  */


//#ifdef RHS_FUNC  //Razi added this part to support RHS functions


#define WTIME  100// 5000   //waiting time in ms used for mypause function throughout the code, 0 results in waiting for a key press
void mypause();
void mypause(float ms);
#ifdef _WIN32
void cls();
#endif
//}//end namespace



#endif /* SETTING_HH_ */
