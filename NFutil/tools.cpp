/*
 * Tools.cpp
 *
 *  Created on: Jun 19, 2016
 *      Author: Abolfazl Razi
 */
#include <iostream>
#include <string>
#include <time.h>
#include "setting.hh"

#ifdef _WIN32
#include <windows.h>
#elif __linux__
#elif __unix__
#elif __APPLE__
#endif


using namespace std;

/*
 void mypause(){
	//std::cin.get();
	 mypause(WTIME);
}
*/


//razi added this function for test purpose, last update 2017-3-29
void mypause(float ms = WTIME){
	//return;

	int cnt=0;
	const char ch[] = "/-\\-";

//	cout<<" pause for "<<ms<<"ms ";
	if (ms<=0){
		//std::cin.ignore();
		std::cin.get();
	}else{
		ms=ms/1000;

		std::cout.flush();
		time_t start, current, prev;
		time(&start); prev=start;
		bool Exitflag=false;
		do{
#ifdef _WIN32
			if (WaitForSingleObject(GetStdHandle( STD_INPUT_HANDLE ), 50) == WAIT_OBJECT_0) {   //wait for 50 ms and check if key is pressed
				Exitflag = true;   //exit
			}
#else
			if (false){}
#endif
			else{   //check timeout
				time(&current);
				if ((current - prev)>0.2){
					prev=current; cnt++; std::cout<<"\b"<<ch[cnt%4];
				}
				if ((current - start) >= ms)Exitflag=true;
			}
		}
		while (!Exitflag);
	}
//	cout<<"Done!  ";
}

#ifdef _WIN32
//razi added this function for test purpose, last update 2017-3-29
void cls()
{
   COORD coordScreen = { 0, 0 };    // home for the cursor
   DWORD cCharsWritten;
   CONSOLE_SCREEN_BUFFER_INFO csbi;
   DWORD dwConSize;

   HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);

   // Get the number of character cells in the current buffer.


   if( !GetConsoleScreenBufferInfo( hStdout, &csbi ))
      return;
   dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

   // Fill the entire screen with blanks.

   if( !FillConsoleOutputCharacter( hStdout, (TCHAR) ' ',
      dwConSize, coordScreen, &cCharsWritten ))
      return;

   // Get the current text attribute.

   if( !GetConsoleScreenBufferInfo( hStdout, &csbi ))
      return;

   // Set the buffer's attributes accordingly.

   if( !FillConsoleOutputAttribute( hStdout, csbi.wAttributes,
      dwConSize, coordScreen, &cCharsWritten ))
      return;

   // Put the cursor at its home coordinates.

   SetConsoleCursorPosition( hStdout, coordScreen );
}
#endif
