/*
 Copyright (C) 2003-2013 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

// File version: 2013-04-01

/*
	tandem.cpp provides a command line interface to the main processing class, mprocess.
	A single command line parameter is accepted, which is the file name for an XML class
	containing the parameters for performing the protein modelling session. The input file 
	also contains the name of an output file, which will contain the output from the session.
*/

#include "stdafx.h"

#ifndef X_P3
#include "tandem.h"

// For defining R_CStackLimit
#ifndef _WIN32
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif



//int main(int argc, char* argv[]) // rTANDEM
// .Call("tandem",RTsexp['param'], RTsexp['peptide'], RTsexp['saps'], RTsexp['mods'], RTsexp['spectrum'])
SEXP tandem(SEXP param, SEXP peptide, SEXP saps, SEXP mods, SEXP spectrum) // rTANDEM
{
    // Disable R c-stack
#ifndef _WIN32
  R_CStackLimit = (uintptr_t) -1;
#endif

	/*
	* Check the argv array for at least one parameter.
	* mprocess checks the validity of the file.
	*/
// rTANDEM : since we don't use argc and argv, usage was disabled
//	if(argc < 2 || argc > 1 && strstr(argv[1],"-L") == argv[1] || argc > 1 && strstr(argv[1],"-h") == argv[1])	{ 
//		cout << "\n\nUSAGE: tandem filename\n\nwhere filename is any valid path to an XML input file.\n\n+-+-+-+-+-+-+\n"; 
// 		cout << "\nX! TANDEM " << VERSION << "\n"; 
//		cout << "\nCopyright (C) 2003-2013 Ronald C Beavis, all rights reserved\n"; 
// 		cout << "This software is a component of the GPM  project.\n"; 
//		cout << "Use of this software governed by the Artistic license.\n"; 
//		cout << "If you do not have this license, you can get a copy at\n"; 
//		cout << "http://www.perl.com/pub/a/language/misc/Artistic.html\n"; 
//		cout << "\n+-+-+-+-+-+-+\n\npress <Enter> to continue ..."; 
//		char *pValue = new char[128]; 
//		cin.getline(pValue,127); 
//		delete pValue; 
//		return -1; 
//	} 
//	cout << "\nX! TANDEM " << VERSION << "\n\n"; 
//

	/*
	* Create an mprocess object array
	*/
	unsigned long lMaxThreads = 256;
	mprocess **pProcess = new mprocess*[lMaxThreads];
	if(pProcess == NULL)	{
//		cout << "An error was detected creating the processing objects.\nPlease contact a GPM administrator.\n";
		Rprintf("An error was detected creating the processing objects.\nPlease contact a GPM administrator.\n");
//		return -2; // rTANDEM
		return R_NilValue; // rTANDEM
	}
#ifdef MSVC
	DWORD *pId = new DWORD[lMaxThreads];
	HANDLE *pHandle = new HANDLE[lMaxThreads];
#else
	int *pId = new int[lMaxThreads];
	int *pHandle = new int[lMaxThreads];
	pthread_t pThreads[lMaxThreads];

#endif
	unsigned long a = 0;
	while(a < lMaxThreads)	{
		pProcess[a] = NULL;
// rTANDEM : There was a warning here during package installation.
//		pHandle[a] = NULL; // rTANDEM
//		pId[a] = NULL; // rTANDEM
		pHandle[a] = 0; // rTANDEM
		pId[a] = 0; // rTANDEM
		a++;
	}
	pProcess[0] = new mprocess;
//	cout << "Loading spectra";
	Rprintf("Loading spectra\n");
	//cout.flush();
	/*
	* Initialize the first mprocess object with the input file name.
	*/
	char *pS = new char[1024];
// rTANDEM: We need to change the code so data is no longer loader from file
//	strcpy(pS,argv[1]); // rTANDEM
//	if(!pProcess[0]->load(pS))	{
	if(!pProcess[0]->load(param, peptide, saps, mods, spectrum));
//		cout << "\n\nAn error was detected while loading the input parameters.\nPlease follow the advice above or contact a GPM administrator to help you."; // rTANDEM
//		delete pProcess[0]; // rTANDEM
//		delete pProcess; // rTANDEM
//		return -4; // rTANDEM
//	} // rTANDEM
//	cout << " loaded.\n";
	Rprintf(" loaded.\n");
	if(pProcess[0]->m_vSpectra.size() == 0)	{
//		cout << "No input spectra met the acceptance criteria.\n";
		Rprintf("No input spectra met the acceptance criteria.\n");
		//cout.flush();
		//delete pProcess[0]; //rTANDEM: Cpp memory management tend to crash the R session. Better leave R do the gc.
		//delete pProcess;    //rTANDEM: Cpp memory management tend to crash the R session. Better leave R do the gc.
//		return 1; // rTANDEM
		return R_NilValue; // rTANDEM
	}
	pProcess[0]->serialize();
//	cout << "Spectra matching criteria = " << (unsigned long)pProcess[0]->m_vSpectra.size() << "\n";
	Rprintf("Spectra matching criteria = %lu\n", (unsigned long)pProcess[0]->m_vSpectra.size());
	//cout.flush();
#ifdef PLUGGABLE_SCORING
// 	cout << "Pluggable scoring enabled.\n";
 	Rprintf("Pluggable scoring enabled.\n");
#endif
      /*

        * Start the mprocess object and wait for it to return.

        */
	unsigned long lThread =	pProcess[0]->get_thread();
	unsigned long lThreads = pProcess[0]->get_threads();
	if(lThreads	> lMaxThreads)	{
		lThreads = lMaxThreads;
	}
	if(pProcess[0]->m_vSpectra.size() <	lThreads)	{
		lThreads = (unsigned long)pProcess[0]->m_vSpectra.size();
		if(lThreads	< 1)		{
			lThreads = 1;
		}
		pProcess[0]->set_threads(lThreads);
	}
#ifdef MSVC
	DWORD dCount = lThreads	- 1;
#else
	int	dCount = lThreads -	1;
#endif
	long lSpectra =	lThreads + (long)pProcess[0]->m_vSpectra.size()/lThreads;
	bool bSpectra =	true;
//	cout <<	"Starting threads .";
	Rprintf("Starting threads .");
	//cout.flush();
    size_t tCount = pProcess[0]->m_vSpectra.size();
	sort(pProcess[0]->m_vSpectra.begin(),pProcess[0]->m_vSpectra.end(),lessThanSpec);
	float fMax = (float)pProcess[0]->m_vSpectra.back().m_dMH;
	float fZ = pProcess[0]->m_vSpectra.back().m_fZ;
	pProcess[0]->m_fMaxMass = fMax;
	pProcess[0]->m_fMaxZ = fZ;
	if(lThread != 0xFFFFFFFF)		{
		while(dCount > 0)		{
			pProcess[dCount] = new mprocess;
			pProcess[dCount]->set_thread(dCount);						 
			/*

			* initialize the new mprocess objects with	the	spectra	already	loaded into	the	first mprocess

			*/
			pProcess[dCount]->m_vSpectra.reserve(lSpectra);
			pProcess[dCount]->m_fMaxMass = fMax;
			pProcess[dCount]->m_fMaxZ = fZ;
			dCount--;
		}

		size_t tProcesses = lThreads;
		size_t tRing = 0;
		vector<mspectrum> vZero;
		vZero.reserve(lSpectra);
		do	{
			if(tRing == 0)	{
				vZero.push_back(pProcess[0]->m_vSpectra.back());
			}
			else	{
				pProcess[tRing]->m_vSpectra.push_back(pProcess[0]->m_vSpectra.back());
			}
			tRing++;
			pProcess[0]->m_vSpectra.pop_back();
			if(tRing == tProcesses)	{
				tRing = 0;

			}
		}	while(pProcess[0]->m_vSpectra.size() != 0);
		pProcess[0]->m_vSpectra.reserve(vZero.size());
		do	{
			pProcess[0]->m_vSpectra.push_back(vZero.back());
			vZero.pop_back();
		}	while(vZero.size() != 0);
		dCount = lThreads - 1;
		while(dCount > 0)		{
		  //if(!pProcess[dCount]->load(pS,pProcess[0]))	{
		  if(!pProcess[dCount]->load(param, peptide, saps, mods, spectrum)){

//				cout <<	"error pProcess->LoadParameters	returned error (main)\r\n";
				Rprintf("error pProcess->LoadParameters	returned error (main)\r\n");
				delete pProcess;
//				return -4; // rTANDEM
				return R_NilValue; // rTANDEM
			}
			dCount--;
//			cout <<	".";
			Rprintf(".");
			//cout.flush();
		}
	}
	delete pS;
	dCount = 0;
#ifdef MSVC
	pHandle[dCount] = CreateThread(NULL,0,ProcessThread,(void *)pProcess[dCount],0,&pId[dCount]);
#else
	pthread_create(&pThreads[dCount],NULL,ProcessThread,(void*)pProcess[dCount]);
#endif
	dCount++;
	/*
	* Initialize more mprocess objects, if lThread is not 0xFFFFFFFF, which signifies default single
	* threaded operation.
	*/
	if(lThread != 0xFFFFFFFF && bSpectra)	{
		while((unsigned long)(dCount) < lThreads)	{
#ifdef MSVC
			pHandle[dCount] = CreateThread(NULL,0,ProcessThread,(void *)pProcess[dCount],0,&pId[dCount]);
#else
			pthread_create(&pThreads[dCount],NULL,ProcessThread,(void*)pProcess[dCount]);
#endif
			dCount++;
		}
	}
//	cout << " started.\n";
	Rprintf(" started.\n");
	//cout.flush();
//	cout << "Computing models:\n";
	Rprintf("Computing models:\n");
	//cout.flush();
	/*
	* wait until all of the mprocess objects return.
	*/
#ifdef MSVC
//	DWORD wait = WaitForMultipleObjects(dCount,pHandle,true,INFINITE);
	a = 0;
	DWORD dwTime = 100000;
	DWORD wait = WAIT_TIMEOUT;
	int iTics = 0;
	while(a < (unsigned long)(dCount))	{
		wait = WaitForSingleObject(pHandle[a],100);
		if(a > 0 && wait == WAIT_TIMEOUT)	{
			if(a == 1)	{
//				cout << "waiting for " << a+1;
				Rprintf("waiting for %lu", a+1);
			}
			else	{
//				cout << a+1;
				Rprintf("%lu",a+1);
			}
			while(wait == WAIT_TIMEOUT)	{
				wait = WaitForSingleObject(pHandle[a],dwTime);
				if(wait == WAIT_TIMEOUT)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					iTics++;
					if(iTics > 50)	{
//						cout << "|\n\t\t";
						Rprintf("|\n\t\t");
						//cout.flush();
						iTics = 0;
					}
				}
			}
		}
		else	{
			while(wait == WAIT_TIMEOUT)	{
				wait = WaitForSingleObject(pHandle[a],dwTime);
				if(wait == WAIT_TIMEOUT)	{
//					cout << ":";
					Rprintf(":");
					//cout.flush();
				}
			}
			if(a == 1)	{
//				cout << "waiting for " << a+1;
				Rprintf("waiting for %lu", a+1);
			}
			else if(a == 0)	{
//				cout << "\n\t";
				Rprintf("\n\t");
				//cout.flush();
			}
			else	{
//				cout << a+1;
				Rprintf("%lu", a+1);
			}
		}
		a++;
	}
	if(dCount > 1)	{
//			cout << " done.\n\n";
			Rprintf(" done.\n\n");
			//cout.flush();
	}
	else	{
//		cout << "\n";
		Rprintf("\n");
		//cout.flush();
	}
#else
	//2003-03-01:note - the declaration below was changed from void **vp;	
	void *vp;
	int x=0;
	int wait;
	for(x=0;x<dCount;x++){
		//2003-03-01:note - the 2nd parameter in the call to pthread_join() was changed from vp
		wait = pthread_join(pThreads[x],&vp);
		wait++; /* fool the compiler, instead, there should be a check wait == 0 */
	}
#endif
//	cout << "\tsequences modelled = "<< (long)(pProcess[0]->get_protein_count()/1000.0 + 0.5) << " ks\n";
	Rprintf("\tsequences modelled = %ld ks\n", (long)(pProcess[0]->get_protein_count()/1000.0 + 0.5));
	//cout.flush();
	pProcess[0]->merge_spectra();
	a = 1;
	/*
	* merge the results into the first object
	*/
	while(a < (unsigned long)(dCount))	{
		pProcess[0]->merge_map(pProcess[a]->m_mapSequences);
		pProcess[0]->merge_spectra(pProcess[a]->m_vSpectra);
		a++;
	}
	a = 1;
	pProcess[0]->load_sequences();
	while(a < (unsigned long)(dCount))	{
		pProcess[a]->merge_map(pProcess[0]->m_mapSequences);
		pProcess[a]->m_vseqBest = pProcess[0]->m_vseqBest;
		pProcess[a]->m_mapCrc.clear();
		a++;
	}

	/*
	* Report the contents of the mprocess objects into an XML file as described
	* in the input file.
	*/
//	cout << "Model refinement:\n";
	Rprintf("Model refinement:\n");
	//cout.flush();
	dCount = 0;
#ifdef MSVC
	pHandle[dCount] = CreateThread(NULL,0,RefineThread,(void *)pProcess[dCount],0,&pId[dCount]);
#else
	pthread_create(&pThreads[dCount],NULL,RefineThread,(void*)pProcess[dCount]);
#endif
	dCount++;
	/*
	* Initialize more mprocess objects, if lThread is not 0xFFFFFFFF, which signifies default single
	* threaded operation.
	*/
	if(lThread != 0xFFFFFFFF)	{
		while((unsigned long)(dCount) < lThreads)	{
#ifdef MSVC
			pHandle[dCount] = CreateThread(NULL,0,RefineThread,(void *)pProcess[dCount],0,&pId[dCount]);
#else
			pthread_create(&pThreads[dCount],NULL,RefineThread,(void*)pProcess[dCount]);
#endif
			dCount++;
		}
	}
	/*
	* wait until all of the mprocess objects return.
	*/
#ifdef MSVC
//	wait = WaitForMultipleObjects(dCount,pHandle,true,INFINITE);
	a = 0;
	iTics = 0;
	while(a < (unsigned long)(dCount))	{
		wait = WaitForSingleObject(pHandle[a],10000);
		if(a > 0 && wait == WAIT_TIMEOUT)	{
			if(a == 1)	{
//				cout << "waiting for " << a+1;
				Rprintf("waiting for %lu", a+1);
			}
			else	{
//				cout << a+1;
				Rprintf("%lu", a+1);
			}
			//cout.flush();
			while(wait == WAIT_TIMEOUT)	{
				wait = WaitForSingleObject(pHandle[a],dwTime);
				if(wait == WAIT_TIMEOUT)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					iTics++;
					if(iTics > 50)	{
//						cout << "|\n\t\t";
						Rprintf("|\n\t\t");
						//cout.flush();
						iTics = 0;
					}
				}
			}
		}
		else	{
			while(wait == WAIT_TIMEOUT)	{
				wait = WaitForSingleObject(pHandle[a],dwTime);
				if(wait == WAIT_TIMEOUT)	{
//					cout << ":";
					Rprintf(":");
					//cout.flush();
				}
			}
			if(a == 1)	{
//				cout << "waiting for " << a+1;
				Rprintf("waiting for %lu", a+1);
			}
			else if(a == 0)	{
//				cout << "\n\t";
				Rprintf("\n\t");
				//cout.flush();
			}
			else	{
//				cout << a+1;
				Rprintf("%lu", a+1);
			}
		}
		a++;
	}
	if(dCount > 1)	{
//		cout << " done.\n\n";
		Rprintf(" done.\n\n");
		//cout.flush();
	}
	else	{
//		cout << "\n";
		Rprintf("\n");
		//cout.flush();
	}
#else
	//2003-03-01:note - the declaration below was changed from void **vp;	
	x=0;
	for(x=0;x<dCount;x++){
		//2003-03-01:note - the 2nd parameter in the call to pthread_join() was changed from vp
		wait = pthread_join(pThreads[x],&vp);
	}
#endif
	a = 1;
	/*
	* merge the results into the first object
	*/
	if(dCount > 1)	{
//		cout << "Merging results:\n";
		Rprintf("Merging results:\n");
		//cout.flush();
	}
	while(a < (unsigned long)(dCount))	{
		if(a == 1)	{
//			cout << "\tfrom " << a+1;
			Rprintf("\tfrom %lu", a+1);
		}
		else	{
//			cout << a+1;
			Rprintf("%lu", a+1);
		}
		//cout.flush();
		if(!pProcess[0]->add_spectra(pProcess[a]->m_vSpectra))	{
//			cout << "adding spectra failed.\n";
			Rprintf("adding spectra failed.\n");
		}
		pProcess[0]->merge_statistics(pProcess[a]);
		pProcess[a]->clear();
		pProcess[a]->m_mapSequences.clear();
		delete pProcess[a];
		pProcess[a] = NULL;
		a++;
	}
	if(dCount > 1)	{
//		cout << "\n\n";
		Rprintf("\n\n");
		//cout.flush();
	}
	//cout.flush();
//	cout << "Creating report:\n";
	Rprintf("Creating report:\n");
	//cout.flush();
	pProcess[0]->report();
	Rcpp::CharacterVector pathName(pProcess[0]->getPathName()); // rTANDEM
	size_t tValid = pProcess[0]->get_valid();
	size_t tUnique = pProcess[0]->get_unique();
	double dE = pProcess[0]->get_error_estimate();
	unsigned long lE = (unsigned long)(0.5+dE);
	unsigned long lEe = (unsigned long)(0.5 + sqrt(dE));
	if(lEe == 0)	{
		lEe = 1;
	}
	if(dE <= 0.0)	{
		dE = 1.0;
	}
//	cout << "\nValid models = " << (unsigned long)tValid << "\n";
	Rprintf("\nValid models = %lu\n", (unsigned long)tValid);
	if(tUnique > 0)	{
//		cout << "Unique models = " << (unsigned long)tUnique << "\n";
		Rprintf("Unique models = %lu\n", (unsigned long)tUnique);
//		cout << "Estimated false positives = " << lE << " +/- ";
		Rprintf("Estimated false positives = %lu +/- ", lE);
//		cout << lEe << "\n";
		Rprintf("%lu\n", lEe);
	}
//	lE = pProcess[0]->get_reversed();
//	cout << "\n\n";
	Rprintf("\n\n");
	/*
	* Delete the mprocess objects and exit
	*/
	a = 0;
	while(a < 16)	{
		if(pProcess[a] != NULL)	{
#ifdef MSVC
			CloseHandle(pHandle[a]);
#endif
			delete pProcess[a];
		}
		a++;
	}
	delete pProcess;
	delete pId;
	delete pHandle;
//	return 0; // rTANDEM
//	return R_NilValue; // rTANDEM
	return pathName;
}
/*
 * Process thread is used to create the worker threads for each mprocess object
 */

#ifdef MSVC
DWORD WINAPI ProcessThread(LPVOID _p)
{
	((mprocess *)_p)->process();
	return 0;
}
DWORD WINAPI RefineThread(LPVOID _p)
{
	((mprocess *)_p)->refine();
	return 0;
}
#else
void* ProcessThread(void *_p){
	((mprocess *)_p)->process();
 	return (void*)0;
}
void* RefineThread(void *_p){
	((mprocess *)_p)->refine();
 	return (void*)0;
}
#endif

#endif
