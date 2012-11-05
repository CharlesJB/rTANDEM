/*
 Copyright (C) 2003 Ronald C Beavis, all rights reserved
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

#ifndef MPROCESS_H
#define MPROCESS_H

// File version: 2003-08-01
// File version: 2004-03-01
// File version: 2004-11-01

typedef map<size_t,string> SEQMAP;

#include <sys/timeb.h>
#include <ctime>
#include "mcleave.h"
#include "mspectrumcondition.h"
#include "mreport.h"
#include "msequenceserver.h"
#include "mplugin.h"
#include "msemistate.h"
#include "mrefine.h" 
#include <set>

#include <Rcpp.h>  // rTANDEM
#include "dataLoader.h" // rTANDEM

/*
 * the process object coordinates the function of tandem. it contains the information
 * loaded from the input XML file in the m_xmlValues object and performance
 * information in the m_xmlPerformance object. The mass spectra to be analyzed are
 * in the m_vSpectra vector container. A set of input parameters are used to
 * initialize constants that are used in processing the mass spectra.
 */
class merrors
{
public:
	merrors(void) { 
		m_bPpm = true; 
		m_bIsotope = false;
		m_fPlus = 100.; 
		m_fMinus = 100.;
	}
	virtual ~merrors(void) { }
	bool m_bPpm;
	bool m_bIsotope;
	float m_fPlus;
	float m_fMinus;
	bool check(double _s,double _m)	{
		float fDelta = (float)(_s - _m);
		float fPlus = m_fPlus;
		float fMinus = m_fMinus;
		if(m_bPpm)	{
			fPlus *= (float)(_m*1.0e-6);
			fMinus*= (float)(_m*1.0e-6);
		}
		if(fDelta < 0.0)	{
			if(fDelta >= fMinus)	{
				return true;
			}
		}
		else	{
			if(fDelta <= fPlus)	{
				return true;
			}
		}
		if(!m_bIsotope)	{
			return false;
		}
		if(_s > 1000.0)	{
			fDelta -= (float)1.00335;
			if(fDelta < 0.0)	{
				if(fDelta >= fMinus)	{
					return true;
				}
			}
			else	{
				if(fDelta <= fPlus)	{
					return true;
				}
			}
		}
		if(_s > 1500.0)	{
			fDelta -= (float)1.00335;
			if(fDelta < 0.0)	{
				if(fDelta >= fMinus)	{
					return true;
				}
			}
			else	{
				if(fDelta <= fPlus)	{
					return true;
				}
			}
		}
		return false;
	}
};

class mprocesslog
{
public:
	virtual ~mprocesslog() { }
	bool open(string &_s)	{
		m_ofLog.open(_s.c_str());
		if(m_ofLog.fail())	{
			return false;
		}
		return true;
	}
	bool log(string &_m)	{
		if(!m_ofLog.is_open())	{
			return false;
		}
		time_t tValue;
		time(&tValue);
		struct tm *tmValue = localtime(&tValue);
		char pLine[256];	
		strftime(pLine, 255,"%Y-%m-%d %H:%M:%S",tmValue);
		m_ofLog << pLine << "\t" << _m.c_str() << "\n";
		m_ofLog.flush();
		return true;
	}
	bool log(const char *_m)	{
		if(!m_ofLog.is_open())	{
			return false;
		}
		string strValue = _m;
		return log(strValue);
	}
	bool close()	{
		if(!m_ofLog.is_open())	{
			return false;
		}
		m_ofLog.close();
		return true;
	}
private:
	ofstream m_ofLog;
};

class mprocess
{
public:
	mprocess(void);
	virtual ~mprocess(void);
	bool serialize(void); //serializes the m_vMI data elements in m_vSpectra onto the disk to save space during calculation
	bool restore(void); //restores the m_vMI data elements from disk
	bool removeMI(void);
	vector<string> m_vstrPaths;
	mprocesslog m_prcLog;
	bool add_spectra(vector<mspectrum> &_v); // adds the spectra contained in _v to m_vSpectra
	bool clear(); // clears a selection of vectors in the mprocess object
	unsigned long get_thread(); // retrieves the number of the current thread (0 - 15)
	unsigned long get_threads(); // retrieves the total number of threads running (1 - 16)
	size_t get_protein_count(); // gets the total number of proteins read
	size_t get_peptide_count(); // gets the total number of peptides used
	size_t get_total_residues(); // gets the total number of residues read
	size_t get_valid(); // gets the number of valid peptide models
	double get_error_estimate() {return m_dEsum;}
	size_t get_unique(); // gets the number of valid peptide models
	long get_reversed(); // gets the number of reversed peptide models
	double get_threshold(); // gets the number of valid peptide models
	bool load(const char *_f,mprocess *_p = NULL); // loads input parameters
	bool load(SEXP param, SEXP taxonomy, SEXP saps, SEXP mods, SEXP spectrum, mprocess *_p = NULL);  // rTANDEM
	bool load_saps(mprocess *_p); // loads sap information, if it exists
	bool load_annotation(mprocess *_p); // loads sequence annotation information, if it exists
	virtual bool merge_spectra(); // adds externally generated mspectrum vector to m_vSpectra
	virtual bool merge_spectra(vector<mspectrum> &_s); // adds externally generated mspectrum vector to m_vSpectra
	bool merge_map(SEQMAP &_s); // adds externally generated mspectrum vector to m_vSpectra
	bool merge_statistics(const mprocess *_p); // adds externally generated spectra to this object
	bool process(void); // performs identifications based on the input parameters
	bool refine(void); // controls the protein model refinement process
	bool report(void); // produces the XML output report, using an mreport class
	bool pyro_check(const char _c);
	bool pyro_reset(void);
	bool score_each_sequence(); // generates a score for all sequences in the m_svrSequences object
	bool set_threads(const unsigned long _t); // sets the object's thread number
	bool set_thread(const unsigned long _t); // sets the object's thread number
	int set_round(const int _r)	{ m_iCurrentRound = _r; return m_iCurrentRound;}
	virtual bool load_sequences();
	XmlParameter m_xmlPerformance; // stores process performance parameters
	XmlParameter m_xmlValues; // store process input parameters
	vector<mspectrum> m_vSpectra; // store spectra to be analyzed
	SEQMAP m_mapSequences; // a map containing all of the protein sequences discovered, indexed by their m_tUid value
	vector<msequence> m_vseqBest; // a vector of msequences used in the model refinement process
	vector<string> m_vstrModifications; //a vector containing the strings defining fixed modifications for a protein
	size_t m_tRefineModels; // total number of models generated by refinement
	size_t m_tRefineInput; // total number of sequences included in a refinement session
	size_t m_tRefinePartial; // the number of models discovered to have partial cleavage
	size_t m_tRefineUnanticipated; // the number of models discovered to have unanticpated cleavage
	size_t m_tRefineNterminal; // the number of models discovered to have modified N-terminii
	size_t m_tRefineCterminal; // the number of models discovered to have modified C-terminii
	size_t m_tRefinePam; // the number of models discovered to have point mutations
	double m_dRefineTime; // the time required to perform a refinement
	size_t m_tActive;	// total number of models remaining after each refinement step
	bool m_bRefineCterm;  //true if processing 'refine, potential C-terminus modifications'. Set in mrefine::refine and 
						//checked in score(), so the start position can be set to the length of the protein sequence
						// minus the value for 'refine, potential N-terminus modification position limit' before performing cleavage
	string getPathName() { return m_pathName; } // rTANDEM
	vector<int> m_viQuality; // contains the data quality scoring vector
	bool m_bReversedOnly;
	bool m_bSaps;
	bool m_bAnnotation;
	bool m_bMinimalAnnotation;
	bool m_bSerialize;
	bool m_bCheckNg;
	double m_dNt;
	double m_dNtAve;
	double m_dNg;
	double m_dNgAve;
	enum	{
		I_Y =	0x01,
		I_B =	0x02,
		I_X =	0x04,
		I_A =	0x08,
		I_C =	0x10,
		I_Z =	0x20,
	} ion_type; // enum for referencing information about specific ion types.

protected:
	/*
	* refinement classes are declared friend classes so they can access private member variables and functions. 
	* If a new refinement class is added to the project, include a declaration here
	*/
	friend class mrefine;
	friend class mpmods;
	friend class mxxcleavage;
	friend class mtermmods;
	friend class mpam;

protected:
	string m_strLastMods;
	int m_iCurrentRound;
	bool m_bPermute;
	bool m_bPermuteHigh;
	bool m_bCrcCheck;
	bool m_bQuickAcetyl;
	bool m_bQuickPyro;
	double m_dEsum;

	set<size_t> m_setRound;
	vector<string> m_vstrSaps;
	vector<string> m_vstrMods;
	map<string,string> m_mapAnnotation;
	msemistate m_semiState; // maintains the state of the semi-enzymatic cleavage state machine
	mpyrostate m_pyroState; // maintains the state of the pyrolidone carboxylic acid detection state machine
	merrors m_errValues;
	double m_dSearchTime; // total time elapsed during a protein modeling session process
	long m_lIonCount; // minimum sum of detected ions that are significant enough to store a sequence 
	unsigned long m_lThread; // thread number of this object
	unsigned long m_lThreads; // the total number of threads current active
	long m_lReversed; // the total number of peptides found where the reversed sequence was better than the forward sequence
	double m_dThreshold; // the current expectation value threshold
	double m_tContrasted; // the number of spectra subtracted using contrast angle redundancy detection
	long m_lStartMax; // set the maximum distance from N-terminus for a peptide
					  // normally set at an impossibly large value = 100,000,000
					  // for ragged N-terminii with potential modifications, set at a low but plausible value = 50
	long m_lCStartMax;
	char *m_pSeq; // a character pointer, used for temporary sequence information
	bool m_bUn; // if true, cleave at all residues. if false, use cleavage specification in data input.
	bool m_bUseHomologManagement; // set to true to use homologue management 
	size_t m_tMinResidues; // the minimum peptide length that will be scored
	size_t m_tMissedCleaves; // the maximum number of cleavage sites that can be missed
	size_t m_tPeptideCount; // the total number of peptide sequences generated during a process
	size_t m_tPeptideScoredCount; // the total number of peptide sequences scored during a process
	size_t m_tProteinCount; // the total number of protein sequences considered during a process
	size_t m_tSpectra; // the total number of spectra being modeled
	size_t m_tSpectraTotal; // the total number of spectra in the input file
	size_t m_tValid; // the number of valid peptide models
	size_t m_tTotalResidues; // the number of residues read
	size_t m_tSeqSize; // current length of the m_pSeq character array
	size_t m_tUnique; // the number of unique peptides found in a result
	string m_strOutputPath; // the path name of the XML output file
	mcleave m_Cleave; // the specification for a cleavage peptide bond
	msequence m_seqCurrent; // the msequence object that is currently being scored
#ifdef X_P3
	p3msequenceServer m_svrSequences; // the msequenceServer object that provides msequences to msequenceCollection
#else
	msequenceServer m_svrSequences; // the msequenceServer object that provides msequences to msequenceCollection
#endif
	mspectrumcondition m_specCondition; // the mspectrumcondition object that cleans up and normalized
										// spectra for further processing
	mscore* m_pScore; // the object that is used to score sequences and spectra
	mrefine* m_pRefine; // the object that is used to refine models

	bool charge(void); // adds additional charge states to the list of spectra, resulting in +1, +2 and +3 all being represented
	virtual bool clean_sequences(void); // remove sequences no longer required from m_mapSequences
	double dot(const size_t _f,const size_t _s,const float _r,const bool _t); // calculated inner product between two spectrum vectors
	bool create_rollback(vector<mspectrum> &_v); // create a temporary rollback vector
	virtual bool create_score(const msequence &_s,const size_t _v,
								const size_t _w, const long _m,const bool _p); // generates scores for a sequence
	double expect_protein(const unsigned long _c,const unsigned long _t,
							const unsigned long _n,const double _d); // assigns an expectation value for a protein
	bool load_best_vector(void); // creates a new m_vseqBest vector and marks assigned spectra inactive
	bool mark_repeats(void); // checks for repeated assignments of the same peptide sequence
	bool modify(void); // sets the initial residue modification parameters in the m_seqUtil object
	bool refine_model(void); // the method that refines protein models
	bool report_all(); // create a report that contains all peptide models, regardless of expectation value
	bool report_valid(const double _d); // create a report that only contains valid peptide models (expect < m_dThreshold)
	bool report_stochastic(const double _d); // create a report that contains only stochastic peptide models (expect > m_dThreshold)
	bool report_expect(const double _m); // calculates expectation values for the report
	bool report_sort(void); // sorts peptide models for use in the report
	bool residues(); // estimates the minimum number of residues that can possibly be contained in a peptide
	bool rollback(vector<mspectrum> &_v,const double _m,const double _f);
	bool score(const msequence &_s); // generates scores of an msequence object
	bool score_single(const msequence &_s); // generates scores of an msequence object
	bool score_terminus(const string &_s); // attempts to find N or C terminal modified peptides
	bool score_terminus_single(const string &_s); // attempts to find N or C terminal modified peptides
	bool spectra(void); // loads the m_vSpectra object using the m_specCondition object
	bool spectra_force(string &_t,string &_v); // forces the spectrum loader to use a specified file type
	bool subtract(void); // remove redundant mass spectra
	virtual bool taxonomy(void); // loads the taxonomy setting into the m_svrSequences object
private:
	string m_pathName; // rTANDEM
};

#include "p3mprocess.h"

#endif
