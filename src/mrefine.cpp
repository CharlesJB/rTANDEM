/*
 Copyright (C) 2005 Ronald C Beavis, all rights reserved
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

// File version: 2005-09-01

#include "stdafx.h"
#include <sys/timeb.h>
#include <ctime>
#include "mrefine.h"
#include "mspectrum.h"
#include "mprocess.h"
#include "mscore.h"
//includes for refinement classes
#include "mpmods.h"
#include "mxxcleavage.h"
#include "mtermmods.h"
#include "mpam.h"


mrefine::mrefine()
{
	m_pProcess=NULL;
}
mrefine::~mrefine(void)
{
}


bool mrefine::set_mprocess(mprocess *_p){
	if(m_pProcess != NULL)
		delete m_pProcess;
	m_pProcess = _p;
	return true;
}

/*
 *	initialize mprocess::m_vSeqBest vector and set refine input spectra count
 */
bool mrefine::initialize(){
	string strKey;
	string strValue;
	size_t a = 0;

	strKey = "refine, saps";
	m_pProcess->m_bSaps = true;
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(strValue == "no")	{
		m_pProcess->m_bSaps = false;
	}
	else	{
		m_pProcess->m_bSaps = true;
	}


/*
 *  load the m_vseqBest vector with assigned sequences
 */
	if(!m_pProcess->load_best_vector())
		return false;
	while(a < m_pProcess->m_vSpectra.size())	{
		if(!m_pProcess->m_vSpectra[a].m_bActive)
			m_pProcess->m_tActive++;
		a++;
	}
	strKey = "refine, modification mass";
	if(m_pProcess->m_xmlValues.get(strKey,strValue) && strValue.size() > 0)	{
		m_pProcess->m_vstrModifications.clear();
		m_pProcess->m_vstrModifications.push_back(strValue);
		int b = 1;
		char *pLine = new char[256];
		sprintf(pLine,"refine, modification mass %i",b);
		strKey = pLine;
		while(m_pProcess->m_xmlValues.get(strKey,strValue) && strValue.size() > 0) {
			m_pProcess->m_vstrModifications.push_back(strValue);
			b++;
			sprintf(pLine,"refine, modification mass %i",b);
			strKey = pLine;
		}
		delete pLine;
	}
	m_pProcess->m_tRefineInput = m_pProcess->m_vSpectra.size() - m_pProcess->m_tActive;
	return true;
}

/*
 *	refine sequentially calls each refinement object using the current mprocess object:
 *
 *   1. comparison with multiple potential modification
 *   2. full [X]|[X] cleavage
 *   3. N-terminal modifications
 *   4. C-terminal modifications
 *   5. point mutations
 *
 *   the refine objects are derived from mrefine and use the same facility to create new objects as mrefine does
 *   (register_factory() and create_plugin() )
 * 
 *   to add a new refinement class follow the steps below:
 *
 *	 1. add a friend class declaration to the mprocess class (see mprocess.h for the existing ones)
 *   2. add a forward class declaration to the mrefine class (see mrefine.h for the existing ones)
 *   3. add a private member variable pointer of the new class type to the mrefine class (see mrefine.h for the existing ones)
 *   4. add an #include statement for the new class header file in mrefine.cpp (see mrefine.cpp for exising ones)
 *   5. add code to mrefine::refine() to a)create an instance of the new class, b)set the mprocess object and c)process the new refinement algorithm
 *
 *		eg:
 			m_pNewClass = mnewclassmanager::create_mnewclass(m_pProcess->m_xmlValues);
			if (m_pNewClass == NULL) {
				cout << "Failed to create NewClass\n";
				return false;
			}
			m_pNewClass->set_mprocess(m_pProcess);
			m_pNewClass->refine();
 *
 *   6. create .cpp and .h files for the new class using an existing refinement plugin as a template
 */		

bool mrefine::refine()
{
	string strKey;
	string strValue;
	size_t a = 0;
	size_t tActiveNow = 0;
	int iRound = 2;
	initialize();
	m_pProcess->set_round(iRound); // round 2

/*
 *  1. check for potential modifications
 */
	m_pPMods = mpmodsmanager::create_mpmods(m_pProcess->m_xmlValues);
	if (m_pPMods == NULL) {
		cout << "Failed to create mpmods\n";
		return false;
	}
	m_pPMods->set_mprocess(m_pProcess);
	m_pPMods->refine();
	iRound = 3;
	m_pProcess->set_round(iRound); // round 3
	
/*
 *  2. perform full [X]|[X] cleavage if required
 */
	strKey = "refine, use potential modifications for full refinement";
	strValue = "no";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(strValue != "yes")	{
		strKey = "residue, potential modification mass";
		m_pProcess->m_xmlValues.get(strKey,strValue);
		m_pProcess->m_pScore->m_seqUtil.modify_maybe(strValue);
		strKey = "residue, potential modification motif";
		m_pProcess->m_xmlValues.get(strKey,strValue);
		m_pProcess->m_pScore->m_seqUtil.modify_motif(strValue);
	}
	strKey = "refine, unanticipated cleavage";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(m_pProcess->m_bSaps)	{
		strKey = "KR";
		m_pProcess->m_pScore->set_allowed_saps(strKey);
	}
	if(strValue == "yes")	{
		m_pXXCleavage = mxxcleavagemanager::create_mxxcleavage(m_pProcess->m_xmlValues);
		if (m_pXXCleavage == NULL) {
			cout << "Failed to create mxxcleavage\n";
			return false;
		}
		m_pXXCleavage->set_mprocess(m_pProcess);
		m_pXXCleavage->refine();
	}

	iRound = 4;
	m_pProcess->set_round(iRound); // round 4
/*
 *  3. check for modified peptide N-terminii
 *
 *  set the maximum distance from the N-terminus for a peptide's N-terminal residue
 *  at 50, unless otherwise specified
 */
	strKey = "refine, potential N-terminus modification position limit";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	long lStartMax = m_pProcess-> m_lStartMax;
	if(strValue.size() > 0)	{
		m_pProcess->m_lStartMax = atoi(strValue.c_str());
	}
	else	{
		m_pProcess->m_lStartMax = 50;
	}
	strKey = "refine, potential N-terminus modifications";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(fabs(atof(strValue.c_str())) > 0.001)	{
		m_pTermMods = mtermmodsmanager::create_mtermmods(m_pProcess->m_xmlValues);
		if (m_pTermMods == NULL) {
			cout << "Failed to create mtermmods\n";
			return false;
		}
		m_pTermMods->set_mprocess(m_pProcess);
		m_pTermMods->refine();
	}
/*
 *  reset the maximum distance from the N-terminus for a peptide's N-terminal residue
 */
	m_pProcess->m_lStartMax = lStartMax;
	m_pProcess->m_pScore->m_seqUtil.m_pdAaMod['['] = 0.0;

	iRound = 5;
	m_pProcess->set_round(iRound); // round 5
/*
 *  4. check for modified peptide C-terminii
 */
	strKey = "refine, potential C-terminus modifications";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(fabs(atof(strValue.c_str())) > 0.001)	{
		m_pProcess->m_bRefineCterm = true;
		m_pTermMods = mtermmodsmanager::create_mtermmods(m_pProcess->m_xmlValues);
		if (m_pTermMods == NULL) {
			cout << "Failed to create mtermmods\n";
			return false;
		}
		m_pTermMods->set_mprocess(m_pProcess);
		m_pTermMods->refine();
	}
	m_pProcess->m_bRefineCterm = false;
	m_pProcess->m_pScore->m_seqUtil.m_pdAaMod[']'] = 0.0;

	iRound = 6;
	m_pProcess->set_round(iRound); // round 6
/*
 * 5. check for point mutations in the sequences
 */
	strKey = "refine, point mutations";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		m_pProcess->m_bSaps = false;
		m_pPam = mpammanager::create_mpam(m_pProcess->m_xmlValues);
		if (m_pPam == NULL) {
			cout << "Failed to create mpam\n";
			return false;
		}
		m_pPam->set_mprocess(m_pProcess);
		m_pPam->refine();
	}

/*
 * 6. new mrefine derived classes here
 */

/*
 * finish up and return
 */
	if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
		cout << "\tfinishing refinement ... ";
		cout.flush();
	}
	m_pProcess->m_tRefineModels = m_pProcess->m_vseqBest.size();
	m_pProcess->m_vseqBest.clear();
	if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
		cout << "done.\n";
		cout.flush();
	}
	return true;
}

/*
 * mrefinemanager contains static short-cuts for dealing with mrefine
 * plug-ins.
 */
const char* mrefinemanager::TYPE = "refinement, algorithm";

/*
 * create_mrefine creates the correct mrefine for the given set of XML
 * parameters.
 */
mrefine* mrefinemanager::create_mrefine(XmlParameter &_x)
{
	string strValue;
	string strKey = TYPE;
	if (!_x.get(strKey,strValue))
		strValue = "tandem";
	return (mrefine*) mpluginmanager::get().create_plugin(TYPE, strValue.data());
}

/*
 * register_factory registers a factory with the mpluginmanager for
 * creating mrefine derived objects.
 */
void mrefinemanager::register_factory(const char* _spec, mpluginfactory* _f)
{
	mpluginmanager::get().register_factory(TYPE, _spec, _f);
}

// Factory instance, registers itself with the mrefinemanager.
static mrefinefactory_tandem factory;
	
mrefinefactory_tandem::mrefinefactory_tandem()
{
	mrefinemanager::register_factory("tandem", this);
}

mplugin* mrefinefactory_tandem::create_plugin()
{
	return new mrefine();
}
