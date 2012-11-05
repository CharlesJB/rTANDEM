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
#include "mxxcleavage.h"
#include <sys/timeb.h>
#include <ctime>
#include "mprocess.h"
#include "mscore.h"

mxxcleavage::mxxcleavage(){
	m_pProcess=NULL;
}
mxxcleavage::~mxxcleavage(){
}

bool mxxcleavage::set_mprocess(mprocess *_p){
	if(m_pProcess != NULL)
		delete m_pProcess;
	m_pProcess = _p;
	return true;
}

/*
 * process cleavage at any site with 50 missed cleavages
 */
bool mxxcleavage::refine(){
	vector<mspectrum> vspRollback;
	size_t a = 0;
	size_t tPips = 0;
	size_t tActiveNow = 0;
	string strTemp;
	string strKey;
	string strValue;

/*
 *  set the maximum expectation value
 */
	strKey = "refine, tic percent";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	double dTicPercent = atof(strValue.c_str());
	if(dTicPercent == 0)	{
		dTicPercent = 20.0;
	}
	size_t tTicMax = (size_t)((double)m_pProcess->m_vseqBest.size()*dTicPercent/100.0);
	if(tTicMax < 1)	{
		tTicMax = 1;
	}
	strKey = "refine, maximum valid expectation value";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(strValue.size() > 0)	{
		m_dMaxExpect = atof(strValue.c_str());
	}
	string strCleave = m_pProcess->m_Cleave.m_strCleave;
	size_t tMissedCleaves = m_pProcess->m_tMissedCleaves;
	bool bUn = m_pProcess->m_bUn;
	bool bSemi = m_pProcess->m_semiState.m_bActive;
	strKey = "refine, full unanticipated cleavage";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		m_pProcess->m_bUn = true;
		strTemp = "[X]|[X]";
		m_pProcess->m_Cleave.load(strTemp);
		m_pProcess->m_tMissedCleaves = 50;
	}
	else	{
		m_pProcess->m_tMissedCleaves = 2;
		m_pProcess->m_semiState.activate(true);
	}
	if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
		cout << "\tunanticipated cleavage ";
		cout.flush();
		m_pProcess->m_prcLog.log("unanticipated cleavage");
	}
	m_pProcess->create_rollback(vspRollback);

/*
 * score the spectra against each sequence
 */
	while(a < m_pProcess->m_vseqBest.size())	{
		m_pProcess->score(m_pProcess->m_vseqBest[a]);
		tPips++;
		if(tPips == tTicMax)	{
			if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
				cout << ".";
				cout.flush();
				m_pProcess->m_prcLog.log(".");
			}
			tPips = 0;
		}
		a++;
	}
/*
 * update active state of spectra
 */
	m_pProcess->load_best_vector();
	a=0;
	while(a < m_pProcess->m_vSpectra.size())	{
		if(!m_pProcess->m_vSpectra[a].m_bActive)
			tActiveNow++;
		a++;
	}
	if(tActiveNow >= m_pProcess->m_tActive)	{
		m_pProcess->m_tRefineUnanticipated = tActiveNow - m_pProcess->m_tActive;
	}
	m_pProcess->m_tActive = tActiveNow;
	m_pProcess->rollback(vspRollback,m_dMaxExpect,0.1);
	m_pProcess->m_bUn = bUn;
	m_pProcess->m_Cleave.load(strCleave);
	m_pProcess->m_tMissedCleaves = tMissedCleaves;
	m_pProcess->m_semiState.activate(bSemi);
	if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
		cout << " done.\n";
		m_pProcess->m_prcLog.log("done");
	}
	cout.flush();
	return true;
}


/*
 * mxxcleavagemanager contains static short-cuts for dealing with mxxcleavage
 * plug-ins.
 */
const char* mxxcleavagemanager::TYPE = "refinement, unanticipated cleavages algorithm";

/*
 * create_mxxcleavage creates the correct mxxcleavage object for the given set of XML
 * parameters.
 */
mxxcleavage* mxxcleavagemanager::create_mxxcleavage(XmlParameter &_x)
{
	string strValue;
	string strKey = TYPE;
	if (!_x.get(strKey,strValue))
		strValue = "tandem";
	return (mxxcleavage*) mpluginmanager::get().create_plugin(TYPE, strValue.data());
}

/*
 * register_factory registers a factory with the mpluginmanager for
 * creating mxxcleavage derived objects.
 */
void mxxcleavagemanager::register_factory(const char* _spec, mpluginfactory* _f)
{
	mpluginmanager::get().register_factory(TYPE, _spec, _f);
}


// Factory instance, registers itself with the mxxcleavagemanager.
static mxxcleavagefactory_tandem factory;
	
mxxcleavagefactory_tandem::mxxcleavagefactory_tandem()
{
	mxxcleavagemanager::register_factory("tandem", this);
}

mplugin* mxxcleavagefactory_tandem::create_plugin()
{
	return new mxxcleavage();
}
