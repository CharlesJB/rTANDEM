/*
 Copyright (C) 2003-2004 Ronald C Beavis, all rights reserved
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

// File version: 2004-02-01
// File version: 2004-03-01
// File version: 2004-03-01
// File version: 2004-07-12
// File version: 2004-10-05
// File version: 2004-11-01
// File version: 2005-01-01

/*
 * the process object coordinates the function of tandem. it contains the information
 * loaded from the input XML file in the m_xmlValues object and performance
 * information in the m_xmlPerformance object. The mass spectra to be analyzed are
 * in the m_vSpectra vector container. A set of input parameters are used to
 * initialize constants that are used in processing the mass spectra.
 * NOTE: see tandem.cpp for an example of how to use an mprocess class
 * NOTE: mprocess uses cout to report errors. This may not be appropriate for
 *       many applications. Feel free to change this to a more appropriate mechanism
 */

#include "stdafx.h"
#include <algorithm>
#include <set>
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "loadmspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"
#include "saxbiomlhandler.h"
#include "saxsaphandler.h"
#include "saxmodhandler.h"
#include "xmltaxonomy.h"
#include "mbiomlreport.h"
#include "mrefine.h"

// #define TANDEM_EXACT 1

/*
 * global less than operators to be used in sort operations
 */
bool lessThanSequence(const msequence &_l,const msequence &_r);
bool lessThanSequenceUid(const msequence &_l,const msequence &_r);
bool lessThanSequenceDes(const msequence &_l,const msequence &_r);
bool lessThanSpectrum(const mspectrum &_l,const mspectrum &_r);
bool lessThanOrder(const mspectrum &l,const mspectrum &r);
bool lessThanMass(const mi &_l,const mi &_r);

bool lessThanSequence(const msequence &_l,const msequence &_r)
{
	return _l.m_dExpect < _r.m_dExpect;
}

bool lessThanSequenceUid(const msequence &_l,const msequence &_r)
{
	return _l.m_tUid < _r.m_tUid;
}

bool lessThanSequenceDes(const msequence &_l,const msequence &_r)
{
	return _l.m_strDes.size() < _r.m_strDes.size();
}

bool lessThanSpectrum(const mspectrum &_l,const mspectrum &_r)
{
	return _l.m_dProteinExpect < _r.m_dProteinExpect;
}

bool lessThanSpectrumSequence(const mspectrum &_l,const mspectrum &_r)
{
	if(_l.m_vseqBest.empty())
		return false;
	if(_r.m_vseqBest.empty())
		return true;
	return _l.m_dExpect < _r.m_dExpect;
}

bool lessThanOrder(const mspectrum &_l,const mspectrum &_r)
{
	if(_l.m_vseqBest.empty())
		return false;
	if(_r.m_vseqBest.empty())
		return true;
	return _l.m_vseqBest[0].m_vDomains[0].m_lS < _r.m_vseqBest[0].m_vDomains[0].m_lS;
}

bool lessThanMass(const mi &_l,const mi &_r)
{
	return _l.m_fM < _r.m_fM;
}

mprocess::mprocess(void)
{

	m_iCurrentRound = 1;
	m_lReversed = -1;
	m_tMissedCleaves = 1;
/*
 * record the process start time
 */
    time_t tValue;
	time(&tValue);
	struct tm *tmValue = localtime(&tValue);
	char pLine[256];	
	strftime(pLine, 255,"%Y:%m:%d:%H:%M:%S",tmValue);
	string strKey;
	string strValue;
	strKey = "process, start time";
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
/*
 * record the version of the software
 */
	strKey = "process, version";
#ifdef X_P3
	strValue = "x! p3 ";
#else
	strValue = "x! tandem ";
#endif
	strValue += VERSION;
	m_xmlPerformance.set(strKey,strValue);
	m_tProteinCount = 0;
	m_tPeptideCount = 0;
	m_tTotalResidues = 0;
	m_tPeptideScoredCount = 0;
	m_tMinResidues = 0;
	m_lThread = 0;
	m_lThreads = 1;
	m_tValid = 0;
	m_tUnique = 0;
	m_tSpectraTotal = 0;
	m_bUn = false;
	m_tSeqSize = 4096*4;
	m_pSeq = new char[m_tSeqSize];
	m_lStartMax = 100000000;
	m_dThreshold = 1000.0;
	m_tRefineInput = 0;
	m_tRefinePartial = 0;
	m_tRefineUnanticipated = 0;
	m_tRefineNterminal = 0;
	m_tRefineCterminal = 0;
	m_tRefinePam = 0;
	m_tRefineModels = 0;
	m_tContrasted = 0;
	m_bUseHomologManagement = false;
	m_pScore = NULL;
	m_lCStartMax = 50;
	m_bCrcCheck = false;
	m_bRefineCterm = false;
	m_semiState.activate(false);
	m_tActive = 0;
	m_bReversedOnly = false;
	m_bSaps = false;
	m_bAnnotation = false;
	m_bMinimalAnnotation = false;
	m_bSerialize = false;
	m_bQuickAcetyl = true;
	m_bQuickPyro = true;
	m_bCheckNg = false;
	m_dNt = 0.0;
	m_dNtAve = 0.0;
	m_dNg = 0.0;
	m_dNgAve = 0.0;

}

mprocess::~mprocess(void)
{
	if(m_pSeq != NULL)
		delete m_pSeq;
	if(m_pScore != NULL)
		delete m_pScore;
	if(m_lThread == 0 || m_lThread == 0xFFFFFFFF)	{
		m_prcLog.log("X! Tandem exiting");
		m_prcLog.close();
	}
}
/*
* add spectra is used to load spectra into the m_vSpectra vector
*/
bool mprocess::add_spectra(vector<mspectrum> &_v)
{
	size_t a = 0;
	// changed from m_vSpectra.resize(m_vSpectra.size() + _v.size()+1) in 2006.02.01
	m_vSpectra.reserve(m_vSpectra.size() + _v.size()+1);
	size_t c = 0;
	while(a < _v.size())	{
		m_vSpectra.push_back(_v[a]);
		if(c == 1000)	{
//			cout << ".";
			Rprintf(".");
			//cout.flush();
			c = 0;
		}
		c++;
		a++;
	}
	return true;
}
/*
 * clean_sequences is used to maintain the m_mapSequences container. If a sequence
 * has been removed from the m_vSpectra list of sequences, then that sequence is
 * subsequently deleted from the m_mapSequences collection
*/
bool mprocess::clean_sequences(void)
{
	map<size_t,long> mapValue;
	map<size_t,long>::iterator itMap;
	size_t a = 0;
	size_t b = 0;
	size_t tLength = m_vSpectra.size();
	size_t tBest = 0;
	while(a < tLength)	{
		b = 0;
		tBest = m_vSpectra[a].m_vseqBest.size();
		while(b < tBest)	{
			mapValue[m_vSpectra[a].m_vseqBest[b].m_tUid] = 1;
			b++;
		}
		a++;
	}
	SEQMAP::iterator itValue = m_mapSequences.begin();
	while(itValue != m_mapSequences.end())	{
		itMap = mapValue.find((*itValue).first);
		if(itMap == mapValue.end())	{
			m_mapSequences.erase(itValue);
			itValue = m_mapSequences.begin();
		}
		else	{
			itValue++;
		}
	}
	return true;
}
/*
* clear is used to reset any vector or map object necessary to reset the mprocess object
*/
bool mprocess::clear(void)
{
	m_vSpectra.clear();
	if (m_pScore != NULL)
		m_pScore->clear();
	return true;
}

/*
* create_rollback is used to create a vector of mspectrum objects that serves as
* the record of values to be used by the rollback method.
*/
bool mprocess::create_rollback(vector<mspectrum>& _v)
{
	_v.clear();
	size_t a = 0;
	const size_t tSize = m_vSpectra.size();
	mspectrum spTemp;
	double dExpect = 0;
	_v.reserve(tSize);
	while(a < tSize)	{
		_v.push_back(spTemp);
		_v.back() *= m_vSpectra[a];
		m_vSpectra[a].m_hHyper.model();
		m_vSpectra[a].m_hHyper.set_protein_factor(1.0);
		dExpect = (double)m_vSpectra[a].m_hHyper.expect_protein(m_pScore->hconvert(m_vSpectra[a].m_fHyper));
		_v.back().m_dExpect = dExpect;
		a++;
	}
	return true;
}
/*
 * create_score takes an msequence object and a start and an end sequence position
 * and saves that msequence, its mdomain and scoring in the spectrum that has just
 * been scored. equivalent msequence objects (based on their hyper score)
 * are stored sequentially in a vector in the mspectrum object
 */
bool mprocess::create_score(const msequence &_s,const size_t _v,const size_t _w,const long _m,bool _p)
{
	long lIonCount = 0;
	float fScore = -1.0;
	float fHyper = -1.0;
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	bool bOk = false;
	bool bDom = false;
	long lCount = 0;
	bool bIonCheck = false;
	bool bMassCheck = false;
/*
 * score each mspectrum identified as a candidate in the m_pScore->m_State object
 */
	set<size_t> setDone;
	while(lCount < m_pScore->m_State.m_lEqualsS)	{
		a = m_pScore->m_State.m_plEqualsS[lCount];
		lCount++;
		lIonCount = 0;
/*
* this check is needed to keep tandem consistent whether running on
* single-threaded, multi-threaded or on a cluster.  otherwise, when
* there are multiple spectra matching a sequence (which is more common
* the fewer mprocess objects there are) one can cause others to score
* more permutation sequences.
*/
		if (!_p && m_vSpectra[a].m_hHyper.m_ulCount >= 400)
			continue;
		fScore = 1.0;
		fHyper = 1.0;
		m_pScore->m_lMaxCharge = (long)(m_vSpectra[a].m_fZ+0.1);
/*
 * in versions prior to 2004.03.01, spectra with m_bActive == false were
 * rejected at this point, to save time & because of a problem with
 * multiple recording of the same sequence. starting with 2004.03.01,
 * the later problem has been corrected, and because of point mutation
 * analysis, it has become important to reexamine all sequences.
 */
		fScore = m_pScore->score(a);
		fHyper = m_pScore->m_fHyper;
/*
 * If the convolution score is greater than 2.0, record information in the ion-type histograms
 */
		if(fScore > 2.0)	{
			m_tPeptideScoredCount++;
			lIonCount = m_pScore->m_plCount[mscore::S_B] + m_pScore->m_plCount[mscore::S_Y];
			lIonCount += m_pScore->m_plCount[mscore::S_C] + m_pScore->m_plCount[mscore::S_Z];
			lIonCount += m_pScore->m_plCount[mscore::S_A] + m_pScore->m_plCount[mscore::S_X];
			m_vSpectra[a].m_hHyper.add(m_pScore->hconvert(fHyper));
			m_vSpectra[a].m_hConvolute.add(m_pScore->hconvert(fScore));
			m_vSpectra[a].m_chBCount.add(m_pScore->m_plCount[mscore::S_A] + 
								m_pScore->m_plCount[mscore::S_B] + 
								m_pScore->m_plCount[mscore::S_C]);
			m_vSpectra[a].m_chYCount.add(m_pScore->m_plCount[mscore::S_X] + 
								m_pScore->m_plCount[mscore::S_Y] + 
								m_pScore->m_plCount[mscore::S_Z]);

		}
/*
* this check is must be outside the above conditional to keep tandem
* consistent whether running on single-threaded, multi-threaded or on
* a cluster.  otherwise, when there are multiple spectra matching a
* sequence (which is more common the fewer mprocess objects there are)
* one can cause others to score differently from how they would score
* alone.
*/
		if (m_vSpectra[a].m_hHyper.m_ulCount < 400)     {
			if(m_bCrcCheck && _p)   {
				m_bPermute = true;
			}
			else if(m_vSpectra[a].m_dMH > 3000.0)     {
				m_bPermuteHigh = true;
			}
		}
		double dRatio = (double)lIonCount/(double)(_w - _v + 1);
		if(m_vSpectra[a].m_fZ < 1.1 && dRatio > 0.2 && lIonCount > m_lIonCount)	{
			bIonCheck = true;
		}
		if(m_vSpectra[a].m_fZ > 1.0 && dRatio > 0.5 && lIonCount > m_lIonCount)	{
			bIonCheck = true;
		}
		if(lIonCount > m_lIonCount)	{
			bIonCheck = true;
		}
		if(m_vSpectra[a].m_fZ > 1.0 && dRatio > 0.5 && lIonCount > m_lIonCount)	{
			bIonCheck = true;
		}
		bMassCheck = m_errValues.check(m_vSpectra[a].m_dMH,m_pScore->seq_mh());
		if(!_p)	{
			bMassCheck = false;
		}
/*
 * if the same score has been recorded for another peptide in the same msequence object,
 * add a domain to that object, but do not update the entire object
 */
		if(bMassCheck && bIonCheck && fHyper == m_vSpectra[a].m_fHyper && m_vSpectra[a].m_tCurrentSequence == _s.m_tUid)	{
			vector<maa> vAa;
			mdomain domValue;
			domValue.m_vAa.clear();
			double dDelta = 0.0;
			if(m_pScore->get_aa(vAa,_v,dDelta))	{
				b = 0;
				while(b < vAa.size())	{
					domValue.m_vAa.push_back(vAa[b]);
					b++;
				}
			}
			if(fabs(dDelta) > 2.5)	{
				domValue.m_fScore = fScore;
				domValue.m_fHyper = fHyper;
				domValue.m_dMH = m_pScore->seq_mh();
				// m_fDelta was changed to m_dDelta in 2006.02.01
				domValue.m_dDelta = m_vSpectra[a].m_dMH - m_pScore->seq_mh();
				domValue.m_lE = (int)_w;
				domValue.m_lS = (int)_v;
				domValue.m_lMissedCleaves = (int)_m;
				domValue.m_bUn = m_bUn;
				unsigned char lType = 1;
				unsigned long sType = 1;
				while(lType < m_pScore->m_lType+1)	{
					m_vSpectra[a].m_mapScore[lType] = m_pScore->m_pfScore[sType];
					m_vSpectra[a].m_mapCount[lType] = (unsigned int)m_pScore->m_plCount[sType];
					lType *= 2;
					sType++;
				}
				b = 0;
				bOk = true;
				while(bOk && b < m_vSpectra[a].m_vseqBest.back().m_vDomains.size())	{
					if(domValue == m_vSpectra[a].m_vseqBest.back().m_vDomains[b]){	
						bOk = false;
					}
					b++;
				}
				if(bOk)	{
					m_vSpectra[a].m_vseqBest.back().m_vDomains.push_back(domValue);
				}
			}
		}
/*
 * if the same hyper score has been recorded for a different msequence object, retain that
 * object and add the new msequence to the back of the m_vseqBest vector
 */
		else if(bMassCheck && bIonCheck && fHyper == m_vSpectra[a].m_fHyper)	{			
			vector<maa> vAa;
			msequence seqValue;
			mdomain domValue;
			domValue.m_vAa.clear();
			double dDelta = 0.0;
			if(m_pScore->get_aa(vAa,_v,dDelta))	{
				b = 0;
				while(b < vAa.size())	{
					domValue.m_vAa.push_back(vAa[b]);
					b++;
				}
			}
			if(fabs(dDelta) > 2.5)	{
				domValue.m_fScore = fScore;
				domValue.m_fHyper = fHyper;
				domValue.m_dMH = m_pScore->seq_mh();
				// m_fDelta was changed to m_dDelta in 2006.02.01
				domValue.m_dDelta = m_vSpectra[a].m_dMH - m_pScore->seq_mh();
				domValue.m_lE = (int)_w;
				domValue.m_lS = (int)_v;
				domValue.m_lMissedCleaves = (unsigned char)_m;
				domValue.m_bUn = m_bUn;
				unsigned char lType = 1;
				unsigned long sType = 1;
				while(lType < m_pScore->m_lType+1)	{
					m_vSpectra[a].m_mapScore[lType] = m_pScore->m_pfScore[sType];
					m_vSpectra[a].m_mapCount[lType] = (unsigned int)m_pScore->m_plCount[sType];
					lType *= 2;
					sType++;
				}
				seqValue = _s;
				seqValue.m_strSeq = " ";
				m_mapSequences.insert(SEQMAP::value_type(_s.m_tUid,_s.m_strSeq));
				seqValue.m_vDomains.clear();
				seqValue.m_vDomains.push_back(domValue);
				seqValue.format_description();
				bOk = true;
				b = 0;
				while(bOk && b < m_vSpectra[a].m_vseqBest.size())	{
					if(m_vSpectra[a].m_vseqBest[b].m_tUid == _s.m_tUid)	{
						c = 0;
						bDom = true;
						while(bDom && c < m_vSpectra[a].m_vseqBest[b].m_vDomains.size())	{
							if(m_vSpectra[a].m_vseqBest[b].m_vDomains[c] == domValue)	{
								bDom = false;
								bOk = false;
							}
							c++;
						}
						if(bDom)	{
							m_vSpectra[a].m_vseqBest[b].m_vDomains.push_back(domValue);
							bOk = false;
						}
					}
					b++;
				}
				if(bOk)	{
					m_vSpectra[a].m_dRatio = dRatio;
					m_vSpectra[a].m_tCurrentSequence = _s.m_tUid;
					seqValue.m_iRound = m_iCurrentRound;
					m_vSpectra[a].m_vseqBest.push_back(seqValue);
				}
			}
		}
/*
 * if the hyper score is the best found so far for the mspectrum, delete the old msequence
 * objects and record this one.
 */
		else if(bMassCheck && bIonCheck && fHyper > m_vSpectra[a].m_fHyper && lIonCount > m_lIonCount)	{
			vector<maa> vAa;
			msequence seqValue;
			mdomain domValue;
			domValue.m_vAa.clear();
			double dDelta = 0.0;
			if(m_pScore->get_aa(vAa,_v,dDelta))	{
				b = 0;
				while(b < vAa.size())	{
					domValue.m_vAa.push_back(vAa[b]);
					b++;
				}
			}
			if(fabs(dDelta) > 2.5)	{
				m_vSpectra[a].m_fScoreNext = m_vSpectra[a].m_fScore;
				m_vSpectra[a].m_fHyperNext = m_vSpectra[a].m_fHyper;
				m_vSpectra[a].m_fScore = fScore;
				m_vSpectra[a].m_fHyper = fHyper;
				domValue.m_fScore = fScore;
				domValue.m_fHyper = fHyper;
				domValue.m_lE = (int)_w;
				domValue.m_lS = (int)_v;
				domValue.m_lMissedCleaves = (unsigned char)_m;
				domValue.m_dMH = m_pScore->seq_mh();
				// m_fDelta was changed to m_dDelta in 2006.02.01
				domValue.m_dDelta = m_vSpectra[a].m_dMH - m_pScore->seq_mh();
				domValue.m_bUn = m_bUn;
				unsigned char lType = 1;
				unsigned long sType = 1;
				while(lType < m_pScore->m_lType+1)	{
					m_vSpectra[a].m_mapScore[lType] = m_pScore->m_pfScore[sType];
					m_vSpectra[a].m_mapCount[lType] = (unsigned int)m_pScore->m_plCount[sType];
					lType *= 2;
					sType++;
				}
				seqValue = _s;
				seqValue.m_strSeq = " ";
				m_mapSequences.insert(SEQMAP::value_type(_s.m_tUid,_s.m_strSeq));
				seqValue.m_vDomains.clear();
				seqValue.m_vDomains.push_back(domValue);
				seqValue.format_description();
				m_vSpectra[a].m_tCurrentSequence = _s.m_tUid;
				m_vSpectra[a].m_vseqBest.clear();
				seqValue.m_iRound = m_iCurrentRound;
				m_vSpectra[a].m_dRatio = dRatio;
				m_vSpectra[a].m_vseqBest.push_back(seqValue);
			}
		}
		else if (fScore > 2.0 && fHyper > m_vSpectra[a].m_fHyperNext)
		{
			m_vSpectra[a].m_fScoreNext = fScore;
			m_vSpectra[a].m_fHyperNext = fHyper;
		}
	}
	return true;
}
//
//	This method is used to perform an inner product between two mass spectra. The fragment
//	ion mass accuracy is used to set the binning for the two vectors
//
__inline__ double mprocess::dot(const size_t _f,const size_t _s,const float _r,const bool _t)
{ 
	float fValue = 0.0;
	vector<mi>::iterator itA = m_vSpectra[_f].m_vMI.begin();
	vector<mi>::iterator itB = m_vSpectra[_s].m_vMI.begin();
	vector<mi>::const_iterator itAEnd = m_vSpectra[_f].m_vMI.end(); 
	vector<mi>::const_iterator itBEnd = m_vSpectra[_s].m_vMI.end();
	// if _t == true, then the mass error, _r, is given in Daltons
	if(_t)	{
		while(itA != itAEnd)	{
			while(itB != itBEnd)	{
				if(fabs(itB->m_fM - itA->m_fM) <= _r)	{
					fValue += itB->m_fI*itA->m_fI;
				}
				if(itB->m_fM > itA->m_fM)	{
					break;
				}
				itB++;
			}
			itA++;
		}
	}
	// deal with the case where _r is in ppm
	else	{
		const float fRes = (float)(_r/1.0e6);
		float fW = fRes;
		while(itA != itAEnd)	{
			fW = itA->m_fM * fRes;
			while(itB != itBEnd)	{
				if(fabs(itB->m_fM - itA->m_fM) <= fW)	{
					fValue += itB->m_fI*itA->m_fI;
				}
				if(itB->m_fM > itA->m_fM)	{
					break;
				}
				itB++;
			}
			itA++;
		}
	}
	return (double)fValue;	
}
/*
 * expect_protein is used to assign the expectation value for a protein, if more
 * than one peptide has been found for that protein. the expectation values for
 * the peptides are combined with a simple Bayesian model for the probability of
 * having two peptides from the same protein having the best score in different
 * spectra.
 */
double mprocess::expect_protein(const unsigned long _c,const unsigned long _t,
								const unsigned long _n,const double _d)
{
	double dValue = _d+log10((double)m_tProteinCount);
	if(_c == 1 && _d < 0.0)	{
		return _d;
	}
	else if(_c == 1)	{
		return 1.0;
	}
	if(_c == 0)	{
		return 1.0;
	}
	double dN = _n;
	double dK = _c;
	double dV = _t;
	unsigned long a = 0;
	while(a < _c)	{	
		dValue += log10((dV - a)/(dK - a));
		a++;
	}
	dValue -= log10(dV);
	dValue -= (dK-1.0)*log10(dN);
	double dP = dN/(double)m_tPeptideCount;
	if(dP >= 1.0)
		dP = 0.9999999;
	double dLog = dK*log10(dP)+(dV-dK)*log10(1.0-dP);
	dValue += dLog;
	char *pLine = new char[256];
	sprintf(pLine,"%.1lf",dValue);
	if(strstr(pLine,"-1.$"))	{	
		dValue = -5999.0;
	}
	if(dValue < -5999.0)	{	
		dValue = -5999.0;
	}
	delete pLine;
	return dValue;
}

/*
 * get_peptide_count returns the total number of peptides that have been scored
 */
size_t mprocess::get_peptide_count()
{
	return m_tPeptideCount;
}

/*
 * get_protein_count returns the total number of proteins that have been scored
 */
size_t mprocess::get_protein_count()
{
	return m_tProteinCount;
}

/*
 * get_reversed returns the number of significant reverse peptide sequences
 */
long mprocess::get_reversed()
{
	return m_lReversed;
}

/*
 * get_thread returns the thread number for the object
 */
unsigned long mprocess::get_thread()
{
	return m_lThread;
}
/*
 * get_threads returns the total number of threads currently in use
 */
unsigned long mprocess::get_threads()
{
	return m_lThreads;
}

/*
 * get_threshold returns the current value of the expectation value threshold
 */
double mprocess::get_threshold()
{
	return m_dThreshold;
}
/*
 * get_total_residues returns the total number of residues that have been processed
 */
size_t mprocess::get_total_residues()
{
	return m_tTotalResidues;
}

/*
 * get_valid returns the number of unique models
 */
size_t mprocess::get_unique()
{
	return m_tUnique;
}

/*
 * get_valid returns the number of valid models
 */
size_t mprocess::get_valid()
{
	return m_tValid;
}
/*
 * load takes a path name to the input XML parameter file and uses that file name
 * to initialize an XmlParameters object
 */

// rTANDEM : The function load is overloaded (no pun intended) to load the data from SEXP variable 
// 	     instead of using the path to the input file
bool mprocess::load(SEXP param, SEXP taxonomy, SEXP saps, SEXP mods, SEXP spectrum, mprocess *_p) {
	// Load the parameters
	dataLoader::convertSEXPToMap(param, &m_xmlValues.m_mapParam);

	// Load the taxonomy
	dataLoader::convertSEXPToVector(taxonomy, &(m_svrSequences.m_vstrFasta));
	dataLoader::convertSEXPToDeque(taxonomy, &(m_svrSequences.m_dstrFasta));

	// Load the scoring object
	bool bReturn = true;
	if (bReturn) {
		m_pScore = mscoremanager::create_mscore(m_xmlValues);
		if (m_pScore != NULL) {
			bReturn = m_pScore->load_param(m_xmlValues);
		} 
		else {
			bReturn = false;
		}
	}

	// Load parameters
	if (bReturn) {
		bReturn = (m_specCondition.load(m_xmlValues));
	}

	// obtain the tandem MS spectra to analyze
	string strValue;
	if(bReturn)	{
		bReturn = spectra();			
		string strKey = "spectrum, check all charges";
		m_xmlValues.get(strKey,strValue);
		if(bReturn && strValue == "yes" && (m_lThread == 0 || m_lThread == 0xFFFFFFFF))	{
			charge();
//			cout << "#";
			Rprintf("#");
		}
	}
	if(bReturn)	{
//		bReturn = load_saps(saps, _p);
		bReturn = load_saps(_p);
	}
	if(bReturn)	{
		bReturn = load_annotation(_p);
	}
/*
 * load the msequenceutilities object in the m_pScore member class with the amino acid
 * modification information
 */
	if(bReturn)	{
		bReturn = modify();
	}
	return bReturn;
}

bool mprocess::load(const char *_f,mprocess *_p)
{
/*
 * check the string
 */
	if(_f == NULL)
		return false;
	string strFile = _f;
/*
 * load the m_xmlValues object
 */
	bool bReturn = m_xmlValues.load(strFile);
	if(!bReturn)	{
//		cout << "The input parameter file \"" << strFile.c_str() << "\" could not be located.\nCheck the file path name and try again.\n";
		Rprintf("The input parameter file \"%s\" could not be located.\nCheck the file path name and try again.\n", strFile.c_str());
		return false;
	}
/*
 * check for the specification of a default parameter list
 */
	string strValue;
	string strKey = "list path, default parameters";
	if(m_xmlValues.get(strKey,strValue))	{
/*
 * if there is a default parameter list, load it and then reload the input list
 * the input list will over ride all settings in the default list
 */
		m_xmlValues.load(strValue);
		m_xmlValues.load(strFile);
		strKey = "list path, default parameters";
		m_xmlValues.get(strKey,strValue);
	}
/*
 * if a parameter list was found, load the msequenceServer object with the taxonomy information
 */
	if(bReturn)	{
		bReturn = taxonomy();			
	}
/*
 * if the msequenceServer object was loaded, create the scoring object
 */
	if (bReturn) {
		m_pScore = mscoremanager::create_mscore(m_xmlValues);
		if (m_pScore != NULL) {
			bReturn = m_pScore->load_param(m_xmlValues);
		} 
		else {
			bReturn = false;
		}
	}
/*
 * if the scoring object was loaded, load parameters
 */
	if (bReturn) {
		bReturn = (m_specCondition.load(m_xmlValues));
	}
/*
 * if  the msequenceServer object was loaded, obtain the tandem MS spectra to analyze
 */
	if(bReturn)	{
		bReturn = spectra();			
		strKey = "spectrum, check all charges";
		m_xmlValues.get(strKey,strValue);
		if(bReturn && strValue == "yes" && (m_lThread == 0 || m_lThread == 0xFFFFFFFF))	{
			charge();
//			cout << "#";
			Rprintf("#");
		}
	}
	if(bReturn)	{
		bReturn = load_saps(_p);
	}
	if(bReturn)	{
		bReturn = load_annotation(_p);
	}
/*
 * load the msequenceutilities object in the m_pScore member class with the amino acid
 * modification information
 */
	if(bReturn)	{
		bReturn = modify();
	}
	return bReturn;
}
/*
 * charge adds additional charge states to the vector of spectra, generating
 * +1, +2 and +3 charge states
 */
bool mprocess::charge(void)
{
	size_t a = 0;
	size_t tLength = m_vSpectra.size();
	while(a < tLength)	{
		if(m_vSpectra[a].m_tId > 100000000)	{
			return true;
		}
		a++;
	}
	a = 0;
	size_t tTest = 0;
	int iZ = 1;
	double dProton = 1.007276;
	double dMH = 0.0;
	while(a < tLength)	{
		tTest = m_vSpectra[a].m_tId + 100000000;
		iZ = (int)(m_vSpectra[a].m_fZ+0.5);
		if(iZ == 2)	{
			m_vSpectra.push_back(m_vSpectra[a]);
			m_vSpectra.back().m_fZ = 3.0;
			dMH = dProton + ((m_vSpectra[a].m_dMH - dProton)/m_vSpectra[a].m_fZ);
			m_vSpectra.back().m_dMH = dProton + ((dMH - dProton)*m_vSpectra.back().m_fZ);
			m_vSpectra.back().m_tId = tTest;
			if(m_vSpectra.back().m_dMH > 4500.0)	{
				m_vSpectra.pop_back();
			}
			m_vSpectra.push_back(m_vSpectra[a]);
			m_vSpectra.back().m_fZ = 1.0;
			dMH = dProton + ((m_vSpectra[a].m_dMH - dProton)/m_vSpectra[a].m_fZ);
			m_vSpectra.back().m_dMH = dProton + ((dMH - dProton)*m_vSpectra.back().m_fZ);
			m_vSpectra.back().m_tId = tTest + 100000000;
			if(m_vSpectra.back().m_dMH > 4500.0)	{
				m_vSpectra.pop_back();
			}
		}
		else if(iZ == 3)	{
			m_vSpectra.push_back(m_vSpectra[a]);
			m_vSpectra.back().m_fZ = 2.0;
			dMH = dProton + ((m_vSpectra[a].m_dMH - dProton)/m_vSpectra[a].m_fZ);
			m_vSpectra.back().m_dMH = dProton + ((dMH - dProton)*m_vSpectra.back().m_fZ);
			m_vSpectra.back().m_tId = tTest;
			if(m_vSpectra.back().m_dMH > 4500.0)	{
				m_vSpectra.pop_back();
			}
			m_vSpectra.push_back(m_vSpectra[a]);
			m_vSpectra.back().m_fZ = 1.0;
			dMH = dProton + ((m_vSpectra[a].m_dMH - dProton)/m_vSpectra[a].m_fZ);
			m_vSpectra.back().m_dMH = dProton + ((dMH - dProton)*m_vSpectra.back().m_fZ);
			m_vSpectra.back().m_tId = tTest + 100000000;
			if(m_vSpectra.back().m_dMH > 4500.0)	{
				m_vSpectra.pop_back();
			}
		}
		else if(iZ == 1)	{
			m_vSpectra.push_back(m_vSpectra[a]);
			m_vSpectra.back().m_fZ = 2.0;
			dMH = dProton + ((m_vSpectra[a].m_dMH - dProton)/m_vSpectra[a].m_fZ);
			m_vSpectra.back().m_dMH = dProton + ((dMH - dProton)*m_vSpectra.back().m_fZ);
			m_vSpectra.back().m_tId = tTest;
			if(m_vSpectra.back().m_dMH > 4500.0)	{
				m_vSpectra.pop_back();
			}
			m_vSpectra.push_back(m_vSpectra[a]);
			m_vSpectra.back().m_fZ = 3.0;
			dMH = dProton + ((m_vSpectra[a].m_dMH - dProton)/m_vSpectra[a].m_fZ);
			m_vSpectra.back().m_dMH = dProton + ((dMH - dProton)*m_vSpectra.back().m_fZ);
			m_vSpectra.back().m_tId = tTest + 100000000;
			if(m_vSpectra.back().m_dMH > 4500.0)	{
				m_vSpectra.pop_back();
			}
		}
		a++;
	}
	return true;
}

/*
 *  load_best_vector loads the m_vseqBest vector with a list of sequences
 * that correspond to assigned, valid peptide models. spectra that have
 * produced valid models are marked as inactive, so that they are
 * not reassigned the same peptides again.
 */
bool mprocess::load_best_vector(void)
{
	string strKey = "refine, maximum valid expectation value";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	double dMaxExpect = 0.01;
	if(strValue.size() > 0)	{
		dMaxExpect = atof(strValue.c_str());
	}
	size_t a = 0;
	while(a < m_vSpectra.size())	{
		m_vSpectra[a].m_hHyper.model();
		m_vSpectra[a].m_hHyper.set_protein_factor(1.0);
		a++;
	}
	a = 0;
	double dExpect = 1.0;
	while(a < m_vSpectra.size())	{
		dExpect = (double)m_vSpectra[a].m_hHyper.expect_protein(m_pScore->hconvert(m_vSpectra[a].m_fHyper));
		if(dExpect <= dMaxExpect)	{
			m_vSpectra[a].m_bActive = false;
		}
		a++;
	}
	return !m_vseqBest.empty();
}
/*
 * mark_repeats is used to determine which models are simply repeats of each other. finding
 * the same model more than once does not help with confirming the model: a stochastic match
 * to the same pattern twice is actually quite likely if the pattern is repeated in
 * subsequent spectra.
 */
bool mprocess::mark_repeats()
{
	size_t a = 0;
	size_t b = 0;
	string strValue;
	string strCurrent;
	size_t tStart = 0;
	size_t tEnd = 0;
	size_t tUid = 0;
	size_t tUidB = 0;
	size_t tStartB = 0;
	size_t tEndB = 0;
	size_t tLength = m_vSpectra.size();
	double dBestExpect = 0.0;
	SEQMAP::iterator itValue;
	size_t tTics = 0;
	size_t tTicLength = (size_t)((double) tLength/5.0);
	while(a < tLength)	{
		dBestExpect = 1.0e32;
		tTics++;
		if(tTics >= tTicLength)	{
//			cout << ".";
			Rprintf(".");
			//cout.flush();
			tTics = 0;
		}
		if(!m_vSpectra[a].m_bRepeat && !m_vSpectra[a].m_vseqBest.empty())	{
			tStart = m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_lS;
			tEnd = m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_lE;
			tUid = m_vSpectra[a].m_vseqBest[0].m_tUid;
			dBestExpect = m_vSpectra[a].m_dExpect;
		}
		if(!m_vSpectra[a].m_bRepeat && !m_vSpectra[a].m_vseqBest.empty())	{
			b = a + 1;
			while(b < tLength)	{
				if(!m_vSpectra[b].m_bRepeat && !m_vSpectra[b].m_vseqBest.empty())	{
					tStartB = m_vSpectra[b].m_vseqBest[0].m_vDomains[0].m_lS;
					tEndB = m_vSpectra[b].m_vseqBest[0].m_vDomains[0].m_lE;
					tUidB = m_vSpectra[b].m_vseqBest[0].m_tUid;
					if(tEndB == tEnd && tStartB == tStart && tUidB == tUid)	{
						if(dBestExpect <= m_vSpectra[b].m_dExpect)	{
							m_vSpectra[b].m_bRepeat = true;
						}
						else	{
							dBestExpect = m_vSpectra[b].m_dExpect;
						}
					}
				}
				b++;
			}
			if(dBestExpect < m_vSpectra[a].m_dExpect)	{
				m_vSpectra[a].m_bRepeat = true;
			}
		}
		a++;
	}
	return true;
}
/*
 * merge_map adds new values to the existing m_mapSequences map.
*/
bool mprocess::merge_map(SEQMAP &_s)
{
	SEQMAP::iterator itValue = _s.begin();
	SEQMAP::iterator itEnd = _s.end();
	while(itValue != itEnd)	{
		if(m_mapSequences.find(itValue->first) == m_mapSequences.end())	{
			m_mapSequences.insert(*itValue);
		}
		itValue++;
	}
	return true;
}

/*
 * merge_spectra takes an mspectrum vector from an external source and merges it
 * with the m_vSpectra vector. this method is used to combine information obtained
 * from other threads with this mprocess object
 */
bool mprocess::merge_spectra()
{
	string strKey = "refine, maximum valid expectation value";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	double dMaxExpect = 0.01;
	if(strValue.size() > 0)	{
		dMaxExpect = atof(strValue.c_str());
	}
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	{
		while(a < m_vSpectra.size())	{
			m_vSpectra[a].m_hHyper.model();
			m_vSpectra[a].m_hHyper.set_protein_factor(1.0);
//			if(m_bUseHomologManagement && m_vSpectra[a].m_vseqBest.size() > 5)	{
//				m_vSpectra[a].m_vseqBest.erase(m_vSpectra[a].m_vseqBest.begin()+5,m_vSpectra[a].m_vseqBest.end());
//			}
			a++;
		}
		a = 0;
		double dExpect = 1.0;
		SEQMAP::iterator itValue;
		while(a < m_vSpectra.size())	{
			b = 0;
			dExpect = (double)m_vSpectra[a].m_hHyper.expect_protein(m_pScore->hconvert(m_vSpectra[a].m_fHyper));
			if(dExpect <= dMaxExpect)	{
				m_vSpectra[a].m_bActive = false;
				while(b < m_vSpectra[a].m_vseqBest.size())	{
					c = 0;
					while(c < m_vseqBest.size())	{
						if(m_vSpectra[a].m_vseqBest[b].m_tUid == m_vseqBest[c].m_tUid)	{
							break;
						}
						c++;
					}
					if(c == m_vseqBest.size())	{
						m_vseqBest.push_back(m_vSpectra[a].m_vseqBest[b]);
						itValue = m_mapSequences.find(m_vseqBest[c].m_tUid);
						m_vseqBest[c].m_strSeq = ((*itValue).second.c_str());
						m_vseqBest[c].m_vDomains.clear();
					}
					b++;
					if(m_bUseHomologManagement && b >= 5)	{
						break;
					}
				}
			}
			a++;
		}
	}
	return true;
}

bool mprocess::merge_spectra(vector<mspectrum> &_s)
{
	string strKey = "refine, maximum valid expectation value";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	double dMaxExpect = 0.01;
	if(strValue.size() > 0)	{
		dMaxExpect = atof(strValue.c_str());
	}
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	{
		while(a < _s.size())	{
			_s[a].m_hHyper.model();
			_s[a].m_hHyper.set_protein_factor(1.0);
			if(m_bUseHomologManagement && _s[a].m_vseqBest.size() > 5)	{
				_s[a].m_vseqBest.erase(_s[a].m_vseqBest.begin()+5,_s[a].m_vseqBest.end());
			}
			a++;
		}
		a = 0;
		double dExpect = 1.0;
		SEQMAP::iterator itValue;
		while(a < _s.size())	{
			b = 0;
			dExpect = (double)_s[a].m_hHyper.expect_protein(m_pScore->hconvert(_s[a].m_fHyper));
			if(dExpect <= dMaxExpect)	{
				while(b < _s[a].m_vseqBest.size())	{
					c = 0;
					while(c < m_vseqBest.size())	{
						if(_s[a].m_vseqBest[b].m_tUid == m_vseqBest[c].m_tUid)	{
							break;
						}
						c++;
					}
					if(c == m_vseqBest.size())	{
						m_vseqBest.push_back(_s[a].m_vseqBest[b]);
						itValue = m_mapSequences.find(m_vseqBest[c].m_tUid);
						m_vseqBest[c].m_strSeq = ((*itValue).second.c_str());
						m_vseqBest[c].m_vDomains.clear();
					}
					b++;
				}
			}
			a++;
		}
	}
	return true;
}
/*
 * merge_statistics updates the processing statistics with external statistics.
 * this method is used to combine information obtained from other threads with 
 * this mprocess object
 */
bool mprocess::merge_statistics(const mprocess *_p)
{ 
	m_tPeptideCount += _p->m_tPeptideCount;
	m_tRefineInput += _p->m_tRefineInput;
	m_tRefinePartial += _p->m_tRefinePartial;
	m_tRefineUnanticipated += _p->m_tRefineUnanticipated;
	m_tRefineNterminal += _p->m_tRefineNterminal;
	m_tRefineCterminal += _p->m_tRefineCterminal;
	m_tRefinePam += _p->m_tRefinePam;
	return true;
}
/*
 * modify checks the input parameters for known parameters that are use to modify 
 * a protein sequence. these parameters are stored in the m_pScore member object's
 * msequenceutilities member object
 */
bool mprocess::modify()
{
	string strKey = "residue, modification mass";
	string strValue;
	m_vstrModifications.clear();
	if(m_xmlValues.get(strKey,strValue) && strValue.size() > 0) {
		m_vstrModifications.push_back(strValue);
	}
	else	{
		strValue = "";
		m_vstrModifications.push_back(strValue);
	};
	int a = 1;
	char *pLine = new char[256];
	sprintf(pLine,"residue, modification mass %i",(int)a);
	strKey = pLine;
	while(m_xmlValues.get(strKey,strValue) && strValue.size() > 0) {
		m_vstrModifications.push_back(strValue);
		a++;
		sprintf(pLine,"residue, modification mass %i",(int)a);
		strKey = pLine;
	}
	delete pLine;
	strKey = "residue, potential modification mass";
	if(m_xmlValues.get(strKey,strValue)) {
		m_pScore->m_seqUtil.modify_maybe(strValue);
		m_pScore->m_seqUtilAvg.modify_maybe(strValue);
	}
	strKey = "residue, potential modification motif";
	if(m_xmlValues.get(strKey,strValue)) {
		m_pScore->m_seqUtil.modify_motif(strValue);
		m_pScore->m_seqUtilAvg.modify_motif(strValue);
	}
	strKey = "protein, N-terminal residue modification mass";
	if(m_xmlValues.get(strKey,strValue)) {
		m_pScore->m_seqUtil.modify_n((float)atof(strValue.c_str()));
		m_pScore->m_seqUtilAvg.modify_n((float)atof(strValue.c_str()));
	}
	strKey = "protein, C-terminal residue modification mass";
	if(m_xmlValues.get(strKey,strValue)) {
		m_pScore->m_seqUtil.modify_c((float)atof(strValue.c_str()));
		m_pScore->m_seqUtilAvg.modify_c((float)atof(strValue.c_str()));
	}
	strKey = "protein, cleavage N-terminal mass change";
	if(m_xmlValues.get(strKey,strValue)) {
		m_pScore->m_seqUtil.m_dCleaveN = atof(strValue.c_str());
		m_pScore->m_seqUtilAvg.m_dCleaveN = atof(strValue.c_str());
	}
	strKey = "protein, cleavage C-terminal mass change";
	if(m_xmlValues.get(strKey,strValue)) {
		m_pScore->m_seqUtil.m_dCleaveC = atof(strValue.c_str());
		m_pScore->m_seqUtilAvg.m_dCleaveC = atof(strValue.c_str());
	}
	strKey = "residue, NG deamidation";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes") {
		m_bCheckNg = true;
	}
	return true;
}
/*
 * process carries out the protein identification
 */
bool mprocess::process(void)
{
	if(m_vSpectra.size() < 1)
		return false;
	string strKey;
	string strValue;
//	m_pScore->set_mini(true);
#ifdef PLUGGABLE_SCORING
	strKey = "scoring, pluggable scoring";
	strValue = "yes";
	m_xmlValues.set(strKey,strValue);
#else
	strKey = "scoring, pluggable scoring";
	strValue = "no";
	m_xmlValues.set(strKey,strValue);
#endif

	strKey = "output, path";
	m_xmlValues.get(strKey,strValue);
	strKey = "output path: ";
	strKey += strValue;
	m_prcLog.log(strKey);
	strKey = "spectrum, path";
	m_xmlValues.get(strKey,strValue);
	strKey = "input path: ";
	strKey += strValue;
	m_prcLog.log(strKey);

#ifndef X_P3
	strKey = "protein, saps";
	m_bSaps = false;
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		m_bSaps = true;
	}
#endif

	strKey = "protein, homolog management";
	m_xmlValues.get(strKey,strValue);
	m_bUseHomologManagement = false;
	if(strValue == "yes")	{
		m_bUseHomologManagement = true;
	}
	strKey = "scoring, cyclic permutation";
	m_xmlValues.get(strKey,strValue);
	m_bCrcCheck = false;
	if(strValue == "yes")	{
		m_bCrcCheck = true;
	}
/*
 * Detect the presence of ion type parameters: default is b + y ions
 */
	unsigned long lType = 0;
	strKey = "scoring, include reverse";
	m_xmlValues.get(strKey,strValue);
	m_lReversed = -1;
	m_bReversedOnly = false;
	if(strValue == "yes")	{
		m_lReversed = 0;
	}
	else if(strValue == "only")		{
		m_bReversedOnly = true;
		m_lReversed = 0;
	}
	strKey = "scoring, a ions";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		lType |= mscore::T_A;
	}
	strKey = "scoring, a ions";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		lType |= mscore::T_A;
	}
	strKey = "scoring, b ions";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		lType |= mscore::T_B;
	}
	strKey = "scoring, c ions";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		lType |= mscore::T_C;
	}
	strKey = "scoring, x ions";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		lType |= mscore::T_X;
	}
	strKey = "scoring, z ions";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		lType |= mscore::T_Z;
	}
	strKey = "scoring, y ions";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		lType |= mscore::T_Y;
	}
	if(lType == 0)	{
		lType = mscore::T_B | mscore::T_Y;
	}
	m_pScore->set_type(lType);
	strKey = "refine, spectrum synthesis";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		m_pScore->m_seqUtil.synthesis(true);
	}
	else	{
		m_pScore->m_seqUtil.synthesis(false);
	}
/*
 * Determine how to process the parent and fragment error numbers:
 * 1. absolute errors in Daltons; or
 * 2. relative errors in ppm.
 */
	strKey = "spectrum, parent monoisotopic mass error units";
	m_xmlValues.get(strKey,strValue);
	lType = 0;
	m_errValues.m_bPpm = false;
	if(strValue == "Daltons")	{
		lType |= mscore::T_PARENT_DALTONS;
		m_errValues.m_bPpm = false;
	}
	else if(strValue == "ppm")	{
		lType |= mscore::T_PARENT_PPM;
		m_errValues.m_bPpm = true;
	}
	strKey = "spectrum, fragment mass error units";
	m_xmlValues.get(strKey,strValue);
	if (strValue.empty()) {
		// try old parameter name
		strKey = "spectrum, fragment monoisotopic mass error units";
		m_xmlValues.get(strKey,strValue);
	}
	if(strValue == "Daltons")	{
		lType |= mscore::T_FRAGMENT_DALTONS;
	}
	else if(strValue == "ppm")	{
		lType |= mscore::T_FRAGMENT_PPM;
	}
	if(lType == 0)	{
		lType = mscore::T_PARENT_DALTONS | mscore::T_FRAGMENT_DALTONS;
	}
	m_pScore->set_error(lType);
/*
 * check the m_xmlValues parameters for a set of known input parameters, and
 * substitute default values if they are not found 
 */
	strKey = "spectrum, fragment mass error";
	m_xmlValues.get(strKey,strValue);
	if (strValue.empty()) {
		// try old parameter name
		strKey = "spectrum, fragment monoisotopic mass error";
		m_xmlValues.get(strKey,strValue);
	}
	float fErrorValue = (float)atof(strValue.c_str());
	if(fErrorValue <= 0.0)
		fErrorValue = (float)0.45;
	m_pScore->set_fragment_error(fErrorValue);
	strKey = "spectrum, parent monoisotopic mass error plus";
	m_xmlValues.get(strKey,strValue);
	fErrorValue = (float)atof(strValue.c_str());
	fErrorValue = (float)fabs(atof(strValue.c_str()));
	m_errValues.m_fPlus = fErrorValue;
	if(m_errValues.m_bPpm)	{
		if(fErrorValue < 95.0)	{
			m_bCrcCheck = true;
		}
		if(fErrorValue < 10.0)	{
			fErrorValue = 10.0;
		}
	}
	else	{
		if(fErrorValue < 0.095)	{
			m_bCrcCheck = true;
		}
		if(fErrorValue < 0.01)	{
			fErrorValue = (float)0.01;
		}
	}
	m_pScore->set_parent_error(fErrorValue,true);
	strKey = "spectrum, homology error";
	m_xmlValues.get(strKey,strValue);
	fErrorValue = (float)atof(strValue.c_str());
	if(fErrorValue <= 0.0)
		fErrorValue = (float)4.5;
	m_pScore->set_homo_error(fErrorValue);
	strKey = "spectrum, parent monoisotopic mass error minus";
	m_xmlValues.get(strKey,strValue);
	fErrorValue = (float)fabs(atof(strValue.c_str()));
	m_errValues.m_fMinus = (float)(-1.0*fErrorValue);
	if(m_errValues.m_bPpm)	{
		if(fErrorValue < 95.0)	{
			m_bCrcCheck = true;
		}
		if(fErrorValue < 10.0)	{
			fErrorValue = 10.0;
		}
	}
	else	{
		if(fErrorValue < 0.095)	{
			m_bCrcCheck = true;
		}
		if(fErrorValue < 0.01)	{
			fErrorValue = (float)0.01;
		}
	}
	m_pScore->set_parent_error(fErrorValue,false);
	strKey = "spectrum, parent monoisotopic mass isotope error";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		m_errValues.m_bIsotope = true;
		m_pScore->set_isotope_error(true);
	}
	else	{
		m_pScore->set_isotope_error(false);
		m_errValues.m_bIsotope = false;
	}
	strKey = "protein, cleavage N-terminal limit";
	if(m_xmlValues.get(strKey,strValue)) {
		if(atoi(strValue.c_str()) > 0)	{
			m_lStartMax = atoi(strValue.c_str());
		}
	}

	strKey = "protein, quick acetyl";
	if(m_xmlValues.get(strKey,strValue)) {
		if(strValue == "no")	{
			m_bQuickAcetyl = false;
		}
	}
	strKey = "protein, quick pyrolidone";
	if(m_xmlValues.get(strKey,strValue)) {
		if(strValue == "no")	{
			m_bQuickPyro = false;
		}
	}
	strKey = "protein, stP bias";
	bool bA = true;
	if(m_xmlValues.get(strKey,strValue)) {
		if(strValue == "no")	{
			bA = false;
		}
	}
	m_pScore->set_phospho_bias(bA);

	strKey = "protein, modified residue mass file";
	if(m_xmlValues.get(strKey,strValue))	{
		m_pScore->m_seqUtil.set_aa_file(strValue);
		m_pScore->m_seqUtilAvg.set_aa_file(strValue);
	}
	strKey = "protein, cleavage site";
	m_xmlValues.get(strKey,strValue);
	m_Cleave.load(strValue);
	strKey = "protein, cleavage semi";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		m_semiState.activate(true);
	}
	strKey = "scoring, minimum ion count";
	m_xmlValues.get(strKey,strValue);
	m_lIonCount = (unsigned long)atoi(strValue.c_str()) - 1;
	strKey = "scoring, maximum missed cleavage sites";
	m_xmlValues.get(strKey,strValue);		
	m_tMissedCleaves = atoi(strValue.c_str());
	if(m_Cleave.m_lType == 0x01 && m_tMissedCleaves < 6)	{
		m_tMissedCleaves = 50;
	}
	strKey = "spectrum, sequence batch size";
	m_xmlValues.get(strKey,strValue);
	size_t tBatch = atoi(strValue.c_str());
	if(tBatch < 1)	{
		tBatch = 1000;
	}
	m_svrSequences.initialize(tBatch);
	strKey = "protein, use annotations";
	m_xmlValues.get(strKey,strValue);
	m_bAnnotation = false;
	m_strLastMods.clear();
	m_bMinimalAnnotation = false;
	if(strValue == "yes")	{
		m_bAnnotation = true;
	}
	else	{
		strKey = "protein, use minimal annotations";
		m_xmlValues.get(strKey,strValue);
		if(strValue == "no")	{
			m_bMinimalAnnotation = false;
			m_xmlValues.set(strKey,strValue);
		}
		else	{
			m_bMinimalAnnotation = true;
			strValue = "yes";
			m_xmlValues.set(strKey,strValue);
		}
	}
	size_t a = 0;
	m_tProteinCount = 0;
	m_tPeptideCount = 0;
	m_tPeptideScoredCount = 0;
	strKey = "output, http";
	m_xmlValues.get(strKey,strValue);
	a = 0;
	long lTics = 0;
/*
 * record the spectrum m/z - intensity pairs into the mscore object. This
 * is done here to speed up processing many spectra, as they are only
 * loaded and processed once.
 */
	m_tSpectra = m_vSpectra.size();
	while(a < m_tSpectra)	{
		m_pScore->add_details(m_vSpectra[a]);
		a++;
	}
	m_pScore->sort_details();
	a =0;
	while(a < m_tSpectra)	{
		m_pScore->add_mi(m_vSpectra[a]);
		a++;
	}
	removeMI();
	a = 0;
/*
 * estimate the smallest number of residues that could possibily correspond to the
 * smallest parent ion in the m_vSpectra vector
 */
	residues();
/*
 * record the start time for the identification process. 
 */
	m_dSearchTime = clock();
	unsigned long lServerCount = 0;
	char *pOut = new char[256];
	sprintf(pOut,"Spectrum-to-sequence matching process in progress");
	strKey = "output, message";
	m_xmlValues.get(strKey,strValue);
	if(strValue.size() > 0)	{
		if(strValue.size() > 255)	{
			delete pOut;
			pOut = new char[strValue.size()+1];
		}
		sprintf(pOut,"%s",strValue.c_str());
	}
	long lOut = 0;
	long lOutLimit = (long)strlen(pOut);
/*
 * use the msequenceServer to get a list of proteins to analyze 
 */
	while(!m_svrSequences.done())	{
		m_svrSequences.next(true);
		lServerCount++;
		score_each_sequence();
/*
 * this section produces the "still alive" messages that are displayed so
 * that a user doesn't get anxious. this implementaion only throws messages
 * from the 0th thread, and uses a character string
 * to generate characters that get flushed to the console.
 */
		if(m_lThread == 0 || m_lThread == 0xFFFFFFFF)	{
			lTics++;
			if(lTics == 50)	{
//				cout << " | " << (unsigned long)m_tProteinCount/1000 << " ks \n";
				Rprintf(" | %l ks \n", (unsigned long)m_tProteinCount/1000l);
				//cout.flush();
				m_prcLog.log(".");
				lTics = 0;
			}
			else	{
				if(lTics == 1)	{
//					cout << "\t";
					Rprintf("\t");
				}
//				cout << pOut[lOut];
				Rprintf("%c", pOut[lOut]);
				//cout.flush();
				m_prcLog.log(".");
				lOut++;
				if(lOut >= lOutLimit)	{
					lOut = 0;
				}
			}
		}
		clean_sequences();
	}
	delete pOut;
/*
 * record the total protein modelling session time 
 */
	m_dSearchTime = (clock() - m_dSearchTime)/(double)CLOCKS_PER_SEC;
	m_vstrPaths = m_svrSequences.m_vstrPaths;
	m_svrSequences.clear();
/*
 * process the scoring histograms in each spectrum, so that expectation values
 * can be calculated. first a survival function is calculated, replacing the
 * original scoring histogram. then, the high scoring tail of the survival function
 * is modeled. 
 * NOTE: only the hyper score histogram is used in the modelling process.
 * NOTE: the protein factor is not used to weigth scores - it is present for future use
 */
	return true;
}

/*
 * refine controls the sequence refinement process.
 */
bool mprocess::refine(void)
{
	size_t tTime = clock();
	m_pScore->set_mini(false);
	// Check for output path for sequences to be loaded to a bioml file
	// Writing only done on the base thread
	string strKey = "output, sequence path";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	if(!strValue.empty() && (m_lThread == 0 || m_lThread == 0xFFFFFFFF))	{
		mbiomlreport rptCurrent;
		rptCurrent.setpath(strValue);
		rptCurrent.write(m_vseqBest,m_vstrPaths);
	}
	strKey = "refine";
	m_xmlValues.get(strKey,strValue);
	bool bReturn = false;
	m_lStartMax = 100000000;
	if(strValue == "yes")	{
		bReturn = refine_model();
		tTime = clock() - tTime;
		m_dRefineTime = (double)tTime/(double)(CLOCKS_PER_SEC);
	}
	return bReturn;
}
/*
 * refine_models takes the model peptides that have been found and tries to
 * find other spectra that can be described by the proteins that these models
 * come from. several layers of modelling are available in the current implementation:
 * 1. comparison with multiple potential modificiation
 * 2. full [X]|[X] cleavage
 * 3. N-terminal modifications
 * 4. C-terminal modifications
 */
bool mprocess::refine_model()
{
	m_pRefine = mrefinemanager::create_mrefine(m_xmlValues);
	if (m_pRefine == NULL) {
//		cout << "Failed to create mrefine\n";
		Rprintf("Failed to create mrefine\n");
		return false;
	}
	m_pRefine->set_mprocess(this);
	m_pRefine->refine();
	return true;
}
/*
 * report outputs the information obtained during the process method, using
 * an mreport object and the reporting parameter settings obtained from the
 * input parameter list. many applications may need to customize this method.
 * two customized output types are given: one which emphasizes protein sequences
 * and one which emphasizes spectra.
 */
bool mprocess::report(void)
{
	m_prcLog.log("creating report");
	vector<mspectrum>::iterator itS = m_vSpectra.begin();
	m_pScore->clear();
	while(itS != m_vSpectra.end())	{
		itS->m_hHyper.model();
		itS->m_hHyper.set_protein_factor(1.0);
		itS++;
	}
	size_t a = 0;
	map<size_t,size_t> mapSpec;
	pair<size_t,size_t> prSpec;
	map<size_t,size_t>::iterator itMap;
	while(a < m_vSpectra.size())	{
		prSpec.first = m_vSpectra[a].m_tId;
		prSpec.second = a;
		mapSpec.insert(prSpec);
		if(!m_vSpectra[a].m_vseqBest.empty() && !m_vSpectra[a].m_vseqBest[0].m_vDomains.empty())	{
			m_vSpectra[a].m_dExpect = m_vSpectra[a].m_hHyper.expect(m_pScore->hconvert(m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_fHyper));
		}
		a++;
	}
	size_t tTest = 0;
	itS = m_vSpectra.begin();
	while(itS != m_vSpectra.end())	{
		if(itS->m_tId < 100000000)	{
			tTest = itS->m_tId + 100000000;
			itMap = mapSpec.find(tTest);
			if(itMap != mapSpec.end())	{
				a = itMap->second;
				if(itS->m_dExpect <= m_vSpectra[a].m_dExpect)	{
					m_vSpectra[a].m_dExpect = 1000;
					m_vSpectra[a].m_fScore = 0.0;
					mapSpec.erase(tTest);
					tTest = tTest + 100000000;
					itMap = mapSpec.find(tTest);
					if(itMap != mapSpec.end())	{
						a = itMap->second;
						if(itS->m_dExpect <= m_vSpectra[a].m_dExpect)	{
							m_vSpectra[a].m_dExpect = 1000;
							m_vSpectra[a].m_fScore = 0.0;
							mapSpec.erase(tTest);
						}
						else	{
							itS->m_dExpect = 1000;
							itS->m_fScore = 0.0;
							mapSpec.erase(tTest);
						}
					}
				}
				else	{
					itS->m_dExpect = 1000;
					itS->m_fScore = 0.0;
					mapSpec.erase(tTest);
					tTest = tTest + 100000000;
					itMap = mapSpec.find(tTest);
					if(itMap != mapSpec.end())	{
						size_t b = itMap->second;
						if(m_vSpectra[a].m_dExpect <= m_vSpectra[b].m_dExpect)	{
							m_vSpectra[b].m_dExpect = 1000;
							m_vSpectra[b].m_fScore = 0.0;
							mapSpec.erase(tTest);
						}
						else	{
							m_vSpectra[a].m_dExpect = 1000;
							m_vSpectra[a].m_fScore = 0.0;
							mapSpec.erase(tTest);
						}
					}
				}
			}
		}
		itS++;
	}
	mapSpec.clear();
	a = 0;
	itS = m_vSpectra.begin();
	m_viQuality.clear();
	while(a < 255)	{
		m_viQuality.push_back(0);
		a++;
	}
	a = 0;
	int iIndex;
	int iQ = 0;
	while(itS != m_vSpectra.end())	{
		a = 0;
		if(itS->m_fScore > 0.0 && itS->m_dExpect < 1.0)	{
			iIndex = int(-1.0*log(itS->m_dExpect));
			if(iIndex >= 0 && iIndex < 255)	{
				m_viQuality[iIndex]++;
				iQ++;
			}
		}
		while(a < itS->m_vMINeutral.size())	{
			itS->m_vMI.push_back(itS->m_vMINeutral[a]);
			a++;
		}
		if(a > 0)	{
			sort(itS->m_vMI.begin(),itS->m_vMI.end(),lessThanMass);
		}
		itS++;
	}
	char *pLine = new char[256];
	string strKey = "quality values";
	string strValue;
	a = 0;
	while(a < 20)	{
		if(a != 0)	{
			strValue += " ";
		}
		sprintf(pLine,"%i",(int)m_viQuality[a]);
		strValue += pLine;
		a++;
	}
	m_xmlPerformance.set(strKey,strValue);
/*
 * store information in the m_smlPerformance object for reporting
 */
	strKey = "timing, initial modelling total (sec)";
	sprintf(pLine,"%.2lf",m_dSearchTime);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	strKey = "timing, initial modelling/spectrum (sec)";
	sprintf(pLine,"%.4lf",m_dSearchTime/(double)m_tSpectraTotal);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	strKey = "timing, load sequence models (sec)";
	sprintf(pLine,"%.2lf",m_svrSequences.get_time());
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	strKey = "modelling, total spectra used";
	sprintf(pLine,"%lu",(unsigned long)m_tSpectraTotal);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	strKey = "modelling, total proteins used";
	sprintf(pLine,"%lu",(unsigned long)m_tProteinCount);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	strKey = "modelling, total peptides used";
	sprintf(pLine,"%lu",(unsigned long)m_tPeptideCount);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	if(m_specCondition.get_noise_suppression())	{
		strKey = "modelling, spectrum noise suppression ratio";
		sprintf(pLine,"%.2lf",(double)(m_tSpectraTotal - m_vSpectra.size())/(double)m_tSpectraTotal);
		strValue = pLine;
		m_xmlPerformance.set(strKey,strValue);
	}
	size_t tSeq = 0;
	while(tSeq < m_svrSequences.m_vstrFasta.size())	{
		sprintf(pLine,"list path, sequence source #%i",(int)(tSeq+1));
		strKey = pLine;
		strValue = m_svrSequences.m_vstrFasta[tSeq];
		m_xmlPerformance.set(strKey,strValue);
		if(tSeq < m_svrSequences.m_vstrDesc.size())	{
			sprintf(pLine,"list path, sequence source description #%i",(int)(tSeq+1));
			strKey = pLine;
			strValue = m_svrSequences.m_vstrDesc[tSeq];
			m_xmlPerformance.set(strKey,strValue);
		}
		tSeq++;
	}
	tSeq = 0;
	while(tSeq < m_vstrSaps.size())	{
		sprintf(pLine,"list path, saps source #%i",(int)(tSeq+1));
		strKey = pLine;
		strValue = m_vstrSaps[tSeq];
		m_xmlPerformance.set(strKey,strValue);
		tSeq++;
	}
	tSeq = 0;
	while(tSeq < m_vstrMods.size())	{
		sprintf(pLine,"list path, mods source #%i",(int)(tSeq+1));
		strKey = pLine;
		strValue = m_vstrMods[tSeq];
		m_xmlPerformance.set(strKey,strValue);
		tSeq++;
	}
	strKey = "output, maximum valid expectation value";
	m_xmlValues.get(strKey,strValue);
	double dMaxExpect = 0.01;
	if(strValue.size() > 0)	{
		dMaxExpect = atof(strValue.c_str());
	}
	m_dThreshold = dMaxExpect;

	strKey = "refining, # input models";
	sprintf(pLine,"%u",(unsigned int)m_tRefineModels);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	
	strKey = "refining, # input spectra";
	sprintf(pLine,"%u",(unsigned int)(m_tRefineInput));
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);

	strKey = "refining, # partial cleavage";
	sprintf(pLine,"%u",(unsigned int)(m_tRefinePartial));
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);

	strKey = "refining, # unanticipated cleavage";
	sprintf(pLine,"%u",(unsigned int)(m_tRefineUnanticipated));
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);

	strKey = "refining, # potential N-terminii";
	sprintf(pLine,"%u",(unsigned int)(m_tRefineNterminal));
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);

	strKey = "refining, # potential C-terminii";
	sprintf(pLine,"%u",(unsigned int)(m_tRefineCterminal));
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);

	strKey = "refining, # point mutations";
	sprintf(pLine,"%u",(unsigned int)(m_tRefinePam));
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	
	strKey = "timing, refinement/spectrum (sec)";
	sprintf(pLine,"%.4lf",m_dRefineTime/(double)m_vSpectra.size());
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);

	strKey = "spectrum, use contrast angle";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		strKey = "modelling, contrast angle rejection ratio";
		sprintf(pLine,"%.2lf",(double)m_tContrasted/((double)m_vSpectra.size()+(double)m_tContrasted));
		strValue = pLine;
		m_xmlPerformance.set(strKey,strValue);
	}
/*
 * calculate the expectation values for the proteins
 */
	m_prcLog.log("calculating expectation values");

	report_expect(dMaxExpect);
/*
 * now, using the expectation values for the proteins, sort the results into proteins or spectra
 */
	m_prcLog.log("sorting peptides");
	report_sort();
/*
 * using the retrieved values, create an mreport and use it to create the output
 */

	dMaxExpect = log10(dMaxExpect);
	strKey = "output, results";
	m_xmlValues.get(strKey,strValue);
	string strResults = strValue;
//	cout << "\twriting results ";
	Rprintf("\twriting results ");
	//cout.flush();
	m_prcLog.log("writing results");
	restore();
	
	m_dEsum = 0.0;

	if(strResults != "all" && strResults != "valid" && strResults != "stochastic")
		strResults = "valid";

	if(strResults == "all")	{
		report_all();
	}
	else if(strResults == "valid")	{
		report_valid(dMaxExpect);
	}
	else if(strResults == "stochastic")	{
		report_stochastic(dMaxExpect);
	}
//	cout << "..... done.\n";
	Rprintf("..... done.\n");
	//cout.flush();
	m_prcLog.log("report complete");
	delete pLine;
	return false;
}

bool mprocess::report_all()
{
/*
 * check the input parameters for known output information
 */
	string strKey = "output, histogram column width";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	long lHistogramColumns = 30;
	if(atoi(strValue.c_str()) > 0)
		lHistogramColumns = atoi(strValue.c_str());

	strKey = "output, spectra";
	m_xmlValues.get(strKey,strValue);
	bool bSpectra = false;
	if(strValue == "yes")
		bSpectra = true;

	strKey = "output, histograms";
	m_xmlValues.get(strKey,strValue);
	bool bHistograms = false;
	if(strValue == "yes")
		bHistograms = true;

	strKey = "output, sequences";
	m_xmlValues.get(strKey,strValue);
	bool bSequences = false;
	if(strValue == "yes")
		bSequences = true;

	strKey = "output, proteins";
	m_xmlValues.get(strKey,strValue);
	bool bProteins = false;
	if(strValue == "yes")
		bProteins = true;

	strKey = "output, parameters";
	m_xmlValues.get(strKey,strValue);
	bool bInput = false;
	if(strValue == "yes")
		bInput = true;

	strKey = "output, performance";
	m_xmlValues.get(strKey,strValue);
	bool bPerf = false;
	if(strValue == "yes")
		bPerf = true;

	strKey = "output, one sequence copy";
	m_xmlValues.get(strKey,strValue);
	bool bCompress = false;
	if(strValue == "yes")
		bCompress = true;

	mreport rptValue(*m_pScore);
	rptValue.set_compression(bCompress);
	rptValue.set_columns(lHistogramColumns);
	rptValue.start(m_xmlValues);
	m_pathName = rptValue.getPathName(); // rTANDEM
	size_t a = 0;	
	size_t tLength = m_vSpectra.size();
	size_t b = 0;
	SEQMAP::iterator itValue;
	map <string,string>::iterator itMod;
	string strMods;

	while(a < tLength)	{
		b = 0;
		if(!m_vSpectra[a].m_vseqBest.empty() && !m_vSpectra[a].m_vseqBest[0].m_vDomains.empty())	{
			m_dEsum += m_vSpectra[a].m_hHyper.expect(m_pScore->hconvert(m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_fHyper));
			while(b < m_vSpectra[a].m_vseqBest.size())	{
				itValue = m_mapSequences.find(m_vSpectra[a].m_vseqBest[b].m_tUid);
				m_vSpectra[a].m_vseqBest[b].m_strSeq = (*itValue).second;
				b++;
			}
			if(bSpectra || bHistograms || bProteins)
				rptValue.group(m_vSpectra[a]);
			if(bProteins)	{
				rptValue.sequence(m_vSpectra[a],bSequences,m_vstrPaths,m_mapAnnotation);
			}
			if(bHistograms)
				rptValue.histogram(m_vSpectra[a]);
			if(bSpectra)
				rptValue.spectrum(m_vSpectra[a]);
			if(bSpectra || bHistograms || bProteins)
				rptValue.endgroup();
		}
		m_vSpectra[a].m_vseqBest.clear();
		a++;
	}
	if(bInput)
		rptValue.info(m_xmlValues);
	if(bPerf)
		rptValue.performance(m_xmlPerformance);
	//CONSIDER(bmaclean): Report average masses too, if using average for fragments?
	if(m_pScore->m_pSeqUtilFrag->is_modified())	{
		rptValue.masses(*(m_pScore->m_pSeqUtilFrag));
	}
	return rptValue.end();
}

bool mprocess::report_valid(const double _d)
{
/*
 * check the input parameters for known output information
 */
	string strKey = "output, histogram column width";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	long lHistogramColumns = 30;
	if(atoi(strValue.c_str()) > 0)
		lHistogramColumns = atoi(strValue.c_str());

	strKey = "output, spectra";
	m_xmlValues.get(strKey,strValue);
	bool bSpectra = false;
	if(strValue == "yes")
		bSpectra = true;

	strKey = "output, histograms";
	m_xmlValues.get(strKey,strValue);
	bool bHistograms = false;
	if(strValue == "yes")
		bHistograms = true;

	strKey = "output, sequences";
	m_xmlValues.get(strKey,strValue);
	bool bSequences = false;
	if(strValue == "yes")
		bSequences = true;

	strKey = "output, proteins";
	m_xmlValues.get(strKey,strValue);
	bool bProteins = false;
	if(strValue == "yes")
		bProteins = true;

	strKey = "output, parameters";
	m_xmlValues.get(strKey,strValue);
	bool bInput = false;
	if(strValue == "yes")
		bInput = true;

	strKey = "output, performance";
	m_xmlValues.get(strKey,strValue);
	bool bPerf = false;
	if(strValue == "yes")
		bPerf = true;

	strKey = "output, one sequence copy";
	m_xmlValues.get(strKey,strValue);
	bool bCompress = false;
	if(strValue == "yes")
		bCompress = true;

	mreport rptValue(*m_pScore);
	rptValue.set_compression(bCompress);
	rptValue.set_columns(lHistogramColumns);
	rptValue.start(m_xmlValues);
	m_pathName = rptValue.getPathName(); // rTANDEM
	size_t a = 0;	
	size_t tLength = m_vSpectra.size();
	size_t tActive = 0;
	double dValue = 0.0;
	size_t b = 0;
	SEQMAP::iterator itValue;
	m_tValid = 0;
	m_tUnique = 1;
	size_t tLast = 0;
	double dProteinMax = pow(10,_d);
	strKey = "output, maximum valid protein expectation value";
	m_xmlValues.get(strKey,strValue);
	if(strValue.size() > 0)	{
		dProteinMax = atof(strValue.c_str());
	}
	dProteinMax = log10(dProteinMax);
	double dProtein = 0.0;
	while(a < tLength)	{
		dValue = 3.0;
		if(m_vSpectra[a].m_fScore > 0.0 && !m_vSpectra[a].m_vseqBest.empty() && !m_vSpectra[a].m_vseqBest[0].m_vDomains.empty())	{
			dValue = m_vSpectra[a].m_hHyper.expect(m_pScore->hconvert(m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_fHyper));
			dProtein = m_vSpectra[a].m_dProteinExpect;
			if(log10(dValue) <= _d && dProtein <= dProteinMax)	{
				m_dEsum += dValue;
			}
			dValue = log10(dValue);
		}
		if(!m_vSpectra[a].m_vseqBest.empty() && !m_vSpectra[a].m_vseqBest[0].m_vDomains.empty() && dValue <= _d && dProtein <= dProteinMax)	{
			b = 0;
			while(b < m_vSpectra[a].m_vseqBest.size())	{
				itValue = m_mapSequences.find(m_vSpectra[a].m_vseqBest[b].m_tUid);
				m_vSpectra[a].m_vseqBest[b].m_strSeq = (*itValue).second;
				b++;
			}
			if(tLast > 0)	{
				if(m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_lS != m_vSpectra[tLast].m_vseqBest[0].m_vDomains[0].m_lS)	{
					if(m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_lE != m_vSpectra[tLast].m_vseqBest[0].m_vDomains[0].m_lE)	{
						m_tUnique++;
						if(m_lReversed != -1 && !m_vSpectra[tLast].m_vseqBest[0].m_bForward)	{
							m_lReversed++;
						}
					}
				}
			}
			tLast = a;
			m_tValid++;
			tActive++;
			if(bSpectra || bHistograms || bProteins)
				rptValue.group(m_vSpectra[a]);
			if(bProteins)
				rptValue.sequence(m_vSpectra[a],bSequences,m_vstrPaths,m_mapAnnotation);
			if(bHistograms)
				rptValue.histogram(m_vSpectra[a]);
			if(bSpectra)
				rptValue.spectrum(m_vSpectra[a]);
			if(bSpectra || bHistograms || bProteins)
				rptValue.endgroup();
		}
		a++;
	}
	if(m_tValid == 0)	{
		m_tUnique = 0;
	}
	strKey = "modelling, total spectra assigned";
	char *pLine = new char[256];
	sprintf(pLine,"%u",(unsigned int)m_tValid);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	strKey = "modelling, total unique assigned";
	sprintf(pLine,"%u",(unsigned int)m_tUnique);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);
	if(m_lReversed != -1)	{
		strKey = "modelling, reversed sequence false positives";
		sprintf(pLine,"%i",(int)m_lReversed);
		strValue = pLine;
		m_xmlPerformance.set(strKey,strValue);
	}
	unsigned long lE = (unsigned long)(0.5+m_dEsum);
	unsigned long lEe = (unsigned long)(0.5 + sqrt(m_dEsum));
	if(lEe == 0)	{
		lEe = 1;
	}
	strKey = "modelling, estimated false positives";
	sprintf(pLine,"%u",(unsigned int)lE);
	strValue = pLine;
	m_xmlPerformance.set(strKey,strValue);

	if(bInput)
		rptValue.info(m_xmlValues);
	if(bPerf)
		rptValue.performance(m_xmlPerformance);
	//TODO(bmaclean): Report average masses too, if using average for fragments?
	if(m_pScore->m_pSeqUtilFrag->is_modified())	{
		rptValue.masses(*(m_pScore->m_pSeqUtilFrag));
	}
	delete pLine;
	return rptValue.end();
}

bool mprocess::report_stochastic(const double _d)
{
/*
 * check the input parameters for known output information
 */
	string strKey = "output, histogram column width";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	long lHistogramColumns = 30;
	if(atoi(strValue.c_str()) > 0)
		lHistogramColumns = atoi(strValue.c_str());

	strKey = "output, spectra";
	m_xmlValues.get(strKey,strValue);
	bool bSpectra = false;
	if(strValue == "yes")
		bSpectra = true;

	strKey = "output, histograms";
	m_xmlValues.get(strKey,strValue);
	bool bHistograms = false;
	if(strValue == "yes")
		bHistograms = true;

	strKey = "output, sequences";
	m_xmlValues.get(strKey,strValue);
	bool bSequences = false;
	if(strValue == "yes")
		bSequences = true;

	strKey = "output, proteins";
	m_xmlValues.get(strKey,strValue);
	bool bProteins = false;
	if(strValue == "yes")
		bProteins = true;

	strKey = "output, parameters";
	m_xmlValues.get(strKey,strValue);
	bool bInput = false;
	if(strValue == "yes")
		bInput = true;

	strKey = "output, performance";
	m_xmlValues.get(strKey,strValue);
	bool bPerf = false;
	if(strValue == "yes")
		bPerf = true;

	strKey = "output, one sequence copy";
	m_xmlValues.get(strKey,strValue);
	bool bCompress = false;
	if(strValue == "yes")
		bCompress = true;

	mreport rptValue(*m_pScore);
	rptValue.set_compression(bCompress);
	rptValue.set_columns(lHistogramColumns);
	rptValue.start(m_xmlValues);
	m_pathName = rptValue.getPathName(); // rTANDEM
	size_t a = 0;	
	size_t b = 0;
	SEQMAP::iterator itValue;
	size_t tLength = m_vSpectra.size();
	double dValue = 0.0;
	while(a < tLength)	{
		dValue = 3.0;
		if(!m_vSpectra[a].m_vseqBest.empty() && !m_vSpectra[a].m_vseqBest[0].m_vDomains.empty())	{
			dValue = m_vSpectra[a].m_hHyper.expect(m_pScore->hconvert(m_vSpectra[a].m_vseqBest[0].m_vDomains[0].m_fHyper));
			if(log10(dValue) > _d)	{
				m_dEsum += dValue;
			}
			dValue = log10(dValue);
		}
		if(m_vSpectra[a].m_vseqBest.empty() || dValue > _d)	{
			b = 0;
			while(b < m_vSpectra[a].m_vseqBest.size())	{
				itValue = m_mapSequences.find(m_vSpectra[a].m_vseqBest[b].m_tUid);
				m_vSpectra[a].m_vseqBest[b].m_strSeq = (*itValue).second;
				b++;
			}
			if(bSpectra || bHistograms || bProteins)
				rptValue.group(m_vSpectra[a]);
			if(bProteins)
				rptValue.sequence(m_vSpectra[a],bSequences,m_vstrPaths,m_mapAnnotation);
			if(bHistograms)
				rptValue.histogram(m_vSpectra[a]);
			if(bSpectra)
				rptValue.spectrum(m_vSpectra[a]);
			if(bSpectra || bHistograms || bProteins)
				rptValue.endgroup();
		}
		m_vSpectra[a].m_vseqBest.clear();
		a++;
	}
	if(bInput)
		rptValue.info(m_xmlValues);
	if(bPerf)
		rptValue.performance(m_xmlPerformance);
	//TODO(bmaclean): Report average masses too, if using average for fragments?
	if(m_pScore->m_pSeqUtilFrag->is_modified() )	{
		rptValue.masses(*(m_pScore->m_pSeqUtilFrag));
	}
	return rptValue.end();
}
/*
 * report_expect calculates the expectation value for proteins
 */
bool mprocess::report_expect(const double _m)
{
//	cout << "\tinitial calculations ";
	Rprintf("\tinitial calculations ");
	//cout.flush();
	size_t tLength = m_vSpectra.size();
	unsigned long a = 0;
	m_tValid = 0;
	unsigned long b = 0;
	size_t tMaxUid = 0;
	SEQMAP::iterator itMap = m_mapSequences.begin();
	while(itMap != m_mapSequences.end())	{
		if(itMap->first > tMaxUid)	{
			tMaxUid = itMap->first;
		}
		itMap++;
	}
	tMaxUid += 10;
	double *pdExpect = new double[tMaxUid];
	unsigned long *plCount = new unsigned long[tMaxUid];
	const double dConstant = 2.0e32;
	while(a < tMaxUid)	{
		pdExpect[a] = dConstant;
		plCount[a] = 0;
		a++;
	}
	double dExpect = 0.0;
	double dSum = 0;
	a = 0;
	while(a < tLength)	{
		dExpect = (double)m_vSpectra[a].m_hHyper.expect_protein(m_pScore->hconvert(m_vSpectra[a].m_fHyper));
		dSum += m_vSpectra[a].m_hHyper.sum();
		if(dExpect <= _m)	{
			m_tValid++;
		}
		m_vSpectra[a].m_dExpect = dExpect;
		m_vSpectra[a].m_dProteinExpect = log10(dExpect);
		a++;
	}
	vector<mspectrum>::iterator itStart = m_vSpectra.begin();
//	cout << " ..... done.\n\tsorting ";
	Rprintf(" ..... done.\n\tsorting ");
	//cout.flush();
	string strKey = "output, results";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	if(strValue == "valid")	{
		sort(m_vSpectra.begin(),m_vSpectra.end(),lessThanSpectrum);
		while(itStart != m_vSpectra.end() && itStart->m_dExpect <= 0.95*_m)	{
			itStart++;
		}
		if(itStart != m_vSpectra.end())	{
			m_vSpectra.erase(itStart,m_vSpectra.end());
		}
	}
	strKey = "output, sort best scores by";
	m_xmlValues.get(strKey,strValue);
	bool bSortScores = (strValue != "sequence");
	unsigned long lSum = (unsigned long)(0.5+dSum/(double)tLength);
//	cout << " ..... done.\n\tfinding repeats ";
	Rprintf(" ..... done.\n\tfinding repeats ");
	//cout.flush();
	mark_repeats();
//	cout << " done.\n\tevaluating results ";
	Rprintf(" done.\n\tevaluating results ");
	//cout.flush();
	a = 0;
	tLength = m_vSpectra.size();
	size_t tUid;
	size_t tBest;
	size_t tTicLength = (size_t)((double)tLength/5.0);
	size_t tTics = 0;
	vector<mdomain>::iterator itDom;
	vector<mdomain>::iterator itDomEnd;
	bool bRound = false;
	while(a < tLength)	{
		tTics++;
		if(tTics >= tTicLength)	{
//			cout << ".";
			Rprintf(".");
			//cout.flush();
			tTics = 0;
		}
		dExpect = log10(m_vSpectra[a].m_dExpect);
		b = 0;
		tBest = m_vSpectra[a].m_vseqBest.size();
		while(b < tBest)	{
			tUid = m_vSpectra[a].m_vseqBest[b].m_tUid;
			bRound = false;
			if(m_vSpectra[a].m_vseqBest[b].m_iRound < 10 || m_lReversed != -1)	{
				bRound = true;
				m_setRound.insert(tUid);
			}
			if(!m_vSpectra[a].m_bRepeat && dExpect <= -1.0)	{
				if(pdExpect[tUid] != dConstant)	{
					if(bRound)	{
						plCount[tUid]++;
					}
					pdExpect[tUid] += dExpect;
				}
				else if(bRound)	{
					plCount[tUid] = 1;
					if(dExpect <= -1.0)	{
						pdExpect[tUid] = dExpect;
					}
					else	{
						pdExpect[tUid] = 0.0;
					}
				}
			}
			b++;
		}
		a++;
	}
	a = 0;
	vector<msequence>::iterator itSeq;
	double dBias = (double)m_tPeptideCount/(double)m_tProteinCount;
	dBias /= (double)m_tTotalResidues/(double)m_tProteinCount;
	SEQMAP::iterator itValue;
//	cout << " done.\n\tcalculating expectations ";
	Rprintf(" done.\n\tcalculating expectations ");
	//cout.flush();
	tTics = 0;
	map<size_t,double> mapExpect;
	pair<size_t,double> pairExpect;
	set<size_t>::iterator itRound = m_setRound.end();
	while(a < tLength)	{
		b = 0;
		tTics++;
		if(tTics >= tTicLength)	{
//			cout << ".";
			Rprintf(".");
			//cout.flush();
			tTics = 0;
		}
		itSeq = m_vSpectra[a].m_vseqBest.begin();
		while(itSeq != m_vSpectra[a].m_vseqBest.end())	{
			tUid = itSeq->m_tUid;
			if(m_setRound.find(tUid) != itRound)	{
				if(pdExpect[tUid] != dConstant)	{
					if(mapExpect.find(tUid) == mapExpect.end())	{
						itSeq->m_dExpect = expect_protein(plCount[tUid],
							(unsigned long)tLength,lSum,
							pdExpect[tUid]);
						pairExpect.first = tUid;
						pairExpect.second = itSeq->m_dExpect;
						mapExpect.insert(pairExpect);
					}
					else	{
						itSeq->m_dExpect = mapExpect.find(tUid)->second;
					}
				}
				else	{
					itSeq->m_dExpect = 0.0;
				}
				if(plCount[tUid] > 1)	{
					itValue= m_mapSequences.find(tUid);
					if((*itValue).second.size()*dBias < 1.0)	{
						itSeq->m_dExpect += plCount[tUid]*log10((*itValue).second.size()*dBias);
					}
				}
				itSeq++;
			}
			else	{
				itSeq = m_vSpectra[a].m_vseqBest.erase(itSeq);
			}
		}
		if(!m_vSpectra[a].m_vseqBest.empty())	{
			if(bSortScores)	{
				sort(m_vSpectra[a].m_vseqBest.begin(),m_vSpectra[a].m_vseqBest.end(),lessThanSequence);
			}
			vector<msequence>::iterator itA = m_vSpectra[a].m_vseqBest.begin();
			vector<msequence>::iterator itB = m_vSpectra[a].m_vseqBest.begin();
			while(itA != m_vSpectra[a].m_vseqBest.end())	{
				itB++;
				if(itB == m_vSpectra[a].m_vseqBest.end() || itB->m_dExpect != itA->m_dExpect)	{
					sort(itA,itB,lessThanSequenceUid);
					itA = itB;
				}
			}
			m_vSpectra[a].m_dProteinExpect = m_vSpectra[a].m_vseqBest[0].m_dExpect;
		}
		a++;
	}
	a = 0;
	b = 0;
	while(a < tMaxUid)	{
		pdExpect[a] = dConstant;
		a++;
	}
	a = 0;
	while(a < tLength)	{
		b = 0;
		tBest = m_vSpectra[a].m_vseqBest.size();
		while(b < tBest)	{
			tUid = m_vSpectra[a].m_vseqBest[b].m_tUid;
			if(pdExpect[tUid] == dConstant)	{
				pdExpect[tUid] = m_vSpectra[a].m_vdStats[0];
			}
			else	{
				pdExpect[tUid] += m_vSpectra[a].m_vdStats[0];
			}
			b++;
		}
		a++;
	}
	a = 0;
	b = 0;
	while(a < tLength)	{
		b = 0;
		tBest = m_vSpectra[a].m_vseqBest.size();
		while(b < tBest)	{
			m_vSpectra[a].m_vseqBest[b].m_fIntensity = (float)(pdExpect[m_vSpectra[a].m_vseqBest[b].m_tUid]);
			b++;
		}
		a++;
	}
	a = 0;
//	cout << " done.\n";
	Rprintf(" done.\n");
	//cout.flush();
	delete plCount;
	delete pdExpect;
	return true;
}
/*
 * report_sort contains the logic for sorting the output results, if sorting
 * is desired.
 */
bool mprocess::report_sort()
{
	string strKey = "output, sort results by";
	string strValue;
	m_xmlValues.get(strKey,strValue);
	vector<msequence>::iterator itSeq;
	vector<msequence>::iterator itSort;
	if(strValue == "protein")	{
		sort(m_vSpectra.begin(),m_vSpectra.end(),lessThanSpectrum);
		vector<mspectrum>::iterator itStart = m_vSpectra.begin();
		vector<mspectrum>::iterator itEnd = itStart;
		while(itStart != m_vSpectra.end() && itEnd != m_vSpectra.end())	{
			while(itStart != m_vSpectra.end() && itStart->m_vseqBest.empty())
				itStart++;
			if(itStart == m_vSpectra.end())
				break;
			itEnd = itStart + 1;
			while(itEnd != m_vSpectra.end())	{
				if(!itEnd->m_vseqBest.empty())	{
					if(itStart->m_vseqBest[0].m_tUid == itEnd->m_vseqBest[0].m_tUid)	{
						itEnd++;
					}
					else	{
						break;
					}
				}
				else	{
					break;
				}
			}
			if(itEnd != itStart + 1)	{
				sort(itStart,itEnd,lessThanOrder);
			}
			itStart = itEnd;
		}
	}
	return true;
}

/*
 * residues is used to estimate the minimum number of residues that are required for a peptide
 * to match the lowest mass spectrum in the m_Spectra vector. this estimate is made to
 * improve performance, by simply rejecting very small sequences from consideration by the
 * calculationally intensive scoring routines.
 */
bool mprocess::residues()
{
/*
 * make a very conservative estimate of the minimum number of residues
 * As of version 2004.03.01, the more complex version of this function has been
 * altered to return simply a constant value. The older method of estimating
 * the minimum number of residues was found to fail badly when there were
 * large moeities added as modifications to the peptide sequences.
 */

	m_tMinResidues = 5; 
	return true;
}

/*
* rollback is a method that is used to reverse changes that have been made to the
* m_vSpectra vector during the refinement process. If a new result, discovered during
* the refinement process, is not sufficiently significant (as determined by _f), then
* the m_vSpectra entry is "rolled back" to the value it had after the initial 
* survey round.
*/
bool mprocess::rollback(vector<mspectrum>& _v,const double _m,const double _f)
{
	if(_v.empty())	{
		return false;
	}
	size_t a = 0;
	const size_t tSize = m_vSpectra.size();
	double dExpect = 1.0;
	double dExpectLast = 1.0;
	bool bDone = false;
	while(a < tSize)	{
		if(!m_vSpectra[a].m_vseqBest.empty() && !_v[a].m_vseqBest.empty())	{
			bDone = false;
			m_vSpectra[a].m_hHyper.model();
			m_vSpectra[a].m_hHyper.set_protein_factor(1.0);
			dExpect = (double)m_vSpectra[a].m_hHyper.expect_protein(m_pScore->hconvert(m_vSpectra[a].m_fHyper));
			dExpectLast = (double)m_vSpectra[a].m_hHyper.expect_protein(m_pScore->hconvert(_v[a].m_fHyper));
			if(dExpect > _m)	{
				m_vSpectra[a] *= _v[a];
				bDone = true;
			}
			else if(dExpect <= _m  && dExpect/dExpectLast > _f)	{
				m_vSpectra[a] *= _v[a];
				bDone = true;
			}
			else if(!bDone && m_vSpectra[a].m_fHyper == _v[a].m_fHyper)	{
				m_vSpectra[a] *= _v[a];
			}
		}
		a++;
	}
	_v.clear();
	return true;
}

/*
 * score takes an msequence object and sequentially creates all possible cleavage peptides
 * from that msequence. each peptide is then tested and scored.
 */
bool mprocess::score(const msequence &_s)
{
	size_t m = 0;
	string strValue;
	if(!m_vstrModifications.empty())	{
		strValue = m_vstrModifications[m];
		m_pScore->m_seqUtil.modify_all(strValue);
		m_pScore->m_seqUtilAvg.modify_all(strValue);
	}
	bool bReturn = score_single(_s);
	m++;
	while(m < m_vstrModifications.size())	{
		strValue = m_vstrModifications[m];
		m_pScore->m_seqUtil.modify_all(strValue);
		m_pScore->m_seqUtilAvg.modify_all(strValue);
		bReturn = score_single(_s);
		m++;
	}
	return bReturn;
}

bool mprocess::score_terminus(const string &_s)
{
	size_t m = 0;
	string strValue;
	if(!m_vstrModifications.empty())	{
		strValue = m_vstrModifications[m];
		m_pScore->m_seqUtil.modify_all(strValue);
		m_pScore->m_seqUtilAvg.modify_all(strValue);
	}
	bool bReturn = score_terminus_single(_s);
	m++;
	while(m < m_vstrModifications.size())	{
		strValue = m_vstrModifications[m];
		m_pScore->m_seqUtil.modify_all(strValue);
		m_pScore->m_seqUtilAvg.modify_all(strValue);
		bReturn = score_terminus_single(_s);
		m++;
	}
	return bReturn;
}

bool mprocess::score_single(const msequence &_s)
{
	m_pScore->m_seqUtil.motif_set(_s);
	if((m_bAnnotation || m_bMinimalAnnotation) && !m_mapAnnotation.empty())	{
		map <string,string>::iterator itMod;
		if(_s.m_strDes.find("IPI") != _s.m_strDes.npos)	{
			size_t tFind = _s.m_strDes.find("IPI");
			size_t tEnd = _s.m_strDes.find(".",tFind);
			itMod = m_mapAnnotation.find(_s.m_strDes.substr(tFind,tEnd-tFind));
		}
		else	{
			itMod = m_mapAnnotation.find(_s.m_strDes);
		}
		string strMods;
		if(itMod != m_mapAnnotation.end())	{
			strMods = itMod->second;
		}
		else	{
			strMods.clear();
		}
		if(m_strLastMods != strMods)	{
			if(!m_bAnnotation && m_bMinimalAnnotation)	{
				if(strMods.find("-1@B") != strMods.npos)	{
					m_pScore->m_seqUtil.modify_annotation(strMods);
				}
				else	{
					strMods.clear();
					m_pScore->m_seqUtil.modify_annotation(strMods);
				}
			}
			else	{
				m_pScore->m_seqUtil.modify_annotation(strMods);
			}
		}
		m_strLastMods = strMods;
	}
	string strDesc = _s.m_strDes;
	m_pScore->set_saps(m_bSaps,strDesc);
	long lLength = (long)_s.m_strSeq.size();
	if(m_tSeqSize < (size_t)(lLength+1))	{
		delete m_pSeq;
		m_tSeqSize = 4096*(size_t)(ceil((double)lLength/4096.0)+1);
		m_pSeq = new char[m_tSeqSize];
	}
	strcpy(m_pSeq,_s.m_strSeq.c_str());
	char cValue;
	m_tTotalResidues += lLength;
	long lStart = 0;
	if(m_bRefineCterm)	{
		lStart = lLength - m_lCStartMax;
		if(lStart < 0)	{
			lStart = 0;
		}
	}
	long lEnd = 0;
	bool bIsFirst = true;
	bool bNg = false;
	long lMissedCleaves = 0;
	long lNextStart = 0;
/*
 * set up variables and make some const local copies of frequently used member variables
 */
	const long lMissedMax = (long)m_tMissedCleaves;
	const long lMinAa = (long)m_tMinResidues;
	long lLastCleave = -1;
	float fMinMax = 0;
	m_semiState.limit(lMinAa);
	if(m_Cleave.m_lType & 0x01)	{
		m_semiState.activate(false);
	}
	const char cAster = '*';
	m_dNt = m_pScore->m_seqUtil.m_pdAaMod['['];
	m_dNtAve = m_pScore->m_seqUtilAvg.m_pdAaMod['['];
	m_dNg = m_pScore->m_seqUtil.m_pdAaMod['n'];
	m_dNgAve = m_pScore->m_seqUtilAvg.m_pdAaMod['n'];
	bool bModsUtil = m_pScore->m_seqUtil.m_bPotential;
	bool bModsUtilAvg = m_pScore->m_seqUtilAvg.m_bPotential;
	const double dAcetyl = 42.010565;
	const double dAcetylAve = 42.0367;
	const double dDeamidate = 0.984016;
	const double dDeamidateAve = 0.9848;
	bool bFirstAcetyl = false;
	unsigned long lRagged = 1;
	bool bMiss = true;
/*
 * continue operations until the start cursor (the beginning of the new peptide) reaches
 * the end of the sequence
 */
	while(lStart < lLength && lStart < m_lStartMax)	{
		while(m_pSeq[lStart] == cAster && lStart < m_lStartMax)	{
			lStart++;
		}
		lEnd = lStart;
		if(lEnd >= lLength)
			lEnd = lLength - 1;
		bIsFirst = true;
		lMissedCleaves = 0;
		lLastCleave = -1;
/*
 * with a given lStart, vary the end cursor (lEnd), using the cleavage rules, until
 * none of the spectra have a parent ion as large as the peptide (tEligible == 0),
 * or the number of missed cleavages is larger than the maximum allowed
 */
		while(lEnd < lLength)	{
/*
 * get the next allowed lEnd
 */
			bMiss = true;
			while(lEnd < lLength)	{
				if(m_pSeq[lEnd] == cAster)	{
					if(lMissedCleaves == 0)
						lNextStart = lEnd+1;
					lMissedCleaves = lMissedMax;
					lEnd--;
					break;
				}
				if(m_pSeq[lEnd] == 'D' && m_pSeq[lEnd+1] == 'P')	{
					bMiss = false;
					if(lMissedCleaves == 0)
						lNextStart = lEnd+1;
					break;
				}
				else if(m_Cleave.m_lType & 0x02)	{
					if(m_pSeq[lEnd+1] != 'P')	{
						if(m_pSeq[lEnd] == 'K' || m_pSeq[lEnd] == 'R')	{
							if(lMissedCleaves == 0)
								lNextStart = lEnd+1;
							break;
						}
					}
				}
				else if(m_Cleave.m_lType & 0x01)	{
					if(lMissedCleaves == 0)
						lNextStart = lEnd+1;
					break;
				}
				else if(m_Cleave.test(m_pSeq[lEnd],m_pSeq[lEnd+1]))	{
					if(lMissedCleaves == 0)
						lNextStart = lEnd+1;
					break;
				}
				lEnd++;
			}
			if(lEnd == lLength && lMissedCleaves == 0)	{
					lNextStart = lEnd+1;
			}
			if(bMiss)	{
				lMissedCleaves++;
			}
/*
 * update the peptide sequence in m_pScore
 */
			if(lEnd >= lLength)
				lEnd = lLength - 1;
			if(lEnd - lStart >= lMinAa && lEnd < lLength)	{
				m_semiState.reset(lStart,lEnd,lLastCleave);
				fMinMax = 0.0;
				if(m_bQuickAcetyl && lStart == 0 && m_pSeq[0] == 'M' && (m_dNt == 0.0 || (int)m_dNt == 42) && m_iCurrentRound < 4)	{
					m_pScore->m_seqUtil.m_pdAaMod['['] = dAcetyl;
					m_pScore->m_seqUtilAvg.m_pdAaMod['['] = dAcetylAve;
					bFirstAcetyl = true;
				}
				do	{
					cValue = m_pSeq[lEnd+1];
					m_pSeq[lEnd+1] = '\0';
					if(m_bQuickAcetyl && lStart <= 2 && m_pSeq[0] == 'M' && (m_dNt == 0.0 || (int)m_dNt == 42) && m_iCurrentRound < 4)	{
						m_pScore->m_seqUtil.m_pdAaMod['['] = dAcetyl;
						m_pScore->m_seqUtilAvg.m_pdAaMod['['] = dAcetylAve;
						m_pScore->m_seqUtil.m_bPotential = true;
						m_pScore->m_seqUtilAvg.m_bPotential = true;
						bFirstAcetyl = true;
					}
					else if(m_bQuickAcetyl && bFirstAcetyl)	{
						m_pScore->m_seqUtil.m_pdAaMod['['] = m_dNt;
						m_pScore->m_seqUtilAvg.m_pdAaMod['['] = m_dNtAve;
						m_pScore->m_seqUtil.m_bPotential = bModsUtil;
						m_pScore->m_seqUtilAvg.m_bPotential = bModsUtilAvg;
					}
					if(m_bCheckNg && m_dNg == 0.0 && strstr(m_pSeq+lStart,"NG"))	{
						m_pScore->m_seqUtil.m_pdAaMod['n'] = dDeamidate;
						m_pScore->m_seqUtilAvg.m_pdAaMod['n'] = dDeamidateAve;
						bNg = true;
					}
					pyro_check(*(m_pSeq+lStart)); // must be done before call to set_seq !!!!
					if(bIsFirst || m_semiState.m_bActive)	{
						m_pScore->set_seq(m_pSeq+lStart,lStart == 0,lEnd == lLength-1,lEnd-lStart+1,lStart);
					}
					else	{
						m_pScore->add_seq(m_pSeq+lStart,lStart == 0,lEnd == lLength-1,lEnd-lStart+1,lStart);
					}
					m_pSeq[lEnd+1] = cValue;
					bIsFirst = false;
	/*
	* use the m_pScore state machine to obtain allowed modified peptides that have
	* the same sequence. then use create_score to score relavent spectra
	*/
					while(m_pScore->load_next())	{
						m_tPeptideCount += m_pScore->m_State.m_lEqualsS;
						m_bPermute = false;
						m_bPermuteHigh = false;
						create_score(_s,lStart,lEnd,lMissedCleaves - 1,true);
/*
 * Perform the permutation if the 
 * 1. permutation must be done and the check is fine.
 *  
 * or 
 *
 * 2. we are in a high permutation zone.
 */
						if((m_bPermute && m_bCrcCheck) || m_bPermuteHigh)	{
							m_pScore->reset_permute();
							while(m_pScore->permute())	{
								create_score(_s,lStart,lEnd,lMissedCleaves - 1,false);
							}
						}
					}
					if(m_pyroState.m_bPyro)	{
						pyro_reset();
					}
					if(bNg)	{
						m_pScore->m_seqUtil.m_pdAaMod['n'] = m_dNg;
						m_pScore->m_seqUtilAvg.m_pdAaMod['n'] = m_dNgAve;
						m_pScore->m_seqUtil.m_bPotential = bModsUtil;
						m_pScore->m_seqUtilAvg.m_bPotential = bModsUtilAvg;
						bNg = false;
					}
	/*
	* as of version 2004.03.01, the test for the number of missed cleavage sites (see the next note)
	* was moved outside of this test.
	*/
				} while(m_semiState.next(lStart,lEnd));
				if(bFirstAcetyl)	{
					m_pScore->m_seqUtil.m_pdAaMod['['] = m_dNt;
					m_pScore->m_seqUtilAvg.m_pdAaMod['['] = m_dNtAve;
						m_pScore->m_seqUtil.m_bPotential = bModsUtil;
						m_pScore->m_seqUtilAvg.m_bPotential = bModsUtilAvg;
					bFirstAcetyl = false;
				}
				if(!m_semiState.m_bActive && m_bQuickAcetyl && lStart == 0 && m_pSeq[0] == 'M' && (m_dNt == 0.0 || (int)m_dNt == 42) && m_iCurrentRound < 4)	{
					lRagged = 1;
					while(lRagged <= 2)	{
						lStart = lRagged;
						if(*(m_pSeq+lStart) != 'Q')	{
							m_pScore->m_seqUtil.m_pdAaMod['['] = dAcetyl;
							m_pScore->m_seqUtilAvg.m_pdAaMod['['] = dAcetylAve;
							m_pScore->m_seqUtil.m_bPotential = true;
							m_pScore->m_seqUtilAvg.m_bPotential = true;
						}
						cValue = m_pSeq[lEnd+1];
						m_pSeq[lEnd+1] = '\0';
						if(m_bCheckNg && m_dNg == 0.0 && strstr(m_pSeq+lStart,"NG"))	{
							m_pScore->m_seqUtil.m_pdAaMod['n'] = dDeamidate;
							m_pScore->m_seqUtilAvg.m_pdAaMod['n'] = dDeamidateAve;
							m_pScore->m_seqUtil.m_bPotential = true;
							m_pScore->m_seqUtilAvg.m_bPotential = true;
						}
						pyro_check(*(m_pSeq+lStart)); // must be done prior to call to set_seq !!!!!
						m_pScore->set_seq(m_pSeq+lStart,lStart == 0,lEnd == lLength-1,lEnd-lStart+1,lStart);
						m_pSeq[lEnd+1] = cValue;
						bIsFirst = true;
		/*
		* use the m_pScore state machine to obtain allowed modified peptides that have
		* the same sequence. then use create_score to score relavent spectra
		*/
						while(m_pScore->load_next())	{
							m_tPeptideCount += m_pScore->m_State.m_lEqualsS;
							m_bPermute = false;
							m_bPermuteHigh = false;
							create_score(_s,lStart,lEnd,lMissedCleaves - 1,true);
						}
						if(m_pyroState.m_bPyro)	{
							pyro_reset();
						}
						lStart = 0;
						m_pScore->m_seqUtil.m_pdAaMod['['] = m_dNt;
						m_pScore->m_seqUtilAvg.m_pdAaMod['['] = m_dNtAve;
						if(bNg)	{
							m_pScore->m_seqUtil.m_pdAaMod['n'] = m_dNg;
							m_pScore->m_seqUtilAvg.m_pdAaMod['n'] = m_dNgAve;
							bNg = false;
						}
						m_pScore->m_seqUtil.m_bPotential = bModsUtil;
						m_pScore->m_seqUtilAvg.m_bPotential = bModsUtilAvg;
						lRagged++;
					}
				}
				if(m_pScore->m_fMinMass > fMinMax)	{
					fMinMax = m_pScore->m_fMinMass;
				}
				if(fMinMax - m_pScore->m_fMaxMass > 100.0)	{
					break;
				}
			}
			lLastCleave = lEnd;
/*
 * as of version 2004.03.01, the test for the number of missed cleavage sites was moved to
 * out of the previous if structure, because it could occasionally result in peptides being
 * considered that contained too many missed cleavage sites.
 */
			if(lMissedCleaves > lMissedMax)	{
				break;
			}
			lEnd++;
		} 
/*
 * get the next allowed value for the start peptide cursor
 */
		if(m_pSeq[lStart] == 'D' && m_pSeq[lStart+1] == 'P')	{
				lStart++;
		}
		else if(m_Cleave.m_lType & 0x02)	{
			if(m_pSeq[lStart+1] != 'P')	{
				if(m_pSeq[lStart] == 'K' || m_pSeq[lStart] == 'R')	{
					lStart++;
				}
				else	{
					lStart = lNextStart;
				}
			}
			else	{
				lStart = lNextStart;
			}
		}
		else if(m_Cleave.test(m_pSeq[lStart],m_pSeq[lStart+1]))	{
			lStart++;
		}
		else	{
			lStart = lNextStart;
		}
	}
	return true;
}

__inline__ bool mprocess::pyro_check(const char _c)
{
	if(m_pScore->m_seqUtil.m_pdAaMod['['] != 0.0 || !m_bQuickPyro)	{
		m_pyroState.m_bPyro = false;
		return false;
	}
	if(_c == 'Q')	{
		m_pyroState.m_dModMass = -1.0*m_pScore->m_seqUtil.m_dAmmonia;
		m_pScore->m_seqUtil.m_pdAaMod['['] = m_pyroState.m_dModMass;
		m_pScore->m_seqUtilAvg.m_pdAaMod['['] = -1.0*m_pScore->m_seqUtilAvg.m_dAmmonia;
		m_pyroState.m_bPotential = m_pScore->m_seqUtil.m_bPotential;
		m_pScore->m_seqUtil.m_bPotential = true;
		m_pScore->m_seqUtilAvg.m_bPotential = true;
		m_pyroState.m_bPyro = true;
		m_pyroState.m_cRes = 'Q';
		return true;
	}
	else if(_c == 'E')	{
		m_pyroState.m_dModMass = -1.0*m_pScore->m_seqUtil.m_dWater;
		m_pScore->m_seqUtil.m_pdAaMod['['] = m_pyroState.m_dModMass;
		m_pScore->m_seqUtilAvg.m_pdAaMod['['] = -1.0*m_pScore->m_seqUtilAvg.m_dWater;
		m_pyroState.m_bPotential = m_pScore->m_seqUtil.m_bPotential;
		m_pScore->m_seqUtil.m_bPotential = true;
		m_pScore->m_seqUtilAvg.m_bPotential = true;
		m_pyroState.m_bPyro = true;
		m_pyroState.m_cRes = 'E';
		return true;
	}
	else if(_c == 'C' && (long)(m_pScore->m_seqUtil.m_pdAaFullMod['C']) == 57)	{
		m_pyroState.m_dModMass = -1.0*m_pScore->m_seqUtil.m_dAmmonia;
		m_pScore->m_seqUtil.m_pdAaMod['['] = m_pyroState.m_dModMass;
		m_pScore->m_seqUtilAvg.m_pdAaMod['['] = -1.0*m_pScore->m_seqUtilAvg.m_dAmmonia;
		m_pyroState.m_bPotential = m_pScore->m_seqUtil.m_bPotential;
		m_pScore->m_seqUtil.m_bPotential = true;
		m_pScore->m_seqUtilAvg.m_bPotential = true;
		m_pyroState.m_bPyro = true;
		m_pyroState.m_cRes = 'C';
		return true;
	}
	m_pyroState.m_bPyro = false;
	return false;
}

__inline__ bool mprocess::pyro_reset()
{
	m_pScore->m_seqUtil.m_bPotential = m_pyroState.m_bPotential;
	m_pScore->m_seqUtilAvg.m_bPotential = m_pyroState.m_bPotential;
	m_pyroState.m_bPyro = false;
	m_pyroState.m_cRes = '\0';
	m_pScore->m_seqUtil.m_pdAaMod['['] = m_dNt;
	m_pScore->m_seqUtilAvg.m_pdAaMod['['] = m_dNtAve;
	m_pyroState.m_dModMass = 0.0;
	return true;
}

/*
 * score_each_sequence runs through the list of sequences in m_svrSequences and
 * generates scores for the peptides in those sequences.
 */
bool mprocess::score_each_sequence()
{
	size_t tLength = m_svrSequences.m_pCol->size();
	size_t a = 0;
//	map <string,string>::iterator itMod;
//	string strMods;
	string strValue;
/*
 * go through the msequenceCollection object and score each msequence
 */
	while(a < tLength)	{
		if(!m_bReversedOnly)	{
			m_svrSequences.m_pCol->m_vASequences[a].m_tUid = m_tProteinCount+1;
			m_svrSequences.m_pCol->m_vASequences[a].m_bForward = true;
			score(m_svrSequences.m_pCol->m_vASequences[a]);
			m_tProteinCount++;
		}
		if(m_lReversed != -1)	{
			m_svrSequences.m_pCol->m_vASequences[a].m_tUid = m_tProteinCount+1;
			m_svrSequences.m_pCol->m_vASequences[a].m_bForward = false;
			string strTemp;
			string::reverse_iterator itS = m_svrSequences.m_pCol->m_vASequences[a].m_strSeq.rbegin();
			string::reverse_iterator itE = m_svrSequences.m_pCol->m_vASequences[a].m_strSeq.rend();
			while(itS != itE)	{
				strTemp += *itS;
				itS++;
			}
			m_svrSequences.m_pCol->m_vASequences[a].m_strSeq = strTemp;
			m_svrSequences.m_pCol->m_vASequences[a].m_strDes += ":reversed";
			score(m_svrSequences.m_pCol->m_vASequences[a]);
			m_tProteinCount++;
		}
		a++;
	}
	return true;
}

/*
 *  score_terminus is used to attempt to find modified terminii. the modifications
 * are entered as a string, e.g.
 *		_s = "42@[,27@["
 * means that the N-terminus may be modified by 42 or 27 daltons. Similarly:
 *		_s = -1@]"
 * means that the C-terminus may be modified by -1 dalton. each potential modification
 * is used separately. the protein sequence is cleaved using [X]|[X].
 */
bool mprocess::score_terminus_single(const string &_s)
{
	if(_s.size() == 0)	{
		return false;
	}
	size_t tStart = 0;
	size_t tAt = 0;
	string strValue = _s.substr(tStart,_s.size()-tStart);
	double dValue = atof(strValue.c_str());
	string strKey = "refine, tic percent";
	m_xmlValues.get(strKey,strValue);
	double dTicPercent = atof(strValue.c_str());
	if(dTicPercent == 0)	{
		dTicPercent = 20.0;
	}
	size_t tTicMax = (size_t)((double)m_vseqBest.size()*dTicPercent/100.0);
	if(tTicMax < 1)	{
		tTicMax = 1;
	}
	size_t a = 0;
	size_t tPips = 0;
	bool bPotential = m_pScore->m_seqUtil.m_bPotential;
	while(fabs(dValue) > 0.001)	{
		tAt = _s.find('@',tStart);
		if(tAt == _s.npos)
			break;
		tAt++;
		m_pScore->m_seqUtil.m_bPotential = true;
		m_pScore->m_seqUtilAvg.m_bPotential = true;
		int index = _s[tAt];
		m_pScore->m_seqUtil.m_pdAaMod[index] = dValue;
		m_pScore->m_seqUtilAvg.m_pdAaMod[index] = dValue;
		a = 0;
		tPips = 0;
		while(a < m_vseqBest.size())	{
			score(m_vseqBest[a]);
			tPips++;
			if(tPips == tTicMax)	{
				if(m_lThread == 0 || m_lThread == 0xFFFFFFFF)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					m_prcLog.log(".");
				}
				tPips = 0;
			}
			a++;
		}
		tStart = _s.find(',',tAt);
		if(tStart == _s.npos)
			break;
//		cout << ". ";
		Rprintf(". ");
		//cout.flush();
		tStart++;
		strValue = _s.substr(tStart,_s.size()-tStart);
		dValue = atof(strValue.c_str());
	}
	m_pScore->m_seqUtil.m_bPotential = bPotential;
	m_pScore->m_seqUtilAvg.m_bPotential = bPotential;
	return true;
}
/*
 * set_thread tells the mprocess it's assigned thread number. this number
 * is used to determine which parts of the sequence list are to be
 * processed by this mprocess object
 */
bool mprocess::set_thread(const unsigned long _t)
{
	m_lThread = _t;
	return true;
}

bool mprocess::set_threads(const unsigned long _t)
{
	m_lThreads = _t;
	if(m_lThreads == 1)	{
		m_lThread = 0xFFFFFFFF;
	}
	return true;
}
/*
 * using the input parameters from the input XML file, retrieve the tandem MS spectra
 */
bool mprocess::spectra()
{
	string strValue;
	string strKey;
	strKey = "spectrum, threads";
	m_xmlValues.get(strKey,strValue);
	unsigned long lThreads = atoi(strValue.c_str());
	if(lThreads > 1)	{
		m_lThreads = lThreads;
	}
	else	{
		m_lThread = 0xFFFFFFFF;
	}
	strKey = "output, log path";
	strValue = "no";
	m_xmlValues.get(strKey,strValue);
	if(m_lThread == 0 || m_lThread == 0xFFFFFFFF)	{
		if(strValue.size() > 0)	{
			m_prcLog.open(strValue);
			strKey = "output, path";
			m_xmlValues.get(strKey,strValue);
			m_prcLog.log("X! Tandem starting");
		}
	}
	if(m_vSpectra.size() > 0)	{
		m_tSpectraTotal = m_vSpectra.size();
		return true;
	}
	m_vSpectra.clear();
	mspectrum spCurrent;
	bool bContinue = true;
	m_tSpectraTotal = 0;
	strKey = "spectrum, path";
	m_xmlValues.get(strKey,strValue);
	FILE *pStream;
	bool bCommon = false;
	pStream = fopen(strValue.c_str(),"r");
	char *pValue = new char[1028];
	memset(pValue,0,1028);
	size_t tRead = 256;
	if(pStream)	{
		tRead = fread(pValue,1,256,pStream);
	}
	if(pStream)	{
		fclose(pStream);

		if((pValue[0] == 1 && pValue[1] == -95)  // TODO: not sure about these parentheses.
			|| (pValue[3] == 'F' && pValue[5] == 'i' && pValue[7] == 'n') )	{

			Rprintf("\nFailed to read spectrum file: %s\n", strValue.c_str());
			Rprintf("Most likely cause: using a Finnigan raw spectrum.\nUse dta, pkl, mgf, mzdata (v.1.05) or mzxml (v.2.0) files ONLY! (1)\n\n");
			m_prcLog.log("error reading spectrum file 1");
			delete pValue;
			return false;
		}
		else if(strstr(pValue,"CMN ") == pValue)	{
			bCommon = true;
		}
		else 	{
			size_t d = 0;
			while(d < tRead)	{
				if(pValue[d] == '\0')	{
//					cout << "\nFailed to read spectrum file: " << strValue.c_str() << "\n";
					Rprintf("\nFailed to read spectrum file: %s\n", strValue.c_str());
//					cout << "Most likely cause: using a binary spectrum file.\nUse dta, pkl, mgf, mzdata (v.1.05) or mzxml (v.2.0) files ONLY! (2)\n\n";
					Rprintf("Most likely cause: using a binary spectrum file.\nUse dta, pkl, mgf, mzdata (v.1.05) or mzxml (v.2.0) files ONLY! (2)\n\n");
					//cout.flush();
					m_prcLog.log("error reading spectrum file 2");
					delete pValue;
					return false;
				}
				d++;
			}
		}
		if(strstr(pValue,"<HTML") != NULL || strstr(pValue,"<!DOCTYPE HTML") != NULL || strstr(pValue,"<html") != NULL)	{
//				cout << "\nFailed to read spectrum file: " << strValue.c_str() << "\n";
				Rprintf("\nFailed to read spectrum file: %s\n", strValue.c_str());
//				cout << "Most likely cause: using an HTML file.\nUse dta, pkl, mgf, mzdata (v.1.05) or mzxml (v.2.0) files ONLY! (2)\n\n";
				Rprintf("Most likely cause: using an HTML file.\nUse dta, pkl, mgf, mzdata (v.1.05) or mzxml (v.2.0) files ONLY! (2)\n\n");
				//cout.flush();
				m_prcLog.log("error reading spectrum file 3");
				delete pValue;
				return false;
		}
	}

	ifstream ifTest;
	ifTest.open(strValue.c_str());
	ifTest.getline(pValue,1024);
	if(strlen(pValue) == 1023)	{
		ifTest.close();
		ifTest.clear();
		ifTest.open(strValue.c_str());
		pValue[0] = '\0';
		ifTest.getline(pValue,1024,'\r');
		ifTest.close();
		if(strlen(pValue) == 1023 && !strchr(pValue,'<'))	{
//			cout << "\nFailed to read spectrum file: " << strValue.c_str() << "\n";
			Rprintf("\nFailed to read spectrum file: %s\n", strValue.c_str());
//			cout << "Most likely: an unsupported data file type:\nUse dta, pkl, mgf, mzdata (v.1.05) or mzxml (v.2.0) files ONLY! (3)\n\n";
			Rprintf("Most likely: an unsupported data file type:\nUse dta, pkl, mgf, mzdata (v.1.05) or mzxml (v.2.0) files ONLY! (3)\n\n");
			//cout.flush();
			m_prcLog.log("error reading spectrum file 3");
			delete pValue;
			return false;
		}
	}
	ifTest.close();
	delete pValue;
	long lLoaded = 0;
	long lLimit = 2000;
/*
 * check the input file for GAML encoded spectra
 */
	m_prcLog.log("loading spectra");
	if(bCommon)	{
		loadcmn ldCmn;
		if(ldCmn.open(strValue))	{
//			cout << " (cmn).";
			Rprintf(" (cmn).");
			while(ldCmn.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					lLoaded = 0;
					m_prcLog.log(".\n");
				}
				if(m_specCondition.condition(spCurrent, *m_pScore))	{
					m_vSpectra.push_back(spCurrent);
				}
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				m_vSpectra.push_back(spCurrent);
			}
			bContinue = false;
		}
	}
	if(bContinue){
		bool bState = m_specCondition.m_bCondition;
//		m_specCondition.use_condition(false); // rTANDEM
//		loadgaml ldGaml(m_vSpectra, m_specCondition, *m_pScore); // rTANDEM
//		string strV; // rTANDEM
//		m_xmlValues.getpath(strV); // rTANDEM
//		if(ldGaml.open(strV))	{ // rTANDEM
//			cout << " (gaml)."; // rTANDEM
//			ldGaml.get(); // rTANDEM
//			m_tSpectraTotal = m_vSpectra.size(); // rTANDEM
//			bContinue = false; // rTANDEM
//		} // rTANDEM
		m_specCondition.use_condition(bState);
	}
/*
 * if the input file did not contain GAML spectra, test for the existence of a
 * "spectrum, path type" parameter. If it exists, use that file format type
 * to force the use of that format.
 * Added in the 2005.08.15 release
 */
	if(bContinue)	{
		strKey = "spectrum, path type";
		string strType;
		m_xmlValues.get(strKey,strType);
		if(strType.size() > 0)	{
			return spectra_force(strType,strValue);
		}
	}
/*
 * if the input file did not contain GAML spectra, test the spectrum path parameter
 * for aGAML spectrum information
 * the ability to use a GAML formated spectrum file was added in v 2004-04-01
 */
	if(bContinue)	{
		bool bState = m_specCondition.m_bCondition;
		m_specCondition.use_condition(false);
		loadgaml ldGaml(m_vSpectra, m_specCondition, *m_pScore);
		if(ldGaml.open(strValue))	{
//			cout << " (gaml).";
			Rprintf(" (gaml).");
			ldGaml.get();
			m_tSpectraTotal = m_vSpectra.size();
			bContinue = false;
		}
		m_specCondition.use_condition(bState);
	}
/*
 * if the spectrum file did not contain GAML spectra, test the spectrum path parameter
 * for a Matrix Science format file
 */
	if(bContinue)	{
		loadmatrix ldMatrix;
		if(ldMatrix.open(strValue))	{
//			cout << " (mgf).";
			Rprintf(" (mgf).");
			while(ldMatrix.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					lLoaded = 0;
					m_prcLog.log(".\n");
				}
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			bContinue = false;
		}
	}
/*
 * if no Matrix Science format data was found, test for PKL format information
 */
	if(bContinue)	{
		loadpkl ldPkl;
		if(ldPkl.open(strValue))	{
//			cout << " (pkl).";
			Rprintf(" (pkl).");
			while(ldPkl.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					lLoaded = 0;
					m_prcLog.log(".\n");

				}
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			bContinue = false;
		}
	}
#ifdef XMLCLASS
  /*
   * if no PKl format data was found, test for mzxml format information
   * Celui ci sera different des autres.
   * La boucle while(ldPkl.get(spCurrent))
   * est remplacee par le travail de la machine xml provenant de mzxml2other
   * ldMzxml.open instantie la machine
   * tandis que ldMzxml.get l'enclenche
   */
  if(bContinue)	{
    loadmzxml ldMzxml( m_vSpectra, m_specCondition, *m_pScore);
    if(ldMzxml.open(strValue))	{
      //Ce .get est different, il va chercher tous les spectres du fichier.
//			cout << " (mzXML).";
			Rprintf(" (mzXML).");
      ldMzxml.get();
	  m_tSpectraTotal = m_vSpectra.size();
      bContinue = false;
    }

  }
  /*
   * if no MzXML format data was found, test for MzData format information
   * Comme pour mzxml
   */
  if(bContinue)	{
   loadmzml ldMzml( m_vSpectra, m_specCondition, *m_pScore);
   if(ldMzml.open(strValue))  {
      //Ce .get est different, il va chercher tous les spectres du fichier.
// 			cout << " (mzML).";
 			Rprintf(" (mzML).");
     ldMzml.get();
	  m_tSpectraTotal = m_vSpectra.size();
      bContinue = false;
    }
  }
  if(bContinue)	{
   loadmzdata ldMzdata( m_vSpectra, m_specCondition, *m_pScore);
    if(ldMzdata.open(strValue))  {
      //Ce .get est different, il va chercher tous les spectres du fichier.
// 			cout << " (mzData).";
 			Rprintf(" (mzData).");
     ldMzdata.get();
	  m_tSpectraTotal = m_vSpectra.size();
      bContinue = false;
    }

  }
#endif
  /*
   * if no mzxml format data was found, test for dta format information
   */
	if(bContinue)	{
		loaddta ldSpec;
		if(!ldSpec.open(strValue))	{
/*
 * report an error if no DTA information was found: there are no more default file types to check
 */
//			cout << "\nFailed to read spectrum file: " << strValue.c_str() << "\n";
			Rprintf("\nFailed to read spectrum file: %s\n", strValue.c_str());
//			cout << "Most likely: an unsupported data file type:\nUse dta, pkl, mgf, mzdata (v.1.05) or mxzml (v.2.0) files ONLY! (4)\n\n";
			Rprintf("Most likely: an unsupported data file type:\nUse dta, pkl, mgf, mzdata (v.1.05) or mxzml (v.2.0) files ONLY! (4)\n\n");
			//cout.flush();
			m_prcLog.log("error loading spectrum file 4");
			return false;
		}
		else	{
//			cout << " (dta).";
			Rprintf(" (dta).");
			while(ldSpec.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					lLoaded = 0;
					m_prcLog.log(".\n");
				}
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
		}
	}
	strKey = "spectrum, use contrast angle";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		// remove spectra that are redundant
		subtract();
	}
	m_prcLog.log("spectra loaded");
	size_t a = 0;
	while(a < m_vSpectra.size())	{
		if(m_vSpectra[a].m_strDescription.find(":ETD:") != m_vSpectra[a].m_strDescription.npos)	{
			m_vSpectra[a].m_uiType = I_C|I_Z;
		}
		else if(m_vSpectra[a].m_strDescription.find(":CID:") != m_vSpectra[a].m_strDescription.npos)	{
			m_vSpectra[a].m_uiType = I_B|I_Y;
		}
		a++;
	}
	return true;
}
/*
 * This method loads sequences from a bioml sequence data file.
 * It was first introduced in version 2006.04.01
*/
bool mprocess::load_sequences(void)
{
	string strKey = "refine, sequence path";
	string strValue;
	size_t a = 0;
	// Check for input sequences for refinement to be loaded from a bioml file
	// Loading done for each thread
	m_xmlValues.get(strKey,strValue);
	if(strValue.size() > 0)	{
		SAXBiomlHandler saxFile;
		saxFile.setFileName(strValue.c_str());
		saxFile.parse();
		a = 0;
		while(a < saxFile.m_vseqBest.size())	{
			if(m_mapSequences.find(saxFile.m_vseqBest[a].m_tUid) == m_mapSequences.end())	{
				short int iPath = saxFile.m_vseqBest[a].m_siPath;
				string strPath = saxFile.m_vstrPaths[iPath];
				size_t b = 0;
				while(b < m_vstrPaths.size())	{
					if(strPath == m_vstrPaths[b])	{
						saxFile.m_vseqBest[a].m_siPath = (short int)b;
						break;
					}
					b++;
				}
				if(b == m_vstrPaths.size())	{
					m_vstrPaths.push_back(strPath);
					saxFile.m_vseqBest[a].m_siPath = (short int)(m_vstrPaths.size() - 1);
				}
				m_vseqBest.push_back(saxFile.m_vseqBest[a]);
				m_mapSequences.insert(SEQMAP::value_type(saxFile.m_vseqBest[a].m_tUid,saxFile.m_vseqBest[a].m_strSeq));
			}
			a++;
		}
	}
	return true;
}


/*
 * This method forces the use of a particular file format, based on the value of the
 * input parameter "spectrum, path type" parameter. This can be useful if
 * for some reason the file type detection routines do not correctly identify your
 * file type. It may cause bad behavior if your file does not match the file type
 * specified by this parameter. 
 * Currently supported values for the "spectrum, path type parameter" are as follows:
 * "dta", "pkl", "mgf", "gaml", "mzdata", "mzxml"
 */
bool mprocess::spectra_force(string &_t,string &_v)
{
	string strValue = _v;
	string strKey;
	mspectrum spCurrent;
	long lLoaded = 0;
	long lLimit = 2000;
//	cout << " (" << _t.c_str() << ").";
	Rprintf(" (%s).", _t.c_str());
	if(_t == "gaml")	{
		bool bState = m_specCondition.m_bCondition;
		m_specCondition.use_condition(false);
		loadgaml ldGaml(m_vSpectra, m_specCondition, *m_pScore);
		if(ldGaml.open_force(strValue))	{
			ldGaml.get();
			m_tSpectraTotal = m_vSpectra.size();
		}
		m_specCondition.use_condition(bState);
	}
	else if(_t == "cmn")	{
		loadcmn ldCmn;
		if(ldCmn.open(strValue))	{
			while(ldCmn.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					lLoaded = 0;
					m_prcLog.log(".\n");
				}
				m_vSpectra.push_back(spCurrent);
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				m_vSpectra.push_back(spCurrent);
			}
		}
	}

/*
 * if the spectrum file did not contain GAML spectra, test the spectrum path parameter
 * for a Matrix Science format file
 */
	else if(_t == "mgf")	{
		loadmatrix ldMatrix;
		if(ldMatrix.open_force(strValue))	{
			while(ldMatrix.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					lLoaded = 0;
					m_prcLog.log(".");
				}
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
		}
	}
/*
 * if no Matrix Science format data was found, test for PKL format information
 */
	else if(_t == "pkl")	{
		loadpkl ldPkl;
		if(ldPkl.open_force(strValue))	{
			while(ldPkl.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					lLoaded = 0;
					m_prcLog.log(".");
				}
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
		}
	}
#ifdef XMLCLASS
	else if(_t == "mzxml")	{
		loadmzxml ldMzxml( m_vSpectra, m_specCondition,	*m_pScore);
		if(ldMzxml.open_force(strValue))	{
			//Ce .get est different, il	va chercher	tous les spectres du fichier.
			ldMzxml.get();
			m_tSpectraTotal	= m_vSpectra.size();
		}

	}
	else if(_t == "mzml")	{
		loadmzxml ldMzml( m_vSpectra, m_specCondition,	*m_pScore);
		if(ldMzml.open_force(strValue))	{
			//Ce .get est different, il	va chercher	tous les spectres du fichier.
			ldMzml.get();
			m_tSpectraTotal	= m_vSpectra.size();
		}

	}
	else if(_t == "mzdata")	{
		loadmzdata ldMzdata( m_vSpectra, m_specCondition, *m_pScore);
		if(ldMzdata.open_force(strValue))	 {
			//Ce .get est different, il	va chercher	tous les spectres du fichier.
			ldMzdata.get();
			m_tSpectraTotal	= m_vSpectra.size();
		}
	}
#endif
	else if(_t == "dta")	{
		loaddta ldSpec;
		if(ldSpec.open_force(strValue))	{
			while(ldSpec.get(spCurrent))	{
				m_tSpectraTotal++;
				lLoaded++;
				if(lLoaded == lLimit)	{
//					cout << ".";
					Rprintf(".");
					//cout.flush();
					m_prcLog.log(".");
					lLoaded = 0;
				}
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
			if(spCurrent.m_vMI.size() > 0)	{
				m_tSpectraTotal++;
				if(m_specCondition.condition(spCurrent, *m_pScore))
					m_vSpectra.push_back(spCurrent);
			}
		}
	}
	else	{
//			cout << "\n" << "The file type \"" << _t.c_str() << " is not supported.\n";
			Rprintf("\nThe file type \"%s\" is not supported.\n", _t.c_str());
//			cout << "Supported values: pkl, dta, mgf, gaml, mzxml, mzdata\n";
			Rprintf("Supported values: pkl, dta, mgf, gaml, mzxml, mzdata\n");
			//cout.flush();
			m_prcLog.log("error loading forced spectrum file 5");
			return false;
	}
	strKey = "spectrum, use contrast angle";
	m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
		// remove spectra that are redundant
		subtract();
	}
	m_prcLog.log("spectra loaded");
	return true;
}
//
// This method is an experiment to decrease the number of redundant spectra that
// enter the calculation. It calculates the inner product (using ::dot) between
// two spectrum vectors and then calculates the cosine of the opening
// angle between the vectors. If the arccos is greater than the parameter
// "spectrum, contrast angle", then the spectrum is considered to be non-redundant and kept.
//
//	The inner product is calculated for parent ion masses within 1000 ppm of each
//  other. If the arccos between two spectrum vectors is below the contrast angle parameter,
//  then the two spectra's summed intensities are compared. The spectrum with the lowest
//  summed intensity is placed in a list of rejected spectra and the process is repeated.
//  Once all of the spectra to be rejected have been marked, a new mspectrum vector is created
//  and copied into m_vSpectra.
//
bool mprocess::subtract(void)
{
	if(m_vSpectra.size() == 0)
		return false;
//	cout << "+";
	Rprintf("+");
	//cout.flush();
	string strValue;
	string strKey = "spectrum, fragment mass error";
	m_xmlValues.get(strKey,strValue);
	if (strValue.empty()) {
		// try old parameter name
		strKey = "spectrum, fragment monoisotopic mass error";
		m_xmlValues.get(strKey,strValue);
	}
	float fRes = (float)atof(strValue.c_str());
	if(fRes <= 0.0)	{
		fRes = 0.5;
	}
	bool bType = true;
	strKey = "spectrum, fragment mass error units";
	m_xmlValues.get(strKey,strValue);
	if (strValue.empty()) {
		// try old parameter name
		strKey = "spectrum, fragment monoisotopic mass error units";
		m_xmlValues.get(strKey,strValue);
	}
	if(strValue != "Daltons")	{
		bType = false;
	}
	// retrieve contrast angle
	strKey = "spectrum, contrast angle";
	m_xmlValues.get(strKey,strValue);
	double dLimit = atof(strValue.c_str());
	// apply limits to the value of the angle
	if(dLimit < 0.0)	{
		dLimit = 0.0;
	}
	if(dLimit > 90.0)	{
		dLimit = 90.0;
	}
	dLimit = cos(dLimit * 3.1415/180.0);
	const size_t tSpectra = m_vSpectra.size();
	size_t a = 0;
	size_t b = 0;
	double dValue = 0.0;
	long lCount = 0;
	double dDot1 = 0.0;
	vector<double> vdDot1;
	// calculate the length of all of the spectrum vectors
	vector<mi>::iterator itM = m_vSpectra[b].m_vMI.begin();
	vector<mi>::iterator itMEnd = m_vSpectra[b].m_vMI.end(); 
	while(b < m_vSpectra.size())	{
		dDot1 = 0.0;
		itM = m_vSpectra[b].m_vMI.begin();
		itMEnd = m_vSpectra[b].m_vMI.end();
		while(itM != itMEnd)	{
			dDot1 += itM->m_fI*itM->m_fI;
			itM++;
		}
		vdDot1.push_back(sqrt(dDot1));
		b++;
	}
	a  = 0;
	long c = 0;
	long lMax = 0;
	set<size_t> vDelete;
	double dMax = 0.0;
	size_t tLastMax;
	double dLastSum;
	size_t tLastA;
	float fMH = 0.0;
	const double dFactor = 0.001;
	double dDelta = 1.0;
	vector<mspectrum>::iterator itA = m_vSpectra.begin();
	vector<mspectrum>::iterator itB = itA;
	// calculated the dot products between the spectrum vector pairs
	while(a < tSpectra)	{
		// calculate only the pairs below the diagonal on the pairs matrix
		// the values on the diagonal are all 1.0
		// the values above the diagonal repeat the below diagonal elements, d(i,j) = d(j,i)
		b = a+1;
		itB = itA;
		itB++;
		lCount = 0;
		dMax  = itA->m_vdStats[0];
		tLastMax = itA->m_tId;
		fMH = (float)itA->m_dMH;
		dDelta = dFactor*fMH;
		dLastSum = itA->m_vdStats[0];
		tLastA = 0;
		while(b < tSpectra)	{
			// only calculate dValue if the parent ion mass is appropriate
			// and if the spectrum has not already been rejected
			if(fabs(fMH - (float)itB->m_dMH) < dDelta && vDelete.find(itB->m_tId) == vDelete.end())	{
				dValue = dot(a,b,fRes,bType)/(vdDot1[a]*vdDot1[b]);
				if(dValue > dLimit)	{
					if(dMax < itB->m_vdStats[0])	{
						dMax = itB->m_vdStats[0];
						tLastA = b;
						vDelete.insert(tLastMax);
						tLastMax = itB->m_tId;
					}
					else	{
						vDelete.insert(itB->m_tId);
						dLastSum += itB->m_vdStats[0];
						tLastA = a;
					}
					lCount++;
				}
			}
			b++;
			itB++;
		}
		if(tLastA != 0)	{
			m_vSpectra[tLastA].m_vdStats[0] += dLastSum;
		}
		if(lCount > lMax)	{	
			lMax = lCount;
		}
		if(c > 1000)	{
			if(m_lThread == 0 || m_lThread == 0xFFFFFFFF)	{
//				cout << "+";
				Rprintf("+");
				//cout.flush();
			}
			c = 0;
		}
		c++;
		a++;
		itA++;
		while(a < tSpectra && vDelete.find(itA->m_tId) != vDelete.end())	{
			a++;
			itA++;
			if(c > 1000)	{
//				cout << "+";
				Rprintf("+");
				//cout.flush();
				c = 0;
			}
			c++;
		}
	}
	vector<mspectrum>::iterator itStart = m_vSpectra.begin();
	vector <mspectrum> vTemp;
	vTemp.reserve(m_vSpectra.size() - vDelete.size()+1);
	m_tContrasted = 0;
	// create a new list of spectra, vTemp and replace m_vSpectra with 
	// those spectra. This is much faster than deleting the spectra
	// from the intact m_vSpectra list
	while(itStart != m_vSpectra.end())	{
		if(vDelete.find(itStart->m_tId) == vDelete.end())	{
			vTemp.push_back(*itStart);
		}
		itStart++;
	}
	m_tContrasted = (double)(m_vSpectra.size() - vTemp.size());
	m_vSpectra.clear();
	m_vSpectra.reserve(vTemp.size() + 1);
	m_vSpectra = vTemp;
	return true;
}

/*
 * taxonomy uses the taxonomy information in the input XML file to load
 * the msequenceServer member object with file path names to the required
 * sequence list files (FASTA format only in the initial release). If these
 */
bool mprocess::taxonomy()
{
	string strValue;
	string strKey = "list path, taxonomy information";
	m_xmlValues.get(strKey,strValue);
	string strTaxonPath = strValue;
	strKey = "protein, taxon";
	m_xmlValues.get(strKey,strValue);
	long lReturn = m_svrSequences.load_file(strTaxonPath,strValue);
/*
 * return false if the load_file method fails
 */
	if(lReturn == 1)	{
//		cout << "\nThe taxonomy parameter file \"" << strTaxonPath.c_str();
		Rprintf("\nThe taxonomy parameter file \"%s", strTaxonPath.c_str());
//		cout << "\" could not be found.\nCheck your settings and try again.\n";
		Rprintf("\" could not be found.\nCheck your settings and try again.\n");
		return false;
	}
	else if(lReturn == 2)	{
//		cout << "\nThe taxonomy parameter file \"" << strTaxonPath.c_str();
		Rprintf("\nThe taxonomy parameter file \"%s", strTaxonPath.c_str());
//		cout << "\" did not contain the value \"" << strValue.c_str() << "\".\nCheck your settings and try again.\n";
		Rprintf("\" did not contain the value \"%s\".\nCheck your settings and try again.\n", strValue.c_str());
		return false;
	}
	else if(lReturn == 3)	{
//		cout << "\nThe taxonomy parameter file \"" << strTaxonPath.c_str();
		Rprintf("\nThe taxonomy parameter file \"%s", strTaxonPath.c_str());
//		cout << "\" contained incorrect entries\nfor the protein sequence files associated with the name: \"" << strValue.c_str() << "\".\nCheck the file names in the taxonomy file and try again.\n";
		Rprintf("\" contained incorrect entries\nfor the protein sequence files associated with the name: \"%s\".\nCheck the file names in the taxonomy file and try again.\n", strValue.c_str());
		return false;
	}
	return true;
}
//
// new code, version 2007/04/01, Ron Beavis
// this section finds the file names associated with the taxonomy specification
// associated with Single Amino acid Polymorphisms (SAPs). The information
// is parsed using SAXSapHandler and stored in the m_mapSap object of mscore.
// This object will be used to control the state machine that generates the
// peptides associated with these polymorphisms
// When the mprocess pointer is NULL, the SAP object is loaded from the object
// in the original file. If the pointer is not NULL, it is assumed that the mprocess
// object contains the necessary SAP information, which is simply copied from memory.
//
bool mprocess::load_saps(mprocess *_p)
{
	string strValue;
	string strKey = "list path, taxonomy information";
	m_xmlValues.get(strKey,strValue);
	string strTaxonPath = strValue;
	strKey = "protein, taxon";
	m_xmlValues.get(strKey,strValue);
	XmlTaxonomy xmlTax;
	string strType = "saps"; //this is the format attribute for the <file> objects in the tax file
	if(!xmlTax.load(strTaxonPath,strValue,strType))
		return true;
	// bail out without error if there aren't any SAP files specified.
	size_t a = 0;
	map<string,multimap<int,prSap> >::iterator itMap;
	map<string,multimap<int,prSap> >::iterator itValue;
	pair<string,multimap<int,prSap> > pairMap;
	multimap<int,prSap>::iterator itMulti;
	// copy from _p if it exists
	if(_p != NULL)	{
		itMap = _p->m_pScore->m_Sap.m_mapSap.begin();
		// load the mscore m_Sap object with the new values
		while(itMap != _p->m_pScore->m_Sap.m_mapSap.end())	{
			pairMap.first = itMap->first;
			pairMap.second.clear();
			itValue = m_pScore->m_Sap.m_mapSap.find(pairMap.first);
			if(itValue == m_pScore->m_Sap.m_mapSap.end())	{
				m_pScore->m_Sap.m_mapSap.insert(pairMap);
				itValue = m_pScore->m_Sap.m_mapSap.find(pairMap.first);
			}
			itMulti = itMap->second.begin();
			while(itMulti != itMap->second.end())	{
				itValue->second.insert(*itMulti);
				itMulti++;
			}
			itMap++;
		}
		return true;
	}
	m_vstrSaps.clear();
	// obtain the lists from files, if no valid mprocess value was passed through _p
	while(a < xmlTax.m_vstrPaths.size())	{
		ifstream ifTest;
		ifTest.open(xmlTax.m_vstrPaths[a].c_str());
		if(!ifTest.fail())	{
			m_vstrSaps.push_back(xmlTax.m_vstrPaths[a]);
			ifTest.close();
		}
		ifTest.clear();
		a++;
	}
	if(!m_vstrSaps.empty())	{
//		cout << " loaded.\nLoading SAPs ";
		Rprintf(" loaded.\nLoading SAPs ");
		//cout.flush();
	}
	a = 0;
	while(a < m_vstrSaps.size())	{
		SAXSapHandler sapXml;
		sapXml.setFileName(xmlTax.m_vstrPaths[a].data());
		sapXml.parse();
		itMap = sapXml.m_mapSap.begin();
		// load the mscore m_Sap object with the new values
		while(itMap != sapXml.m_mapSap.end())	{
			pairMap.first = itMap->first;
			pairMap.second.clear();
			itValue = m_pScore->m_Sap.m_mapSap.find(pairMap.first);
			if(itValue == m_pScore->m_Sap.m_mapSap.end())	{
				m_pScore->m_Sap.m_mapSap.insert(pairMap);
				itValue = m_pScore->m_Sap.m_mapSap.find(pairMap.first);
			}
			itMulti = itMap->second.begin();
			while(itMulti != itMap->second.end())	{
				itValue->second.insert(*itMulti);
				itMulti++;
			}
			itMap++;
		}
//		cout << ".";
		Rprintf(".");
		//cout.flush();
		a++;
	}
	return true;
}

//
// new code, version 2007/04/01, Ron Beavis
// this section finds the file names associated with the taxonomy specification
// associated with Single Amino acid Polymorphisms (SAPs). The information
// is parsed using SAXSapHandler and stored in the m_mapSap object of mscore.
// This object will be used to control the state machine that generates the
// peptides associated with these polymorphisms
// When the mprocess pointer is NULL, the SAP object is loaded from the object
// in the original file. If the pointer is not NULL, it is assumed that the mprocess
// object contains the necessary SAP information, which is simply copied from memory.
//
bool mprocess::load_annotation(mprocess *_p)
{
	string strValue;
	string strKey = "list path, taxonomy information";
	m_xmlValues.get(strKey,strValue);
	string strTaxonPath = strValue;
	strKey = "protein, taxon";
	m_xmlValues.get(strKey,strValue);
	XmlTaxonomy xmlTax;
	string strType = "mods"; //this is the format attribute for the <file> objects in the tax file
	if(!xmlTax.load(strTaxonPath,strValue,strType))
		return true;
	// bail out without error if there aren't any mods files specified.
	size_t a = 0;
	map<string,string>::iterator itMods;
	if(_p != NULL)	{
		itMods = _p->m_mapAnnotation.begin();
		// load the m_mapAnnotation object with the new values
		while(itMods != _p->m_mapAnnotation.end())	{
			m_mapAnnotation[itMods->first] = itMods->second;
			itMods++;
		}
		return true;
	}
	m_vstrMods.clear();
	// obtain the lists from files, if no valid mprocess value was passed through _p
	while(a < xmlTax.m_vstrPaths.size())	{
		ifstream ifTest;
		ifTest.open(xmlTax.m_vstrPaths[a].c_str());
		if(!ifTest.fail())	{
			m_vstrMods.push_back(xmlTax.m_vstrPaths[a]);
			ifTest.close();
		}
		ifTest.clear();
		a++;
	}
	if(!m_vstrMods.empty())	{
//		cout << " loaded.\nLoading annotation ";
		Rprintf(" loaded.\nLoading annotation ");
		//cout.flush();
	}
	a = 0;
	while(a < m_vstrMods.size())	{
		SAXModHandler modXml;
		modXml.setFileName(xmlTax.m_vstrPaths[a].data());
		modXml.parse();
		itMods = modXml.m_mapMod.begin();
		// load the mscore m_Sap object with the new values
		while(itMods != modXml.m_mapMod.end())	{
			m_mapAnnotation[itMods->first] = itMods->second;
			itMods++;
		}
//		cout << ".";
		Rprintf(".");
		//cout.flush();
		a++;
	}
	return true;
}

bool mprocess::serialize(void)
{	
	if(!m_bSerialize)	{
		return true;
	}
	string strKey = "output, path";
	string strValue;
	if(!m_xmlValues.get(strKey,strValue))	{
		return false;
	}
	FILE *pFile = fopen(strValue.c_str(),"wb");
	if(!pFile)	{
//		cout << "Warning: serialization did not occur.\n";
		Rprintf("Warning: serialization did not occur.\n");
		//cout.flush();
		return false;
	}
	vector<mspectrum>::iterator itS = m_vSpectra.begin();
	vector<mspectrum>::iterator itE = m_vSpectra.end();
	vector<mi>::iterator itMI;
	vector<mi>::iterator itMIE;
	size_t tLength = m_vSpectra.size();
	size_t elements = fwrite((const void*)&tLength,sizeof(size_t),1,pFile);
	while(itS != itE)	{
		tLength = itS->m_vMI.size();
		elements = fwrite((const void*)&(itS->m_tId),sizeof(size_t),1,pFile);
		elements = fwrite((const void*)&tLength,sizeof(size_t),1,pFile);
		itMI = itS->m_vMI.begin();
		itMIE = itS->m_vMI.end();
		while(itMI != itMIE)	{
			elements = fwrite((const void*)&(itMI->m_fM),sizeof(float),1,pFile);
			elements = fwrite((const void*)&(itMI->m_fI),sizeof(float),1,pFile);

			elements++; /* fool the compiler */

			itMI++;
		}
		itS++;
	}
	fclose(pFile);
	return true;
}

bool mprocess::restore(void)
{
	if(!m_bSerialize)	{
		return true;
	}
	string strKey = "output, path";
	string strValue;
	if(!m_xmlValues.get(strKey,strValue))	{
		return false;
	}
	FILE *pFile = fopen(strValue.c_str(),"rb");
	if(!pFile || feof(pFile))	{
		Rprintf("Warning: could not find serialization file \"%s\", spectrum restoration not performed.\n", strValue.c_str());
		return false;
	}
	vector<mspectrum>::iterator itS = m_vSpectra.begin();
	vector<mspectrum>::iterator itE = m_vSpectra.end();
	size_t tLength = 0;
	size_t elements = fread((void *)&tLength,sizeof(size_t),1,pFile);

	if(tLength == 0 || feof(pFile))	{
		Rprintf("Warning: could not find serialization file \"%s\" appears to be corrupt.\n", strValue.c_str());
		fclose(pFile);
		return false;
	}
	size_t tL = 0;
	size_t tId = 0;
	size_t a = 0;
	size_t b = 0;
	map<size_t,size_t> mapId;
	pair<size_t,size_t> pairId;
	while(itS < itE)	{
		pairId.first = itS->m_tId;
		pairId.second = a;
		mapId.insert(pairId);
		a++;
		itS++;
	}
	vector<mi> vMI;
	mi miV;
	float fV;
	a = 0;
	while(a < tLength && !feof(pFile))	{
		vMI.clear();
		elements = fread((void *)&tId,sizeof(size_t),1,pFile);
		elements = fread((void *)&tL,sizeof(size_t),1,pFile);

		b = 0;
		while(b < tL && !feof(pFile))	{
			elements = fread((void *)&fV,sizeof(float),1,pFile);
			miV.m_fM = fV;
			elements = fread((void *)&fV,sizeof(float),1,pFile);

			elements++; /* use it */

			miV.m_fI = fV;
			vMI.push_back(miV);
			b++;
		}
		if(mapId.find(tId) != mapId.end())	{
			m_vSpectra[mapId.find(tId)->second].m_vMI = vMI;
		}
		a++;
	}
	fclose(pFile);
	return true;
}

bool mprocess::removeMI(void)
{
	if(!m_bSerialize)	{
		return true;
	}
	size_t a = 0;
	while(a < m_vSpectra.size())	{
		m_vSpectra[a].clear_intensity_values();
		a++;
	}
	return true;
}
