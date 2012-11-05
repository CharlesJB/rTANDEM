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
#include "xmltaxonomy.h"
#include "mscore.h"
#include "mprocess.h"
#include "saxbiomlhandler.h"
#include "saxsaphandler.h"
#include "mbiomlreport.h"
#include "mrefine.h"

// #define TANDEM_EXACT 1

/*
 * global less than operators to be used in sort operations
 */
p3mprocess::p3mprocess(void)
{
	string strKey = "process, version";
	string strValue = "X! P3 ";
	strValue += VERSION;
	m_xmlPerformance.set(strKey,strValue);
}

p3mprocess::~p3mprocess(void)
{
}
/*
 * taxonomy uses the taxonomy information in the input XML file to load
 * the msequenceServer member object with file path names to the required
 * sequence list files (FASTA format only in the initial release). If these
 */
bool p3mprocess::taxonomy()
{
#ifdef X_P3
	string strValue;
	string strKey = "list path, taxonomy information";
	m_xmlValues.get(strKey,strValue);
	string strTaxonPath = strValue;
	strKey = "protein, taxon";
	m_xmlValues.get(strKey,strValue);
	long lReturn = m_svrSequences.u_load_file(strTaxonPath,strValue);
/*
 * return false if the load_file method fails
 */
	if(lReturn == 2)	{
		cout << "\nThe taxonomy parameter file \"" << strTaxonPath.c_str();
		cout << "\" did not contain the value \"" << strValue.c_str() << "\".\nCheck your settings and try again.\n";
		return false;
	}
	if(lReturn == 1)	{
		cout << "\nThe taxonomy parameter file \"" << strTaxonPath.c_str();
		cout << "\" could not be found.\nCheck your settings and try again.\n";
		return false;
	}
	return true;
#endif
	return true;
}
/*
 * merge_spectra takes an mspectrum vector from an external source and merges it
 * with the m_vSpectra vector. this method is used to combine information obtained
 * from other threads with this mprocess object
 */
bool p3mprocess::merge_spectra()
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
			m_vSpectra[a].reset();
			a++;
		}
	}
	return true;
}

bool p3mprocess::merge_spectra(vector<mspectrum> &_s)
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
			a++;
			if(m_bUseHomologManagement && _s[a].m_vseqBest.size() > 5)	{
				_s[a].m_vseqBest.erase(_s[a].m_vseqBest.begin()+5,_s[a].m_vseqBest.end());
			}
		}
		a = 0;
		double dExpect = 1.0;
		SEQMAP::iterator itValue;
		while(a < _s.size())	{
			b = 0;
			dExpect = (double)_s[a].m_hHyper.expect_protein(_s[a].m_fHyper);
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
			_s[a].reset();
			a++;
		}
	}
	return true;
}
/*
 * merge_map adds new values to the existing m_mapSequences map.
*/
bool p3mprocess::merge_map(SEQMAP &_s,DESMAP &_d)
{
	SEQMAP::iterator itValue = _s.begin();
	SEQMAP::iterator itEnd = _s.end();
	while(itValue != itEnd)	{
		if(m_mapSequences.find(itValue->first) == m_mapSequences.end())	{
			m_mapSequences.insert(*itValue);
		}
		itValue++;
	}
	DESMAP::iterator itV = _d.begin();
	DESMAP::iterator itE = _d.end();
	while(itV != itE)	{
		if(m_mapDesc.find(itV->first) == m_mapDesc.end())	{
			m_mapDesc.insert(*itV);
		}
		itV++;
	}
	return true;
}

/*
 * merge_map adds new values to the existing m_mapSequences map.
*/
bool p3mprocess::merge_map(SEQMAP &_s)
{
	return mprocess::merge_map(_s);
}

bool p3mprocess::load_sequences()
{
#ifdef X_P3
	vector<msequence> vSeq;
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
	if(lReturn == 2)	{
		cout << "\nThe taxonomy parameter file \"" << strTaxonPath.c_str();
		cout << "\" did not contain the value \"" << strValue.c_str() << "\".\nCheck your settings and try again.\n";
		return false;
	}
	if(lReturn == 1)	{
		cout << "\nThe taxonomy parameter file \"" << strTaxonPath.c_str();
		cout << "\" could not be found.\nCheck your settings and try again.\n";
		return false;
	}
	m_svrSequences.u_maps(m_mapDesc,vSeq);
	size_t a = 0;
	DESMAP mapDes;
	m_mapSequences.clear();
	while(a < vSeq.size())	{
		mapDes[vSeq[a].m_strDes] = a;
		m_mapSequences[vSeq[a].m_tUid] = vSeq[a].m_strSeq;
		a++;
	}
	m_vstrPaths = m_svrSequences.m_vstrPaths;
	a = 0;
	while(a < m_vseqBest.size())	{
		m_vseqBest[a] = vSeq[mapDes[m_vseqBest[a].m_strDes]];
		a++;
	}
	strKey = "refine, sequence path";
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
#endif
	return true;
}

bool p3mprocess::clean_sequences(void)
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
	DESMAP::iterator itV = m_mapDesc.begin();
	while(itV != m_mapDesc.end())	{
		itMap = mapValue.find((*itV).second);
		if(itMap == mapValue.end())	{
			m_mapDesc.erase(itV);
			itV = m_mapDesc.begin();
		}
		else	{
			itV++;
		}
	}
	return true;
}

/*
 * create_score takes an msequence object and a start and an end sequence position
 * and saves that msequence, its mdomain and scoring in the spectrum that has just
 * been scored. equivalent msequence objects (based on their hyper score)
 * are stored sequentially in a vector in the mspectrum object
 */
bool p3mprocess::create_score(const msequence &_s,const size_t _v,const size_t _w,const long _m,bool _p)
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
	while(lCount < m_pScore->m_State.m_lEqualsS)	{
		a = m_pScore->m_State.m_plEqualsS[lCount];
		lCount++;
		lIonCount = 0;
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
			if(m_bCrcCheck && _p && m_vSpectra[a].m_hHyper.m_ulCount < 400)	{
				m_bPermute = true;
			}
			else if(m_vSpectra[a].m_dMH > 3000.0 && m_vSpectra[a].m_hHyper.m_ulCount < 400)	{
				m_bPermuteHigh = true;
			}
		}
		bIonCheck = false;
		if(lIonCount > m_lIonCount)	{
			if( (m_pScore->m_plCount[mscore::S_A] || m_pScore->m_plCount[mscore::S_B] || m_pScore->m_plCount[mscore::S_C]) && 
				(m_pScore->m_plCount[mscore::S_X] || m_pScore->m_plCount[mscore::S_Y] || m_pScore->m_plCount[mscore::S_Z]))	{
				bIonCheck = true;
			}
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
				domValue.m_lMissedCleaves = (unsigned char)_m;
				domValue.m_bUn = m_bUn;
				unsigned char lType = 1;
				unsigned long sType = 1;
				while(lType < m_pScore->m_lType+1)	{
					m_vSpectra[a].m_mapScore[lType] = m_pScore->m_pfScore[sType];
					m_vSpectra[a].m_mapCount[lType] = (m_pScore->m_plCount[sType]);
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
					m_vSpectra[a].m_mapCount[lType] = (m_pScore->m_plCount[sType]);
					lType *= 2;
					sType++;
				}
				seqValue = _s;
				seqValue.m_strSeq = " ";
				m_mapSequences.insert(SEQMAP::value_type(_s.m_tUid,_s.m_strSeq));
				m_mapDesc.insert(DESMAP::value_type(_s.m_strDes,_s.m_tUid));
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
					m_vSpectra[a].m_mapCount[lType] = m_pScore->m_plCount[sType];
					lType *= 2;
					sType++;
				}
				seqValue = _s;
				seqValue.m_strSeq = " ";
				m_mapSequences.insert(SEQMAP::value_type(_s.m_tUid,_s.m_strSeq));
				m_mapDesc.insert(DESMAP::value_type(_s.m_strDes,_s.m_tUid));
				seqValue.m_vDomains.clear();
				seqValue.m_vDomains.push_back(domValue);
				seqValue.format_description();
				m_vSpectra[a].m_tCurrentSequence = _s.m_tUid;
				seqValue.m_iRound = m_iCurrentRound;
				m_vSpectra[a].m_vseqBest.clear();
				m_vSpectra[a].m_vseqBest.push_back(seqValue);
			}
		}
		else if (_p && bIonCheck && fHyper > m_vSpectra[a].m_fHyperNext)
		{
			m_vSpectra[a].m_fScoreNext = fScore;
			m_vSpectra[a].m_fHyperNext = fHyper;
		}
	}
	return true;
}
