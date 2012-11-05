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

// File version: 2003-07-01

/*
 * msequtilities contains a set of constants and arrays for scoring and setting the molecular
 * mass of a peptide sequence. Biemann notation is used for fragment ions.
 */

#include "stdafx.h"
#include "msequence.h"
#include "msequtilities.h"

/*
 * sets the default values for the necessary constants. see msequtilities.h for descriptions 
 */

msequtilities::msequtilities(masscalc::massType _t)
	: m_calc(_t)
{
	m_dCleaveCdefault = m_calc.calcMass("OH");
	m_dCleaveNdefault = m_calc.calcMass("H");
	m_dCleaveC = m_dCleaveCdefault;
	m_dCleaveN = m_dCleaveNdefault;
	const unsigned long lValue = 128;
	m_pfAaMass = new float[lValue];
	memset(m_pfAaMass,0,lValue*sizeof(float));
	m_pdAaMass = new double[lValue];
	memset(m_pdAaMass,0,lValue*sizeof(double));
	m_pdAaPrompt = new double[lValue];
	m_pdAaMod = new double[lValue];
	m_pdAaFullMod = new double[lValue];
	m_pfAScore = new float[lValue];
	m_pfBScore = new float[lValue];
	m_pfYScore = new float[lValue];
	m_pfXScore = new float[lValue];
	m_pfA18Score = new float[lValue];
	m_pfB18Score = new float[lValue];
	m_pfY18Score = new float[lValue];
	m_pfX18Score = new float[lValue];
	m_pfA17Score = new float[lValue];
	m_pfB17Score = new float[lValue];
	m_pfY17Score = new float[lValue];
	m_pfX17Score = new float[lValue];
	long a = 0;
	while(a < lValue)	{
		m_pdAaMod[a] = 0.0;
		m_pdAaPrompt[a] = 0.0;
		m_pdAaFullMod[a] = 0.0;
		m_pfXScore[a] = 1.0;
		m_pfYScore[a] = 1.0;
		m_pfAScore[a] = 1.0;
		m_pfBScore[a] = 1.0;
		m_pfX18Score[a] = 1.0;
		m_pfY18Score[a] = 1.0;
		m_pfA18Score[a] = 1.0;
		m_pfB18Score[a] = 1.0;
		m_pfX17Score[a] = 1.0;
		m_pfY17Score[a] = 1.0;
		m_pfA17Score[a] = 1.0;
		m_pfB17Score[a] = 1.0;
		a++;
	}
	set_aa();
	m_fNT = 0.0;
	m_fCT = 0.0;
	// Mass of a proton, i.e. +1 positive charge, hydrogen atom without its electron.
	// See http://antoine.frostburg.edu/chem/senese/101/atoms/index.shtml
	m_dProton = 1.007276;
	m_dHydrogen = 1.007825035;
	m_dWater = m_calc.calcMass("H2O");
	m_dAmmonia = m_calc.calcMass("NH3");
	m_dA = -m_calc.calcMass("CO");
	m_dB = 0.0;
	m_dC = m_calc.calcMass("NH3");
	m_dY = m_calc.calcMass("H2O");;
	m_dX = m_calc.calcMass("CO2");
	// changed from m_dZ = m_dY - m_calc.calcMass("NH3") in 2006.02.01
	m_dZ = m_dY - m_calc.calcMass("NH2");
	m_bComplete = false;
	m_bPotential = false;
	m_vMotifs.clear();
	m_mapMotifs.clear();
	m_bPotentialMotif = false;
	m_bSequenceMods = false;
	m_mapMods.clear();
	m_bIsModified = false;
	m_bPrompt = false;
	m_bPhosphoSerine = false;
	m_bPhosphoThreonine = false;
}
/*
 * clean up the arrays created using the new operator
 */

msequtilities::~msequtilities(void)
{
	if(m_pfAaMass != NULL)
		delete m_pfAaMass;
	if(m_pdAaMass != NULL)
		delete m_pdAaMass;
	if(m_pdAaMod != NULL)
		delete m_pdAaMod;
	if(m_pdAaPrompt != NULL)
		delete m_pdAaPrompt;
	if(m_pdAaFullMod != NULL)
		delete m_pdAaFullMod;
	if(m_pfAScore != NULL)
		delete m_pfAScore;
	if(m_pfBScore != NULL)
		delete m_pfBScore;
	if(m_pfYScore != NULL)
		delete m_pfYScore;
	if(m_pfXScore != NULL)
		delete m_pfXScore;
	if(m_pfA18Score != NULL)
		delete m_pfA18Score;
	if(m_pfB18Score != NULL)
		delete m_pfB18Score;
	if(m_pfY18Score != NULL)
		delete m_pfY18Score;
	if(m_pfX18Score != NULL)
		delete m_pfX18Score;
	if(m_pfA17Score != NULL)
		delete m_pfA17Score;
	if(m_pfB17Score != NULL)
		delete m_pfB17Score;
	if(m_pfY17Score != NULL)
		delete m_pfY17Score;
	if(m_pfX17Score != NULL)
		delete m_pfX17Score;
}

bool msequtilities::add_mod(const char _c,const size_t _v)
{
	m_mapMotifMods[_c] = _v;
	return true;
}

bool msequtilities::clear_motifs(const bool _b)	
{
	map<char,size_t>::iterator itValue = m_mapMotifMods.begin();
	while(itValue != m_mapMotifMods.end())	{
		m_pdAaMod[itValue->first+32] = 0.0;
		m_pdAaPrompt[itValue->first+32] = 0.0;
		itValue++;
	}
	if(_b)	{
		m_mapMotifMods.clear();
	}
	return true;
}

bool msequtilities::set_motifs(void)	
{
	map<char,size_t>::iterator itValue = m_mapMotifMods.begin();
	while(itValue != m_mapMotifMods.end())	{
		m_pdAaMod[itValue->first+32] = m_vMotifs[itValue->second].m_fMass;
		m_pdAaPrompt[itValue->first+32] = m_vMotifs[itValue->second].m_fMassPrompt;
		itValue++;
	}
	return true;
}

bool msequtilities::synthesis(const bool _b)
{
	if(_b)	{
		set_y('P',5.0);
		set_b('D',5.0);
		set_b('N',2.0);
		set_b('V',3.0);
		set_b('E',3.0);
		set_b('Q',2.0);
		set_b('I',3.0);
		set_b('L',3.0);
	}
	else	{
		char a = 0;
		while(a < 127)	{
			set_y(a,1.0);
			set_b(a,1.0);
			a++;
		}
	}
	return true;
}

/*
 * imports a list of complete residue modifications, using a string with the following format:
 * m1@r1,m2@r2,m3@r3, ...
 * where mX is the mass of the modification, in Daltons, and rX is the single letter code for the
 * modified residue. there is no limit to the length of the string (but there are only 20 amino acids)
 * NOTE: lower case index values are reserved for potential modifications and may not be used here
 */

bool msequtilities::modify_all(string &_s)
{
	size_t tStart = 'A';
	while(tStart < 'Z'+1)	{
		m_pdAaFullMod[tStart] = 0.0;
		tStart++;
	}
	tStart = 'a';
	while(tStart < 'z'+1)	{
		m_pdAaFullMod[tStart] = 0.0;
		tStart++;
	}
	m_pdAaFullMod[']'] = 0.0;
	m_pdAaFullMod['['] = 0.0;
	tStart = 0;
	if(_s.empty())
		return false;
	size_t tAt = 0;
	double fValue;
	string strValue = _s.substr(tStart,_s.size()-tStart);
	fValue = atof(strValue.c_str());
	while(fValue != 0.0)	{
		m_bComplete = true;
		tAt = _s.find('@',tStart);
		if(tAt == _s.npos)
			break;
		tAt++;
		if(isalpha(_s[tAt]))	{
			m_pdAaFullMod[_s[tAt]] = fValue;
			m_pdAaFullMod[_s[tAt]+32] = fValue;
		}
		else	{
			m_pdAaFullMod[_s[tAt]] = fValue;
		}
		tStart = _s.find(',',tAt);
		if(tStart == _s.npos)
			break;
		tStart++;
		strValue = _s.substr(tStart,_s.size()-tStart);
		fValue = atof(strValue.c_str());
	}
	return true;
}
/*
 * imports a list of potential residue modifications, using a string with the following format:
 * m1@r1,m2@r2,m3@r3, ...
 * where mX is the mass of the modification, in Daltons, and rX is the single letter code for the
 * modified residue. there is no limit to the length of the string (but there are only 20 amino acids)
 * NOTE: lower case index values are reserved for potential modifications
 */

bool msequtilities::modify_maybe(string &_s)
{
	size_t tStart = 'a';
	while(tStart < 'z'+1)	{
		m_pdAaMod[tStart] = 0.0;
		m_pdAaPrompt[tStart] = 0.0;
		tStart++;
	}
	m_pdAaMod['['] = 0.0;
	m_pdAaMod[']'] = 0.0;

	m_bPotential = false;
	if(_s.empty())
		return false;
	tStart = 0;
	size_t tAt = 0;
	size_t tColon = 0;
	double fValue = 0.0;
	double fPrompt = 0.0;
	string strValue = _s.substr(tStart,_s.size()-tStart);
	m_strDefaultMaybe = strValue;
	fValue = atof(strValue.c_str());
	char cRes = '\0';
	double fSum = 0.0;
	while(fValue != 0.0)	{
		m_bPotential = true;
		tAt = _s.find('@',tStart);
		if(tAt == _s.npos)
			break;
		tColon = _s.find(':',tStart);
		fPrompt = 0.0;
		if(tColon != _s.npos && tColon < tAt)	{
			fPrompt = atof(_s.substr(tColon+1,tAt - tColon).c_str());
		}
		tAt++;
		cRes = _s[tAt];
		if(cRes >= 'A' && cRes <= 'Z')	{
			cRes += 32;
		}
		m_pdAaMod[cRes] = fValue;
		m_pdAaPrompt[cRes] = fPrompt;
		fSum += fPrompt;
		tStart = _s.find(',',tAt);
		if(tStart == _s.npos)
			break;
		tStart++;
		strValue = _s.substr(tStart,_s.size()-tStart);
		fValue = atof(strValue.c_str());
	}
	if(fPrompt == 0.0)	{
		m_bPrompt = false;
	}
	else	{
		m_bPrompt = true;
	}
	return true;
}

bool msequtilities::modify_annotation(string &_s)
{
	size_t tStart = 'a';
	const size_t tEnd = 'z'+1;
	while(tStart < tEnd)	{
		m_pdAaMod[tStart] = 0.0;
		m_pdAaPrompt[tStart] = 0.0;
		tStart++;
	}
	m_pdAaMod['['] = 0.0;
	m_pdAaMod[']'] = 0.0;

	m_bPotential = false;
	tStart = 0;
	size_t tAt = 0;
	size_t tColon = 0;
	double fValue = 0.0;
	double fPrompt = 0.0;
	string strValue = _s;
	if(!_s.empty())	{
		strValue += ",";
	}
	strValue += m_strDefaultMaybe;
	if(strValue.empty())	{
		return false;
	}
	string strV = strValue;
	fValue = atof(strV.c_str());
	char cRes = '\0';
	while(fValue != 0.0)	{
		m_bPotential = true;
		tAt = strValue.find('@',tStart);
		if(tAt == strValue.npos)
			break;
		tColon = strValue.find(':',tStart);
		fPrompt = 0.0;
		if(tColon != strValue.npos && tColon < tAt)	{
			fPrompt = atof(strValue.substr(tColon+1,tAt - tColon).c_str());
		}
		tAt++;
		cRes = strValue[tAt];
		if(cRes >= 'A' && cRes <= 'Z')	{
			cRes += 32;
		}
		m_pdAaMod[cRes] = fValue;
		m_pdAaPrompt[cRes] = fPrompt;
		tStart = strValue.find(',',tAt);
		if(tStart == strValue.npos)
			break;
		tStart++;
		strV = strValue.substr(tStart,_s.size()-tStart);
		fValue = atof(strV.c_str());
	}
	m_bPhosphoThreonine = false;
	if(fabs(m_pdAaMod['t']-79.966331) < 0.1)	{
		m_bPhosphoThreonine = true;
	}
	m_bPhosphoSerine = false;
	if(fabs(m_pdAaMod['s']-79.966331) < 0.1)	{
		m_bPhosphoSerine = true;
	}
	m_bPhosphoTyrosine = false;
	if(fabs(m_pdAaMod['y']-79.966331) < 0.1)	{
		m_bPhosphoTyrosine = true;
	}
	return true;
}


bool msequtilities::motif_set(const msequence &_s)
{
	m_mapMods.clear();
	m_mapMods = _s.m_mapMods;
	if(!m_mapMods.empty())	{
		m_bSequenceMods = true;
	}
	else	{
		m_bSequenceMods = false;
	}
	if(m_vMotifs.empty())	{
		return false;
	}
	m_mapMotifs.clear();
	char *pString = new char[_s.m_strSeq.size()+1];
	strcpy(pString,_s.m_strSeq.c_str());
	char *pValue = pString;
	size_t a = 0;
	size_t b = 0;
	size_t tPos;
	const size_t tMotifs = m_vMotifs.size();
	while(*pValue != '\0')	{
		a = 0;
		while(a < tMotifs)	{
			if(m_vMotifs[a].check(pValue,tPos))	{
				m_mapMotifs[pValue-pString+tPos] = a;
			}
			a++;
		}
		pValue++;
	}
	delete pString;
	return true;
}

bool msequtilities::modify_motif(const string &_s)
{
	m_vMotifs.clear();
	m_bPotentialMotif = false;
	if(_s.empty())	{
		return false;
	}
	size_t tStart = 0;
	size_t tAt = 0;
	float fValue;
	string strValue = _s.substr(tStart,_s.size()-tStart);
	fValue = (float)atof(strValue.c_str());
	char *pOut = new char[1024];
	mmotif motValue;
	while(fValue != 0.0)	{
		tAt = _s.find('@',tStart);
		if(tAt == _s.npos)
			break;
		tAt = tStart;
		tStart = _s.find(',',tAt);
		if(tStart == _s.npos)	{
			strValue = _s.substr(tAt,tStart-tAt);
			strcpy(pOut,strValue.c_str());
			motValue.initialize();
			if(motValue.set(pOut))	{
				m_vMotifs.push_back(motValue);
			}
			break;
		}
		else	{
			strValue = _s.substr(tAt,_s.size()-tStart);
			strcpy(pOut,strValue.c_str());
			motValue.initialize();
			if(motValue.set(pOut))	{
				m_vMotifs.push_back(motValue);
			}
		}
		tStart++;
		strValue = _s.substr(tStart,_s.size()-tStart);
		fValue = (float)atof(strValue.c_str());
	}
	if(!m_vMotifs.empty())	{
		m_bPotential = true;
		m_bPotentialMotif = true;
	}
	return true;
}

/*
 * sets the mass of the c-terminal residue modification for a protein
 */
bool msequtilities::modify_c(const float _f)
{
	m_fCT = _f;
	return true;
}
/*
 * sets the mass of the n-terminal residue modification for a protein
 */
bool msequtilities::modify_n(const float _f)
{
	m_fNT = _f;
	return true;
}
/*
 * initializes the mass of amino acid residues (less water)
 *     chemical formulas found at http://haven.isb-sib.ch/tools/isotopident/htdocs/aa-list.html
 */
bool msequtilities::set_aa()
{
	float *pValue = m_pfAaMass;
	if(pValue == NULL)
		return false;
	double *pdValue = m_pdAaMass;
	if(pdValue == NULL)
		return false;

	if (m_calc.getMassType() == masscalc::monoisotopic)	{
		pdValue['a'] = pdValue['A'] = m_calc.calcMass("C3H5ON");
		pValue['a'] = pValue['A'] = (float)pdValue['A'];

		pdValue['b'] = pdValue['B'] = m_calc.calcMass("C4H6O2N2");	// Same as N
		pValue['b'] = pValue['B'] = (float)pdValue['B'];

		pdValue['c'] = pdValue['C'] = m_calc.calcMass("C3H5ONS");
		pValue['c'] = pValue['C'] = (float)pdValue['C'];
		
		pdValue['d'] = pdValue['D'] = m_calc.calcMass("C4H5O3N");
		pValue['d'] = pValue['D'] = (float)pdValue['D'];
		
		pdValue['e'] = pdValue['E'] = m_calc.calcMass("C5H7O3N");
		pValue['e'] = pValue['E'] = (float)pdValue['E'];
		
		pdValue['f'] = pdValue['F'] = m_calc.calcMass("C9H9ON");
		pValue['f'] = pValue['F'] = (float)pdValue['F'];
		
		pdValue['g'] = pdValue['G'] = m_calc.calcMass("C2H3ON");
		pValue['g'] = pValue['G'] = (float)pdValue['G'];
		
		pdValue['h'] = pdValue['H'] = m_calc.calcMass("C6H7ON3");
		pValue['h'] = pValue['H'] = (float)pdValue['H'];
		
		pdValue['i'] = pdValue['I'] = m_calc.calcMass("C6H11ON");
		pValue['i'] = pValue['I'] = (float)pdValue['I'];
		
		pdValue['j'] = pdValue['J'] = 0.0;
		pValue['j'] = pValue['J'] = (float)pdValue['J'];
		
		pdValue['k'] = pdValue['K'] = m_calc.calcMass("C6H12ON2");
		pValue['k'] = pValue['K'] = (float)pdValue['K'];
		
		pdValue['l'] = pdValue['L'] = m_calc.calcMass("C6H11ON");
		pValue['l'] = pValue['L'] = (float)pdValue['L'];
		
		pdValue['m'] = pdValue['M'] = m_calc.calcMass("C5H9ONS");
		pValue['m'] = pValue['M'] = (float)pdValue['M'];
		
		pdValue['n'] = pdValue['N'] = m_calc.calcMass("C4H6O2N2");
		pValue['n'] = pValue['N'] = (float)pdValue['N'];
		
		pdValue['o'] = pdValue['O'] = m_calc.calcMass("C4H6O2N2");	// Same as N
		pValue['o'] = pValue['O'] = (float)pdValue['O'];
		
		pdValue['p'] = pdValue['P'] = m_calc.calcMass("C5H7ON");
		pValue['p'] = pValue['P'] = (float)pdValue['P'];
		
		pdValue['q'] = pdValue['Q'] = m_calc.calcMass("C5H8O2N2");
		pValue['q'] = pValue['Q'] = (float)pdValue['Q'];
		
		pdValue['r'] = pdValue['R'] = m_calc.calcMass("C6H12ON4");
		pValue['r'] = pValue['R'] = (float)pdValue['R'];
		
		pdValue['s'] = pdValue['S'] = m_calc.calcMass("C3H5O2N");
		pValue['s'] = pValue['S'] = (float)pdValue['S'];
		
		pdValue['t'] = pdValue['T'] = m_calc.calcMass("C4H7O2N");
		pValue['t'] = pValue['T'] = (float)pdValue['T'];
		
		pdValue['u'] = pdValue['U'] = 150.953640;	// Why?
		pValue['u'] = pValue['U'] = (float)pdValue['U'];
		
		pdValue['v'] = pdValue['V'] = m_calc.calcMass("C5H9ON");
		pValue['v'] = pValue['V'] = (float)pdValue['V'];
		
		pdValue['w'] = pdValue['W'] = m_calc.calcMass("C11H10ON2");
		pValue['w'] = pValue['W'] = (float)pdValue['W'];
		
		pdValue['x'] = pdValue['X'] = 111.060000;	// Why?
		pValue['x'] = pValue['X'] = (float)pdValue['X'];
		
		pdValue['y'] = pdValue['Y'] = m_calc.calcMass("C9H9O2N");
		pValue['y'] = pValue['Y'] = (float)pdValue['Y'];
		
		pdValue['z'] = pdValue['Z'] = m_calc.calcMass("C5H8O2N2");	// Same as Q
		pValue['z'] = pValue['Z'] = (float)pdValue['Z'];
	}
	else {
		/*
		 * unfortunately, the average masses for amino acids do not
		 * seem to be straight sums of the average masses for the atoms
		 * they contain.
		 *
		 * instead of using the mass calculator, these numbers are taken
		 * as constants from the web page referenced above.
		 */
		pdValue['a'] = pdValue['A'] = 71.0788;
		pValue['a'] = pValue['A'] = (float)pdValue['A'];

		pdValue['b'] = pdValue['B'] = 114.1038;	// Same as N
		pValue['b'] = pValue['B'] = (float)pdValue['B'];

		pdValue['c'] = pdValue['C'] = 103.1388;
		pValue['c'] = pValue['C'] = (float)pdValue['C'];
		
		pdValue['d'] = pdValue['D'] = 115.0886;
		pValue['d'] = pValue['D'] = (float)pdValue['D'];
		
		pdValue['e'] = pdValue['E'] = 129.1155;
		pValue['e'] = pValue['E'] = (float)pdValue['E'];
		
		pdValue['f'] = pdValue['F'] = 147.1766;
		pValue['f'] = pValue['F'] = (float)pdValue['F'];
		
		pdValue['g'] = pdValue['G'] = 57.0519;
		pValue['g'] = pValue['G'] = (float)pdValue['G'];
		
		pdValue['h'] = pdValue['H'] = 137.1411;
		pValue['h'] = pValue['H'] = (float)pdValue['H'];
		
		pdValue['i'] = pdValue['I'] = 113.1594;
		pValue['i'] = pValue['I'] = (float)pdValue['I'];
		
		pdValue['j'] = pdValue['J'] = 0.0;
		pValue['j'] = pValue['J'] = (float)pdValue['J'];
		
		pdValue['k'] = pdValue['K'] = 128.1741;
		pValue['k'] = pValue['K'] = (float)pdValue['K'];
		
		pdValue['l'] = pdValue['L'] = 113.1594;
		pValue['l'] = pValue['L'] = (float)pdValue['L'];
		
		pdValue['m'] = pdValue['M'] = 131.1926;
		pValue['m'] = pValue['M'] = (float)pdValue['M'];
		
		pdValue['n'] = pdValue['N'] = 114.1038;
		pValue['n'] = pValue['N'] = (float)pdValue['N'];
		
		pdValue['o'] = pdValue['O'] = 114.1038;	// Same as N
		pValue['o'] = pValue['O'] = (float)pdValue['O'];
		
		pdValue['p'] = pdValue['P'] = 97.1167;
		pValue['p'] = pValue['P'] = (float)pdValue['P'];
		
		pdValue['q'] = pdValue['Q'] = 128.1307;
		pValue['q'] = pValue['Q'] = (float)pdValue['Q'];
		
		pdValue['r'] = pdValue['R'] = 156.1875;
		pValue['r'] = pValue['R'] = (float)pdValue['R'];
		
		pdValue['s'] = pdValue['S'] = 87.0782;
		pValue['s'] = pValue['S'] = (float)pdValue['S'];
		
		pdValue['t'] = pdValue['T'] = 101.1051;
		pValue['t'] = pValue['T'] = (float)pdValue['T'];
		
		pdValue['u'] = pdValue['U'] = 0.0;	// Why?
		pValue['u'] = pValue['U'] = (float)pdValue['U'];
		
		pdValue['v'] = pdValue['V'] = 99.1326;
		pValue['v'] = pValue['V'] = (float)pdValue['V'];
		
		pdValue['w'] = pdValue['W'] = 186.2132;
		pValue['w'] = pValue['W'] = (float)pdValue['W'];
		
		pdValue['x'] = pdValue['X'] = 113.1594;	// Why?
		pValue['x'] = pValue['X'] = (float)pdValue['X'];
		
		pdValue['y'] = pdValue['Y'] = 163.1760;
		pValue['y'] = pValue['Y'] = (float)pdValue['Y'];
		
		pdValue['z'] = pdValue['Z'] = 128.1307;	// Same as Q
		pValue['z'] = pValue['Z'] = (float)pdValue['Z'];
	}

	return true;
}
/*
 * sets the score for a particular a ion
 */
bool msequtilities::set_a(const char _c,const float _f)
{
	m_pfAScore[_c] = _f;
	m_pfA18Score[_c] = _f;
	m_pfA17Score[_c] = _f;
	return true;
}
/*
 * sets the score for a particular b ion
 */
bool msequtilities::set_b(const char _c,const float _f)
{
	m_pfBScore[_c] = _f;
	m_pfB18Score[_c] = _f;
	m_pfB17Score[_c] = _f;
	return true;
}
/*
 * sets the score for a particular x ion
 */
bool msequtilities::set_x(const char _c,const float _f)
{
	m_pfXScore[_c] = _f;
	m_pfX18Score[_c] = _f;
	m_pfX17Score[_c] = _f;
	return true;
}
/*
 * sets the score for a particular y ion
 */
bool msequtilities::set_y(const char _c,const float _f)
{
	m_pfYScore[_c] = _f;
	m_pfY18Score[_c] = _f;
	m_pfY17Score[_c] = _f;
	return true;
}

bool msequtilities::set_aa_file(string &_p)
{
	m_bIsModified = false;
	ifstream ifIn;
	ifIn.open(_p.c_str());
	if(ifIn.fail())	{
		return false;
	}
	vector<string> vstrValues;
	string strValue;
	char *pLine = new char[1024];
	while(ifIn.good() && !ifIn.eof())	{
		ifIn.getline(pLine,1023);
		strValue = pLine;
		vstrValues.push_back(strValue);
	}
	ifIn.close();
	size_t a = 0;
	size_t tStart;
	size_t tEnd;
	double dMass = 0.0;
	char cAa = '\0';
	const size_t tLength = vstrValues.size();
	string strType;
	string strLine;
	while(a < tLength)	{
		cAa = '\0';
		strType = " ";
		dMass = -1.0;
		strLine = vstrValues[a];
		if(strLine.find("<aa ") != strLine.npos)	{
			tStart = strLine.find("mass=\"");
			if(tStart != strLine.npos)	{
				tStart = strLine.find("\"",tStart);
				tStart++;
				dMass = atof(strLine.substr(tStart,strLine.size() - tStart).c_str());
			}
			tStart = strLine.find("type=\"");
			if(tStart != strLine.npos)	{
				tStart = strLine.find("\"",tStart);
				tStart++;
				cAa = strLine[tStart];
			}
			if(cAa != '\0' && dMass >= 0.0)	{
				if(cAa > 'Z')	{
					cAa -= 32;
				}
				if(isalpha(cAa))	{
					m_bIsModified = true;
					m_pdAaMass[cAa+32] = m_pdAaMass[cAa] = dMass;
					m_pfAaMass[cAa+32] = m_pfAaMass[cAa] = (float)dMass;
				}
			}
		}
		if(strLine.find("<molecule ") != strLine.npos)	{
			tStart = strLine.find("mass=\"");
			if(tStart != strLine.npos)	{
				tStart = strLine.find("\"",tStart);
				tStart++;
				dMass = atof(strLine.substr(tStart,strLine.size() - tStart).c_str());
			}
			tStart = strLine.find("type=\"");
			if(tStart != strLine.npos)	{
				tStart = strLine.find("\"",tStart);
				tStart++;
				tEnd = strLine.find("\"",tStart);
				if(tEnd != strLine.npos)	{
					strType = strLine.substr(tStart,tEnd - tStart);
				}
			}
			if(dMass >= 0.0 && strType != " ")	{
				if(strType == "NH3")	{
					m_bIsModified = true;
					m_dAmmonia = dMass;
				}
				if(strType == "H2O")	{
					m_bIsModified = true;
					m_dWater = dMass;
				}
			}
		}
		a++;
	}
	m_dC = m_dAmmonia;
	m_dY = m_dWater;
	return true;
}

bool msequtilities::is_modified()
{
	return m_bIsModified;
}

bool msequtilities::set_modified(const bool _b)
{
	m_bIsModified = _b;
	return m_bIsModified;
}
