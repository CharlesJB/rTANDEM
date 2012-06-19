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

#ifndef MCLEAVE_H
#define MCLEAVE_H

// File version: 2006-03-10
/*
 * mcleave is a specialty class meant to store information about protein cleavage specificity
 * and rapidly test a peptide sequence to see if it is cleavable.
 * NOTE: mcleave.h does not have a corresponding .cpp file
 */
class mcleave_single
{
public:
	mcleave_single(void)	{
		strcpy(m_pNCleave,"KR");
		strcpy(m_pCCleave,"P");
		m_bN = true;
		m_bC = false;
		m_bNX = false;
		m_bCX = false;
		m_lType = 0;
	}
	virtual ~mcleave_single(void) { }

	char m_pNCleave[32]; // residues that are valid cleavage sites (or invalid if m_bN = false)
	char m_pCCleave[32]; // residues that are valid cleavage sites (or invalid if m_bC = false)
	bool m_bN; // if true, all residues in m_pNCleave can be N-temrinal to a cleavable bond
	bool m_bC; // if true, all residues in m_pNCleave can be C-temrinal to a cleavable bond
	bool m_bCX;
	bool m_bNX;
	unsigned long m_lType;
	string m_strCleave;
	mcleave_single& operator=(const mcleave_single &rhs)	{
		strcpy(m_pNCleave,rhs.m_pNCleave);
		strcpy(m_pCCleave,rhs.m_pCCleave);
		m_bN = rhs.m_bN;
		m_bC = rhs.m_bC;
		m_bNX = rhs.m_bNX;
		m_bCX = rhs.m_bCX;
		m_lType = rhs.m_lType;
		return *this;
	}
	/*
 * load takes a string containing cleavage information and parses it
 * these strings take the format:
 * A|B, where A or B are expressions of the form [xyz] or {xyz}. if square brackets, only the single letter
 * abbreviations in the brackets are valid at that position for a cleavage, and if french brackets, those
 * single letter abbreviations are the only non-valid residues at that position. For example, trypsin 
 * can be represented by [KR]|{P}, i.e. cleave C-terminal to a K or R residue, except when followed by
 * a P. The expression [X]|[X] means cleave at all residues.
 * NOTE: use of upper case characters is required for amino acid abbreviations
 */
	bool load(string &_s)	{
		m_strCleave = _s;
		if(_s == "[X]|[X]")	{
			m_lType = 0x01;
			return true;
		}
		else if(_s == "[KR]|{P}" || _s == "[RK]|{P}")	{
			m_lType = 0x02;
			return true;
		}
		m_lType = 0x04;
		size_t a = 0;
		size_t b = 0;
		if(_s[a] == '[')	{
			m_bN = true;
			a++;
			b = 0;
			while(a < _s.size() && _s[a] != ']')	{
				m_pNCleave[b] = _s[a];
				a++;
				b++;
			}
			m_pNCleave[b] = '\0';
			a = _s.find('|');
			if(a == _s.npos)
				return false;
			a++;
			if(_s[a] == '{')	{
				m_bC = false;
				a++;
				b = 0;
				while(a < _s.size() && _s[a] != '}')	{
					m_pCCleave[b] = _s[a];
					a++;
					b++;
				}
				m_pCCleave[b] = '\0';
			}
			else if(_s[a] == '[')	{
				m_bC = true;
				a++;
				b = 0;
				while(a < _s.size() && _s[a] != ']')	{
					m_pCCleave[b] = _s[a];
					a++;
					b++;
				}
				m_pCCleave[b] = '\0';
			}
		}
		else if(_s[a] == '{')	{
			m_bN = false;
			a++;
			b = 0;
			while(a < _s.size() && _s[a] != '}')	{
				m_pNCleave[b] = _s[a];
				a++;
				b++;
			}
			m_pNCleave[b] = '\0';
			a = _s.find('|');
			if(a == _s.npos)
				return false;
			a++;
			if(_s[a] == '{')	{
				m_bC = false;
				a++;
				b = 0;
				while(a < _s.size() && _s[a] != '}')	{
					m_pCCleave[b] = _s[a];
					a++;
					b++;
				}
				m_pCCleave[b] = '\0';
			}
			else if(_s[a] == '[')	{
				m_bC = true;
				a++;
				b = 0;
				while(a < _s.size() && _s[a] != ']')	{
					m_pCCleave[b] = _s[a];
					a++;
					b++;
				}
				m_pCCleave[b] = '\0';
			}
		}
		if(m_pNCleave[0] == 'X')	{
			m_bNX = true;
		}
		if(m_pCCleave[0] == 'X')	{
			m_bCX = true;
		}
		return true;
	}
/*
 * test takes the abbreviations for the residue N-terminal to a potentially cleaved bond
 * and the residue C-terminal to the bond and checks to see if the bond can be cleaved
 * according to the rules stored in the load method. load must always be called at least once
 * prior to using test, or the results may be unpredictable
 */
	bool test(const char _n,const char _c)	{
		if(m_lType & 0x01)
			return true;
		if(m_lType & 0x02)	{
			if(_n == 'K' || _n == 'R')	{
				if(_c != 'P')	{
					return true;
				}
			}
			return false;
		}
		bool bReturn = false;
		bool bN = false;
		bool bC = false;
		if(m_bNX)	{
			bN = true;
		}
		else	{
			if(strchr(m_pNCleave,_n))
				bN = true;
		}
		if(m_bN)	{
			bReturn = bN;
		}
		else	{
			bReturn = !bN;
		}
		if(!bReturn)
			return false;
		if(m_bCX)	{
			bC = true;
		}
		else	{
			if(strchr(m_pCCleave,_c))
				bC = true;
		}
		if(m_bC && bC)	{
			return true;
		}
		else if(!m_bC && !bC)	{
			return true;
		}
		return false;
	}
};

class mcleave
{
public:
	mcleave(void)	{
		m_lType = 0;
		m_vCleaves.clear();
		m_itStart = m_vCleaves.begin();
		m_itEnd = m_vCleaves.end();
	}
	virtual ~mcleave(void) { }

	vector<mcleave_single> m_vCleaves;
	vector<mcleave_single>::iterator m_itStart;
	vector<mcleave_single>::iterator m_itEnd;
	vector<mcleave_single>::iterator m_itValue;
	string m_strCleave;
	unsigned long m_lType;
	/*
 * load takes a string containing cleavage information and parses it
 * these strings take the format:
 * A|B, where A or B are expressions of the form [xyz] or {xyz}. if square brackets, only the single letter
 * abbreviations in the brackets are valid at that position for a cleavage, and if french brackets, those
 * single letter abbreviations are the only non-valid residues at that position. For example, trypsin 
 * can be represented by [KR]|{P}, i.e. cleave C-terminal to a K or R residue, except when followed by
 * a P. The expression [X]|[X] means cleave at all residues.
 * NOTE: use of upper case characters is required for amino acid abbreviations
 */
	bool load(string &_s)	{
		m_strCleave = _s;
		m_lType = 0x04;
		size_t a = 0;
		size_t tEnd = _s.size();
		string strTemp;
		mcleave_single clvTemp;
		m_vCleaves.clear();
		while(a < tEnd)	{
			if(_s[a] == ',')	{
				if(clvTemp.load(strTemp))	{
					m_vCleaves.push_back(clvTemp);
				}
				strTemp.erase(0,strTemp.size());
			}
			else if(strchr("ABCDEFGHIJKLMNOPQRSTUVWXYZ[]{}|",_s[a]))	{
				strTemp += _s[a];
			}
			else if(_s[a] >= 'a' && _s[a] <= 'z')	{
				strTemp += (char)(_s[a]-32);
			}
			a++;
		}
		if(!strTemp.empty())	{
			if(clvTemp.load(strTemp))	{
				m_vCleaves.push_back(clvTemp);
			}
		}
		m_itStart = m_vCleaves.begin();
		m_itEnd = m_vCleaves.end();
		if(m_vCleaves.size() == 1)	{
			m_lType = m_vCleaves[0].m_lType;
		}
		return !m_vCleaves.empty();
	}
/*
 * test takes the abbreviations for the residue N-terminal to a potentially cleaved bond
 * and the residue C-terminal to the bond and checks to see if the bond can be cleaved
 * according to the rules stored in the load method. load must always be called at least once
 * prior to using test, or the results may be unpredictable
 */
	bool test(const char _n,const char _c)	{
		if(m_lType & 0x01)
			return true;
		if(m_lType & 0x02)	{
			if(_n == 'K' || _n == 'R')	{
				if(_c != 'P')	{
					return true;
				}
			}
			return false;
		}
		m_itValue = m_itStart;
		while(m_itValue != m_itEnd)	{
			if(m_itValue->test(_n,_c))	{
				return true;
			}
			m_itValue++;
		}
		return false;
	}
};
#endif //ifdef MCLEAVE_H
