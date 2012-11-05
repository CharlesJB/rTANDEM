/*
Copyright (C) 2004 Ronald C Beavis, all rights reserved
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

#ifndef MSCOREPAM_H
#define MSCOREPAM_H

#include <set>

// File version: 2004-02-01
/*
 * mscorepam is a specialty class used to store information necessary for the state machine that
 * mscore uses to generate all potential point mutations
 */

class mscorepam
{
public:

	mscorepam(void)
	{
		m_pAa = new char[26];
		strcpy(m_pAa,"ACDEFGHKLMNPQRSTVWY");
		m_tAaTotal = strlen(m_pAa);
		m_pSeqTrue = new char[256];
		m_pSeqTrue[0] = '\0';
		m_tLength = 0;
		m_tPos = 0;
		m_tAa = 0;
		m_tCount = 0;
		m_cCurrent = '\0';
		m_tSeqTrue = 255;
		m_fSeqTrue = 0.0;
	}

	virtual ~mscorepam(void)
	{
		if(m_pAa != NULL)	{
			delete m_pAa;
		}
		if(m_pSeqTrue != NULL)	{
			delete m_pSeqTrue;
		}
	}

	bool initialize(const char *_p,const size_t _s,const float _f)
	{
		if(m_tSeqTrue < _s)	{
			if(m_pSeqTrue != NULL)	{
				delete m_pSeqTrue;
			}
			m_tSeqTrue = _s;
			m_pSeqTrue = new char[m_tSeqTrue + 1];
		}
		strcpy(m_pSeqTrue,_p);
		m_fSeqTrue = _f;
		m_tLength = strlen(_p);
		m_tPos = 0;
		m_tAa = 0;
		m_tCount = 0;
		m_cCurrent = '\0';
		return true;
	}

	size_t m_tLength;
	size_t m_tPos;
	size_t m_tAa;
	size_t m_tCount;
	char m_cCurrent;
	size_t m_tAaTotal;
	char *m_pAa;
	char *m_pSeqTrue;
	size_t m_tSeqTrue;
	float m_fSeqTrue;
};

class mscoresap
{
public:
	mscoresap(void)	{
		m_bOk = false;
		m_pSeqTrue = new char[256];
		m_pSeqTrue[0] = '\0';
		m_cCurrent = '\0';
		m_tSeqTrue = 255;
		m_fSeqTrue = 0.0;
		m_tCount = 0;
		m_tPos = 0;
		m_iPos = 0;
		m_strId.clear();
		m_bEnd = true;
	}
	virtual ~mscoresap(void) { 
		delete m_pSeqTrue;
	}
	int m_iStart;
	int m_iEnd;
	int m_iPos;
	string m_strId;
	map<string,multimap<int,prSap> > m_mapSap;
	multimap<int,prSap>::iterator m_itSap;
	multimap<int,prSap>::iterator m_itSapEnd;
	map<string,multimap<int,prSap> >::iterator m_itAcc;
	char *m_pSeqTrue;
	set<char> m_setAllowed;
	size_t m_tSeqTrue;
	size_t m_tLength;
	float m_fSeqTrue;
	size_t m_tCount;
	size_t m_tPos;
	char m_cCurrent;
	bool m_bOk;
	bool m_bEnd;
	bool reset_value(string &_s,const bool _b)	{
		m_bOk = false;
		if(!_b)	{
			return m_bOk;
		}
		if(m_mapSap.empty())	{
			return m_bOk;
		}
		m_tCount = 0;
		m_itAcc = m_mapSap.find(_s);
		m_bOk = true;
		if(m_itAcc == m_mapSap.end())	{
			m_bOk = false;
		}
		m_bEnd = true;
		m_iStart = 0;
		m_iEnd = 0;
		cout.flush();
		return m_bOk;
	}
	bool allowed(string &_s)
	{
		m_setAllowed.clear();
		if(_s.empty())	{
			return false;
		}
		size_t a = 0;
		while(a < _s.length())	{
			m_setAllowed.insert(_s[a]);
			a++;
		}
		return true;
	}
	bool initialize(const char *_p,const size_t _s,const float _f,const int _i)
	{
		if(!m_bOk)	{
			return false;
		}
		if(m_tSeqTrue < _s)	{
			if(m_pSeqTrue != NULL)	{
				delete m_pSeqTrue;
			}
			m_tSeqTrue = _s;
			m_pSeqTrue = new char[m_tSeqTrue + 1];
		}
		strcpy(m_pSeqTrue,_p);
		m_fSeqTrue = _f;
		m_tLength = strlen(_p);
		m_cCurrent = '\0';
		m_iStart = _i + 1;
		m_iEnd = m_iStart + (int)m_tLength - 1;
		m_tCount = 0;
		m_tPos = 0;
		m_strId.clear();
		m_bEnd = false;
		m_iPos = m_iStart;
		m_itSap = m_itAcc->second.lower_bound(m_iPos);
		m_itSapEnd = m_itAcc->second.end();
		if(m_itSap == m_itSapEnd)	{
			m_bEnd = true;
			return false;
		}
#if 0
		//WCH: wrong pair attribute checking
		if(m_itSap->second.first > m_iEnd)	{
#else
		if(m_itSap->first > m_iEnd)	{
#endif
			m_bEnd = true;
			return false;
		}
		return true;
	}
	bool initialize(const char *_p,const size_t _s,const float _f)
	{
		if(!m_bOk)	{
			return false;
		}
		if(m_tSeqTrue < _s)	{
			if(m_pSeqTrue != NULL)	{
				delete m_pSeqTrue;
			}
			m_tSeqTrue = _s;
			m_pSeqTrue = new char[m_tSeqTrue + 1];
		}
		strcpy(m_pSeqTrue,_p);
		m_fSeqTrue = _f;
		m_tLength = strlen(_p);
		m_cCurrent = '\0';
		m_iEnd = m_iStart + (int)m_tLength - 1;
		m_tCount = 0;
		m_tPos = 0;
		m_strId.clear();
		m_bEnd = false;
		m_iPos = m_iStart;
		m_itSap = m_itAcc->second.lower_bound(m_iPos);
		m_itSapEnd = m_itAcc->second.end();
		if(m_itSap == m_itSapEnd)	{
			m_bEnd = true;
			return false;
		}
#if 0
		//WCH: wrong pair attribute checking
		if(m_itSap->second.first > m_iEnd)	{
#else
		if(m_itSap->first > m_iEnd)	{
#endif
			m_bEnd = true;
			return false;
		}
		return true;
	}
	bool next(void)
	{
		if(m_bEnd || !m_bOk)	{
			return false;
		}
		m_iPos = m_itSap->first;
		if(m_iPos > m_iEnd)	{
			m_bEnd = true;
			return false;
		}
		m_cCurrent = m_itSap->second.first;
		if(!m_setAllowed.empty())	{
			while(m_setAllowed.find(m_cCurrent) == m_setAllowed.end())	{
				m_tCount++;
				m_itSap++;
				if(m_itSap == m_itSapEnd)	{
					m_bEnd = true;
					return false;
				}
				m_iPos = m_itSap->first;
				if(m_iPos > m_iEnd)	{
					m_bEnd = true;
					return false;
				}
				m_cCurrent = m_itSap->second.first;
			}
		}
		m_tPos = m_iPos - m_iStart;
		m_strId = m_itSap->second.second;
		m_tCount++;
		m_itSap++;
		if(m_itSap == m_itSapEnd || m_itSap->first > m_iEnd)	{
			m_bEnd = true;
		}
		return true;
	}
};
#endif
