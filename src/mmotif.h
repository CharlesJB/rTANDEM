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

#ifndef MMOTIF_H
#define MMOTIF_H

#include <string.h> // rTANDEM

// File version: 2004-06-01
// File version: 2005-01-01
/*
 * mmotif is a specialty class meant to store information about protein sequence motifs.
 * NOTE: mmotif.h does not have a corresponding .cpp file
 */

class mmotifres
{
public:

	mmotifres(void)
	{
		m_bRes = false;
		m_bPos = true;
		m_bIsX = false;
		strcpy(m_pMotif,"");
	}

	virtual ~mmotifres(void)
	{
	}
	char m_pMotif[32];
	bool m_bRes;
	bool m_bPos;
	bool m_bIsX;

	bool set(const char *_p)	{
		if(_p == NULL)	{
			return false;
		}
		if(strlen(_p) == 0)	{
			return false;
		}
		strcpy(m_pMotif,"");
		m_bIsX = false;
		m_bRes = true;
		m_bPos = true;
		if(strchr(_p,'!') == NULL)	{
			m_bRes = false;
		}
		if(strchr(_p,'X') != NULL)	{
			m_bIsX = true;
			m_bPos = true;
			strcpy(m_pMotif,"X");
			return true;
		}
		char *pValue = new char[strlen(_p)+1];
		strcpy(pValue,_p);
		char *pStart = pValue;
		size_t a = 0;
		if(strchr(pValue,'[') != NULL)	{
			pStart = strchr(pValue,'[');
			pStart++;
			while(*pStart != ']' && *pStart != '\0')	{
				if(isalpha(*pStart))	{
					m_pMotif[a] = *pStart;
				}
				pStart++;
				a++;
			}
			m_pMotif[a] = '0';
			m_bPos = true;
			return true;
		}
		else if(strchr(pValue,'{') != NULL)	{
			pStart = strchr(pValue,'{');
			pStart++;
			while(*pStart != ']' && *pStart != '\0')	{
				if(isalpha(*pStart))	{
					m_pMotif[a] = *pStart;
				}
				pStart++;
				a++;
			}
			m_pMotif[a] = '0';
			m_bPos = false;
			return true;
		}
		else	{
			pStart = pValue;
			while(*pStart != '0')	{
				if(isalpha(*pStart))	{
					m_pMotif[a] = *pStart;
					a++;
					break;
				}
				pStart++;
			}
			m_pMotif[a] = '\0';
			m_bPos = true;
		}
		return false;
	}

	bool check(const char _c)	{
		if(m_bIsX)	{
			return true;
		}
		if(m_bPos)	{
			if(strchr(m_pMotif,_c) != NULL)	{
				return true;
			}
			return false;
		}
		if(strchr(m_pMotif,_c) == NULL)	{
			return true;
		}
		return false;
	}
	mmotifres& operator=(const mmotifres &rhs)	{
		strcpy(m_pMotif,rhs.m_pMotif);
		m_bRes = rhs.m_bRes;
		m_bPos = rhs.m_bPos;
		m_bIsX = rhs.m_bIsX;
		return *this;
	}

};

class mmotif
{
public:

	mmotif(void)
	{
		m_fMass = 0.0;
		m_fMassPrompt = 0;
		m_tPos = 0;
	}

	virtual ~mmotif(void)
	{
	}
	vector<mmotifres> m_vresMotifs;
	float m_fMass;
	float m_fMassPrompt;
	size_t m_tPos;
	bool set(const char *_p)	{
		if(_p == NULL)	{
			return false;
		}
		if(strlen(_p) == 0)	{
			return false;
		}
		const char *pStart = (char *)strchr(_p,'@');
		const char *pColon = (char *)strchr(_p,':');
		if(pStart == NULL)	{
			return false;
		}
		m_fMassPrompt = 0.0;
		if(pColon != NULL && pColon < pStart)	{
			m_fMassPrompt = (float)atof(pColon+1);
		}
		m_fMass = (float)(atof(_p));
		pStart++;
		size_t a = 0;
		size_t b = 0;
		size_t tEnd = strlen(pStart);
		char *pOut = new char[256];
		mmotifres motRes;
		m_vresMotifs.clear();
		while(a < tEnd)	{
			if(pStart[a] == '[')	{
				a++;
				b = 0;
				pOut[b] = '[';
				b++;
				while(pStart[a] != '\0' && pStart[a] != ']')	{
					pOut[b] = pStart[a];
					b++;
					a++;
				}
				pOut[b] = ']';
				b++;
				pOut[b] = '\0';
				motRes.set(pOut);
				if(motRes.m_bRes)	{
					m_tPos = m_vresMotifs.size();
				}
				m_vresMotifs.push_back(motRes);
			}
			else if(pStart[a] == '{')	{
				a++;
				b = 0;
				pOut[b] = '{';
				b++;
				while(pStart[a] != '\0' && pStart[a] != '}')	{
					pOut[b] = pStart[a];
					b++;
					a++;
				}
				pOut[b] = '}';
				b++;
				pOut[b] = '\0';
				motRes.set(pOut);
				if(motRes.m_bRes)	{
					m_tPos = m_vresMotifs.size();
				}
				m_vresMotifs.push_back(motRes);
			}
			else if(pStart[a] == '(')	{
				a++;
				unsigned long c = (unsigned long)atoi(pStart+a);
				if(c > 0)
					c--;
				while(pStart[a] != '\0' && pStart[a] != ')')	{
					a++;
				}
				if(!m_vresMotifs.empty())	{
					b = 0;
					while(b < c)	{
						m_vresMotifs.push_back(motRes);
						b++;
					}
				}
			}
			else if(isalpha(pStart[a]))	{
				b = 0;
				pOut[b] = pStart[a];
				b++;
				if(pStart[a+1] == '!')	{
					pOut[b] = pStart[a+1];
					b++;
				}
				pOut[b] = '\0';
				motRes.set(pOut);
				if(motRes.m_bRes)	{
					m_tPos = m_vresMotifs.size();
				}
				m_vresMotifs.push_back(motRes);
			}
			a++;
		}
		delete pOut;
		return true;
	}

	bool check(const char *_p,size_t &_o)	{
		size_t a = 0;
		const size_t tSize = m_vresMotifs.size();
		while(*(_p+a) != '\0' && a < tSize)	{
			if(!m_vresMotifs[a].check(*(_p+a)))	{
				return false;
			}
			a++;
		}
		if(a != tSize)	{
			return false;
		}
		_o = m_tPos;
		return true;
	}

	bool initialize()	{
		m_vresMotifs.clear();
		m_fMass = 0.0;
		m_tPos = 0;
		return true;
	}

	mmotif& operator=(const mmotif &rhs)	{
		m_fMass = rhs.m_fMass;
		m_tPos = rhs.m_tPos;
		m_vresMotifs.clear();
		size_t a = 0;
		while(a < rhs.m_vresMotifs.size())	{
			m_vresMotifs.push_back(rhs.m_vresMotifs[a]);
			a++;
		}
		return *this;
	}
};
#endif 
