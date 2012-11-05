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
// File version: 2004-03-01

/*
 * msequenceserver takes a sequence list file (in FASTA format) and uses it to load an
 * msequencecontainer object with a set of amino acid sequences and descriptions. the container
 * can be loaded repeatedly until the end-of-file is reached. the next list file path in a 
 * deque of file names is then extracted and the process continues until the last
 * sequence from the last file in the deque is read into the msequencecontainer.
 */

#include "stdafx.h"
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "xmltaxonomy.h"
#include <ctime>

msequenceServer::msequenceServer(void)
{
	m_pCol = new msequenceCollection;
	m_bStarted = false;
	m_bDone = false;
	m_bError = false;
	m_strStatus = "msequenceServer initialized<BR>\r\n";
	m_tColMax = 1000;
	m_tStartAt = 1;
	m_dTime = 0.0;
	m_lFileType = FASTA;
	// 2006.11.21 - increased the size from 32*4096 to 512*1024 because of very long lines in nr FASTA files
	m_lSize = 512*1024-1;
	m_pLine = new char[m_lSize+1];
}

msequenceServer::~msequenceServer(void)
{
	if(m_pCol != NULL)	{
		delete m_pCol;
	}
	delete m_pLine;
}
/*
 * clears the sequence collection object	
 */
bool msequenceServer::clear(void)
{
	if(m_pCol == NULL)
		return false;
	m_pCol->m_vASequences.clear();
	return true;
}
/*
 * returns true if the server has retrieved all possible sequences	
 */
bool msequenceServer::done(void)
{
	return m_bDone;
}
/*
 * returns true if the server is in an error conditions	
 */
bool msequenceServer::error(void)
{
	return m_bError;
}
/*
 * called to finish the sequence loading process	
 */
bool msequenceServer::finish(void)
{
	m_bDone = true;
	fclose(m_pInput);
	m_strStatus += "Server finished properly\n";
	return m_bDone;
}
/*
 * returns the time elapsed retrieving sequences in seconds	
 */
double msequenceServer::get_time(void)
{
	return m_dTime/(double)CLOCKS_PER_SEC;
}
bool msequenceServer::initialize(const size_t _t)
{
	m_pCol->initialize(_t);
	m_tColMax = _t;
	return true;
}
/*
 * retreives the list of FASTA sequence list files using an XmlTaxonomy object
 * and a taxon string. the list of file paths is stored in the m_dstrFasta deque and
 * the m_vstrFasta vector
 */
long msequenceServer::load_file(const string &_p,const string &_t)
{
	m_strTaxonPath = _p;
	m_strTaxon = _t;
	XmlTaxonomy xmlTax;
	string strType = "peptide";
	if(!xmlTax.load(m_strTaxonPath,m_strTaxon,strType))
		return 1;
	size_t a = 0;
	ifstream ifTest;
	m_vstrFasta.clear();
	long lFailed = 0;
	while(a < xmlTax.m_vstrPaths.size())	{
		ifTest.open(xmlTax.m_vstrPaths[a].c_str());
		if(!ifTest.fail())	{
			m_dstrFasta.push_back(xmlTax.m_vstrPaths[a]);
			m_vstrFasta.push_back(xmlTax.m_vstrPaths[a]);
			ifTest.close();
		}
		else	{
			lFailed++;
		}
		ifTest.clear();
		a++;
	}
	if(m_dstrFasta.empty())	{
		if(lFailed != 0)	{
			return 3;
		}
		return 2;
	}
	return 0;
}
/*
 * refill the msequencecontainer object with the next set of sequences and descriptions	
 */
unsigned long msequenceServer::next(const bool _f)
{
/*
 * exit on completion	
 */
	if(done())
		return 0;
/*
 * start the reading process, if it hasn't been started	
 */
	if(!started())	{
		if(!start())	{
			m_bDone = true;
			m_bError = true;
			m_strStatus += "Server would not start.\r\n";
			return 0;
		}
	}
	if(m_lFileType == XBANG)
		return next_pro(_f);
	if(!_f)
		return next_l();
/*
 * initialized the time	
 */
	double dStart = clock();
	unsigned long iLength = 0;
	msequence seqTemp;
	char cValue = '\0';
	char *pValue = NULL;
	char *pEol = NULL;
	m_pCol->clear();
	while(!feof(m_pInput) && iLength < m_pCol->m_tMax)	{
/*
 * store the description in a temporary msequence object, obtained in the previous read	
 */
		m_pCol->m_vASequences[iLength].m_strDes = m_strFirst;
/*
 * strip whitespace characters from the sequence line
 */
		pValue = m_pLine;
		fgets(pValue,m_lSize,m_pInput);
/*
 * clear the sequence in a temporary msequence object	
 */
		while(pValue[0] != '>' && !feof(m_pInput))	{
/*
 * store initial sequence line, and repeat until the next description line is encountered	
 */
			pValue += strlen(pValue);
			pValue--;
			if(pValue > m_pLine)	{
				while(pValue > m_pLine && isspace(*pValue))	{
					pValue--;
				}
				if(!isspace(*pValue) && *pValue != '\0')	{
					pValue++;
					*pValue = '\0';
				}
			}
			fgets(pValue,m_lSize,m_pInput);
		}
		cValue = *pValue;
		*pValue = '\0';
		bz(m_pLine);
		m_pCol->m_vASequences[iLength].m_strSeq = m_pLine;
		m_pCol->m_vASequences[iLength].m_siPath = (short int)(m_vstrPaths.size() - 1);
		*pValue = cValue;
/*
 *	store the next description line
 */
		if(pValue[0] == '>')	{
			if(strchr(pValue,0x01))	{
				pEol = strchr(pValue,0x01);
				*pEol = '\0';
			}
			else	{
				pEol = pValue + strlen(pValue) - 1;
				while(pEol > pValue && isspace(*pEol))	{
					*pEol = '\0';
					pEol--;
				}
			}
			pEol = strchr(pValue,'\r');
			if(pEol)	{
				*pEol = '\0';
			}
			pEol = strchr(pValue,'\n');
			if(pEol)	{
				*pEol = '\0';
			}
			m_strFirst = pValue + 1;
		}
		m_pCol->m_tLength++;
		iLength++;
	}
/*
 * if the current sequence list file is finished, close it and get the next one, otherwise finish	
 */
	if(feof(m_pInput))	{
		if(m_dstrFasta.empty())	{
			finish();
		}
		else	{
			fclose(m_pInput);
			start();
		}
	}
/*
 * store the time required to load the msequencecontainer	
 */
	m_dTime += (double)(clock() - dStart);
	return iLength;
}
bool msequenceServer::bz(char *_p)
{
	if(_p == NULL)	{
		return false;
	}
	char *pFind = strchr(_p,'B');
	while(pFind != NULL)	{
		*pFind = 'N';
		pFind = strchr(_p,'B');
	}
	pFind = strchr(_p,'Z');
	while(pFind != NULL)	{
		*pFind = 'Q';
		pFind = strchr(_p,'Z');
	}
	return true;
}
/*
 * refill the msequencecontainer object with the next set of sequences and descriptions	
 */

unsigned long msequenceServer::next_pro(const bool _f)
{
/*
 * exit on completion	
 */
	if(done())
		return 0;
/*
 * start the reading process, if it hasn't been started	
 */
	if(!started())	{
		if(!start())	{
			m_bDone = true;
			m_bError = true;
			m_strStatus += "Server would not start.\r\n";
			return 0;
		}
	}
/*
 * initialized the time	
 */
	double dStart = clock();
	unsigned long iLength = 0;
	msequence seqTemp;
	register char cValue = '\0';
	char *pValue = NULL;
	char *pEol = NULL;
	m_pCol->clear();
	unsigned long lLength = 0;
	seqTemp.m_strDes = " ";
	seqTemp.m_strSeq = " ";
	size_t tS = 0;
	while(!feof(m_pInput) && iLength < m_tColMax)	{
/*
 * store the description in a temporary msequence object, obtained in the previous read	
 */
		tS = fread(&lLength,4,1,m_pInput);
		if(feof(m_pInput))	{
			break;
		}
#ifdef OSX
		lLength = mac_rev(lLength);
#endif
		tS = fread(m_pLine,lLength,1,m_pInput);
		if(feof(m_pInput))	{
			break;
		}
		if(_f)
			m_pCol->m_vASequences[iLength].m_strDes = m_pLine;
		tS = fread(&lLength,4,1,m_pInput);
#ifdef OSX
		lLength = mac_rev(lLength);
#endif
		tS = fread(m_pLine,lLength,1,m_pInput);
		if(feof(m_pInput))	{
			break;
		}
		if(_f)	{
			bz(m_pLine);
			m_pCol->m_vASequences[iLength].m_strSeq = m_pLine;
			m_pCol->m_vASequences[iLength].m_siPath = (short int)(m_vstrPaths.size() - 1);
		}
		m_pCol->m_vASequences[iLength].m_mapMods.clear();
		m_pCol->m_tLength++;
		iLength++;
	}
/*
 * if the current sequence list file is finished, close it and get the next one, otherwise finish	
 */
	if(feof(m_pInput))	{
		if(m_dstrFasta.empty())	{
			finish();
		}
		else	{
			fclose(m_pInput);
			start();
		}
	}
/*
 * store the time required to load the msequencecontainer	
 */
	m_dTime += (double)(clock() - dStart);
	return iLength;
}
/* 
 * mac_rev was added in version 2004.04.01 so that the OSX version can read files
 * in .pro format that were compiled on a windows or linux box. mac OSX uses the
 * reverse format for reading integers from the disk, so the bytes in the
 * integer have to be reversed
*/
unsigned long msequenceServer::mac_rev(const unsigned long _l)
{
	union sValue	{
		unsigned long ul;
		unsigned char cl[4];
	} lValue;
	lValue.ul = _l;
	union oValue {
		unsigned long ul;
		unsigned char cl[4];
	} olValue;
	olValue.cl[3] = lValue.cl[0];
	olValue.cl[2] = lValue.cl[1];
	olValue.cl[1] = lValue.cl[2];
	olValue.cl[0] = lValue.cl[3];
	return olValue.ul;
}
/*
 * refill the msequencecontainer object with the next set of sequences and descriptions	
 */
unsigned long msequenceServer::next_l(void)
{
/*
 * exit on completion	
 */
	if(done())
		return 0;
/*
 * start the reading process, if it hasn't been started	
 */
	if(!started())	{
		if(!start())	{
			m_bDone = true;
			m_bError = true;
			m_strStatus += "Server would not start.\r\n";
			return 0;
		}
	}
/*
 * initialized the time	
 */
	double dStart = clock();
	unsigned long iLength = 0;
	msequence seqTemp;
	long lSize = (10*4096)-1;
	register char cValue = '\0';
	char *pValue = NULL;
	char *pEol = NULL;
	char *pLine = new char[lSize+1];
	size_t a = 0;
	while(!feof(m_pInput) && iLength < m_pCol->m_tMax)	{
/*
 * store the description in a temporary msequence object, obtained in the previous read	
 */
		fgets(pLine,lSize,m_pInput);
/*
 * clear the sequence in a temporary msequence object	
 */
		while(pLine[0] != '>' && !feof(m_pInput))	{
/*
 * store initial sequence line, and repeat until the next description line is encountered	
 */
			fgets(pLine,lSize,m_pInput);
		}
/*
 *	store the next description line
 */
		if(pLine[0] == '>')	{
			if(strchr(pLine,0x01))	{
				pEol = strchr(pLine,0x01);
				*pEol = '\0';
			}
			else	{
				pEol = pLine + strlen(pLine) - 1;
				while(pEol > pLine && isspace(*pEol))	{
					*pEol = '\0';
					pEol--;
				}
			}
			pEol = strchr(pLine,'\r');
			if(pEol)	{
				*pEol = '\0';
			}
			pEol = strchr(pLine,'\n');
			if(pEol)	{
				*pEol = '\0';
			}
			m_strFirst = pLine + 1;
		}
		iLength++;
	}
	delete pLine;
/*
 * if the current sequence list file is finished, close it and get the next one, otherwise finish	
 */
	if(feof(m_pInput))	{
		if(m_dstrFasta.empty())	{
			finish();
		}
		else	{
			fclose(m_pInput);
			start();
		}
	}
/*
 * store the time required to load the msequencecontainer	
 */
	m_dTime += (double)(clock() - dStart);
	return iLength;
}
/*
 * start the process of loading a sequence list file	
 */
bool msequenceServer::start(void)
{
	m_bStarted = false;
/*
 * return false if there are no more sequence list files	
 */
	if(m_dstrFasta.empty())	{
		return false;
	}
	m_strPath = m_dstrFasta.front();
	m_dstrFasta.pop_front();
	m_vstrPaths.push_back(m_strPath);
/*
 * open the file	
 */
	m_pInput = fopen(m_strPath.c_str(),"rb");
	if(m_pInput == NULL)	{
		m_bError = true;
		m_strStatus = "\n*********\nWarning:\n  Sequence list path '";
		m_strStatus += m_strPath;
		m_strStatus += "'\n  could not be opened and was skipped.\n*********\n\n";
		cout << m_strStatus.c_str();
		return m_bStarted;
	}
	size_t tS = 0;
	char *pS = NULL;
	tS = fread(m_pLine,256,1,m_pInput);
	string strDesc = "no description";
	if(strstr(m_pLine,"xbang-pro-fasta-format") != NULL)	{
		m_lFileType = XBANG;
		char *pV = m_pLine+64;
		if(strlen(pV) > 0)	{
			strDesc = pV;
		}
	}
	else if(m_pLine[0] == '>')	{
		fclose(m_pInput);
		m_lFileType = FASTA;
		m_pInput = fopen(m_strPath.c_str(),"r");
	}
	else	{
		m_lFileType = UNKNOWN;
		m_bError = true;
		m_strStatus = "\n*********\nWarning:\n  Sequence list path '";
		m_strStatus += m_strPath;
		m_strStatus += "'\n  was not in a recognized file format and was skipped.\n*********\n\n";
		cout << m_strStatus.c_str();
		return m_bStarted;
	}
	m_vstrDesc.push_back(strDesc);
	m_bStarted = true;
	m_strStatus += "Path '";
	m_strStatus += m_strPath;
	m_strStatus += "' was opened.\n";
/*
 * read down to the first valid FASTA description line	
 */
	if(m_lFileType == XBANG)
		return m_bStarted;
	pS = fgets(m_pLine,m_lSize,m_pInput);
	while(m_pLine[0] != '>' && !feof(m_pInput))	{
		pS = fgets(m_pLine,m_lSize,m_pInput);
	}
	if(m_pLine[0] == '>')	{
		char *pEol = NULL;
		if(strchr(m_pLine,0x01))	{
			pEol = strchr(m_pLine,0x01);
			*pEol = '\0';
		}
		else	{
			pEol = m_pLine + strlen(m_pLine) - 1;
			while(pEol > m_pLine && isspace(*pEol))	{
				*pEol = '\0';
				pEol--;
			}
		}
		pEol = strchr(m_pLine,'\r');
		if(pEol)	{
			*pEol = '\0';
		}
		pEol = strchr(m_pLine,'\n');
		if(pEol)	{
			*pEol = '\0';
		}
		m_strFirst = m_pLine+1;
	}
/*
 * create the msequencecollection object, if necessary	
 */
	return m_bStarted;
}
/*
 * return true if the server is started	
 */
bool msequenceServer::started(void)
{
	return m_bStarted;
}

/*
 * return true if started, but not finished	
 */
bool msequenceServer::working(void)
{
	return m_bStarted && !m_bDone;
}
