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

#ifndef LOADSPECTRUM_H
#define LOADSPECTRUM_H


// File version: 2004-01-07
/*
loadmspectrum is the base class for deriving classes that take the data from
a file containing one or more ms/ms spectra and load that data into an
mspectrum class for modeling. These derived classes consist of overrides of
the open and get methods that perform file format testing and provide 
compatibility with non-standard file formats.
*/

class loadmspectrum
{
public:
	loadmspectrum(void) { m_tId = 0; m_cEol = 0x0A; m_tSize = 4096*4096;}
	virtual ~loadmspectrum(void) { }

	size_t m_tId; // the id number of a spectrum
	std::streamsize m_tSize;
	string m_strPath; // the path information for the data file
	string m_strTest;
	virtual bool get() {return true; }
	virtual bool get(mspectrum &_m) {return true; } // retrieves a single spectrum from a data file
	virtual bool open(string &_s) {return true; } // attaches the data file to m_ifIn and checks formats
	virtual bool open_force(string &_s) {return true; } // attaches the data file to m_ifIn and checks formats
	int load_test(const char *_p)	{
		m_ifIn.open(m_strPath.c_str());
		if(m_ifIn.fail())	{
//			cout << "<br>Fatal error: input file could not be opened.<BR>"; // rTANDEM
			return 0;
		}
		/*
		 * test for the file extension .mzdata
		 */
		string strTest = m_strPath;
		int (*pf)(int) = tolower; 
		transform(strTest.begin(), strTest.end(), strTest.begin(), pf); 
		if(strTest.find(_p) != strTest.npos)	{
			m_ifIn.close();
			return 2;
		}
		m_strTest.clear();
		size_t tBuffer = 128*1024;
		char *pLine = new char[tBuffer];
		memset((void *)pLine,'\0',tBuffer);
		m_ifIn.getline(pLine,tBuffer,'\n');
		m_strTest += pLine;
		while(m_ifIn.good() && !m_ifIn.eof() && m_strTest.size() < tBuffer)	{
			memset((void *)pLine,'\0',tBuffer);
			m_ifIn.getline(pLine,tBuffer-1,'\n');
			m_strTest += pLine;
		}
		delete pLine;
		m_ifIn.close();
//		cout.flush();
		return 1;
	}
	enum	{
		I_Y =	0x01,
		I_B =	0x02,
		I_X =	0x04,
		I_A =	0x08,
		I_C =	0x10,
		I_Z =	0x20,
	} ion_type; // enum for referencing information about specific ion types.

protected:
	char m_cEol; // the character chosen to mark the end-of-line, 0X0A or 0X0D
	ifstream m_ifIn; // the input file stream
};
/*
loaddta is publicly derived from loadmspectrum. It overrides get and open to work with
DTA files. Multiple ms/ms spectra can be placed in a single file by separating the spectra
with at least one blank line. The DTA file format is as follows:

parent charge
mz1 I1
mz2 I2
...

where: 
parent = the M+H mass of the parent ion in units
charge = the parent ion charge (1,2,3, etc.)
mzX = the mass-to-charge ratio of the Xth fragment ion
IX = the intensity of the Xth fragment ion
*/
class loaddta : public loadmspectrum
{
public:
	loaddta(void);
	virtual ~loaddta(void);

	virtual bool get(mspectrum &_m);
	virtual bool open(string &_s);
	virtual bool open_force(string &_s);
};
/*
loadpkl is publicly derived from loadmspectrum. It overrides get and open to work with
PKL files. Multiple ms/ms spectra can be placed in a single file by separating the spectra
with at least one blank line. The PKL file format is as follows:

parent intensity charge
mz1 I1
mz2 I2
...

where: 
parent = the mass-to-charge of the parent ion
intensity = the intensity of the parent ion
charge = the parent ion charge (1,2,3, etc.)
mzX = the mass-to-charge ratio of the Xth fragment ion
IX = the intensity of the Xth fragment ion
*/

class loadpkl : public loadmspectrum
{
public:
	loadpkl(void);
	virtual ~loadpkl(void);

	virtual bool get(mspectrum &_m);
	virtual bool open(string &_s);
	virtual bool open_force(string &_s);
};
/*
loadmatrix is publicly derived from loadmspectrum. It overrides get and open to work with
Matrix Science files. Please see www.matrixscience.com for a definition of this format.
*/

class loadmatrix : public loadmspectrum
{
public:
	loadmatrix(void);
	virtual ~loadmatrix(void);

	virtual bool get(mspectrum &_m);
	virtual bool open(string &_s);
	virtual bool open_force(string &_s);
};
/*
loadgaml is publicly derived from loadmspectrum. It overrides get and open to work with
the gaml-named nodes in tandem output files. Please see www.gaml.org for a definition of this
format. You can also examine tandem output files to see how this format works.
*/
#include "saxgamlhandler.h"

class loadgaml : public loadmspectrum
{
public:
	loadgaml(vector<mspectrum>& _pv, mspectrumcondition& _p, mscore& _m);
	virtual ~loadgaml(void);

	virtual bool get();
	virtual bool open(string &_s);
	virtual bool open_force(string &_s);

 private:

  mspectrum specCurrent;
  SAXGamlHandler handler;
};

#ifdef XMLCLASS
#include "saxmzxmlhandler.h"
#include "saxmzdatahandler.h"
#include "saxmzmlhandler.h"

/*
loadmzxml est derive publiquement de loadspectrum. Il redefinit get et open. Voi http://sashimi.sourceforge.net/ pour connaitre le mzxml.
Ajouter par Patrick Lacasse en decembre 2004.
 */
class loadmzxml : public loadmspectrum
{
 public:
  loadmzxml(vector<mspectrum>& _pv, mspectrumcondition& _p, mscore& _m);
  virtual ~loadmzxml(void );
  
  virtual bool get();
  virtual bool open(string &_s);
	virtual bool open_force(string &_s);

 private:

  mspectrum specCurrent;
  SAXMzxmlHandler handler;

  string strLine; // Usage temporaire
};

/*
loadmzdata est derive publiquement de loadspectrum. Il redefinit get et open. Voi http://sashimi.sourceforge.net/ pour connaitre le mzxml.
Ajouter par Patrick Lacasse en janvier 2005.
 */
class loadmzdata : public loadmspectrum
{
 public:
  loadmzdata(vector<mspectrum>& _pv, mspectrumcondition& _p, mscore& _m);
  virtual ~loadmzdata(void );
  
  virtual bool get();
  virtual bool open(string &_s);
	virtual bool open_force(string &_s);

 private:

  mspectrum specCurrent;
  SAXMzdataHandler handler;

  string strLine; // Usage temporaire
};
/*
loadmzml was added as an experimental interface to EBI's mzML format
The initial version was adapted from loadmzdata by Ron Beavis in August 2008.
 */
class loadmzml : public loadmspectrum
{
 public:
  loadmzml(vector<mspectrum>& _pv, mspectrumcondition& _p, mscore& _m);
  virtual ~loadmzml(void );
  
  virtual bool get();
  virtual bool open(string &_s);
	virtual bool open_force(string &_s);

 private:

  mspectrum specCurrent;
  SAXMzmlHandler handler;

  string strLine; // Usage temporaire
};

class loadcmn : public loadmspectrum
{
 public:
  loadcmn(void);
  virtual ~loadcmn(void );
  
	virtual bool get(mspectrum &_m);
	virtual bool open(string &_s);
	virtual bool open_force(string &_s);
	int m_iVersion;

 private:

    FILE *m_pFile;
	mspectrum specCurrent;
};
#endif //ifdef XMLCLASS
#endif //ifdef LOADSPECTRUM_H
