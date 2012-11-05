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

#ifndef MSCORESTATE_H
#define MSCORESTATE_H

// File version: 2003-07-01
// File version: 2004-03-01
/*
 * mscorestate is a specialty class used to store information necessary for the state machine that
 * mscore uses to generate all potentially modified peptide sequences and masses
 */


class mscorestate
{
public:
      mscorestate(void) {    
            m_lSizeS = 128;
            m_pSeqS = new char[m_lSizeS];
            m_ppModsS = new char*[m_lSizeS];
            m_piMods = new unsigned long[m_lSizeS];
			m_lSizeEqualsS = 256;
            m_plEqualsS = new long[m_lSizeEqualsS];
            m_lEqualsS = 0;
            m_lEligibleS = 0;
			m_bIsPossible = true;
      }
      virtual ~mscorestate(void) {
            if(m_plEqualsS != NULL)
                  delete m_plEqualsS;
            if(m_pSeqS != NULL)
                  delete m_pSeqS;
            if(m_ppModsS != NULL)
                  delete m_ppModsS;
      }
      bool m_bStateS; // true if there are more potential sequences
	  bool m_bIsPossible; // true if a modified state is a potential solution
      double m_dSeqMHS; // M+H mass of the unmodified peptide
      double m_dSeqMHFailedS; // the value of M+H that failed the check_parents test last time through
      long m_lEligibleS; // number of spectra with M+H <= the current modified peptide M+H
      long m_lEqualsS; // number of spectra with M+H within error of the current modified peptide M+H
      unsigned long m_lCursorS; // cursor counting the position of the N-terminal modified residue
      unsigned long m_lFilledS; // cursor index to the most C-terminal residue modified,
                                            //where all residues between m_lFirstS & this one are modified
      unsigned long m_lLastS; // cursor index to the most C-terminal residue modified in a peptide,
                                        // effectively the length of m_ppModsS
      unsigned long m_lStates;
      long m_lSizeS; // maximum peptide length
      long m_lSizeEqualsS; // maximum length of EqualS array
      long *m_plEqualsS; // indexes of spectra with M+H within error of the current modified peptide M+H
      unsigned long *m_piMods;
      char **m_ppModsS; // pointers to the potentially modified peptide residues
      char *m_pSeqS; // the unmodified peptide sequence
/*
 * make sure that m_lSizeS is big enough for a peptide, length _s. if not, update array sizes
 */   bool check_size(const long _s)     
      {
            if(_s > m_lSizeS) {
                  m_lSizeS = _s + 1;
                  delete m_pSeqS;
                  delete m_ppModsS;
                  m_pSeqS = new char[m_lSizeS];
                  m_ppModsS = new char*[m_lSizeS];
            }
            return true;
      }
/*
 * create the appropriate m_plEquals array
 */
      bool create_equals(const long _s)
      {
		  if(_s < m_lSizeEqualsS)
			  return true;
          if(m_plEqualsS != NULL)
               delete m_plEqualsS;
		  m_lSizeEqualsS = _s + 1;
          m_plEqualsS = new long[m_lSizeEqualsS];
          if(m_plEqualsS == NULL)
				return false;
          return true;
      }
/*
 * set up the state machine for the next sequence
 */
      bool initialize(const char *_p,const long _s)
      {
            check_size(_s);
            strcpy(m_pSeqS,_p);
            m_lStates = 0;
            m_lFilledS = 0;
            m_lCursorS = 0;
            m_bStateS = true;
            m_lEligibleS = 0;
            m_dSeqMHFailedS = 0.0;
			m_bIsPossible = true;
            return true;
     }
/*
* add a pointer to a potentially modified residue to the m_ppModsS list
*/
      long add_mod(char *_p)
      {
            m_ppModsS[m_lLastS] = _p;
            m_lLastS++;
            return m_lLastS;
      }

      static unsigned long M_lMaxModStates;
};


class mscoreterm
{
public:
	mscoreterm(void) { 	
		m_lN = 0;
		m_lC = 0;
		m_bN = false;
		m_bC = false;
	}
	virtual ~mscoreterm(void) { 
	}
	bool m_bN; // true if N terminal modification allowed
	bool m_bC; // true if C terminal modification allowed
	long m_lC; // 0 if C-terminal not modified, 1 if C-terminal modified
	long m_lN; // 0 if N-terminal not modified, 1 if N-terminal modified
	long m_lState; // Tracks the current state of the machine: value is the binary number formed by
				   // assuming m_lC and m_lN are ordered bits in a two bit number, e.g.
				   // state => m_lC = 1 & m_lN = 0 => m_lState = binary number '10' = 2.
/*
 * set up the state machine for the next sequence
 */
	bool initialize(const double _n,const double _c)
	{
		m_bN = false;
		m_bC = false;
		if(fabs(_n) > 0.001)	{
			m_bN = true;
		}
		if(fabs(_c) > 0.001)	{
			m_bC = true;
		}
		m_lState = 0;
		m_lN = 0;
		m_lC = 0;
		return true;
	}
};
#endif
