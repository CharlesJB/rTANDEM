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

#ifndef MSPECTRUMCONDITION_H
#define MSPECTRUMCONDITION_H

// File version: 2003-07-01
// File version: 2005-01-01

/*
 * the mspectrumcondition object analyzes tandem mass spectra and alters them to
 * be suitable for searching. it loads a set of parameters from an XmlParameters
 * object and uses those values for the spectrum analysis. this object contains
 * more "get" and "set" type methods than most of the other objects in this
 * project: there was considerable experimentation necessary to refine the spectrum
 * analysis methods and these additional methods were used for these experiments.
 *
 * mspectrumcondition performs a series of operations on a spectrum list to make
 * compatible with the scoring methods. the spectrum can be normalized, adjacent
 * peaks removed, parent ion masses excluded, low mass parents excluded, highly
 * charged parents excluded, low mass fragments excluded, and low intensity
 * fragment ions excluded.
 */

#include "xmlparameter.h"
#include "mspectrum.h"

class mscore;

class mspectrumcondition
{
public:
	mspectrumcondition(void);
	virtual ~mspectrumcondition(void);

	bool m_bCondition; // enables the use of all conditioning methods

	bool condition(mspectrum &_s, mscore &_score); // method to apply current conditioning to a spectrum

	bool get_noise_suppression(void); // retrieves the value of the m_bUseNoiseSuppression member

	bool load(XmlParameter &_x); // sets the values for member values using an XmlParameter object

	bool set_max_peaks(const long _p); // sets the maximum number of peaks in a spectrum
	bool set_dynamic_range(const float _p); // sets the maximum dynamic range of a spectrum
	bool set_parent_exclusion(const float _l,const float _u); // sets an upper and lower bound for excluding neutral losses around a parent ion
	bool set_lowest_mass(const float _m); // sets the lowest m/z value in a spectrum
	bool set_max_charge(const long _m); // sets the maximum charge allowable for a spectrum
	bool set_min_size(const long _m); // sets the minimum number of ions in a spectrum
	bool set_min_mass(const float _m); // sets the minimum parent ion mass allowed

	bool use_dynamic_range(const bool _f); // enables using the dynamic range
	bool use_charge_suppression(const bool _f); // enables the rejection of highly charged parents
	bool use_condition(const bool _f); // enables the use of all spectrum conditioning (must be true)
	bool use_lowest_mass(const bool _f); // enables the removal of very low m/z peaks from a spectrum
	bool use_max_peaks(const bool _f); // enables removing low intensity peaks
	bool use_min_size(const bool _f); // enables using the minimum number of peaks to exclude spectra
	bool use_min_mass(const bool _f); // enables using a minimum parent ion mass to exclude spectra
	bool use_noise_suppression(const bool _f); // enables the rejection of noise spectra
	bool use_parent_exclusion(const bool _f); // enables the exclusion of neutrals near the parent ion mass
private:
	bool m_bUseChargeSuppression; // enables the rejection of highly charge parent ions
	bool m_bUseDynamicRange; // enables using the dynamic range
	bool m_bUseLowestMass; // enables the removal of very low m/z peaks from a spectrum
	bool m_bUseMaxPeaks; // enables removing low intensity peaks
	bool m_bUseMinMass; // sets the minimum parent ion mass allowed
	bool m_bUseMinSize; // enables using the minimum number of peaks to exclude spectra
	bool m_bUseNoiseSuppression; // enables the detection of purely noise spectra
	bool m_bUseParent; // enables the exclusion of spectra by the parent ion mass
	bool m_bUseNeutralLoss;
	bool m_bUsePhosphoDetection;
	bool m_bUseAllowedNeutralLosses;
	size_t m_tMaxPeaks; // the maximum number of peaks in a spectrum
	float m_fDynamicRange; // the normalized intensity of the most intense peak in a spectrum
	float m_fLowestMass; // the lowest m/z in a spectrum
	long m_lMinSize; // the minimum number of peaks in a spectrum
	float m_fMinMass; // the minimum parent ion mass in a spectrum
	float m_fParentLower; // the low end of the mass window for excluding parent ion neutral losses
	float m_fParentUpper; // the high end of the mass window for excluding parent ion neutral losses (just passed the last C13 isotope peak)
	long m_lMaxCharge; // the maximum parent ion charge allowed
	float m_fNeutralLoss;
	float m_fNeutralLossWidth;
	float m_fFactor;
	float m_fMaxZ;
	vector<double> m_vdNeutralLosses;

	bool find_loss(mspectrum &_s,const float _d,const float _t,const float _p = 0.10);
	bool check_neutral(mspectrum &_s); // determines if there are neutral loss peaks in a spectrum
	bool clean_isotopes(mspectrum &_s); // removes isotope peaks
	bool dynamic_range(mspectrum &_s); // normalizes the spectrum using the dynamic range
	bool is_noise(mspectrum &_s); // determines if a spectrum is purely noise
	bool remove_isotopes(mspectrum &_s); // removes peaks too close to be real isotope peaks
	bool remove_neutral(mspectrum &_s); // removes neutral loss ions
	bool remove_parent(mspectrum &_s); // removes parent ions & neutrals
	bool remove_small(mspectrum &_s); // removes low intensity peaks
	bool remove_low_masses(mspectrum &_s); // removes low m/z peaks
};
#endif
