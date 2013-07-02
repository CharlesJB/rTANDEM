/*
 Copyright (C) 2003 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project
 
 link to  libguide.lib if using the Intel Compiler

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
// File version: 2004-04-01
// 2005-12-01

/*
	stdafx.h contains the include files used in common by the files within the project. The
	includes are mostly parts of the Standard Template Library for C++. STL is used frequently
	throughout the project, to provide containers and algorithms.
	This file also contains the preprocessor directives necessary to compile the project on
	different platforms.
*/

//#define OSX // to compile with cc on Macintosh OSX
//#define OSX_TIGER // use with OSX if OSX Tiger is being used
//#define OSX_INTEL // compile with cc on Macintosh OSX on Intel-style processors
//#define GCC // to compile with gcc (or g++) on LINUX
//#define MSVC // to compile with Microsoft Studio

//#define X_P3 // when this is defined, the proteotypic peptide version of X! Tandem, X! P3 is built

/*
	NOTE: Change with version 2007.04.01 in handling the OS define statements.
		From this version forward, these defines will be set either in the
		makefile or the preprocessor definitions with MSVC.
		The past practice of changing the defines in this file will be
		through commenting out unwanted definitions will be discontinued.
		The pre-existing define statements will be left in place as comments.
*/

#define XMLCLASS // to compile with MzData and MzXML classes
#define XML_STATIC // to statically link the expat libraries

/* CROSS-PLATFORM COMPILATION NOTES
	1. 32 and 64 bit integers do not have any standard type definitions that can
	be relied on to work properly across different platforms. Therefore
	there are specific #define statements made to compensate for platform
	specific difference. If your platform is complaining about
	any of the data types defined below, it may be necessary to tinker
	with the #defines to get them to work properly.
	2. Extensive use is made of the STL data type size_t. This type may be defined
	as different types on different platforms.
	3. The file "base64.cpp" tends to cause problems with this type of compatibility:
	if you change the defines below, you may have to alter base64.cpp as well.
*/
/* // rTANDEM
#ifdef MSVC
	#define __inline__ _inline
	#define __int64_t _int64
	#define uint32_t unsigned long
	#define uint64_t unsigned _int64
	// this define was made necessary in VC++ 2005 to avoid 
	// unnecessary warnings associated with the new "secure" versions
	// of strcpy and related routines.
	#define _CRT_SECURE_NO_DEPRECATE
	#include <string>
#endif
*/

#include <stdint.h>

#ifdef OSX
	#define __inline__ inline
	#ifndef OSX_TIGER
		#define __int64_t int64_t
	#endif
	#include <string>
#endif
// this define for the INTEL-based OSX platform is untested and may not work
#ifdef OSX_INTEL
	#define __inline__ inline
	#include <string>
#endif
// this define should work for most LINUX and UNIX platforms
#ifdef GCC
	#include <stdio.h>
#endif

#ifdef GCC4
	#include <stdio.h>
	// this switch was added for compatibility with GCC v. 4
#endif
#ifdef GCC4_3
	#include <cstring>
	// this switch was added for compatibility with GCC v. 4.3
#endif

// common includes for the standard template library
	// this test was suggested by Steve Wiley to correct a problem
	// associated with compiling using the 64-bit version of Redhat Linux

// rTANDEM This was all replaced by #include <stdint.h>
// #ifndef __APPLE__ 
// #if defined(__x86_64__) && defined(__linux__) 
//#ifndef uint32_t
//  #define uint32_t unsigned int
//#endif
//#ifndef uint64_t
//  #define uint64_t unsigned long long
//#endif

// #endif 
#if defined(WIN32) || defined(_WIN32)
	#define __inline__ inline
	#define __int64_t long long
  //	#define uint32_t unsigned long
  //    #define uint64_t unsigned long long
#endif
// #endif
// rTANDEM - 2 last lines

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <deque>
#include <cmath>
#include <vector>
#include <map>
using namespace std;


// These typedefs are used in several places in the code and are included here
// for simplicity.
typedef map<size_t,double> SMap;
typedef map<string,string,less<string> > xMap;
typedef map<string,bool> bMap;
typedef pair<char,string> prSap;

#define VERSION "Jackhammer (2013.06.15)"

