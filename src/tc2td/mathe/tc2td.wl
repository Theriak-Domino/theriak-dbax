(* ::Package:: *)


(*
    Main code of program tc2td that converts Thermocalc formatted 
    datasets and ax files to Theriak-Domino format.

    Copyright (C) 2022  Doug Tinkham

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*)

(*
   Note: this package is currently dependent on functions from other
   packages that are not included. These functions are currently being
   added to this package and will eventually appear here, or in the 
   cpp version of this package.
*)

(* ::Title::Regular:: *)
(*Package ngen`tc2td`*)


(* ::Text:: *)
(*Doug Tinkham, March 2022*)
(**)
(*2022-03-29: *)
(*Initial version of new TC->TD ds and ax file converter. It reuses parts of existing*)
(*ngen`dset`, `thermo` and `utils` but should be free of ngen`dsaxcnvrt`; Most *)
(*functions started as copies from ngen`dsaxcnvrt` and then modified. Eventually *)
(*move into ngen2`, when parts of the dataset structures will likely need updated. *)
(*Some new functions were added to ngen`utils` and need copied over to ngen2`. *)
(*Developed in file  TC2TD-metabasite-hpxeos-2022-03-20.nb.*)
(**)
(*Package to read a TC ds and ax file, and convert them to a TD formatted ds and ax file.*)


(* ::Text:: *)
(*......................................................................................................................................................................*)


(* ::Section::Closed:: *)
(*Begin*)


axrefs[_]:={};

BeginPackage["ngen`tc2td`", {"ngen`utils`", "ngen`dset`", "ngen`thermo`",
"ngen`axdata`","ngen`axmodel`"}];


(* ::Section:: *)
(*Usage*)


(* usage *)

GetHP11TDHeader::usage=
"";

GetBufferStuff::usage=
"";

GetTC2TDAXReferences::usage=
""

MatchPatt::usage=
"MatchPatt[arg_?StringQ]:=
  Returns a string inidating a match to a custom StringPattern for
  the arg arg_. The matches returned are useful for determining if
  the string can be converted to a number. 
  Return strings:
    MultExpPatt: string of form 1.23*10^-07.
    ScientificPatt: string is of form -1.23E-07,
      +1.23D-07, 1.23e-07, -1.23d+07.
    RealPatt: string is of form +1.23, etc.
    IntegerPatt: string is of form -123.
    FlagPatt: string is of form %<AString> that is used as 
      a key flag in tc files I convert to td format.
    String: string is not in the form of a number.
    ERROR: string not matched by any of the patterns above.
";

GetEx::usage=
"Works well to convert all strings representing single "<>
"number to a number, including simple rational pattern "<>
"like 1/3 or 1.0/3.0";

Options[GetEx]:={"ConvertFortranNum"->False,
                 "DoArithmetic"->False};

GenExListFromString::usage=
"";

Options[GenExListFromString]:={"ConvertFortranNum"->False,
                               "DoArithmetic"->False};
                               
GenExListFromStringList::usage=
"";

Options[GenExListFromStringList]:={"ConvertFortranNum"->False,
                                   "DoArithmetic"->False};                                                              
                               
GetCleanTCAxLines::usage=
"";

ParseCleanTCAxLines::usage=
"Returns a nested master ax data list, where each sublist is all\n"<>
" the data for each a-x model extracted from the tc-ax file. The \n"<>
" last sublist is actually the pure phases listed at bottom of tc \n"<>
" a-x file (so is extracted in LoadTCAXModes). Each sublist has the\n"<>
" following 19 element data structure:\n"<>
" {ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,\n"<>
" nsitex,sitemult,effsitemult,sitexls,aidlns,makelns,dqflns,checklns,info}\n"<>
" the last element, info, is a nested list with the following structure:\n"<>
" {{skip,False},{name,CAMP},{ref,{HGP18,HPx-eos}},{note,Some text.}}";

GenMinusSpeciesData::usage=
"";
                              
ProcessMakeVecTC350::usage=
"";

NewAnalSites::usage=
"";       

LoadTCAXModels::usage=
"Returns {models,dsabbls,pures,mades,madefroms} \n"<>
"where models is a 22 element list with this data:\n"<>
"{ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,\n"<>
"nsitex,sitemult,esitemult,sitexls,aidlns,mdlmklns,dqflns,\n"<>
"checklns,mdlmades,mdlmadefroms,verbatims[[i]],info} ";

Options[LoadTCAXModels]={"LoadExtras"->False,"ExtrasFile"->"Unknown"};

WriteTDAXFile::usage=
"";

Options[WriteTDAXFile]={"UseF3"->False,"TCDSFile"->"HPxeos ds website file",
  "TCAXFile"->"HPxeos ax website file","TCDSver"->"unknown dsver"};


(* ::Section::Closed:: *)
(*Unprotect*)


(* unprotect *)
Unprotect[GetHP11TDHeader,GetTC2TDAXReferences,GetBufferStuff,MatchPatt,
GetEx,GenExListFromString,GenExListFromStringList,GetCleanTCAxLines,
ParseCleanTCAxLines,GenMinusSpeciesData,ProcessMakeVecTC350,
NewAnalSites,LoadTCAXModels,WriteTDAXFile];

Begin["Private`"];


(* ::Section:: *)
(*Code*)


(* ::Text:: *)
(*......................................................................................................................................................................*)


(* ::Subsubsection::Closed:: *)
(*GetHP11TDHeader, GetTC2TDAXReferences, GetBufferStuff, GetPCNames*)


GetHP11TDHeader[]:=Block[{},
  Return[
    "19  8.31441  0\n"<>
    "     O         SI        TI        AL        FE        MG        MN\n"<>
    "     CA        NA        K         H         C         CL        NI\n"<>
    "     ZR        S         CU        CR        E\n"<>
    "  15.99940  28.08550  47.86700  26.98154  55.84500  24.30500  54.93805\n"<>
    "  40.07800  22.98977  39.09830   1.00794  12.01070  35.45300  58.69340\n"<>
    "  91.22400  32.06500  63.54600  51.9961    1.00000\n"<>
    "   0.0       2.0       2.0       1.5       1.0       1.0       1.0\n"<>
    "   1.0       0.5       0.5       0.5       2.0       0.0       1.0\n"<>
    "   2.0       0.0       1.0       1.5       0.0"
  ];
];

GetTC2TDAXReferences[]:=Block[{refstr,tstr,strls1,strls,refs={},crec={},ts,lnum=0},
  refstr=
    "  REFERENCES
    ref HP96:
        Holland, TJB & Powell, R (1996). Thermodynamics of order-disorder in 
        minerals: II. Symmetric formalism applied to solide solutions. American
        Mineralogist, 81, 1425-1437.
    ref HP98:
        Holland, TJB & Powell, R (1998). An internally consistent thermodynamic 
        data set for phases of petrological interest. Journal of Metamorphic
        Geology, 16, 309-343. 
    ref WPH00: 
        White, RW, Powell, R, Holland, TJB & Worley, BA (2000). The effect of 
        TiO2 and Fe2O3 on metapelitic assemblages at greenschist and amphibolite 
        facies conditions: mineral equilibria calculations in the system K2O-FeO-
        MgO-Al2O3-SiO2-H2O-TiO2-Fe2O3. Journal of Metamorphic Geology, 18, 
        497-511.
    ref WPC02: 
        White, RW, Powell, R & Clarke, GL (2002). The interpretation of reaction 
        textures in Fe-rich metapelitic granulites of the Musgrave Block, central 
        Australia: constraints from mineral equilibria calculations in the system 
        K2O-FeO-MgO-Al2O3-SiO2-H2O-TiO2-Fe2O3. Journal of Metamorphic Geology,  
        20, 41-55.
    ref HP03: 
        Holland, T & Powell, R (2003). Activity-composition relations for phases
        in petrological calculations; an asymmetric multicomponent formulation.
        Contributions to Mineralogy and Petrology 145, 492-501.
        doi: 10.1007/s00410-003-0464-z
    ref DPW07: 
        Diener, JFA, Powell, R, White, RW & Holland, TJB (2007). A new thermo-
        dynamic model for clino- and orthoamphiboles in Na2O-CaO-FeO-MgO-Al2O3
        -SiO2-H2O-O. Journal of Metamorphic Geology, 25, 631-656.
    ref DWP08:
        Diener, JFA, White, RW, & Powell, R (2008). Granulite facies metamorphism
        and subsolidus fluid-absent reworking, Strangways Range, Arunta Block,
        central Australia. Journal of Metamorphic Geology, 26, 603-622.
    ref EPH10: 
        Evans KA, Powell, R & Holland, TJB (2010). Internally consistent data
        for sulphur-bearing phases and application to the construction of
        pseudosections for mafic greenschist facies rocks in Na2O-CaO-K2O-FeO-
        MgO-Al2O3-SiO2-CO2-O-S-H2O. Journal of Metamorphic Geology, 28, 667-687.
    ref HP11:  
        Holland, TJB & Powell, R (2011). An improved and extended internally
        consistent thermodynamic dataset for phases of petrological interest,
        involving a new equation of state for solids. Journal of Metamorphic
        Geology, 29, 333-383.
    ref GHP12:  
        Green, ECR, Holland, TJB, Powell, R & White, RW (2012). Garnet and spinel
        lherzolite assemblages in MgO-Al2O3-SiO2 and CaO-MgO-Al2O3-SiO2: thermo-
        dynamic models and an experimental conflict. Journal of Metamorphic
        Geology, 30, 561-577.
    ref EPF13: 
        Evans KA, Powell, R & Frost, BR (2013). Using equilibrium thermodynamics
        in the study of metasomatic alteration, illustrated by an application to
        serpentinites. Lithos 168-169, 67-84.
    ref HHP13:  
        Holland, TJB, Hudson, NFC, Powell, R & Harte, B (2013). New thermodynamic
        models and calculated phase equilibria in NCFMAS for basic and ultrabasic
        compositions through the transition zone into the uppermost lower mantle.
        Journal of Petrology, 54, 1901-1920.
    ref PWG14:  
        Powell, R, White, RW, Green, ECR, Holland, TJB, & Diener, JFA (2014).
        On parameterising thermodynamic descriptions of minerals for petrological
        calculations. Journal of Metamorphic Geology, 32 (3), 245-260.
    ref WP14: 
        Wheller, CJ & Powell, R (2014). A new thermodynamic model for sapphirine: 
        calculated phase equilibria in K2O-FeO-MgO-Al2O3-SiO2-H2O-TiO2-Fe2O3.
        Journal of Metamorphic Geology, 32, 287-299.
    ref WPH14:  
        White, RW, Powell, R, Holland, TJB, Johnson, TE & Green, ECR (2014). New
        mineral activity-composition relations for thermodynamic calculations in
        metapelitic systems. Journal of Metamorphic Geology, 32, 261-286. 
        doi: 10.1111/jmg.12071
    ref WPJ14:  
        White, RW, Powell, R & Johnson, TE (2014). The effect of Mn on mineral
        stability in metapelites revisited: new a-x relations for manganese-
        bearing minerals. Journal of Metamorphic Geology, 32, 809-828. 
        doi: 10.1111/jmg.12095
    ref EP15: 
        Evans, KA & Powell, R (2015). The effect of subduction on the sulphur
        carbon and redox budget of lithospheric mantle. Journal of Metamorphic
        Geology, 33, 649-670.
    ref JH15:  
        Jennings, ES & Holland, TJB (2015). A simple thermodynamic model for
        melting of peridotite in the system NCFMASOCr. Journal of Petrology,
        v. 56, no. 5, 869-892. doi: 10.1093/petrology/egv020
    ref GWD16:  
        Green, ECR, White, RW, Diener, JFA, Powell, R, Holland, TJB & Palin, RM
        (2016). Activity-composition relations for the calculation of partial
        melting equilibria in metabasic rocks.  Journal of Metamorphic Geology,
        v. 34, 845-869. doi:10.1111/jmg.12211
    ref PWG16:  
        Palin, RM, White, RW, Green, ECR, Diener, JFA., Powell, R & Holland, TJB
        (2016). High-grade metamorphism and partial melting of basic and inter-
        mediate rocks. Journal of Metamorphic Geology, v. 34, 871-892.
        doi:10.1111/jmg.12212
    ref HGP18:  
        Holland, TJB, Green, ECR & Powell, R (2018). Melting of peridotites
        through to granites: a simple thermodynamic model in the system 
        KNCFMASHTOCr. Journal of Petrology v. 59, 881-900. 
        doi: 10.1093/petrology/egy048.
    ref EF21: 
        Evans, KA & Frost, BR (2021). Deserpentinization in subduction zones as 
        a source of oxidation in arcs: A reality check. Journal of Petrology,
        1-32. doi: 10.1093/petrology/egab016
    ref HGP21: 
        Holland, TJB, Green, ECR & Powell, R (2021). A thermodynamic model for  
        feldspars in KAlSi3O8-NaAlSi3O8-CaAl2Si2O8 for mineral equilibrium   
        calculations. Journal of Metamorphic Geology, 1-14. doi:10.1111/jmg.12639
    ref TH21:
        Tomlinson, EL & Holland, TJB (2021). A thermodynamic model for the sub-
        solidus evolution and melting of peridotite. Journal of Petrology, 1-23.
        doi: 10.1093/petrology/egab012
    ref HPx-eos:
        Model contains modifications posted in HPx-eos website data files since
        the time of publication by cited authors.
        https://hpxeosandthermocalc.org/downloads/#thermocalc
  "; (* end refstr *)
  
  strls1=ReadList[StringToStream[refstr],String];
  (*Print[strls1//InputForm];*)
  (*Print[Length[strls1]];*)
  strls={};
  Do[
    tstr=StringTrim[strls1[[i]]];
    If[StringLength[tstr]>0, AppendTo[strls,tstr]];
    ,{i,Length[strls1]}
  ];
  (*Print["ok1"];*)
  (*Print[{"strls = ",strls//InputForm}];*)
  lnum=0;
  Do[
    If[ts=StringTrim[strls[[++lnum]]];StringStartsQ[ts,"ref "]&&StringEndsQ[ts,":"],
      If[Length[crec]>0,AppendTo[refs,crec]; crec={StringDrop[StringDrop[ts,-1],4]};];
      ,
      AppendTo[crec,ts]
    ];
    ,{i,Length[strls]}
  ];
  AppendTo[refs,crec];
  If[StringStartsQ[refs[[1,1]],"REFERENCES"],refs=refs[[2;;Length[refs]]]];
  
  (*Print[{"refs[[-1]] = ",refs[[-1]]//InputForm}];*)
  
  Print["Done"];
  Return[refs];
];

GetBufferStuff[]:=Module[{bstr},
  bstr=
"! BUFFER STUFF
! HHM-BUFFER: H(2); requires: magnetite, hematite, H2O
! OHM-BUFFER: O(2); requires: magnetite, hematite
!
! HQFM-BUFFER:H(2); requires: quartz, fayalite, magnetite, H2O
! QFM-BUFFER: O(2); requires: quartz, fayalite, magnetite
!
! NNO-BUFFER: O(2); requires: Ni, NiO
! GC-BUFFER:  O(2); requires; graphite, CO2
! PPS-BUFFER: S(2); requires: pyrite, trot (FeS)
! Ox-Buffer:  O(2); fixed activity
!
! Replaces species in COM statements with those you are using in
! your calculations (can/should be minus species, H2O or CO2)
!----------------------------------------------------------------
!
***** GAS DATA *****
!----------------------------------------------------------------
!                   H2 = 2 Fe3O4 + 1 H2O - 3 Fe2O3
HHM-BUFFER     H(2)                            HHM           *HHM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          mtD-[2]H2O[1]hemD-[-3]
!--------------
-HHM-BUFFER    H(-2)                          -HHM           *HHM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          HHM-BUFFER[-1]
!----------------------------------------------------------------
!                      O2 = 6 Fe2O3 - 4 Fe3O4
OHM-BUFFER     O(2)                            OHM           *OHM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          hemD-[6]mtD-[-4]      0  hemD-
!--------------
-OHM-BUFFER    O(-2)                          -OHM           *OHM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          OHM-BUFFER[-1] 
!----------------------------------------------------------------
!                   O2 = 3 SiO2 + 2 Fe3O4 - 3 Fe2SiO4 
QFM-BUFFER     O(2)                            QFM           *QFM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          fa[-3]mtD-[2]q[3]     0  mtD-
!--------------
-QFM-BUFFER    O(-2)                          -QFM           *QFM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          QFM-BUFFER[-1]
!----------------------------------------------------------------
!           H2 = 3/2 Fe2SiO4 + 1 H2O - 1 Fe3O4 - 3/2 SiO2
HQFM-BUFFER     H(2)                          HQFM          *HQFM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          fa[3/2]H2O[1]mtD-[1]q[-3/2]     
!--------------
-HQFM-BUFFER    H(-2)                        -HQFM          *HQFM
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          HQFM-BUFFER[-1]
!----------------------------------------------------------------
!                         O2 = 2 NiO - 2 Ni
NNO-Buffer     O(2)                            NNO           *NNO
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          NiO[2]Ni[-2]          0         NiO
!--------------
-NNO-Buffer    O(-2)                          -NNO           *NNO
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          NNO-Buffer[-1]
!----------------------------------------------------------------
!                         O2 = 1 CO2 - 1 C
GC-BUFFER      O(2)                             GC            *GC
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          CO2[1]gph[-1]          0        CO2
!--------------
-GC-BUFFER     O(-2)                           -GC            *GC
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          GC-BUFFER[-1]
!----------------------------------------------------------------
!                        S2 = 2 FeS2 - 2 FeS
PPS-BUFFER     S(2)                            PPS           *PPS
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          pyr[2]trot[-2]         0        pyr
!--------------
-PPS-BUFFER    S(-2)                          -PPS           *PPS
  ST           0.000          0.000         0.0000          0.000
  CP1          0.000          0.000         0.0000          0.000
  COM          PPS-BUFFER[-1]
!----------------------------------------------------------------
!
";
  Return[bstr];
];

(* These are not used yet. Need completed. *)
GetPCNames[]:=Module[{arr,ds634dat},
  arr={
    {"fo", "forsterite", "olivine"},
    {"fa", "fayalite", "olivine"},
    {"teph", "tephroite", "olivine"},
    {"lrn", "larnite", "group"},
    {"mont", "monticellite", "group"},
    {"chum", "clinohumite", "group"},
    {"chdr", "chondrodite", "group"},
    {"mwd", "name", "group"},
    {"fwd", "name", "group"},
    {"mrw", "name", "group"},
    {"frw", "name", "group"},
    {"mpv", "name", "group"},
    {"fpv", "name", "group"},
    {"apv", "name", "group"},
    {"npv", "name", "group"},
    {"ppv", "name", "group"},
    {"cpv", "name", "group"},
    {"mak", "name", "group"},
    {"fak", "name", "group"},
    {"maj", "majorite", "garnet"},
    {"nagt", "name", "garnet"},
    {"py", "pyrope", "garnet"},
    {"alm", "almandine", "garnet"},
    {"spss", "spessartine", "garnet"},
    {"gr", "grossular", "garnet"},
    {"andr", "andradite", "garnet"},
    {"ski", "skiagite", "garnet"},
    {"knor", "knorringite", "garnet"},
    {"uv", "uvarovite", "garnet"},
    {"osma", "name", "group"},
    {"osmm", "name", "group"},
    {"osfa", "name", "group"},
    {"vsv", "name", "group"},
    {"and", "andalusite", "group"},
    {"ky", "kyanite", "group"},
    {"sill", "sillimanite", "group"},
    {"smul", "name", "group"},
    {"amul", "name", "group"},
    {"tpz", "topaz", "group"},
    {"mst", "Mg-staurolite", "group"},
    {"fst", "Fe-staurolite", "group"},
    {"mnst", "Mn-staurolite", "group"},
    {"mctd", "Mg-chloritoid", "group"},
    {"fctd", "Fe-chloritoid", "group"},
    {"mnctd", "Mn-chloritoid", "group"},
    {"merw", "name", "group"},
    {"spu", "name", "group"},
    {"zo", "zoisite", "group"},
    {"cz", "clinozoisite", "group"},
    {"ep", "epidote", "group"},
    {"fep", "name", "group"},
    {"pmt", "name", "group"},
    {"law", "lawsonite", "group"},
    {"mpm", "name", "group"},
    {"fpm", "name", "group"},
    {"jgd", "name", "group"},
    {"geh", "name", "group"},
    {"ak", "name", "group"},
    {"rnk", "name", "group"},
    {"ty", "name", "group"},
    {"crd", "Mg-cordierite", "cordierite"},
    {"hcrd", "hyd-Mg-cordierite", "cordierite"},
    {"fcrd", "Fe-cordierite", "cordierite"},
    {"mncrd", "Mn-cordierite", "cordierite"},
    {"phA", "name", "group"},
    {"phD", "name", "group"},
    {"phE", "name", "group"},
    {"shB", "name", "group"},
    {"sph", "name", "group"},
    {"cstn", "name", "group"},
    {"zrc", "name", "group"},
    {"zrt", "name", "group"},
    {"tcn", "name", "group"},
    {"en", "name", "pyroxene"},
    {"pren", "name", "pyroxene"},
    {"cen", "name", "pyroxene"},
    {"hen", "name", "pyroxene"},
    {"hfs", "name", "pyroxene"},
    {"fs", "name", "pyroxene"},
    {"mgts", "name", "pyroxene"},
    {"di", "name", "pyroxene"},
    {"hed", "name", "pyroxene"},
    {"jd", "name", "pyroxene"},
    {"kjd", "name", "pyroxene"},
    {"acm", "name", "pyroxene"},
    {"kos", "name", "pyroxene"},
    {"cats", "name", "pyroxene"},
    {"caes", "name", "group"},
    {"rhod", "name", "group"},
    {"pxmn", "name", "group"},
    {"wo", "wollastonite", "group"},
    {"pswo", "name", "group"},
    {"wal", "name", "group"},
    {"tr", "name", "amphibole"},
    {"fact", "name", "amphibole"},
    {"ts", "name", "amphibole"},
    {"parg", "name", "amphibole"},
    {"gl", "name", "amphibole"},
    {"fgl", "name", "amphibole"},
    {"nyb", "name", "amphibole"},
    {"rieb", "name", "amphibole"},
    {"anth", "name", "amphibole"},
    {"fanth", "name", "amphibole"},
    {"cumm", "name", "amphibole"},
    {"grun", "name", "amphibole"},
    {"ged", "name", "amphibole"},
    {"spr4", "name", "group"},
    {"spr5", "name", "group"},
    {"fspr", "name", "group"},
    {"mcar", "name", "group"},
    {"fcar", "name", "group"},
    {"deer", "name", "group"},
    {"mu", "muscovite", "group"},
    {"cel", "celadonite", "group"},
    {"fcel", "Fe-celadonite", "group"},
    {"pa", "paragonite", "group"},
    {"ma", "margarite", "group"},
    {"phl", "phlogopite", "group"},
    {"ann", "annite", "group"},
    {"fbi", "name", "group"},
    {"obi", "name", "group"},
    {"mnbi", "Mn-biotite", "group"},
    {"east", "eastonite", "group"},
    {"naph", "Na-phlogopite", "group"},
    {"tan", "name", "group"},
    {"clin", "name", "group"},
    {"ames", "name", "group"},
    {"afchl", "Al-free-chlorite", "group"},
    {"daph", "daphnite", "group"},
    {"mnchl", "Mn-chlorite", "group"},
    {"sud", "sudoite", "group"},
    {"fsud", "Fe-sudoite", "group"},
    {"prl", "pyrophyllite", "group"},
    {"ta", "talc", "group"},
    {"fta", "name", "group"},
    {"tats", "name", "group"},
    {"tap", "name", "group"},
    {"nta", "name", "group"},
    {"minn", "minnesotaite", "group"},
    {"minm", "name", "group"},
    {"kao", "kaolinite", "group"},
    {"pre", "prehnite", "group"},
    {"fpre", "Fe-prehnite", "group"},
    {"chr", "name", "group"},
    {"liz", "name", "group"},
    {"glt", "name", "group"},
    {"fstp", "Fe-stilpnolmelane", "group"},
    {"mstp", "Mg-stilpnomelane", "group"},
    {"atg", "antigorite", "group"},
    {"ab", "albite", "group"},
    {"abh", "albite-high", "group"},
    {"mic", "microcline", "group"},
    {"san", "sanidine", "group"},
    {"an", "anorthite", "group"},
    {"kcm", "name", "group"},
    {"wa", "name", "group"},
    {"hol", "name", "group"},
    {"q", "quartz", "group"},
    {"trd", "tridymite", "group"},
    {"crst", "cristobalite", "group"},
    {"coe", "coesite", "group"},
    {"stv", "stishovite", "group"},
    {"ne", "nepheline", "group"},
    {"cg", "name", "group"},
    {"cgh", "name", "group"},
    {"macf", "name", "group"},
    {"mscf", "name", "group"},
    {"fscf", "name", "group"},
    {"nacf", "name", "group"},
    {"cacf", "name", "group"},
    {"manal", "name", "group"},
    {"nanal", "name", "group"},
    {"msnal", "name", "group"},
    {"fsnal", "name", "group"},
    {"canal", "name", "group"},
    {"sdl", "name", "group"},
    {"kls", "name", "group"},
    {"lc", "name", "group"},
    {"me", "name", "group"},
    {"wrk", "name", "group"},
    {"lmt", "name", "group"},
    {"heu", "name", "group"},
    {"stlb", "name", "group"},
    {"anl", "name", "group"},
    {"lime", "name", "group"},
    {"ru", "name", "group"},
    {"per", "name", "group"},
    {"fper", "name", "group"},
    {"wu", "name", "group"},
    {"mang", "name", "group"},
    {"cor", "name", "group"},
    {"mcor", "name", "group"},
    {"hem", "name", "group"},
    {"esk", "name", "group"},
    {"bix", "name", "group"},
    {"NiO", "name", "group"},
    {"pnt", "name", "group"},
    {"geik", "name", "group"},
    {"ilm", "name", "group"},
    {"bdy", "name", "group"},
    {"bdyT", "name", "group"},
    {"bdyC", "name", "group"},
    {"ten", "name", "group"},
    {"cup", "name", "group"},
    {"sp", "name", "group"},
    {"herc", "name", "group"},
    {"mt", "name", "group"},
    {"mft", "name", "group"},
    {"qnd", "name", "group"},
    {"usp", "name", "group"},
    {"picr", "name", "group"},
    {"br", "name", "group"},
    {"dsp", "name", "group"},
    {"gth", "name", "group"},
    {"cc", "name", "group"},
    {"arag", "name", "group"},
    {"mag", "name", "group"},
    {"sid", "name", "group"},
    {"rhc", "name", "group"},
    {"dol", "name", "group"},
    {"ank", "name", "group"},
    {"syv", "name", "group"},
    {"hlt", "name", "group"},
    {"pyr", "name", "group"},
    {"trot", "name", "group"},
    {"tro", "name", "group"},
    {"lot", "name", "group"},
    {"trov", "name", "group"},
    {"any", "name", "group"},
    {"iron", "name", "group"},
    {"Ni", "name", "group"},
    {"Cu", "name", "group"},
    {"gph", "name", "group"},
    {"diam", "name", "group"},
    {"S", "name", "group"},
    {"H2O", "name", "group"},
    {"CO2", "name", "group"},
    {"CO", "name", "group"},
    {"CH4", "name", "group"},
    {"O2", "name", "group"},
    {"H2", "name", "group"},
    {"S2", "name", "group"},
    {"H2S", "name", "group"},
    {"syvL", "name", "group"},
    {"hltL", "name", "group"},
    {"perL", "name", "group"},
    {"limL", "name", "group"},
    {"corL", "name", "group"},
    {"eskL", "name", "group"},
    {"hemL", "name", "group"},
    {"qL", "name", "group"},
    {"h2oL", "name", "group"},
    {"foL", "name", "group"},
    {"faL", "name", "group"},
    {"woL", "name", "group"},
    {"enL", "name", "group"},
    {"diL", "name", "group"},
    {"silL", "name", "group"},
    {"anL", "name", "group"},
    {"kspL", "name", "group"},
    {"abL", "name", "group"},
    {"neL", "name", "group"},
    {"lcL", "name", "group"},
    {"ruL", "name", "group"},
    {"bdyL", "name", "group"},
    {"H+", "name", "group"},
    {"Cl-", "name", "group"},
    {"OH-", "name", "group"},
    {"Na+", "name", "group"},
    {"K+", "name", "group"},
    {"Ca++", "name", "group"},
    {"Mg++", "name", "group"},
    {"Fe++", "name", "group"},
    {"Al+++", "name", "group"},
    {"CO3--", "name", "group"},
    {"AlOH3", "name", "group"},
    {"AlOH4-", "name", "group"},
    {"KOH", "name", "group"},
    {"HCl", "name", "group"},
    {"KCl", "name", "group"},
    {"NaCl", "name", "group"},
    {"CaCl2", "name", "group"},
    {"CaCl+", "name", "group"},
    {"MgCl2", "name", "group"},
    {"MgCl+", "name", "group"},
    {"FeCl2", "name", "group"},
    {"aqSi", "name", "group"},
    {"HS-", "name", "group"},
    {"HSO3-", "name", "group"},
    {"SO42-", "name", "group"},
    {"HSO4-", "name", "group"}
  };
  ds634dat={
    {"abb",   "name","group","formula"},
    {"fo",    "forsterite", "group", "norm", "Mg(2)Si(1)O(4)"},
    {"fa",    "fayalite", "group", "norm", "Fe(2)Si(1)O(4)"},
    {"teph",  "tephroite", "group", "norm", "Mn(2)Si(1)O(4)"},
    {"lrn",   "larnite", "group", "land", "Ca(2)Si(1)O(4)"},
    {"mont",  "monticellite", "group", "norm", "Ca(1)Mg(1)Si(1)O(4)"},
    {"chum",  "clinohumite", "group", "norm", "Mg(9)Si(4)H(2)O(18)"},
    {"chdr",  "name", "group", "norm", "Mg(5)Si(2)H(2)O(10)"},
    {"mwd",   "name", "group", "norm", "Mg(2)Si(1)O(4)"},
    {"fwd",   "name", "group", "norm", "Fe(2)Si(1)O(4)"},
    {"mrw",   "name", "group", "norm", "Mg(2)Si(1)O(4)"},
    {"frw",   "name", "group", "norm", "Fe(2)Si(1)O(4)"},
    {"mpv",   "name", "group", "norm", "Mg(1)Si(1)O(3)"},
    {"fpv",   "name", "group", "norm", "Fe(1)Si(1)O(3)"},
    {"apv",   "name", "group", "norm", "Al(2)O(3)"},
    {"npv",   "name", "group", "norm", "Na(0.5)Al(0.5)Si(1)O(3)"},
    {"ppv",   "name", "group", "norm", "Mg(1)Si(1)O(3)"},
    {"cpv",   "name", "group", "norm", "Ca(1)Si(1)O(3)"},
    {"mak",   "name", "group", "norm", "Mg(1)Si(1)O(3)"},
    {"fak",   "name", "group", "norm", "Fe(1)Si(1)O(3)"},
    {"maj",   "majorite", "group", "norm", "Mg(4)Si(4)O(12)"},
    {"nagt",  "name", "group", "norm", "Na(1)Mg(2)Al(1)Si(4)O(12)"},
    {"py",    "pyrope", "group", "norm", "Mg(3)Al(2)Si(3)O(12)"},
    {"alm",   "almandine", "group", "norm", "Fe(3)Al(2)Si(3)O(12)"},
    {"spss",  "spessartine", "group", "norm", "Mn(3)Al(2)Si(3)O(12)"},
    {"gr",    "grossular", "group", "norm", "Ca(3)Al(2)Si(3)O(12)"},
    {"andr",  "andradite", "group", "norm", "Ca(3)Fe(2)Si(3)O(12)"},
    {"ski",   "skiagite", "group", "norm", "Fe(5)Si(3)O(12)"},
    {"knor",  "knorringite", "group", "norm", "Mg(3)Cr(2)Si(3)O(12)"},
    {"uv",    "uvarovite", "group", "norm", "Ca(3)Cr(2)Si(3)O(12)"},
    {"osma",  "name", "group", "norm", "K(1)Mg(2)Al(5)Si(10)O(30)"},
    {"osmm",  "name", "group", "norm", "K(1)Mg(3)Al(3)Si(11)O(30)"},
    {"osfa",  "name", "group", "norm", "K(1)Fe(2)Al(5)Si(10)O(30)"},
    {"vsv",   "name", "group", "norm", "Ca(19)Mg(2)Al(11)Si(18)H(9)O(78)"},
    {"and",   "andalusite", "group", "norm", "Al(2)Si(1)O(5)"},
    {"ky",    "kyanite", "group", "norm", "Al(2)Si(1)O(5)"},
    {"sill",  "sillimanite", "group", "bw", "Al(2)Si(1)O(5)"},
    {"smul",  "name", "group", "norm", "Al(2)Si(1)O(5)"},
    {"amul",  "name", "group", "norm", "Al(2.5)Si(0.5)O(4.75)"},
    {"tpz",   "topaz", "group", "norm", "Al(2)Si(1)H(2)O(6)"},
    {"mst",   "Mg-staurolite", "group", "norm", "Mg(4)Al(18)Si(7.5)H(4)O(48)"},
    {"fst",   "Fe-staurolite", "group", "norm", "Fe(4)Al(18)Si(7.5)H(4)O(48)"},
    {"mnst",  "Mn-staurolite", "group", "norm", "Mn(4)Al(18)Si(7.5)H(4)O(48)"},
    {"mctd",  "Mg-chloritoid", "group", "norm", "Mg(1)Al(2)Si(1)H(2)O(7)"},
    {"fctd",  "Fe-chloritoid", "group", "norm", "Fe(1)Al(2)Si(1)H(2)O(7)"},
    {"mnctd", "Mn-chloritoid", "group", "norm", "Mn(1)Al(2)Si(1)H(2)O(7)"},
    {"merw",  "name", "group", "norm", "Ca(3)Mg(1)Si(2)O(8)"},
    {"spu",   "name", "group", "norm", "Ca(5)Si(2)C(1)O(11)"},
    {"zo",    "zoisite", "group", "norm", "Ca(2)Al(3)Si(3)H(1)O(13)"},
    {"cz",    "clinozoisite", "group", "norm", "Ca(2)Al(3)Si(3)H(1)O(13)"},
    {"ep",    "epidote", "group", "norm", "Ca(2)Fe(1)Al(2)Si(3)H(1)O(13)"},
    {"fep",   "name", "group", "norm", "Ca(2)Fe(2)Al(1)Si(3)H(1)O(13)"},
    {"pmt",   "name", "group", "norm", "Ca(2)Mn(1)Al(2)Si(3)H(1)O(13)"},
    {"law",   "lawsonite", "group", "norm", "Ca(1)Al(2)Si(2)H(4)O(10)"},
    {"mpm",   "name", "group", "norm", "Ca(4)Mg(1)Al(5)Si(6)H(7)O(28)"},
    {"fpm",   "name", "group", "norm", "Ca(4)Fe(1)Al(5)Si(6)H(7)O(28)"},
    {"jgd",   "name", "group", "norm", "Ca(4)Fe(6)Si(6)H(7)O(28)"},
    {"geh",   "name", "group", "bw", "Ca(2)Al(2)Si(1)O(7)"},
    {"ak",    "name", "group", "norm", "Ca(2)Mg(1)Si(2)O(7)"},
    {"rnk",   "name", "group", "norm", "Ca(3)Si(2)O(7)"},
    {"ty",    "name", "group", "norm", "Ca(5)Si(2)C(2)O(13)"},
    {"crd",   "Mg-cordierite", "group", "bw", "Mg(2)Al(4)Si(5)O(18)"},
    {"hcrd",  "Hyd-cordierite", "group", "bw", "Mg(2)Al(4)Si(5)H(2)O(19)"},
    {"fcrd",  "Fe-cordierite", "group", "bw", "Fe(2)Al(4)Si(5)O(18)"},
    {"mncrd", "Mn-cordierite", "group", "bw", "Mn(2)Al(4)Si(5)O(18)"},
    {"phA",   "name", "group", "norm", "Mg(7)Si(2)H(6)O(14)"},
    {"phD",   "name", "group", "norm", "Mg(1)Si(2)H(2)O(6)"},
    {"phE",   "name", "group", "norm", "Mg(2.4)Si(1.2)H(2.4)O(6)"},
    {"shB",   "name", "group", "norm", "Mg(10)Si(3)H(4)O(18)"},
    {"sph",   "sphene", "group", "land", "Ca(1)Ti(1)Si(1)O(5)"},
    {"cstn",  "name", "group", "norm", "Ca(1)Si(2)O(5)"},
    {"zrc",   "zircon", "group", "norm", "Zr(1)Si(1)O(4)"},
    {"zrt",   "name", "group", "norm", "Zr(1)Ti(1)O(4)"},
    {"tcn",   "name", "group", "norm", "Ti(1)Si(1)O(4)"},
    {"en",    "enstatite", "group", "norm", "Mg(2)Si(2)O(6)"},
    {"pren",  "name", "group", "norm", "Mg(2)Si(2)O(6)"},
    {"cen",   "clinoenstatite", "group", "norm", "Mg(2)Si(2)O(6)"},
    {"hen",   "name", "group", "norm", "Mg(2)Si(2)O(6)"},
    {"hfs",   "name", "group", "norm", "Fe(2)Si(2)O(6)"},
    {"fs",    "name", "group", "norm", "Fe(2)Si(2)O(6)"},
    {"mgts",  "name", "group", "norm", "Mg(1)Al(2)Si(1)O(6)"},
    {"di",    "name", "group", "norm", "Ca(1)Mg(1)Si(2)O(6)"},
    {"hed",   "name", "group", "norm", "Ca(1)Fe(1)Si(2)O(6)"},
    {"jd",    "name", "group", "norm", "Na(1)Al(1)Si(2)O(6)"},
    {"kjd",   "name", "group", "norm", "K(1)Al(1)Si(2)O(6)"},
    {"acm",   "name", "group", "norm", "Na(1)Fe(1)Si(2)O(6)"},
    {"kos",   "name", "group", "norm", "Na(1)Cr(1)Si(2)O(6)"},
    {"cats",  "name", "group", "bw", "Ca(1)Al(2)Si(1)O(6)"},
    {"caes",  "name", "group", "norm", "Ca(0.5)Al(1)Si(2)O(6)"},
    {"rhod",  "name", "group", "norm", "Mn(1)Si(1)O(3)"},
    {"pxmn",  "name", "group", "norm", "Mn(1)Si(1)O(3)"},
    {"wo",    "wollastonite", "group", "norm", "Ca(1)Si(1)O(3)"},
    {"pswo",  "name", "group", "norm", "Ca(1)Si(1)O(3)"},
    {"wal",   "name", "group", "norm", "Ca(1)Si(1)O(3)"},
    {"tr",    "tremolite", "group", "norm", "Ca(2)Mg(5)Si(8)H(2)O(24)"},
    {"fact",  "Fe-actinolite", "group", "norm", "Ca(2)Fe(5)Si(8)H(2)O(24)"},
    {"ts",    "tschermakite", "group", "norm", "Ca(2)Mg(3)Al(4)Si(6)H(2)O(24)"},
    {"parg",  "pargasite", "group", "norm", "Na(1)Ca(2)Mg(4)Al(3)Si(6)H(2)O(24)"},
    {"gl",    "glaucophane", "group", "norm", "Na(2)Mg(3)Al(2)Si(8)H(2)O(24)"},
    {"fgl",   "name", "group", "norm", "Na(2)Fe(3)Al(2)Si(8)H(2)O(24)"},
    {"nyb",   "name", "group", "norm", "Na(3)Mg(3)Al(3)Si(7)H(2)O(24)"},
    {"rieb",  "riebeckite", "group", "norm", "Na(2)Fe(5)Si(8)H(2)O(24)"},
    {"anth",  "anthophyllite", "group", "norm", "Mg(7)Si(8)H(2)O(24)"},
    {"fanth", "Fe-anthophyllite", "group", "norm", "Fe(7)Si(8)H(2)O(24)"},
    {"cumm",  "cummingtonite", "group", "norm", "Mg(7)Si(8)H(2)O(24)"},
    {"grun",  "grunerite", "group", "norm", "Fe(7)Si(8)H(2)O(24)"},
    {"ged",   "gedrite", "group", "norm", "Mg(5)Al(4)Si(6)H(2)O(24)"},
    {"spr4",  "name", "group", "norm", "Mg(4)Al(8)Si(2)O(20)"},
    {"spr5",  "name", "group", "norm", "Mg(3)Al(10)Si(1)O(20)"},
    {"fspr",  "name", "group", "norm", "Fe(4)Al(8)Si(2)O(20)"},
    {"mcar",  "name", "group", "norm", "Mg(1)Al(2)Si(2)H(4)O(10)"},
    {"fcar",  "name", "group", "norm", "Fe(1)Al(2)Si(2)H(4)O(10)"},
    {"deer",  "name", "group", "norm", "Fe(18)Si(12)H(10)O(50)"},
    {"mu",    "muscovite", "group", "norm", "K(1)Al(3)Si(3)H(2)O(12)"},
    {"cel",   "celadonite", "group", "norm", "K(1)Mg(1)Al(1)Si(4)H(2)O(12)"},
    {"fcel",  "Fe-celadonite", "group", "norm", "K(1)Fe(1)Al(1)Si(4)H(2)O(12)"},
    {"pa",    "paragonite", "group", "norm", "Na(1)Al(3)Si(3)H(2)O(12)"},
    {"ma",    "margarite", "group", "norm", "Ca(1)Al(4)Si(2)H(2)O(12)"},
    {"phl",   "phlogopite", "group", "norm", "K(1)Mg(3)Al(1)Si(3)H(2)O(12)"},
    {"ann",   "annite", "group", "norm", "K(1)Fe(3)Al(1)Si(3)H(2)O(12)"},
    {"fbi",   "ferri-biotite", "group", "norm", "K(1)Fe(1)Mg(2)Al(2)Si(2)H(2)O(12)"},
    {"obi",   "ord-biotite", "group", "norm", "K(1)Fe(1)Mg(2)Al(1)Si(3)H(2)O(12)"},
    {"mnbi",  "Mn-biotite", "group", "norm", "K(1)Mn(3)Al(1)Si(3)H(2)O(12)"},
    {"east",  "eastonite", "group", "norm", "K(1)Mg(2)Al(3)Si(2)H(2)O(12)"},
    {"naph",  "Na-phlogopite", "group", "norm", "Na(1)Mg(3)Al(1)Si(3)H(2)O(12)"},
    {"tan",   "name", "group", "norm", "Mg(3)Si(4)H(2)O(12)"},
    {"clin",  "name", "group", "norm", "Mg(5)Al(2)Si(3)H(8)O(18)"},
    {"ames",  "name", "group", "norm", "Mg(4)Al(4)Si(2)H(8)O(18)"},
    {"afchl", "name", "group", "norm", "Mg(6)Si(4)H(8)O(18)"},
    {"daph",  "name", "group", "norm", "Fe(5)Al(2)Si(3)H(8)O(18)"},
    {"mnchl", "name", "group", "norm", "Mn(5)Al(2)Si(3)H(8)O(18)"},
    {"sud",   "name", "group", "norm", "Mg(2)Al(4)Si(3)H(8)O(18)"},
    {"fsud",  "name", "group", "norm", "Fe(2)Al(4)Si(3)H(8)O(18)"},
    {"prl",   "pyrophyllite", "group", "norm", "Al(2)Si(4)H(2)O(12)"},
    {"ta",    "name", "group", "norm", "Mg(3)Si(4)H(2)O(12)"},
    {"fta",   "name", "group", "norm", "Fe(3)Si(4)H(2)O(12)"},
    {"tats",  "name", "group", "norm", "Mg(2)Al(2)Si(3)H(2)O(12)"},
    {"tap",   "name", "group", "norm", "Al(2)Si(4)H(2)O(12)"},
    {"nta",   "name", "group", "norm", "Na(1)Mg(3)Al(1)Si(3)H(2)O(12)"},
    {"minn",  "name", "group", "norm", "Fe(3)Si(4)H(2)O(12)"},
    {"minm",  "name", "group", "norm", "Mg(3)Si(4)H(2)O(12)"},
    {"kao",   "kaolitite", "group", "norm", "Al(2)Si(2)H(4)O(9)"},
    {"pre",   "prehnite", "group", "norm", "Ca(2)Al(2)Si(3)H(2)O(12)"},
    {"fpre",  "ferri-prehnite", "group", "norm", "Ca(2)Fe(1)Al(1)Si(3)H(2)O(12)"},
    {"chr",   "name", "group", "norm", "Mg(3)Si(2)H(4)O(9)"},
    {"liz",   "name", "group", "norm", "Mg(3)Si(2)H(4)O(9)"},
    {"glt",   "name", "group", "norm", "Fe(3)Si(2)H(4)O(9)"},
    {"fstp",  "Fe-stilpnomelane", "group", "norm", "K(0.5)Fe(5)Al(2)Si(8)H(12.5)O(30.5)"},
    {"mstp",  "Mg-stilpnomelane", "group", "norm", "K(0.5)Mg(5)Al(2)Si(8)H(12.5)O(30.5)"},
    {"atg",   "name", "group", "norm", "Mg(48)Si(34)H(62)O(147)"},
    {"ab",    "albite", "group", "bw", "Na(1)Al(1)Si(3)O(8)"},
    {"abh",   "albite-high", "group", "norm", "Na(1)Al(1)Si(3)O(8)"},
    {"mic",   "microcline", "group", "norm", "K(1)Al(1)Si(3)O(8)"},
    {"san",   "sanidine", "group", "bw", "K(1)Al(1)Si(3)O(8)"},
    {"an",    "anorthite", "group", "bw", "Ca(1)Al(2)Si(2)O(8)"},
    {"kcm",   "name", "group", "norm", "K(1)Al(1)Si(3)H(2)O(9)"},
    {"wa",    "name", "group", "norm", "K(2)Si(4)O(9)"},
    {"hol",   "name", "group", "norm", "K(1)Al(1)Si(3)O(8)"},
    {"q",     "quartz", "group", "land", "Si(1)O(2)"},
    {"trd",   "tridymite", "group", "norm", "Si(1)O(2)"},
    {"crst",  "cristobalite", "group", "norm", "Si(1)O(2)"},
    {"coe",   "coesite", "group", "norm", "Si(1)O(2)"},
    {"stv",   "stishovite", "group", "norm", "Si(1)O(2)"},
    {"ne",    "nepheline", "group", "land", "Na(1)Al(1)Si(1)O(4)"},
    {"cg",    "name", "group", "norm", "Na(1)Al(1)Si(1)O(4)"},
    {"cgh",   "name", "group", "norm", "Na(1)Al(1)Si(1)O(4)"},
    {"macf",  "name", "group", "norm", "Mg(1)Al(2)O(4)"},
    {"mscf",  "name", "group", "norm", "Mg(2)Si(1)O(4)"},
    {"fscf",  "name", "group", "norm", "Fe(2)Si(1)O(4)"},
    {"nacf",  "name", "group", "norm", "Na(1)Al(1)Si(1)O(4)"},
    {"cacf",  "name", "group", "norm", "Ca(1)Al(2)O(4)"},
    {"manal", "name", "group", "norm", "Mg(3)Al(6)O(12)"},
    {"nanal", "name", "group", "norm", "Na(1)Mg(2)Al(5)Si(1)O(12)"},
    {"msnal", "name", "group", "norm", "Mg(6)Si(3)O(12)"},
    {"fsnal", "name", "group", "norm", "Fe(6)Si(3)O(12)"},
    {"canal", "name", "group", "norm", "Ca(1)Mg(2)Al(6)O(12)"},
    {"sdl",   "name", "group", "norm", "Na(8)Al(6)Si(6)Cl(2)O(24)"},
    {"kls",   "name", "group", "norm", "K(1)Al(1)Si(1)O(4)"},
    {"lc",    "name", "group", "bw", "K(1)Al(1)Si(2)O(6)"},
    {"me",    "name", "group", "norm", "Ca(4)Al(6)Si(6)C(1)O(27)"},
    {"wrk",   "name", "group", "norm", "Ca(1)Al(2)Si(4)H(4)O(14)"},
    {"lmt",   "name", "group", "norm", "Ca(1)Al(2)Si(4)H(8)O(16)"},
    {"heu",   "name", "group", "norm", "Ca(1)Al(2)Si(7)H(12)O(24)"},
    {"stlb",  "name", "group", "norm", "Ca(1)Al(2)Si(7)H(14)O(25)"},
    {"anl",   "name", "group", "norm", "Na(1)Al(1)Si(2)H(2)O(7)"},
    {"lime",  "name", "group", "norm", "Ca(1)O(1)"},
    {"ru",    "rutile", "group", "norm", "Ti(1)O(2)"},
    {"per",   "periclase", "group", "norm", "Mg(1)O(1)"},
    {"fper",  "Fe-periclase", "group", "norm", "Fe(1)O(1)"},
    {"wu",    "name", "group", "norm", "Fe(1)O(1)"},
    {"mang",  "name", "group", "norm", "Mn(1)O(1)"},
    {"cor",   "corundum", "group", "norm", "Al(2)O(3)"},
    {"mcor",  "Mg-corundum", "group", "norm", "Mg(1)Si(1)O(3)"},
    {"hem",   "hematite", "group", "land", "Fe(2)O(3)"},
    {"esk",   "name", "group", "norm", "Cr(2)O(3)"},
    {"bix",   "name", "group", "norm", "Mn(2)O(3)"},
    {"NiO",   "nickel-oxide", "group", "land", "Ni(1)O(1)"},
    {"pnt",   "name", "group", "norm", "Mn(1)Ti(1)O(3)"},
    {"geik",  "name", "group", "bw", "Mg(1)Ti(1)O(3)"},
    {"ilm",   "ilmenite", "group", "bw", "Fe(1)Ti(1)O(3)"},
    {"bdy",   "name", "group", "norm", "Zr(1)O(2)"},
    {"bdyT",  "name", "group", "norm", "Zr(1)O(2)"},
    {"bdyC",  "name", "group", "norm", "Zr(1)O(2)"},
    {"ten",   "name", "group", "norm", "Cu(1)O(1)"},
    {"cup",   "cuprite", "group", "norm", "Cu(2)O(1)"},
    {"sp",    "spinel", "group", "bw", "Mg(1)Al(2)O(4)"},
    {"herc",  "hercynite", "group", "bw", "Fe(1)Al(2)O(4)"},
    {"mt",    "magnetite", "group", "land", "Fe(3)O(4)"},
    {"mft",   "name", "group", "land", "Fe(2)Mg(1)O(4)"},
    {"qnd",   "name", "group", "norm", "Mg(2)Ti(1)O(4)"},
    {"usp",   "ulvospinel", "group", "bw", "Fe(2)Ti(1)O(4)"},
    {"picr",  "name", "group", "norm", "Mg(1)Cr(2)O(4)"},
    {"br",    "name", "group", "norm", "Mg(1)H(2)O(2)"},
    {"dsp",   "name", "group", "norm", "Al(1)H(1)O(2)"},
    {"gth",   "name", "group", "norm", "Fe(1)H(1)O(2)"},
    {"cc",    "calcite", "group", "land", "Ca(1)C(1)O(3)"},
    {"arag",  "aragonite", "group", "norm", "Ca(1)C(1)O(3)"},
    {"mag",   "magnesite", "group", "norm", "Mg(1)C(1)O(3)"},
    {"sid",   "siderite", "group", "norm", "Fe(1)C(1)O(3)"},
    {"rhc",   "rhodochrosite", "group", "norm", "Mn(1)C(1)O(3)"},
    {"dol",   "dolomite", "group", "bw", "Ca(1)Mg(1)C(2)O(6)"},
    {"ank",   "ankerite", "group", "bw", "Ca(1)Fe(1)C(2)O(6)"},
    {"syv",   "sylvite", "group", "norm", "K(1)Cl(1)"},
    {"hlt",   "halite", "group", "norm", "Na(1)Cl(1)"},
    {"pyr",   "pyrite", "group", "norm", "Fe(1)S(2)"},
    {"trot",  "name", "group", "land", "Fe(1)S(1)"},
    {"tro",   "name", "group", "land", "Fe(1)S(1)"},
    {"lot",   "name", "group", "land", "Fe(1)S(1)"},
    {"trov",  "name", "group", "land", "Fe(0.875)S(1)"},
    {"any",   "anhydrite", "group", "norm", "Ca(1)S(1)O(4)"},
    {"iron",  "iron", "group", "land", "Fe(1)"},
    {"Ni",    "nickel", "group", "land", "Ni(1)"},
    {"Cu",    "copper", "group", "norm", "Cu(1)"},
    {"gph",   "graphite", "group", "norm", "C(1)"},
    {"diam",  "diamond", "group", "norm", "C(1)"},
    {"S",     "sulfur", "group", "norm", "S(1)"},
    {"H2O",   "h2o", "group", "h2o", "H(2)O(1)"},
    {"CO2",   "carbon-dioxide", "group", "co2", "C(1)O(2)"},
    {"CO",    "carbon-monoxide", "group", "cs", "C(1)O(1)"},
    {"CH4",   "methane", "group", "cs", "C(1)H(4)"},
    {"O2",    "name", "group", "cs", "O(2)"},
    {"H2",    "name", "group", "cs", "H(2)"},
    {"S2",    "name", "group", "cs", "S(2)"},
    {"H2S",   "name", "group", "cs", "S(1)H(2)"},
    {"syvL",  "liq-sylvite", "group", "liq", "K(1)Cl(1)"},
    {"hltL",  "liq-halite", "group", "liq", "Na(1)Cl(1)"},
    {"perL",  "liq-periclase", "group", "liq", "Mg(1)O(1)"},
    {"limL",  "liq-lime", "group", "liq", "Ca(1)O(1)"},
    {"corL",  "liq-corundum", "group", "liq", "Al(2)O(3)"},
    {"eskL",  "liq-eskolaite", "group", "liq", "Cr(2)O(3)"},
    {"hemL",  "liq-hematite", "group", "liq", "Fe(2)O(3)"},
    {"qL",    "liq-quartz", "group", "liq", "Si(1)O(2)"},
    {"h2oL",  "liq-h2o", "group", "liq", "H(2)O(1)"},
    {"foL",   "liq-forsterite", "group", "liq", "Mg(2)Si(1)O(4)"},
    {"faL",   "liq-fayalite", "group", "liq", "Fe(2)Si(1)O(4)"},
    {"woL",   "liq-wollastonite", "group", "liq", "Ca(1)Si(1)O(3)"},
    {"enL",   "liq-enstatite", "group", "liq", "Mg(2)Si(2)O(6)"},
    {"diL",   "liq-diopside", "group", "liq", "Ca(1)Mg(1)Si(2)O(6)"},
    {"silL",  "liq-sillimanite", "group", "liq", "Al(2)Si(1)O(5)"},
    {"anL",   "liq-anorthite", "group", "liq", "Ca(1)Al(2)Si(2)O(8)"},
    {"kspL",  "liq-k-feldspar", "group", "liq", "K(1)Al(1)Si(3)O(8)"},
    {"abL",   "liq-albite", "group", "liq", "Na(1)Al(1)Si(3)O(8)"},
    {"neL",   "liq-nepheline", "group", "liq", "Na(1)Al(1)Si(1)O(4)"},
    {"lcL",   "liq-leucite", "group", "liq", "K(1)Al(1)Si(2)O(6)"},
    {"ruL",   "liq-rutile", "group", "liq", "Ti(1)O(2)"},
    {"bdyL",  "liq-badellyite", "group", "liq", "Zr(1)O(2)"},
    {"H+",    "name", "group", "aq", "H(1)el(-1)"},
    {"Cl-",   "name", "group", "aq", "Cl(1)el(1)"},
    {"OH-",   "name", "group", "aq", "H(1)O(1)el(1)"},
    {"Na+",   "name", "group", "aq", "Na(1)el(-1)"},
    {"K+",    "name", "group", "aq", "K(1)el(-1)"},
    {"Ca++",  "name", "group", "aq", "Ca(1)el(-2)"},
    {"Mg++",  "name", "group", "aq", "Mg(1)el(-2)"},
    {"Fe++",  "name", "group", "aq", "Fe(1)el(-2)"},
    {"Al+++", "name", "group", "aq", "Al(1)el(-3)"},
    {"CO3--", "name", "group", "aq", "C(1)O(3)el(2)"},
    {"AlOH3", "name", "group", "aq", "Al(1)H(3)O(3)"},
    {"AlOH4-","name", "group", "aq", "Al(1)H(4)O(4)el(1)"},
    {"KOH",   "name", "group", "aq", "K(1)H(1)O(1)"},
    {"HCl",   "name", "group", "aq", "Cl(1)H(1)"},
    {"KCl",   "name", "group", "aq", "K(1)Cl(1)"},
    {"NaCl",  "name", "group", "aq", "Na(1)Cl(1)"},
    {"CaCl2", "name", "group", "aq", "Ca(1)Cl(2)"},
    {"CaCl+", "name", "group", "aq", "Ca(1)Cl(1)el(-1)"},
    {"MgCl2", "name", "group", "aq", "Mg(1)Cl(2)"},
    {"MgCl+", "name", "group", "aq", "Mg(1)Cl(1)el(-1)"},
    {"FeCl2", "name", "group", "aq", "Fe(1)Cl(2)"},
    {"aqSi",  "name", "group", "aq", "Si(1)O(2)"},
    {"HS-",   "name", "group", "aq", "S(1)H(1)el(1)"},
    {"HSO3-", "name", "group", "aq", "S(1)H(1)O(3)el(1)"},
    {"SO42-", "name", "group", "aq", "S(1)O(4)el(2)"},
    {"HSO4-", "name", "group", "aq", "S(1)H(1)O(4)el(1)"}
  };  
  Return[ds634dat];
]


(* ::Subsubsection::Closed:: *)
(*MatchPatt*)


(*Don't think this is used anywhere. Move to utils*)
RationalToString[rat_Rational]:=ToString[rat,InputForm];


MatchPatt[arg_?StringQ]:=Module[
  {pmstr,multexppatt,rationalpatt,ftnformpatt,realpatt,
   intpatt,flagpatt},
  multexppatt=NumberString~~"*10^"~~{"+","-",DigitCharacter..};
  rationalpatt=NumberString~~"/"~~NumberString;
  ftnformpatt=NumberString~~{"d","e","D","E"}~~{"+","-",DigitCharacter..};
  realpatt={{"-","+"}|DigitCharacter}..~~"."~~{EndOfString,DigitCharacter}..;
  intpatt={"+","-",StartOfString}~~DigitCharacter..;
  flagpatt="%<"~~Except[">"]..~~">";
  pmstr="None";
  
  Which[
    StringMatchQ[arg,multexppatt],
      (*Print["Matches multexppatt"];*)
      pmstr="MultExpPatt",
    StringMatchQ[arg,ftnformpatt],
      (*Print["Matches ftnformpatt"];*)
      pmstr="ScientificPatt",
    StringMatchQ[arg,rationalpatt],
      (*Print["Matches rationalpatt"];*)
      pmstr="RationalPatt",
    StringMatchQ[arg,realpatt],
      (*Print["Matches realpatt"];*)
      pmstr="RealPatt",
    StringMatchQ[arg,intpatt],
      (*Print["Matches intpatt"];*)
      pmstr="IntegerPatt",
    StringMatchQ[arg,flagpatt],
      (*Print["Matches flagpatt"];*)
      pmstr="FlagPatt",
    True,
      (*Print["Doesn't match my patterns, so a String."];*)
      pmstr="String",
    _,
      pmstr="ERROR"
      (*Print["ERROR"]*)
  ];

  Return[pmstr];
];


(* ::Subsubsection::Closed:: *)
(*GetEx, GenExListFromString, GenExListFromStringList*)


GetEx[arg_?StringQ, opts:OptionsPattern[]]:=Module[{cfnum,
  doarith,mt,isconv,res},
  cfnum=OptionValue[GetEx,"ConvertFortranNum"];
  doarith=OptionValue[GetEx,"DoArithmetic"];
  res=StringTrim[arg];
  isconv=False;
  (*Print[res];*)
  
  mt=MatchPatt[arg];
  
  Which[
    StringMatchQ[mt,"MultExpPatt"]||StringMatchQ[mt,"RealPatt"]
      ||StringMatchQ[mt,"RationalPatt"]||StringMatchQ[mt,"IntegerPatt"],
      res=ToEx[res],
    StringMatchQ[mt,"ScientificPatt"],
      res=System`Convert`TableDump`ParseTable[{{res}},{{{},{}},
        {"-","+"},"."},False][[1,1]],
    StringMatchQ[arg,"FlagPatt"],
      res=StringDrop[res,2];
      res=StringDrop[res,-1],
    True,
     (* nothing to do *)
     ,
    _,
      Print["ERROR in GetEx trying to classify "<>res];
      Abort[];
  ];
  
  Return[res];
];



GenExListFromString[str_, opts:OptionsPattern[]]:=Module[{cfnum,
  doarith,spl,res},
  cfnum=OptionValue[GenExListFromString,"ConvertFortranNum"];
  doarith=OptionValue[GenExListFromString,"DoArithmetic"];
  spl=StringSplit[str];
  res=spl;
  Do[res[[i]]=GetEx[spl[[i]],"ConvertFortranNum"->cfnum]
    ,{i,Length[spl]}
  ];
  Return[res];
];



GenExListFromStringList[slis_, opts:OptionsPattern[]]:=
  Module[{cfnum,doarith,res},
  cfnum=OptionValue[GenExListFromStringList,"ConvertFortranNum"];
  doarith=OptionValue[GenExListFromStringList,"DoArithmetic"];
  res=slis;
  Do[res[[i]]=GetEx[slis[[i]],"ConvertFortranNum"->cfnum]
    ,{i,Length[slis]}
  ];
  Return[res];
];


(* ::Subsubsection::Closed:: *)
(*GetCleanTCAxLines*)


GetCleanTCAxLines[tcfile_?StringQ]:=Module[
  {dat,ndat,blpos,hdrpos,vrbpos,hpos,vpos,verbatims,rvpos,clndat,
  header,tstr},
  dat=GenListOfLinesFromFile[tcfile];
  ndat=Flatten[dat];
  (*Remove blank lines*)
  blpos={};
  Do[
    If[StringLength[StringTrim[ndat[[i]]]]==0,AppendTo[blpos,i]];
    ,{i,1,Length[ndat]}
  ];
  blpos=Reverse[blpos];
  Do[ndat=Drop[ndat,{blpos[[i]]}],{i,Length[blpos]}];
  (*Find header and verbatim positions*)
  hdrpos={}; vrbpos={};
  Do[
    If[StringStartsQ[ndat[[i]],"header"],AppendTo[hdrpos,i]];
    If[StringStartsQ[ndat[[i]],"verbatim"],AppendTo[vrbpos,i]];
    ,{i,1,Length[ndat]}
  ];
  (*{hpos,vpos}={hdrpos,Partition[vrbpos,2]};*)
  vpos=Partition[vrbpos,2];
  (*Store header*)
  header=Take[ndat,hdrpos];
  (*Store verbatims*)
  verbatims={};
  Do[AppendTo[verbatims,Take[ndat,vpos[[i]]]],{i,Length[vpos]}];
  rvpos=Reverse[vpos];
  (*Remove header and verbatim sections from ndat, in reverse*)
  Do[ndat=Drop[ndat,rvpos[[i]]],{i,Length[rvpos]}];
  ndat=Drop[ndat,hdrpos];
  (*Remove Comment Lines and Trailing Comments, put in clndat*)
  clndat={};
  Block[{cont},
    Do[
      tstr=StringTrim[ndat[[i]]];
    
      (* 2022-04-05: need to check for * as first char of line in case
         storing text below the star. *)
      If[StringStartsQ[tstr,"*"], AppendTo[clndat,"*"]; Break[]; ];
    
      (* Grab special lines I insert that are needed to process. For now site substitutions *)
      (*If[StringStartsQ[tstr,"%<sub>"~~"%<flag>"~~"%<name>"~~"%<note>"],AppendTo[clndat,tstr]; Continue[]];*)
      cont=False;
      Map[
        If[StringStartsQ[tstr,#],AppendTo[clndat,tstr];cont=True;]&
        ,{"%<sub>","%<flag>","%<name>","%<ref>","%<note>","%<sm>","%<esm>"}
      ];
      If[TrueQ[cont],Continue[]];
    
    
      If[StringStartsQ[tstr,"%"],Continue[],
        If[StringContainsQ[tstr,"%"],tstr=StringTake[tstr,{1,StringPosition[tstr,"%",1][[1,1]]-1}]];
        AppendTo[clndat,tstr]
      ];
      ,{i,Length[ndat]}
    ];
  ];
  Print[clndat];
  Return[{clndat,header,verbatims,dat}];
];


(* ::Subsubsection::Closed:: *)
(*ParseCleanTCAxLines*)


(*
ToDo=Done: use GetSIPs to convert sitemult and effsitemult to nested numeric lists:
sitemult2=Map[GenNumListFromStringList[GetSIPs[#]]&,sitemult];
esitemult2=Map[GenNumListFromStringList[GetSIPs[#]]&,esitemult];
mult=MapThread[#1[[2]]/#2[[2]]&,{sitemult2,esitemult2}]  If any el not==1, then has eff site mult.
*)

ParseCleanTCAxLines[lns_List,dsv_String]:=Module[
  {dat,info,EOD=False,CLLS,nlns,i,j,k,l,m,n,ln,cln,ph,npc,ver,tls2,tstr,pcnames={},
  vars,guesses,ppnls,xcstype,wlist,alist,nsitex,sitexrep,dositexrep,sitemult,
  effsitemult,sitexls,mklns,aidlns,checklns,makelns,dqflns,MSTRAXDATA={},
  purephases={},GetCLLS,IsAXModelStart,StandardizeMakeLine,skipax},
  (* 
     CLLS = current line, as list 
     EOD = end of data flag indicating on last line read with
           last line of * read and clls set to {*} 
     lns = clean lines; checked on entry to add * as last line
     ln  = index indicating last line # read, and line put in CLLS
  *)
  (*lns passed in is clean with no blank lines. Each line is a single string in top level list*)
  (*
    This puts each string line into a list, and strings in number form are turned into numbers.
    Is there any case where we need string rationals left as a string?
  *)
  dat=Map[GenExListFromString[#,"ConvertFortranNum"->False]&,lns];
  (*Print["ok 1.1"];*)
  
  nlns=Length[dat];
  If[StringMatchQ[dat[[nlns]][[1]],"\\*"],
    dat[[nlns]]={"ENDOFFILE"};
    ,
    AppendTo[dat,{"ENDOFFILE"}];
    nlns=Length[dat];
  ];
  (*Print["ok 1.2"];*)
  
  
  (* Function to grab next line ln if it is available. *)
  GetCLLS[]:=If[ln<nlns-1,
      CLLS=dat[[++ln]];
      Return[CLLS];
      ,
      (*Print[{ln,nlns}];*)
      EOD=True;
      CLLS={"ENDOFFILE"};
      Return[CLLS];
    ];
  IsAXModelStart[ls_List]:=
    If[Length[ls]==3,
      If[StringQ[ls[[1]]] && IntegerQ[ls[[2]]] && IntegerQ[ls[[3]]],
        Return[True]
        ,
        Return[False]
      ]
      ,
      Return[False]
    ];
   StandardizeMakeLine[ln_?ListQ]:=Block[{idx=1,nc,nln={}},
     (*{"make",1,"ordered","cats",1}*)
     (*Print[{"make line on entry: ",ln}];*)
     nc=ln[[2]];
     AppendTo[nln,ln[[1]]];
     AppendTo[nln,ln[[2]]];
     idx=3;
     Do[
       If[StringMatchQ[ln[[idx]],{"ordered","disordered","equilibrium"}],
         AppendTo[nln,ln[[idx++]]];(*flag*)
         AppendTo[nln,ln[[idx++]]];(*pc*)
         AppendTo[nln,ln[[idx++]]];(*rxn coeff*)
         ,
         AppendTo[nln,"norm"];(*flag*)
         If[!StringMatchQ[tdtype[ln[[idx]],dsv],"norm"]
            &&!StringMatchQ[tdtype[ln[[idx]],dsv],"liq"],
           Print["In make line, a non-normal and non-liq species does not have a flag: "<>ln[[idx]]];
           Abort[];
         ];
         AppendTo[nln,ln[[idx++]]];(*pc*)
         AppendTo[nln,ln[[idx++]]];(*rxn coeff*)
         ];
       ,{i,1,nc}
     ];
     (*Print[{"make line on exit: ",nln}];*)
     Return[nln];
   ];
    
  (* Start Processing *)
  ln=0;
  CLLS=GetCLLS[];
  If[IsAXModelStart[CLLS],
    {ph,npc,ver}=CLLS;
    Print[{ph,npc,ver}];
    ,
    Print["First line is not start of an ax model; Aborting."];
    Abort[];
  ];
  
 
  While[ (StringQ[ph]&&IntegerQ[npc]&&ver==1),
    skipax=False;
	(* Check first if next line is any of the %<****> special lines *)
	(* and process them. When done, backup line by 1 to continue. *)
    info={{"skip",False},{"name",ph},{"ref",""},{"note",""},{"ext",False}};
    Print["Parsing "<>ph];
    CLLS=GetCLLS[];
    (*Print[CLLS];*)
    While[StringStartsQ[CLLS[[1]],"%<"],
      Which[
        StringContainsQ[CLLS[[1]],"<flag>"],
          Do[
            If[StringMatchQ[CLLS[[2]],"SKIP",IgnoreCase->True],
              info[[1,2]]=True;
              skipax=True;
            ];
            If[StringMatchQ[CLLS[[2]],"EXT",IgnoreCase->True],
              info[[5,2]]=True;
              (*skipax=True;*)
            ];
            (* add check for other flag types here *)
            ,{i,2,Length[CLLS]}
          ],
        StringContainsQ[CLLS[[1]],"<name>"],
          info[[2,2]]=CLLS[[2]];
          (*can't change ph yet b/c it is used as check in var names in SIP*)
          (*ph=CLLS[[2]]*)
          ,
        StringContainsQ[CLLS[[1]],"<ref>"],
          
          tls2={};
          Do[AppendTo[tls2,CLLS[[i]]],{i,2,Length[CLLS]}];
          (*Print["CLLS tls2 = ",tls2//InputForm];*)
          info[[3,2]]=tls2;
          tls2={};
          ,
        StringContainsQ[CLLS[[1]],"<note>"],
          (* combine remaining vals to 1 string, append to note *)
          tstr="";
          If[StringLength[info[[4,2]]]>0,
            tstr=" \n"<>"!          ",
            tstr=dsv<>".  "
          ];
          Do[
            tstr=tstr<>ToString[CLLS[[i]]]<>" ";
            ,{i,2,Length[CLLS]}
          ];
          (* join, but don't string trim, to leave space if another line *)
          info[[4,2]]=info[[4,2]]<>tstr;
          ,
        True,
          Print["ERROR in ParseCleanTCAxLines in info. Unknown flag."];
          Abort[];
      ];
      CLLS=GetCLLS[];
      If[ !(Head[CLLS[[1]]]===String),Break[]];
    ];
    ln=ln-1;
    
    
    
    vars={};
    guesses={};
    Block[{var,val},
      (*Get var and guesses*)
      (*Print["   Parsing: Gettings vars and guesses: "<>ph];*)
      Do[
        CLLS=GetCLLS[];
        {var,val}={CLLS[[1]],CLLS[[2]]};
        (*Make str in parens match ph name, unless str is L, which is
        var they seem to always use for a liquid regardless of ph name*)
        If[GetSIP[var]!=ph && GetSIP[var]!="L",Print["SIP != ph WARNING in variable name: "<>var]];
        var=StringTake[var,StringPosition[var,"(",1][[1,1]]-1];
        AppendTo[vars,var];
        AppendTo[guesses,val];
        
        (*ToDo: Check for range, isQ*)
        ,{i,npc-1}
      ];
      (*Print["ok 2.1"];*)
      (*Get ppn expressions*)
      ppnls={};
      pcnames={};
      i=1;
      (*Need to read-ahead of while loop and at end of each If case in 
        order to catch multi-line ppn spec for last pc. This means last 
        read will contain xcstype, so don't need to read xcstyp outside
        of while loop *)
      (*Print["   Parsing: Gettings ppns: "<>ph];*)
      CLLS=GetCLLS[];
      While[i<=npc||IntegerQ[CLLS[[1]]],
        If[IntegerQ[CLLS[[1]]],
          (*Add tls2 contents onto tls contents, don't increment i*)
          ppnls[[i-1]]=Join[ppnls[[i-1]],CLLS];
          CLLS=GetCLLS[];
          ,
          AppendTo[ppnls,CLLS];
          AppendTo[pcnames,GetSIP[CLLS[[1]]]];
          ++i;
          CLLS=GetCLLS[];
        ];
      ];
    ];(*end Block*)
    (*phase name, npc, vars & vals, pcnames, and ppnls are now set. pc  
      names in pcnames & ppnls structures may be replaced later after 
      makes processed. *)
      
    (*get excess type*)
    (* xcstype read in last while loop above, so grab it, don't
       need to read again to get it. *)
    (*Print["ok 2.2"];*)
    xcstype=CLLS[[1]];
    (*Print["ok 2.3"];*)
    
    (*Construct wlist*)
    wlist={};
    (*Print["   Parsing: Checking for W's: "<>ph];*)
    If[StringMatchQ[xcstype,"sf"]||StringMatchQ[xcstype,"asf"],
      Do[
        (*Print["i="<>ToString[i]];*)
        CLLS=GetCLLS[];
        (*Print[{"tls",tls}];*)
        If[Length[CLLS] != 4,
          Print[{"ERROR: w line does not have length of 4 for phase "<>ph<>": ",CLLS}];
          Abort[];
        ];
        Assert[Length[CLLS]==4];
        tls2=GetSIPs[CLLS[[1]]];
        (*Print[{"tls2",tls2}];*)
        (*Print[tls2];*)
        Assert[Length[tls2]==2];
        AppendTo[wlist,Join[tls2,CLLS[[2;;4]]]];
        
        ,{i,1,1/2*npc*(npc-1)}
      ],
      If[!StringMatchQ[xcstype,"ideal"],
        Print["xcstype StringMatchQ failed"];
      ];
      wlist={};
    ];
    alist={};
    (*Print["   Parsing: Checking for a's: "<>ph];*)
    If[StringMatchQ[xcstype,"asf"],
      Do[
        CLLS=GetCLLS[];
        Assert[Length[CLLS]==4];
        (*Assert[StringStartsQ[CLLS[[1]],"a("]];*)
        
        (*]a terms not required to be in format a(pc), as shown by pl4TR, so
        can't use GetSIP without checking for ( and ) first. *)
        tstr=CLLS[[1]];
        If[StringContainsQ[tstr,"("] 
           && StringContainsQ[tstr,")"],
          tstr=GetSIP[CLLS[[1]]];
          Assert[StringQ[tstr]]
          ,
          tstr=CLLS[[1]]
        ];
        (* Check tstr is a valid pcname *)
        If[ !StringMatchQ[tstr,pcnames], 
          Print["a term string does not match a pc name: "<>tstr]; 
          Abort[];
        ];
        AppendTo[alist,Join[{tstr},CLLS[[2;;4]]]];
        ,{i,1,npc}
      ];
    ];
    nsitex=ToEx[lns[[++ln]]];
    Assert[IntegerQ[nsitex]];
    
    (* Check for %<sub> lines and append into a sub list. *)
    (*
    sitexrep={};
    dositexrep=False;
    While[StringMatchQ[GetCLLS[][[1]],"%<sub>"],
      CLLS=Delete[CLLS,1];
      dositexrep=True;
      AppendTo[sitexrep,CLLS];
    ];
    sitexrep=Flatten[sitexrep];
    ln=ln-1;
    *)
    
    (*TESTING NEW FLAGS %<sm> and %<esm>. If broken
    reactivate above code and remove these flags*)
    sitexrep={};
    dositexrep=False;
    sitemult={};
    effsitemult={};
    (*Print["   Parsing: Reading site X flags: "<>ph];*)
    While[StringStartsQ[GetCLLS[][[1]],"%<"],
      (*Print[{"CLLS: ",CLLS//InputForm}];*)
      (*Print[{"   ",CLLS}];*)
      Which[
        StringMatchQ[CLLS[[1]],"%<sub>"],
          CLLS=Delete[CLLS,1];
          dositexrep=True;
          AppendTo[sitexrep,CLLS];
          ,
        StringMatchQ[CLLS[[1]],"%<sm>"],
          CLLS=Delete[CLLS,1];
          AppendTo[sitemult,CLLS];
          ,
        StringMatchQ[CLLS[[1]],"%<esm>"],
          CLLS=Delete[CLLS,1];
          AppendTo[effsitemult,CLLS];
          ,
        True,
          Print["ERROR in ParseCleanTCAxLines in %<..> parsing at sites. Unknown flag."];
          Abort[];
      ];
    ];
    (*Print["ok 2.4"];*)
    sitexrep=Flatten[sitexrep];
    sitemult=Flatten[sitemult];
    effsitemult=Flatten[effsitemult];
    If[Len[effsitemult]==0 && Len[sitemult]>0,
      effsitemult=sitemult;
    ];
    sitemult=Map[GenNumListFromStringList[GetSIPs[#]]&,sitemult];
    effsitemult=Map[GenNumListFromStringList[GetSIPs[#]]&,effsitemult];
    (*Print[{"site mult: ",sitemult}//InputForm];*)
    
    ln=ln-1;
    
    
    (* Read sitex lines, replace sitex names as go along *)
    sitexls={};
    i=1;
(*    
    While[i<=nsitex,
      CLLS=GetCLLS[];
      Print[{"sitex: ",i,CLLS}];
      If[IntegerQ[CLLS[[1]]],
        (*Add tls2 contents onto tls contents, don't increment i*)
        sitexls[[i-1]]=Join[sitexls[[i-1]],CLLS];
        ,
        (*build sitex replacement name list to use on aid lines below. 
          Replace in sitexls now while here*)
        If[dositexrep,
          sitexrep[[i]]={CLLS[[1]],sitexrep[[i]]};
          (*Print[" before "<>"CLLS[[1]]=sitexrep[[i,2]];"];*)
          CLLS[[1]]=sitexrep[[i,2]];
        ];
        Print[{"  appending to sitexls: ",CLLS}];
        AppendTo[sitexls,CLLS];
        ++i;
      ];
    ];
    *)
    (*Above is broken. If last sitex is a multiline sitex, then 
    it is not appended to sitexls. You need to write a real parser
    of TC site x lines, parse el by el. Or do a diff check and then
    do a ln=ln-1, like in ppn section.*)
    
    CLLS=GetCLLS[];
    While[i<=nsitex||IntegerQ[CLLS[[1]]],
      
      (*Print[{"sitex: ",i,CLLS}];*)
      If[IntegerQ[CLLS[[1]]],
        (*Add tls2 contents onto tls contents, don't increment i*)
        sitexls[[i-1]]=Join[sitexls[[i-1]],CLLS];
        CLLS=GetCLLS[];
        ,
        (*build sitex replacement name list to use on aid lines below. 
          Replace in sitexls now while here*)
        If[dositexrep,
          sitexrep[[i]]={CLLS[[1]],sitexrep[[i]]};
          (*Print[" before "<>"CLLS[[1]]=sitexrep[[i,2]];"];*)
          CLLS[[1]]=sitexrep[[i,2]];
        ];
        (*Print[{"  appending to sitexls: ",CLLS}];*)
        AppendTo[sitexls,CLLS];
        ++i;
        CLLS=GetCLLS[];
      ];
    ];
    ln=ln-1;
    
    
    
    (*Print[{"initial sitexrep",sitexrep}];*)
    
    aidlns={};
    checklns={};
    makelns={};
    dqflns={};
    Do[
      AppendTo[aidlns,{}];
      AppendTo[checklns,{}];
      AppendTo[makelns,{}];
      AppendTo[dqflns,{}];
      ,{i,npc}
    ];
    
    i=1;
    CLLS=GetCLLS[];
    (*Do[If[StringQ[CLLS[[j]]],CLLS[[j]]=StringTrim[CLLS[[j]]]],{j,Length[CLLS]}];*)
    (*Print["ok 3.1"];*)
    (*Print[{"   Parsing: aids section. CLLS =  ",CLLS}];*)
    While[ StringMatchQ[CLLS[[1]],"check"]
        || StringMatchQ[CLLS[[1]],"make"]
        || StringStartsQ[CLLS[[1]],"delG"]
        || StringStartsQ[CLLS[[1]],"DQF"]
        || ElPos[CLLS[[1]],pcnames]>0 ,
      (*Print[{"CLLS: ",CLLS}];*)
      Which[
        StringMatchQ[CLLS[[1]],"check"],
          checklns[[i-1]]=CLLS,
        StringMatchQ[CLLS[[1]],"make"],
          (*Add "norm" flag to make from where needed*)
          makelns[[i-1]]=StandardizeMakeLine[CLLS];
          (*Print[makelns[[i-1]]];*)
          (*makelns[[i-1]]=CLLS*)
          ,
        StringStartsQ[CLLS[[1]],"delG"],
          dqflns[[i-1]]=CLLS,
        StringStartsQ[CLLS[[1]],"DQF"],
          dqflns[[i-1]]=CLLS,
        ElPos[CLLS[[1]],pcnames]>0,
          (* make sitex replacements in aidlns if there are replacements *)
          If[Length[sitexrep]>0,
            Do[CLLS[[k]]=GetReplace[CLLS[[k]],sitexrep],{k,4,Length[CLLS]-1,2}];
          ];
          aidlns[[i]]=CLLS;
          ++i;
      ];
      (*Print["ok 3.2"];*)
      
      
      
      CLLS=GetCLLS[];
      (*ToDo: Need a check here if * or no more lines, before doing loop
              Move Gettls up, use it everywhere, and add ability to flag
              if at end of file or hit *. Set EOData=True, then !TrueQ
              on it to continue.
      *)
      Do[If[StringQ[CLLS[[j]]], CLLS[[j]]=StringTrim[CLLS[[j]]]], {j,Length[CLLS]}];
      (*Print["..."];*)
    ];
    (*Print[aidlns];*)
    (*Last read of lns in tls does not match aid,make,check nor dqf, so is 
    either the first line of a new model, or is end of solutions and is pure phases or * *)
    If[!TrueQ[skipax],
      Print["   Parsing: Appending to MSTRAXDATA "<>ph];
      If[ !StringMatchQ[info[[2,2]],ph], ph=info[[2,2]] ];
      AppendTo[MSTRAXDATA,{ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,
               nsitex,sitemult,effsitemult,sitexls,aidlns,makelns,dqflns,checklns,info}];
    ];
    
    (* Check if last read data is compatible with start of ax model; *)
    (* If not, Break out of loop and process pure phases below *)
    If[IsAXModelStart[CLLS],
      {ph,npc,ver}=CLLS,
      Break[]
    ];
  ];
  (*Print[{"ok1",CLLS}];*)
    
  If[!StringMatchQ[CLLS[[1]],"ENDOFFILE"],
    (*Print[{"ok2",CLLS}];*)
    AppendTo[purephases,CLLS];
    CLLS=GetCLLS[];
    (*Print[{"ok3",CLLS}];*)
    While[!StringMatchQ[CLLS[[1]],"ENDOFFILE"],
      (*Print[{"ok4",CLLS}];*)
      AppendTo[purephases,CLLS];
      CLLS=GetCLLS[];
      (*Print[{"ok5",CLLS}];*)
    ]
  ];
  (*Print[{"ok6",CLLS}];*)
  purephases=Flatten[purephases];
  (*If[Length[purephases]>0,*)
    AppendTo[MSTRAXDATA,purephases];
  (*];*)
  Print["RETURNING from ParseCleanTCAXLines"];
  Return[MSTRAXDATA];
  (*Return[dat];*)
];


(* ::Subsubsection::Closed:: *)
(*GenMinusSpeciesData*)


(* ::Text:: *)
(*For TD +- phase construction only, tc350 fmt. Internal use only. *)
(*2021-02-31, updated 2022-03-25.*)
(*ToDo: Integrate behavior into a master ConstructSpecies func to reduce*)
(*multiple functions that do same/similar thing of adding species to td dicts.*)
(*Initially copied from dsaxncvrt`GenTDODEMinusSpecidesDatFromDQFtc350fmt... *)
(*but null changed to norm plus some other minor changes. *)
(*Should potentially reside directly in dset.wl*)


GenMinusSpeciesData[{newsp_String,fromsp_String,dsv_String,dqfls_List,plusminus_String,dropod_:True,otype_:"norm"}]:=Module[
  {
    pc,
    hadj,sadj,vadj
  },
  pc=newsp; 
  {hadj,sadj,vadj}=dqfls;
  Which[
    otype=="norm",
      pc=pc<>plusminus,
    otype=="equilibrium",
      pc=pc<>"E"<>plusminus,
    otype=="ordered",
      pc=pc<>"O"<>plusminus,
    otype=="disordered",
      pc=pc<>"D"<>plusminus
  ];
    
  
  tdabb[pc,dsv]=pc;
  tdver[pc,dsv]=tdver[fromsp,dsv];
  tdfmt[pc,dsv]=tdfmt[fromsp,dsv];
  (* **** *)
  AppendTo[tddsetabbs,{pc,dsv}];
  (* **** *)
                         
  If[ TrueQ[dropod], tdtype[pc,dsv]="norm", tdtype[pc,dsv]=tdtype[fromsp,dsv]];
              
  (*ToDo: Put a check for non-zero length tdreact and tddqf*)
  If[Length[tdreact[fromsp,dsv]]>0,
    Print["ERROR in GenMinusSpeciesData: made from species "<>tdabb[fromsp,dsv]<>" has non-zero length of tdreact, and"<>
          "this is not fully checked yet. tdreact = "<>ToString[tdreact[fromsp,dsv]]];
    Abort[];
  ];
  
  (*Should be ||, not &&; 2022-05-26*)
  If[tddqf[fromsp,dsv][[1]]!=0 || tddqf[fromsp,dsv][[2]]!=0 || tddqf[fromsp,dsv][[3]]!=0,
    Print["WARNING in GenMinusSpeciesData: made from species "<>tdabb[fromsp,dsv]<>" has non-zero tddqf, and "<>
          "this is not implemented yet. "<>
          "tddqf = "<>ToString[tddqf[fromsp,dsv]]];
    Abort[];
  ];
  tdreact[pc,dsv]=tdreact[fromsp,dsv];(* This should be {}. If not, need a more elaborate processing *)
  tddqf[pc,dsv]=tddqf[fromsp,dsv];(* This should be {}. If not, need a more elaborate processing *)
  
  tdcomp[pc,dsv]=tdcomp[fromsp,dsv];
  tdformula[pc,dsv]=tdformula[fromsp,dsv];
  tdhf[pc,dsv]=tdhf[fromsp,dsv]+hadj;(* hadj passed in from calling module. Should just be o/d mods, not tc ax file dqf value *)
  tdsr[pc,dsv]=tdsr[fromsp,dsv]+sadj;(* sadj passed in from calling module. Should just be o/d mods, not tc ax file dqf value *)
  tdvr[pc,dsv]=tdvr[fromsp,dsv]+vadj;(* vadj passed in from calling module. Should just be o/d mods, not tc ax file dqf value *)
  tdcpa[pc,dsv]=tdcpa[fromsp,dsv];
  tdcpb[pc,dsv]=tdcpb[fromsp,dsv];
  tdcpc[pc,dsv]=tdcpc[fromsp,dsv];
  tdcpd[pc,dsv]=tdcpd[fromsp,dsv];
  tdaugb[pc,dsv]=tdaugb[fromsp,dsv];
  tdao[pc,dsv]=tdao[fromsp,dsv];
  tdaoc[pc,dsv]=tdaoc[fromsp,dsv];
  tdk0[pc,dsv]=tdk0[fromsp,dsv];
  tdk0p[pc,dsv]=tdk0p[fromsp,dsv];
  tdk0pp[pc,dsv]=tdk0pp[fromsp,dsv];
  tdtheta[pc,dsv]=tdtheta[fromsp,dsv];
  tddkdt[pc,dsv]=tddkdt[fromsp,dsv];
  If[TrueQ[dropod],
    tdtcr[pc,dsv]=0;
    tdpcr[pc,dsv]=0;
    tdsmax[pc,dsv]=0;
    tdvmax[pc,dsv]=0;
    tdsfdh[pc,dsv]=0;
    tdsfdhv[pc,dsv]=0;
    tdsfw[pc,dsv]=0;
    tdsfwv[pc,dsv]=0;
    tdsfn[pc,dsv]=0;
    tdsffac[pc,dsv]=0
    ,
    tdtcr[pc,dsv]=tdtcr[fromsp,dsv];
    tdpcr[pc,dsv]=tdpcr[fromsp,dsv];
    tdsmax[pc,dsv]=tdsmax[fromsp,dsv];
    tdvmax[pc,dsv]=tdvmax[fromsp,dsv];
    tdsfdh[pc,dsv]=tdsfdh[fromsp,dsv];
    tdsfdhv[pc,dsv]=tdsfdhv[fromsp,dsv];
    tdsfw[pc,dsv]=tdsfw[fromsp,dsv];
    tdsfwv[pc,dsv]=tdsfwv[fromsp,dsv];
    tdsfn[pc,dsv]=tdsfn[fromsp,dsv];
    tdsffac[pc,dsv]=tdsffac[fromsp,dsv];
  ];
  (* just return modified name; (-,E-,D-,O-) or (+,E+,D+,O+) *)
  Return[pc];
];


(* ::Subsubsection::Closed:: *)
(*ProcessMakeVecTC350*)


(* ::Text:: *)
(*Copied from dsaxcnvrt.wl ProcessMakeVecTC350fmt, which takes a list of only strings. *)
(*New TCAx processors GenExListFromString... have already converted string numbers to *)
(*actual numbers, so dsaxcnvrt`ProcessMakeVecTC350fmt[..] does not work directly.*)
(*Expects an imvec list of strings, in following form:  *)
(*{"Mmt", "ds633", "make", 1, "equilibrium", "mt", 1, "delG(od)", 10.0, 0.001, 0.1}*)


(*Clear[ProcessMakeVec];*)
(*Options[ProcessMakeVecTC350fmt]={"AdjFromSpeciesHSV"->True};*)
ProcessMakeVecTC350[imvec_List,plusminus_String,opts_:OptionsPattern[]]:=Module[
  {adjFromSpHSV, mvec,     nspname="null", dsv="null", res={},      specres={},
   dqfpos=0,    dqfflags,  odpos={},       odflags,    hasod=False, nfrom=0,
   hdqf=0.0,    sdqf=0.0,  vdqf=0.0, pos, epos, newmfspname={},
   newcomp,     newtdreac, mainmakepcname, dumb
  },
  (*Print["=== MakeVecTC350:"<>imvec[[1]]<>" ===========================\n"];*)
  (*Print[imvec];*)
  (*adjFromSpHSV=OptionValue[ProcessMakeVec,opts,"AdjFromSpeciesHSV"];*)
  
  (*For now, always force True; old behavior of directly adjusting made 
    from species if need to adjust h/s/v od props*)
  adjFromSpHSV=True;
  
  mvec=imvec;
  If[mvec[[3]]!="make",Print["ERROR ProcessMakeVec: mvec[[3]] != make."]; Abort[]];
  (*nfrom=ToExpression[mvec[[4]]];*)
  nfrom=mvec[[4]];
  If[Head[nfrom]!=Integer,Print["ERROR ProcessMakeVec: mvec[[4]] != Integer."]; Abort[]];
  
  nspname=mvec[[1]];
  dsv=mvec[[2]];
  
  dqfflags={"DQF", "delG(make)", "delG(mod)", "delG(od)", "delG(tran)", "delG(rcal)"};
  
  (*This gets the dqf of the main make phase as listed in the tc-ax file.*)
  (*Sets hdqf, sdqf, vdqf, which will be dqf's on the COM phase*)
  If[ ContainsAny[mvec,dqfflags], 
    Do[
      If[MemberQ[dqfflags,mvec[[i]]],dqfpos=i;Break[]];
      ,{i,5,Length[mvec]}
    ];
    hdqf=mvec[[dqfpos+1]];
    sdqf=mvec[[dqfpos+2]];
    vdqf=mvec[[dqfpos+3]];
  ];
  
  (*New ParseCleanTCAxLines requires a flag for being made from ds species, 
    even for norm or liq pc's*)
  odflags={"norm","equilibrium","ordered","disordered"};
  
  (*Sets idx of every odflag in mvec; so can now have multiple odflags!*)
  If[ ContainsAny[mvec,odflags], 
    Do[
      If[MemberQ[odflags,mvec[[i]]],AppendTo[odpos,i];hasod=True];
      ,{i,5,Length[mvec]}
    ];
  ];
  
  epos=If[dqfpos>0,dqfpos-1,Length[mvec]];
  
  (*Build up res, a nested list of pc to make from*)
  (*New ParseCleanTCAxLines now adds norm for pc's w/o flags,
    so once checked can remove appending of null below*)
  pos=4;
  While[pos+1<epos,
    specres={};
    pos+=1;
    If[MemberQ[odpos,pos], 
      AppendTo[specres,mvec[[pos]]]; 
      pos+=1
      , 
      Print["ERROR in specres loop in ProcessMakeVecTC350.. aborting."];
      Abort[];
      (*AppendTo[specres,"null"]*)
    ];
    AppendTo[specres,mvec[[pos]]];
    pos+=1;
    AppendTo[specres,mvec[[pos]]];
    
    AppendTo[res,specres];
  ];
  
  If[Length[res]!=nfrom, Print["ERROR ProcessMakeVec: Length[res] != nfrom. Abort[]."]; Abort[] ];
  
  Block[{dqfcont,hadj,sadj,vadj,cdat,odtype,spname,sptype,adjspname,dropod},
    dqfcont={0,0,0};
    
    (*2022-03-25; Moved this array initialization to before do loop so
      I can update the make from sp name if I adjust the sp name being
      made from within the do loop.*)
    newtdreac=Table[{res[[j,2]],res[[j,3]]},{j,1,Length[res]}];(*  Check format of tdreac for pc *)
    
    (*Loop over the existing pc to make from*)
    Do[
      
      {hadj,sadj,vadj}={0,0,0};
      cdat=res[[i]];(*curr pc making from*)
      odtype=cdat[[1]];(*norm,ordered,disordered or equilibrium*)
      spname=tdabb[cdat[[2]],dsv];(*kind of dumd as they are the same name, but will fail if not loaded*)
      (*Print["--- Making from ds species: "<>spname<>"/n"];*)
      adjspname=spname;
      sptype=tdtype[cdat[[2]],dsv];(*norm,liq,land,bw,h2o,co2,cs,aq*)
      
      
      
      
 
      
      (*If[odtype != "norm",*)
      If[True,
      
        (*If[sptype=="land" || sptype=="bw",*)
        If[True,
          (*Print[{"Making minus species for make from species of type: ",odtype}];*)
          Switch[odtype,
            
            "equilibrium",
              (
                If[sptype!="land" && sptype!="bw",
                  Print["ERROR in TCax file: a non-od phase has equilibrium flag."];
                  Abort[];
                ];
                (*This gives default adjusted species name, before checking and appending
                a num to it in GenTDODE... if it already exists*)
                adjspname=adjspname<>"E"<>plusminus;
                If[TrueQ[adjFromSpHSV],
                  dropod=False; (*False so don't trop od data from td dicts in E case*)
                  adjspname=GenMinusSpeciesData[{spname,spname,dsv,{0,0,0},plusminus,dropod,"equilibrium"}];
                  (*2022-03-25 addition of adjusting name in tdreact*)
                  newtdreac[[i,1]]=adjspname;(* This is saved for the make/com phase tdreact, not madefrom tdreact dict*)
                  Print[{spname,adjspname,spname->adjspname}];
                  AppendTo[newmfspname,{spname->adjspname}];(*2022-05-26: was getting pa-rg in tdreact of kprg*)
                  (*AppendTo[newmfspname,{adjspname}];*)
                ];
              ),
              
            "disordered", (* checked and good (compared dilm) *)
              (
                dropod=True; (*True so trop od data from td dicts in D case*)
                If[sptype!="land" && sptype!="bw",
                  Print["ERROR in TCax file: a non-od sp made from has disordered flag."];
                  Abort[];
                ];
                adjspname=adjspname<>"D"<>plusminus;
                (*Print["Before GenTDODEMinusSpeciesDataFromDQF, adjspname = "<>adjspname];*)
                (*2022-03-25: Taking ilm with make 1 disordered ilm 1 for ds62, tc350 
                  output shows adjustment that is equiv to {hf,sr,vr}_ilm + LandRefHSV + DQFHSV 
                  However, tc-project.csf output shows the phase constructed from as ilm with land
                  smax and vmax properties just dropped but keeps hf,sr,vr exactly same as real ilm 
                  (so it does not strip lnd contribution from hsf,sr,vr of the base ilm); It then 
                  incorporates LandRefHSV contributions to the dqf listed in tc-ax file to use a 
                  diff dqf internally than listed in tc-ax file; the end result is the same as  
                  my ilmD- + reported dqf on COM approach here.  *)
                If[sptype=="land", 
                  (* Landau and LandRefHSV give exact same results for dis case *)
                  {hadj,sadj,vadj}=Landau[spname,dsv,0.001,298.15,"dis","RetVal"->{"Hlnddis","Slnddis","Vlnddis"}];
                  (*{hadj,sadj,vadj}=LandRefHSV[spname,dsv];*)
                  , 
                  {dumb,hadj,sadj,vadj}=SFBWeq[spname,dsv,0.001,298.15,"dis"];
                  (*{hadj,sadj,vadj}=BWRefHSV[spname,dsv]*)
                ];
                If[TrueQ[adjFromSpHSV],                                       
                  adjspname=GenMinusSpeciesData[{spname,spname,dsv,{hadj,sadj,vadj},plusminus,dropod,"disordered"}];
                  (*2022-03-25 addition of adjusting name in tdreact*)
                  newtdreac[[i,1]]=adjspname;
                  Print[{spname,adjspname,spname->adjspname}];
                  AppendTo[newmfspname,{spname->adjspname}];(*2022-05-26: was getting pa-rg in tdreact of kprg*)
                  (*AppendTo[newmfspname,{adjspname}];*)
                ];
              ),
              
            "ordered", (* possible error in mp50-03 mp50-04 and mp50-05 td ax files; checked against dilm *) 
              (
                dropod=True; (*True so trop od data from td dicts in O case*)
                If[sptype!="land" && sptype!="bw",
                  Print["ERROR in TCax file: a non-od sp made from has ordered flag."];
                  Abort[];
                ];
                adjspname=adjspname<>"O"<>plusminus;
                If[sptype=="land", 
                  (* 2022-04-29: Switched to use Landau as LandRefHSV does not give correct vals for ord case *)
                  (*do a test with qtz to see if need 0.001 or 0.000*)
                  {hadj,sadj,vadj}=Landau[spname,dsv,0.001,298.15,"ord","RetVal"->{"Hlndord","Slndord","Vlndord"}];
                  
                  (*2022-03-25: POSSIBLE ERROR, TO FIX IN DSAXCNVRT` also *)
                  (*This hadj using Landds6 is in error as index [[1]] is not correct. Need to just use
                    simple LandRefHSV modification based on ilm->dilm tc350. *)
                  (*hadj = Landds6[spname,dsv,0.001,0.7][[1]];*)(* ERROR? Landds6[][[1]] is Gval, not h. [[1]]->[[3]] 2022-03-25*)
                  (*sadj = -tdsmax[spname,dsv] - SLnd[spname,dsv,0.001,298.15];*)
                  (* ToDo: vadj *)
                  (*Print["NOTE: In OD ordered case when making from a land species. adjust sdqf by Slnd_0 - Slnd_298.. to check ******."];*)
                  ,
                  (*Print["Note: In OD ordered case when making from a BW species. Not making any adjustments.."];*)
                  {dumb,hadj,sadj,vadj}=SFBWeq[spname,dsv,0.001,298.15,"ord"];
                ];
                If[TrueQ[adjFromSpHSV],
                  adjspname=GenMinusSpeciesData[{spname,spname,dsv,{hadj,sadj,vadj},plusminus,dropod,"ordered"}];
                  (*2022-03-25 addition of adjusting name in tdreact*)
                  newtdreac[[i,1]]=adjspname;
                  Print[{spname,adjspname,spname->adjspname}];
                  AppendTo[newmfspname,{spname->adjspname}];(*2022-05-26: was getting pa-rg in tdreact of kprg*)
                  (*AppendTo[newmfspname,{adjspname}];*)
                ];
              ),
            "norm",
              (
                dropod=False; (* False there is no od data in td dicts in norm case. ERROR if it does have od data! *)
                If[sptype=="land" || sptype=="bw", Print["ERROR in TCax file: a land or bw sp made from has norm flag."]];
                If[TrueQ[adjFromSpHSV],
                  adjspname=GenMinusSpeciesData[{spname,spname,dsv,{hadj,sadj,vadj},plusminus,False,"norm"}];
                  (*2022-03-25 addition of adjusting name in tdreact*)
                  newtdreac[[i,1]]=adjspname;
                  Print[{spname,adjspname,spname->adjspname}];
                  AppendTo[newmfspname,{spname->adjspname}];(*2022-05-26: was getting pa-rg in tdreact of kprg*)
                  (*AppendTo[newmfspname,{adjspname}];*)
                ];
              ),  
            _,
            Print["ERROR ProcessMakeVec: Unknown od type in Switch. Abort[]."; Abort[] ]
          ];
          ,
          Print["ERROR ProcessMakeVec: type is not an odtype. Abort[]."; Abort[] ];
        ];(*END If[sptype=="land" || sptype=="bw",*)
        
        
        (* Still in If[odtype != "norm" but pc made from is not land nor bw, so still
           Need to make a straight +/- phase of normal type
         *)
        
        
        
        (*2022-03-25 My logic here seems to be to make from a pc for an od case w/o ever adjusting 
          the thermo properties of that minus phase, when adjFromSpHSV=False. This would mean the 
          ord modifications for ordered and disordered must go to the COM phase, which appears what
          tc350 are doing based on it's output thermo table. For now I have not
          tested this carefully, so keep adjFromSpHSF True until sure all the rest of code is 
          working correctly so the difference is isolated to this. *)
        (*Check dsaxcnvrt` b/c null was misspelled. However, should not really be here b/c we are in 
          an If case of !null, so odtype=="null" will never occure here. 2022-03-25*)
        If[!TrueQ[adjFromSpHSV],
          (*Print[{"dbg 013:odtyp=null",spname}];*)
          Print["DANGER: In If[!TrueQ[adjFromSpHSV] and we don't want that."];
          dropod=False;(* Should this be true?? *)
          adjspname=GenMinusSpeciesData[{spname,spname,dsv,{0,0,0},plusminus,dropod,"null"}];
          (*Print[{"dbg 014:odtyp=null",spname}];*)
          
          (*2022-03-25 addition of adjusting name in tdreact*)
          newtdreac[[i,1]]=adjspname;
          AppendTo[newmfspname,{spname->adjspname}];
          (*,
          Print["Likely ERROR: Made it through w/o creating a make from pc with odtype of null!"];*)
        ];
        
      ];(*END If odtype!=norm*)
      
      
      
      
      
      
      
      
      (*When adjFromSpHSV==False, store the adjustments for current make from
      in dqfcont to be added to the COM phase further below (through h/s/vdqf)*)
      If[adjFromSpHSV==False,
        Print["Storing adjustment to add to COM phase later"];
        dqfcont=dqfcont+{hadj,sadj,vadj}*cdat[[3]];
      ];
      
      ,{i,1,nfrom}
    ];(*END If[odtype != "null",*)
    
    (*Note: dqfcont will always be zero when adjFromSpHSV is True.*)
    hdqf=hdqf+dqfcont[[1]];
    sdqf=sdqf+dqfcont[[2]];
    vdqf=vdqf+dqfcont[[3]];
    (*Print[{"{hdqf,sdqf,vdqf} for COM phase = ",{hdqf,sdqf,vdqf}}];*)
    
  ];
  (*Print[{"newmfspname: ",newmfspname}];*)
  
  (* ok, at this stage all the made from species exist, with their new names in newmfspname list. However,
     it appears if od was set to null due to no ordered or disordered keyword in tc make statement, the made 
     from is/was not made b/c it is a straight HP ds pc w/o any modification; in that case it needs added to
     the minus phase list unless it is in the same ax model that the COM phase belongs to; check this 
     happens (2022-03-25).
     Now need to make actual pc species for ax model as a COM phase. Requires checking if exists, increment 
     num on name, make tdreact, tddqf, comp, for COM behavior, and pass back new pc name.  Get this code from 
     part of MakePCds6 in thermo.*)
     
  (* MARK - adjust name if needed (already exists) *)
  Block[{exists,pc,suf,strsuf},
    pc=mvec[[1]];
    exists=pcLoaded[pc,dsv];
    (*Print[{"phase to make already exist?",exists,pc}];*)
    If[exists,
      suf=1;
      strsuf=TextString[suf];
      pc=StringInsert[pc,strsuf,-1];
      While[pcLoaded[pc,dsv],
        pc=StringTake[pc,{1,StringLength[pc]-StringLength[strsuf]}];
        ++suf;
        strsuf=TextString[suf];
        pc=StringInsert[pc,strsuf,-1];
      ];
      Print["ax pc "<>mvec[[1]]<>" changed to "<>pc];
    ];
    
    (*newcomp=Table[0.,{Length[tdcomp[mvec[[5]],dsv]]}];*)(*Active until 2022-03-24*)
    newcomp=Table[0.,{Length[tdcomp[res[[1,2]],dsv]]}];
    (*Print["dbg 020, newcomp= ",newcomp," mvec[[5]]= ",mvec[[5]]];*)
  
    (*ToDo: 2022-03-24: Fix the following do loop in dsaxcnvrt`, from using newmfspname to res[[1,2]]. For
    od=null case, newmfspname not getting set. For now the od pc made should have same
    comp as original, so this should work. TO CHECK MANUAL OUTPUT to make sure this is
    always the case.*)
    Do[
      (*AppendTo[newmvec,{mvec[[i]],ToEx[mvec[[i+1]]]}];*)(*Dont think this is needed*)
      (*Print[{newmfspname[[i]]}];*)
      (*Print[{res[[i,2]]}];*)
      (*newcomp+=tdcomp[newmfspname[[i]],dsv]*res[[i,3]];*)(*Active until 2022-03-24*)
      (*Assumes Rationals like 3/7 are not in string format, so works with new MatchPatt and GetEx.*)
      newcomp+=tdcomp[res[[i,2]],dsv]*res[[i,3]];
      ,{i,1,Length[res]}
    ];
    (*Print[{"newcomp = ",newcomp}];*)
    (*Print["dbg 021"];*)
    
    (*2022-03-25; This needs to pick up new sp name being made from, b/c that can be adjusted above
    when making new sp made from. Needs fixed in dsaxcnvrt` also.
    Moved to above do loop so can change newtdreact as sp names are adjusted*)
    (*newtdreac=Table[{res[[i,2]],res[[i,3]]},{i,1,Length[res]}];*)(*  Check format of tdreac for pc *)
    
    (*Print[{"newtdreac = ",newtdreac}];*)
    
    (*Print["dbg 022"];*)
    
    tdabb[pc,dsv]=pc;
    tdver[pc,dsv]=dsv;
    (* **** *)
    AppendTo[tddsetabbs,{pc,dsv}];
    (* **** *)
    tdfmt[pc,dsv]=StringTake[dsv,3];
    tdtype[pc,dsv]="make";
    tdreact[pc,dsv]=newtdreac;
    tddqf[pc,dsv]={hdqf,sdqf,vdqf};
    tdcomp[pc,dsv]=newcomp;
    (*Print["Calling tdformula with TDCompOrder: ",TDCompOrder];*)
    tdformula[pc,dsv]=GenFormulaString[newcomp,TDCompOrder];
    (*Print[{"tdformula[pc,dsv]= ",tdformula[pc,dsv]}];*)
    (*Print["Done Calling tdformula with TDCompOrder: ",TDCompOrder];*)
    (*Print[tdformula[pc,dsv]];*)
    (* set all the rest to zero, completeness *)
    tdhf[pc,dsv]=0.;
    tdsr[pc,dsv]=0.;
    tdvr[pc,dsv]=0.;
    tdcpa[pc,dsv]=0.;
    tdcpb[pc,dsv]=0.;
    tdcpc[pc,dsv]=0.;
    tdcpd[pc,dsv]=0.;
    tdaugb[pc,dsv]=0.;
    tdao[pc,dsv]=0.;
    tdaoc[pc,dsv]=0.;
    tdk0[pc,dsv]=0.;
    tdk0p[pc,dsv]=0.;
    tdk0pp[pc,dsv]=0.;
    tdtheta[pc,dsv]=0.;
    tddkdt[pc,dsv]=0.;
    tdtcr[pc,dsv]=0.;
    tdpcr[pc,dsv]=0.;
    tdsmax[pc,dsv]=0.;
    tdvmax[pc,dsv]=0.;
    tdsfdh[pc,dsv]=0.;
    tdsfdhv[pc,dsv]=0.;
    tdsfw[pc,dsv]=0.;
    tdsfwv[pc,dsv]=0.;
    tdsfn[pc,dsv]=0.;
    tdsffac[pc,dsv]=0.;
    
    mainmakepcname=pc;
  ];
  
  
     
  (*tpcabb=MakePCds6[mvec];*)
  Return[{mainmakepcname,newmfspname}];
];


(* ::Subsubsection::Closed:: *)
(*NCompTbl  -  Modified from ngen`axmodel`GenNCompTbl;*)


(*Clear[NCompTbl];*)
(*
Copied from GenNCompTbl in axmodel`, 2022-04-04. Modified by passing in 
 1. pcsitcomp instead of having func use global ncmp arrays.
 2. passing in a list of just sites, sites_: e.g. {M1,M2,T}. This is 
    constructed in func in GenNCompTbl, but calling code already 
    has this list prepared, so pass in and remove sitesin from func..
    
 3. sitesin_ arg changed to ellsin_, which is nested list of elements
    mixing on each site. This is constructed in func in GenNCompTbl, but
    calling code already has this list prepared, so pass in and remove
    elsites from func.
*)
NCompTbl[pcs_, sites_, ellsin_, pcsitcomp_, opts:OptionsPattern[]]:=Module[
  {(*pcs={},*) cts,tbl={},hdr={"pc's"},(*sites={},*)(*elsites={},*) tval,tls,cmptbl},
  (*
  examples of formats:
  sites = {"M1","M2","T"}
  elsites={ {"Fe","Mg","Mn"},{"Fe","Mg","Mn","Ca"},{"Si","Al"} }
  *)
  cts=OptionValue[GenNCompTbl,CheckTableStructure];
  (*Map[(AppendTo[sites,#[[1]]];AppendTo[elsites,Drop[#,1]])&,sitesin];*)
  Do[
    tval=sites[[i]];
    tls=ellsin[[i]];
    Map[(AppendTo[hdr,"n"<>#<>tval])&,tls];
    ,{i,Length[sites]}
  ];
  AppendTo[tbl,hdr];
  Block[{tdatentry,csite,cel,tentry},
  (* 1. loop through pc's *)
  Do[
    tdatentry={pcs[[i]]};
    (*tls=Global`ncmp[pcs[[i]]];*)
    tls=pcsitcomp[[i]];
    If[!Head[tls]===List,
      Print["Head[tls]!=List so pc "<>pcs[[i]]
        <>" must not be defined in your session. Please define site n comps for all pc's.\nAborting calculation."];
      Abort[];
    ];
    (* 2. Loop through sites *)
    Do[
      csite=sites[[j]];
      (* 3. Loop through elsites *)
      Do[  
        cel=ellsin[[j,k]];
        (* 4. now loop through pc ncmp (tls) structure to see if 
           it has an entry for cel in csite. If so, 
           append it's value to data tbl. If not, append
           a 0 to data tbl.
        *)
        tentry=0;
        Do[
          If[cts==True,
          If[ ElPos[tls[[l,1]],sites]==0,
            Print["ERROR: Site "<>tls[[l,1]]<>" for pc "<>pcs[[i]]<>" is not in your table."];
            Print["Aborting calculation!"];
            Abort[];
            ,
            If[tls[[l,1]]==csite && ElPos[tls[[l,2]],ellsin[[j]]]==0,
              Print["ERROR: Element "<>tls[[l,2]]<>" on site "<>csite<>" for pc "<>pcs[[i]]<>" is missing from your table."];
              Print["Aborting calculation!"];
              Abort[];
            ];
          ];
          ];
          If[tls[[l,1]]==csite&&tls[[l,2]]==cel,tentry=tls[[l,3]];Break;];
          ,{l,Length[tls]}
        ];(* end loop 4, over ncmp structure of curr pc *)
        AppendTo[tdatentry,tentry];
        ,{k,Length[ellsin[[j]]]}
      ];(* end loop 3, over elsites *)
      ,{j,Length[sites]}
    ];(* end loop 2, over sites *)
    AppendTo[tbl,tdatentry];
    ,{i,Length[pcs]}
  ];(* end loop 1, over pcs *)
  ];(* end Block *)

(*Return[{hdr,sites,elsites}];*)
cmptbl=tbl[[2;;Length[tbl],2;;Length[tbl[[1]]]]];
(*Print[Grid[tbl,Dividers->{{2->True},{2->True}}]];*)
Return[{cmptbl,Grid[tbl,Dividers->{{2->True},{2->True}}]}];
];


(* ::Subsubsection::Closed:: *)
(*TblToTextTbl:= write tbl as string*)


(*Clear[TblToTextTbl];*)
(*To Copy to ngen`utils` after make output string more flexible (additonal header arg, spacing arg, etc.*)
TblToTextTbl[itbl_?ArrayQ, hdr_?StringQ]:=Module[{res="",stbl="",nc=0,nr=0,maxw={},nsep=1,sl,ts},
  (* set table dimension vals *)
  If[Len[Dimensions[itbl]]!=2,
    Print["Length of itbl passed to TblToTextTbl !=2. Aborting."];
    Abort[];
  ];
  nr=Len[itbl];
  nc=Len[itbl[[1]]];
  (* --- Build stbl where each cell is a trimmed string *)
  stbl=Table[".",{nr},{nc}];
  maxw=Table[0,nc];
  
  Do[
    Do[
      If[Head[itbl[[i,j]]]===String,
        stbl[[i,j]]=itbl[[i,j]];
        ,
        stbl[[i,j]]=StringTrim[ToString[(itbl[[i,j]])//InputForm]];
      ];
      If[j==1,
        stbl[[i,j]]="!   "<>stbl[[i,j]];
        (*stbl[[i,j]]=stbl[[i,j]]<>" |";*)
      ];
      sl=StringLength[stbl[[i,j]]];
      If[sl>maxw[[j]], maxw[[j]]=sl]
      ,{j,nc}
    ];
    ,{i,nr}
  ];
  maxw=maxw+nsep;(*array add*)
  Do[
    If[i==1,
      (*res=res<>"!   Site mixing table (with effective atom #'s given by ax model).\n";*)
      res=res<>"!   "<>hdr<>"\n";
    ];
    Do[
      
      If[j==1,
        ts=SPadR[stbl[[i,j]],maxw[[j]]];
        ts=ts<>"|  ";
        ,
        ts=SPadL[stbl[[i,j]],maxw[[j]]];
      ];
      res=res<>ts;
      ,{j,nc}
    ];
    If[i<nr,res=res<>"\n"];
    ,{i,nr}
  ];
  
  Return[res];
];


(* ::Subsubsection::Closed:: *)
(*TblToTextTbl2:= write tbl as string, neater format than TblToTextTbl. NotFinished*)


(*Clear[TblToTextTbl];*)
(*To Copy to ngen`utils` after make output string more flexible (additonal header arg, spacing arg, etc.*)
TblToTextTbl2[sitnames_, ells_, itbl_?ArrayQ, hdr_?StringQ]:=Module[{res="",stbl="",nc=0,nr=0,
  maxw={},nsep=1,sl,ts},
  (* set table dimension vals *)
  If[Len[Dimensions[itbl]]!=2,
    Print["Length of itbl passed to TblToTextTbl !=2. Aborting."];
    Abort[];
  ];
  
  nr=Len[itbl];
  nc=Len[itbl[[1]]];
  maxw=Table[0,nc];
  Print[{"itble: ",nr,nc}];
  (* --- Build stbl where each cell is a trimmed string *)
  stbl=Table[".",{nr+1},{nc}];
  
  Block[{cc},
    stbl[[1,1]]="!   ";
    stbl[[2,1]]="!   ";
    cc=1;
    Do[
      Do[
        cc=cc+1;
        stbl[[1,cc]]=sitnames[[i]];
        stbl[[2,cc]]=ells[[i,j]];
        If[StringLength[stbl[[1,cc]]]>maxw[[cc]],maxw[[cc]]=StringLength[stbl[[1,cc]]]];
        If[StringLength[stbl[[2,cc]]]>maxw[[cc]],maxw[[cc]]=StringLength[stbl[[2,cc]]]];
        ,{j,1,Len[ells[[i]]]}
      ];
      ,{i,1,Len[sitnames]}
    ];
  ];
  
  Do[
    Do[
      If[Head[itbl[[i,j]]]===String,
        stbl[[i+1,j]]=itbl[[i,j]];
        ,
        stbl[[i+1,j]]=StringTrim[ToString[(itbl[[i,j]])//InputForm]];
      ];
      If[j==1,
        stbl[[i+1,j]]="!   "<>stbl[[i+1,j]];
        (*stbl[[i+1,j]]=stbl[[i+1,j]]<>" |";*)
      ];
      sl=StringLength[stbl[[i+1,j]]];
      If[sl>maxw[[j]], maxw[[j]]=sl]
      ,{j,nc}
    ];
    ,{i,2,nr}
  ];
  
  maxw=maxw+nsep;(*array add*)
  (*Loops over entire stbl, which has 1 more row than itbl. So increment nr*)
  nr=nr+1;
  Block[{cc},
    Do[
      If[i==1,
        (*res=res<>"!   Site mixing table (with effective atom #'s given by ax model).\n";*)
        res=res<>"!   "<>hdr<>"\n";
      ];
      cc=1;
      ts=SPadR[stbl[[i,1]],maxw[[1]]];
      ts=ts<>"   ";
      res=res<>ts;
      
      Do[
        Do[
          cc=cc+1;
          ts=SPadL[stbl[[i,cc]],maxw[[cc]]];
          res=res<>ts;
          ,{k,Len[ells[[j]]]}
        ];
        res=res<>"  ";
        ,{j,Len[sitnames]}
      ];
      
      If[i<nr,res=res<>"\n"];
      ,{i,nr}
    ];
  ];
  Print[{"res = ",res}];
  Return[res];
];


(* ::Subsubsection:: *)
(*NewAnalSites*)


(*ToDo: Would be good to get data in order as needed by GenNCompTbl and GenXCompTbl
in ngen`axmodel`
Currently, sm and esm passed in are not actually used... to do. *)

NewAnalSites[ph_?StringQ, sxs_?ListQ, aids_?ListQ, sm_?ListQ, esm_?ListQ]:=Module[{
  smesmfac, smtblstr, sfnames, sitnames, sitels, 
  ells={}, nells={}, sitmult={}, pcsitcomp={}, pcs={},
  aidnorm, tnells, tpcsitcomp, ncmptbl, tvec, tbls1, tbls2, tbls3, ti, 
  redincstr },
  (* Args *)
  (* sxs: have site name at first pos.*)
  (* aids: have pc name at first pos.*)
  (* sm:  site mult (atoms on site).     { {M1,1}, {M2,3}, {T,2} } *)
  (* esm: eff site mult (atoms on site). { {M1,1}, {M2,3}, {T,1} } *)
  
  (*Determine which if any sites are effective sites, and mult.*)
  (*Print["nas1"];*)
  smesmfac=MapThread[#1[[2]]/#2[[2]]&,{esm,sm}];
  (*Print["nas1.1"];*)
  smtblstr="Site mixing table. ";
  redincstr="reduced";
  ti=0;
  Do[
    If[smesmfac[[i]]!=1,
      ti+=1;
      If[smesmfac<1,redincstr="reduced",redincstr="increased"];
      If[ti>1,smtblstr=smtblstr<>"\n"];
      smtblstr=smtblstr<>"!   Site "<>sm[[i,1]]<>" is "<>redincstr<>" to ";
      smtblstr=smtblstr<>ToString[smesmfac[[i]],InputForm]<>" * the true multiplicity. ";
    ];
    ,{i,Length[smesmfac]}
  ];
  ti=0;
  (*Print["nas2"];*)
  (* Code *)
  (* sfnames:  {x(Mg,M1),x(Fe,M1),x(Fe3,M1),x(Al,M1),x(Mg,M2),x(Fe,M2),x(Ca,M2),x(Al,T),x(Si,T)}*)
  (* sitnames: {M1,M2,T} *)
  (* sitels :  {Mg,Fe,Fe3,Al},{Mg,Fe,Ca},{Al,Si}} *)
  sfnames = Map[#[[1]]&,sxs];
  sitnames = DeleteDuplicates[Map[(GetSIPs[#[[1]]][[2]])&,sxs]];
  (*Print["nas3"];*)
  sitels = {};
  (* Build sitels: technically not needed 
     unless reorder results below to match *)
  Block[{cs,cels},
    Do[
      cs = sitnames[[i]];
      cels = {};
      Map[ If[StringMatchQ[cs, (GetSIPs[#])[[2]]],
           AppendTo[cels,(GetSIPs[#])[[1]] ]]&
        , sfnames
      ];
      AppendTo[sitels,cels];
      ,{i,Length[sitnames]}
    ];
  ];
  (*Print["nas4"];*)
  (* Loop through sits, aids and get site mixing *)
  (* old analsites does a 3-nested loop, over len(sitenames), len(aids), len(curr aid) *) 
  (* ToDo: This can end up with sites in diff order than in sitels, so adjust code below
           to keep els in same order as siteels. *)
  Block[{tls,tsit,tel,tnel,tpos,sielidx},
    (* sielidx is idx of current read el in sitels; to implement *)
    
    Do[ (*loop on each site, i = curr site *)
    
      AppendTo[ells, {}];
      AppendTo[nells, {}];
      (* sitmult is sum of # atoms on each sites, calculated as parsed 
         based on actual tc a-x model (not as passed in by sm and esm),
         and based on the sum of each site in 1st pc only. *)
      AppendTo[sitmult, 0];
      (*tpos={};*)
      (* surely can do in tigher loop code *)
      
      Do[(* loop over j = curr pc, looking at curr site i (M1, etc) *)
      
        tls=aids[[j]];
        If[i == 1, AppendTo[pcs, tls[[1]]]];
        If[i == 1, AppendTo[pcsitcomp, {}]];
        
        Do[ (* loop on k = curr pos in aid ln of curr pc j *) 
        
          (* {en, 1, 3, x(Mg,M1), 1, x(Mg,M2), 1, x(Si,T), 1/2} *)
          (*Print["nas5"];*)
          tsit = GetSIPs[tls[[k]]][[2]];
          (*Print["nas6"];*)
          If[StringMatchQ[sitnames[[i]],tsit],
            tel = GetSIPs[tls[[k]]][[1]];
            tnel = tls[[k+1]];
            AppendTo[pcsitcomp[[j]], {tsit, tel, tnel}];
            If[j == 1, sitmult[[i]] = sitmult[[i]] + tnel];
            If[FirstPosition[ells[[i]], tel, {0}][[1]] == 0, 
              AppendTo[ells[[i]], tel]; AppendTo[nells[[i]], tnel]
            ];
            tpos = FirstPosition[ells[[i]], tel, {0}][[1]];
            If[tnel > nells[[i, tpos]], nells[[i, tpos]] = tnel];
          ];
          ,{k,4,Length[tls]-1,2}
        ];
        ,{j,Length[aids]}
      ];
      ,{i,Length[sitnames]}
    ](* end loop through aids to get pc site mixing *)
  ];(* end block *)
  (*Print["nas7"];*)
  (* GCD requires ints or rationals....  Add a check for IntegerQ or RationalQ *)
  (* For rational numbers Subscript[r, i], GCD[Subscript[r, 1],Subscript[r, 2],\[Ellipsis]] 
     gives the greatest rational number r for which all the Subscript[r, i]/r are 
     integers. *)
  aidnorm = Apply[GCD, Flatten[nells]];
  (*Print["nells before =",Flatten[nells]," with aidnorm = ",aidnorm];*)
  If[aidnorm>1,
    (* ToDo: Add test if no rationals in this case. *)
    Print[Style["ERROR in NewAnalSites for phase "<>ph<>
      ": Trying to recast sites with GCD > 1! aidnorm = "<>
      ToString[aidnorm,InputForm]<>
      "....  Continuing, but setting aidnorm=1, so check output."<>
      "  nells="<>ToString[nells,InputForm],12,Red]];
    aidnorm=1; (*2022-04-10 to have OL work normally. A single mix site with sm>1 and no fractions.*)
    (*Print[nells//InputForm];*)
  ];
  tnells = nells;
  (*Print["nas8"];*)
  nells = Map[(#/aidnorm) &, tnells, {2}];
  (*Print["aidnorm = ",aidnorm];*)
  (*Print["nells after =",Flatten[nells]];*)
  (*Print["nas9"];*)
  (* 
    Now normalize #atoms by aidnorm to convert tcax site
    comp (with effective site logic) to tdax site format
    (with increased # atom to make whole numbers for all
    atoms, in the cases were eff reduces to non-whole numbers).
    Interestingly, neither case is true #atoms in cases
    where eff site mult reduced. tpcsitcomp will be tcax
    fmt, and pcsitcomp will hold tdax fmt. So, after this
    normalization, should implement setting ncmp and xcmp
    arrays to gen site mixing tables.
  *)
  tpcsitcomp = pcsitcomp; 
  pcsitcomp = 
    Map[({#[[1]], #[[2]], #[[3]]/aidnorm}) &, tpcsitcomp, {2}];
  
  tvec = sitmult;
  sitmult = Map[(#/aidnorm) &, tvec];
  
  (* TODO: Generate ncmp and xcmp table data arrays from the following *)
  (* tpcsitcomp: { {{"M1","Mg",1},{"M2","Mg",1},{"T","Si",1/2}},..other pc's..} *)
  (* pcsitcomp: {{{"M1","Mg",4},{"M2","Mg",4},{"T","Si",2}},..other pc's... } *)
  (* can use new sm and esm arrays args now passed in to get correct ncmp atoms *)
  
  ncmptbl=NCompTbl[pcs,sitnames,sitels,tpcsitcomp][[2]];
  (*tbls1=TblToTextTbl[ncmptbl[[1]],smtblstr];*)
  tbls1=TblToTextTbl2[sitnames,sitels,ncmptbl[[1]],smtblstr];
  ncmptbl=NCompTbl[pcs,sitnames,sitels,pcsitcomp][[2]];
  (*tbls2=TblToTextTbl[ncmptbl[[1]],
    "Normalized sites for TD logic, compatible with ax model via (X^n)^("<>ToString[aidnorm,InputForm]<>")."
    ];*)
  tbls2=TblToTextTbl2[sitnames,sitels,ncmptbl[[1]],
    "Normalized sites for TD logic, compatible with ax model via (X^n)^("<>ToString[aidnorm,InputForm]<>")."
    ];
  
  (*Block[{writelogoutput, outstyle},
    writelogoutput = True;
    If[writelogoutput,
      outstyle = InputForm;
      Print["Site Analysis"];
      Print["site names: "];
      Print[sitnames // outstyle];
      Print["normalized for td sitmult: "];
      Print[sitmult // outstyle];
      Print["aid normalizer: "];
      Print[aidnorm // outstyle];
      Print["ells: "];
      Print[ells // outstyle];
      Print["nells: "];
      Print[nells // outstyle];
      Print["pcs:"];
      Print[pcs // outstyle];
      Print["pcsitcomp"];
      Print[pcsitcomp //TableForm];
    ];
  ];  *)
  
  Return[{sitnames, sitmult, aidnorm, ells, nells, pcs, pcsitcomp, {tbls1,tbls2}}];
];


(* ::Subsubsection::Closed:: *)
(*SetAXDicts*)


SetAXDicts[mdl_?ListQ]:=Module[{},

Return[0];
];


(* ::Subsubsection::Closed:: *)
(*LoadTCAXModels*)


LoadTCAXModels[tcaxfile_?StringQ, tcdsfile_?StringQ, dsv_?StringQ, 
               dsfmt_?StringQ, skipval_:Null, opts:OptionsPattern[]]:=
Module[
  {doxtra, xtrasfile, exclnlines, exhdr, exverbatims, exoriginal, 
   expures, exaxdata, taxdata,
   clnlines, hdr, verbatims, original, nmdls, models={}, axdata={}, 
   pures={}, mades={}, madefroms={}, mdlnames={}, dsabbls, dstypels
  },
  (* SETUP OPTIONS *)
  doxtra=OptionValue[LoadTCAXModels,"LoadExtras"];
  xtrasfile=OptionValue[LoadTCAXModels,"ExtrasFile"];
  
  (* READ DATABASE *)
  Block[{mstrdat, infarr, abbls, spcnumls, scomparr, sdatarr,
    scarr2, rcarr, formls},
    {mstrdat, infarr, dsabbls, dstypels, spcnumls, scomparr, sdatarr,
     scarr2, rcarr, formls} = ReadLocalTCDS[tcdsfile,dsv,dsfmt,0.0];
  ];
  
  (* READ EXTRA MODELS FILE *)
  If[TrueQ[doxtra],
    {exclnlines,exhdr,exverbatims,exoriginal}
      = GetCleanTCAxLines[xtrasfile];
    exaxdata = ParseCleanTCAxLines[exclnlines,dsv];
    expures  = exaxdata[[-1]];
    exaxdata = Delete[exaxdata,-1];
    Assert[Length[exaxdata]==Length[exverbatims]];
  ];
  
  (* READ AX FILE *)
  (* each axdata list has this: {ph,npc,ver,pcnames,vars,guesses,ppnls,
     xcstype,wlist,alist,nsitex,sitexls,aidlns,makelns,dqflns,checklns},
     except last list which is just list of pure phases on 1 line *)
  {clnlines,hdr,verbatims,original}=GetCleanTCAxLines[tcaxfile];
  axdata=ParseCleanTCAxLines[clnlines,dsv];
  pures=axdata[[-1]];
  axdata=Delete[axdata,-1];
  Assert[Length[axdata]==Length[verbatims]];
  
  (* COMBINE EXTRA AX AND AX FILE *)
  If[TrueQ[doxtra],
    If[Length[exaxdata]>0,
      Do[AppendTo[axdata,exaxdata[[i]]],{i,Length[exaxdata]}];
      Do[AppendTo[verbatims,exverbatims[[i]]],{i,Length[exverbatims]}];
      pures=DeleteDuplicates[Flatten[PrependTo[pures,expures]]];
    ];
  ];
  nmdls=Length[axdata];
  Print["Parsed "<>ToString[nmdls]<>" models."];
  
  (* LOOP OVER MODELS, IN A BLOCK *)
  Block[{ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,nsitex,sitemult,esitemult,
    sitexls,aidlns,makelns,dqflns,checklns,info,i,j,domake,mkln,smplmkln,
    npcname,mfromrepls,renamels,tstr,mdlmklns,mdlmades,mdlmadefroms},
    
    (* Loop over all models; Process make lines and  do string replacements on names *)
    Do[
    
      mdlmklns={};
      mdlmades={};
      mdlmadefroms={};
      
      
      {ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,nsitex,sitemult,
       esitemult,sitexls,aidlns,makelns,dqflns,checklns,info} = axdata[[i]];
       
       
       (*info has form: {{"skip",False},{"name",ph},{"ref",""},{"note",""},{"ext",False};*)
       (* all models that make it here from ParseCleanTCAxLines have skip=False. *)
       
       
      renamels={};
      (*Print["PHASE: "<>ph<>" ================================================================"];*)
      Print["PHASE: "<>ph];
      
      (*Loop over each pc, process its make lines if it has one, make replacements*)
      Do[
        If[Length[makelns[[j]]]>1, domake=True, domake=False];
        
        If[TrueQ[domake],
          (*
          makelns have norm inserted in front of sp lacking a key,
          add pc being made to start of mkln. Add dqf to end of
          mkln. I also make a simple make line w/o the dqf to simplify
          later TD printouts.
          *)
          smplmkln=Join[{pcnames[[j]]},Insert[makelns[[j]],dsv,1]];
          mkln=Join[{pcnames[[j]]},Insert[makelns[[j]],dsv,1],dqflns[[j]]];
          (* get new pc name if it needs increment, and names of 
             make froms in form of a replacement list. *)
          {npcname,mfromrepls}=ProcessMakeVecTC350[mkln,"-"];
          mfromrepls=Flatten[mfromrepls];
          AppendTo[mdlmades,npcname];
          
          Do[
            AppendTo[mdlmadefroms,Extract[mfromrepls[[k]],2]],{k,Length[mfromrepls]}
          ];
          
          (*If npcname diff than pcnames[[j]], then add to a replacement list 
            to run through appropriate data lists*)
          If[ !StringMatchQ[npcname,pcnames[[j]]],
            AppendTo[renamels,pcnames[[j]]->npcname];
            
            (*Make replacements in ppnls, alist, aidlns, makelns; do wlist outside 
              this loop after all replace names known.*)
              
            (*StringReplace is ok to use in this loop because comparing specific
            strings in pcnames with specific string that is to replace it; so,
            should not get substring matches within original name being replaced. 
            However, could replace StringReplace func here without using a do loop
            to do so, so it is an easy change. Got rid of StringReplace on 
            2022-05-27*)
            
            (*format of ppnls[[j,1]] is p(pcname), so if get rid of StringReplace
            for ppnls, need to use GetSIP to extract pc name, then put back
            together. *)
            Print["Replacement: "<>ppnls[[j,1]]<>"->"<>"p("<>npcname<>")"];
            ppnls[[j,1]]="p("<>npcname<>")";
            (*ppnls[[j,1]]=StringReplace[ppnls[[j,1]],pcnames[[j]]->npcname]; *)
            (*Print[{ppnls[[j,1]],pcnames[[j]],npcname}];*)
           
            If[Length[alist]>0,
              Print["Replacement: "<>alist[[j,1]]<>"->"<>npcname];
              alist[[j,1]]=npcname;
              (*alist[[j,1]]=StringReplace[alist[[j,1]],pcnames[[j]]->npcname]; *)
            ];
            Print["Replacement: "<>aidlns[[j,1]]<>"->"<>npcname];
            aidlns[[j,1]]=npcname;
            (*aidlns[[j,1]]=StringReplace[aidlns[[j,1]],pcnames[[j]]->npcname]; *)
            (*Print[{aidlns[[j,1]],pcnames[[j]],npcname}];*)
            
            Print["Replacement: "<>mkln[[1]]<>"->"<>npcname];
            mkln[[1]]=npcname;
            (*mkln[[1]]=StringReplace[mkln[[1]],pcnames[[j]]->npcname];*)
            (*Print[{mkln[[1]],pcnames[[j]],npcname}];*)
            
            Print["Replacement: "<>smplmkln[[1]]<>"->"<>npcname];
            smplmkln[[1]]=npcname;
            (*smplmkln[[1]]=StringReplace[smplmkln[[1]],pcnames[[j]]->npcname];*)
            (*Print[{smplmkln[[1]],pcnames[[j]],npcname}];*)
            
            (* do last as pcnames[[j]] is used above *)
            pcnames[[j]]=npcname;
          ];
          
          (* Deal with mkln and smplmkln by making replacements
             of made from names *)
          If[ Length[mfromrepls]>0,
            Do[
              (* Don't use string replace here, as it can match on partial
                 string. E.g. kprg: has both pa->pa- and parg->parg-, which
                 resulted in parg being replaced by pa-rg 
                 Do loop 2022-05-26 *)
              Do[
                If[StringMatchQ[mkln[[k]],mfromrepls[[m,1]]],
                  mkln[[k]]=mfromrepls[[m,2]];
                  smplmkln[[k]]=mfromrepls[[m,2]];
                ];
                ,{m,Length[mfromrepls]}
              ];
              (*
              mkln[[k]]=StringReplace[mkln[[k]],mfromrepls];
              smplmkln[[k]]=StringReplace[smplmkln[[k]],mfromrepls];
              *)
              
              ,{k,6,6+3*(mkln[[4]]-1),3}
            ];
          ];
          
          (*AppendTo[mdlmklns,mkln];*)(*change to smplmkln*)
          AppendTo[mdlmklns,smplmkln];
          ,
          AppendTo[mdlmklns,{}];  
        ]; (*If[TrueQ[domake]*)
        
        ,{j,npc}
      ]; (*end loop over each pc of model*)
      
      (*Make replacements for pc names in wlist that need changed *)
      If[Length[renamels]>0,
        (*wlist*)
        Do[ 
          (*2022-05-27: Replaced StringReplace with do loop*)
          Do[
            If[StringMatchQ[wlist[[j,1]],renamels[[k,1]]],wlist[[j,1]]=renamels[[k,2]]];
            If[StringMatchQ[wlist[[j,2]],renamels[[k,1]]],wlist[[j,2]]=renamels[[k,2]]];
            ,{k,Length[renamels]}
          ];
          (*
          wlist[[j,1]]=StringReplace[wlist[[j,1]],renamels];  
          wlist[[j,2]]=StringReplace[wlist[[j,2]],renamels]; 
          *)
          ,{j,Length[wlist]}
        ];
      ];
      
      (* ToDo: want to store mklns or smplmklns for later use, as makelns is just
         slight modification of tc-ax make lines.*)
      
      AppendTo[models,{ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,
        nsitex,sitemult,esitemult,sitexls,aidlns,mdlmklns,dqflns,checklns,mdlmades,mdlmadefroms,verbatims[[i]],info}];
      
      AppendTo[mades,mdlmades];
      AppendTo[madefroms,mdlmadefroms];
      
      (*Print["Starting ax dict construction for phase "<>ph];*)
      (* Generate ax dicts *)
      axmdl[ph]=ph;
      (*In future, should be able to skip NewAnalSites if info[[5,2]]=True (an EXT model).
      For now it is evidently working with EXT models*)
      (*axsiteanal[ph]=NewAnalSites[models[[i,14]],models[[i,15]]];*)
      axsiteanal[ph]=NewAnalSites[ph,sitexls,aidlns,sitemult,esitemult];
      (*Print["ok1"];*)
      axref[ph]=info[[3,2]];
      (*Print["ok2"];*)
      axnpc[ph]=npc;
      axipc[ph]=pcnames;
      axopc[ph]={};
      axdpc[ph]={};
      axnsx[ph]=nsitex;
      axsxnames[ph]=Map[#[[1]]&,sitexls];
      (*Print["ok3"];*)
      axtype[ph]=xcstype;
      axvar[ph]=vars;
      axvv[ph]=Map[vars[[#]]->guesses[[#]]&,Range[Length[vars]]];
      (*Print["ok4"];*)
      axovar[ph]={};
      axpcrepls[ph]={}; (*ToDo: Either delete dict or figure out how to populate*)
      axmakes[ph]=mdlmklns;
      axdqfs[ph]=dqflns;
      
      (*Below are still to implement*)
      axsm[ph]={}; (*axsm["G","M2"] = 2; ToDo: Figure out how to get. From aidlns would give eff site mult.*)
      axesm[ph]={};
      axsiteels[ph]={};
      axnmtx[ph]={};
      axxmtx[ph]={};
      Do[
        axncmp[ph,pcnames[[ii]]]={};
        axxcmp[ph,pcnames[[ii]]]={};
        ,{ii,npc}
      ];
      (*Print["ok5"];*)
      (* --- *)
      
      Do[
        axwb[ph,wlist[[ii,1]],wlist[[ii,2]]]=wlist[[ii,3;;]];
        ,{ii,1,Length[wlist]}
      ];
      (*Print["ok6"];*)
      (*no axwt to build for tc ax models*)
     
      Do[
        If[Length[alist]>=ii,
          axaterm[ph,pcnames[[ii]]]=alist[[ii,2;;]],
          axaterm[ph,pcnames[[ii]]]={}
        ];
        ,{ii,npc}
      ];
      (*Print["ok7"];*)
      
      (*taxis*)
      axppnstx[ph]=Map[TCExpLoop[#[[2;;]]]&,ppnls]; 
      (*Print["ok8"];*)
      axsxstx[ph]=Map[TCExpLoop[#[[2;;]]]&,sitexls];
      (*Print["ok9"];*)
      axaidlns[ph]=aidlns;
      axaidstx[ph]=Map[(TCAidExpLoop[#[[2;;]]//.Thread[axsxnames[ph]->axsxstx[ph]]])&,aidlns];
      (*Print["ok10"];*)
      Do[  
        axppntx[ph,pcnames[[ii]]]=axppnstx[ph][[ii]];
        axaidtx[ph,pcnames[[ii]]]=axaidstx[ph][[ii]];
        axaidln[ph,pcnames[[ii]]]=aidlns[[ii]];
        ,{ii,npc}
      ];
      (*Print["ok11"];*)
      
      Do[axsxtx[ph,sitexls[[ii]][[1]]]=axsxstx[ph][[ii]],{ii,Length[axsxstx[ph]]}];
      (*Print["ok12"];*)
      
      ,{i,nmdls}
    ];
    (*END DO loop over ax models*)
    Print["Done looping over ax models. ax dicts created."];
  ];(*End Block*)
  
  (*
    models: nested list of ax model data, with elements in this order:
      {ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,
        nsitex,sitexls,aidlns,mdlmklns,dqflns,checklns,verbatims[[i]],info}
    dsabbls = list of dataset species in tc-dsxx.txt file.
    pures = list of pure phases from the tc-axmodelfile.txt
    mades = pc members of ax models that were made; will have tdtype=make
    madefrom = species made from for make phases; minus phases in my
      td- files.
  *)
  (*Do[
    NewAnalSites[models[[i,14]],models[[i,15]]];
    ,{i,nmdls}
  ];*)
  
  Return[{models,dsabbls,pures,mades,madefroms}];
];


(* ::Subsubsection::Closed:: *)
(*PredetermineRefs*)


PredetermineRefs[mdls_?ListQ]:=Module[{ph,npc,ver,pcnames,vars,guesses,ppnls,
  xcstype,wlist,alist,nsitex,sitemult,esitemult,sitexls,aidlns,mdlmklns,
  dqflns,checklns,mdlmades,mdlmadefroms,verbatims,info,strls,skipstr,namestr,
  refstr,notestr,ts,refsusedbyax={}},
  
  Do[
    {ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,nsitex,sitemult,esitemult,sitexls,
     aidlns,mdlmklns,dqflns,checklns,mdlmades,mdlmadefroms,verbatims,info} = mdls[[i]];
  
    (*info has form: {{"skip",False},{"name",ph},{"ref",""},{"note",""}};*)
    (* all models that make it here have skip=False. *)
    (*Print["info = ",info];*)
    (*strls=ReadList[StringToStream[info[[3,2]]],String];*)
    If[(*Length[strls]>0*)Length[info[[3,2]]]>0,
      AppendTo[refsusedbyax,{ToUpperCase[info[[2,2]]],info[[3,2]]}];
      (*in ngen2`, would have an axrefs dict to add to; don't have that in ngen,
        so I added axrefs[_]:={} to tc2td.wl specifically for this! *)
      (*Global`axrefs[ToUpperCase[info[[2,2]]]]=strls;*)
      axrefs[ToUpperCase[info[[2,2]]]]=info[[3,2]];
    ];
    ,{i,Length[mdls]}
  ];
  Return[refsusedbyax];
];


(* ::Subsubsection::Closed:: *)
(*WriteTDAXFile*)


WriteTDAXFile[ofile_?StringQ, mdls_?ListQ, dsabbls_?ListQ, pures_?ListQ, 
  mades_?ListQ, madefroms_?ListQ, opts:OptionsPattern[] ]:=
Module[{usef3,tcax,tcds,tcdsver,ATmstr,GetTypeKey,WritePCs,WriteDivider,mstr,dsmstr,
        mmstr,mfmstr,hdrstr,refstr,timestr,mstrmadefroms,ostrm,tstr,tls,tls2,
        tkeys,allrefs={},allusedrefs={},refsusedbyax={},TYPEKEY="MINERAL",
        solnmades={}, solnmadefroms={}, writtenpcs={}, remdbspcs={}},
  usef3=OptionValue[WriteTDAXFile,"UseF3"];
  tcds=OptionValue[WriteTDAXFile,"TCDSFile"];
  tcax=OptionValue[WriteTDAXFile,"TCAXFile"];
  tcdsver=OptionValue[WriteTDAXFile,"TCDSver"];
  
  (* Local Funcs To Append To mstr *)
  ATmstr[str_?StringQ]:=AppendTo[mstr,str];
  ATmstr[ls_?ListQ]:=Map[AppendTo[mstr,#]&,ls];
  
  (*
    MINERAL: make, norm, land, bw, liq
    GAS: h2o, co2, cs
    AQU: aq
  *)
  GetTypeKey[pct_?StringQ]:=Block[{},
    If[StringMatchQ[pct,{"make","norm","land","bw","liq"}],
      Return["MINERAL"];
    ];
    If[StringMatchQ[pct,{"h2o","co2","cs"}],
      Return["GAS"];
    ];
    If[StringMatchQ[pct,{"aq"}],
      Return["AQU"];
    ];
  ];
  WritePCs[ls_?ListQ,tdsv_?StringQ]:=Block[{i,tmstr,pctk},
    Do[
      If[pcLoaded[ls[[i]],tdsv],
        tmstr=GetPCInTDFormat[ls[[i]],tdsv];
        pctk=GetTypeKey[tdtype[ls[[i]],tdsv]];
        If[!StringMatchQ[TYPEKEY,pctk],
          If[StringMatchQ[pctk,"MINERAL"],
            ATmstr["*** MINERAL DATA ***"];
            TYPEKEY="MINERAL";
          ];
          If[StringMatchQ[pctk,"GAS"],
            ATmstr["*** GAS DATA ***"];
            TYPEKEY="GAS";
          ];
          If[StringMatchQ[pctk,"AQU"],
            ATmstr["*** AQU DATA ***"];
            TYPEKEY="AQU";
          ];
        ];
        ,
        tmstr="! MISSING pc "<>ls[[i]];
        Print["ERROR - missing pc when writing dset: "<>ls[[i]]<>"  "<>tdsv];
      ];
      Map[ATmstr[#]&, tmstr];
      ,{i,Length[ls]}
    ]
  ];
  (* startstr should be 2 char; usually '! ' *)
  WriteDivider[start2char_?StringQ, char_?StringQ, n_?IntegerQ]:=Block[{s},
    s = start2char;
    Do[s = s <> char, {Range[Ceiling[(96 - 2)/2]]}];
    Do[ATmstr[s],{i,1,n}];
  ];
  (* ***** *)
  
  If[Global`fctdbg==True, Print["Entering WriteTDAXFile.."]];
  mstr=dsmstr=mmstr=mfmstr=hdrstr=refstr={};
  mstrmadefroms={};
  (* ***** *)
  
  (*
     ToDo:
   - Write tc-dsxx file name with footer
   - implement adding corresponding states parameters automatically
   - reorder to locate database species in this order
     1. pures listed as pures in tc-ax, H2O, CO2 first
     2. minus phases
     3. soln models (add comments to ferric-bearing species?)
     4. buffers
     4. all remaining pure min species in dataset that are not already
        added through one of the groups above
     6. aqueous species
   - Write comments & modifications
  *)
  
  (*Step 1. Gather all pc's to be written if diff categories,
  determine which members of hp dbs are not set to be written.*)
  Block[{ph,pcnames,mdlmades,mdlmadefroms,info,allmdlpcs2write={},dbsnotwritten={}},
    Do[
      (* sets 22 variables to mdls[[i]]; *)
      (*{ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,
       nsitex,sitemult,esitemult,sitexls,aidlns,mdlmklns,dqflns,checklns,mdlmades,mdlmadefroms,
       verbatims,info} = mdls[[i]];*)
      ph=mdls[[i,1]];
      pcnames=mdls[[i,4]];
      mdlmades=mdls[[i,19]];
      mdlmadefroms=mdls[[i,20]];
      AppendTo[allmdlpcs2write,Flatten[{mdlmades,mdlmadefroms,pcnames}]];
      AppendTo[solnmades,{ph,mdlmades}]; 
      AppendTo[solnmadefroms,{ph,mdlmadefroms}];
      AppendTo[mstrmadefroms,mdlmadefroms];
      ,{i,Length[mdls]}
    ];
    allmdlpcs2write=Flatten[allmdlpcs2write];
    mstrmadefroms=Flatten[mstrmadefroms];
    tls=dsabbls;
    Print[{"tls",tls}];
    Print[{"allmdlpcs2write",allmdlpcs2write}];
    Do[tls=DeleteCases[tls,allmdlpcs2write[[i]]],{i,Length[allmdlpcs2write]}];
    Do[If[!StringMatchQ[tdtype[tls[[i]],tcdsver],"aq"],AppendTo[remdbspcs,tls[[i]]]],{i,Length[tls]}];
    Print[{"Remaining dbs pcs to write (non-aq)",remdbspcs}];
  ];  
  
  
  
  
  
  
  (* header for tcds6 *)
  ATmstr[GetHP11TDHeader[]];
  WriteDivider["  ","==",1];
  
  ATmstr[
    "  File "<>ofile<>"\n"<>
    "  Converted by D.K. Tinkham (tc2td.wl) to Theriak-Domino format from files "
    <>FileNameTake[tcds]<>" and \n  "<>
    FileNameTake[tcax]<>" on "<>DateString[]<>".\n\n\n\n"<>
    "  FILE UPDATES/ADDITIONS:\n"<>
    "  ver 01: "<>DateString[]<>". DKT\n"<>
    "  ver 02: \n\n\n"
  ];
  
  
  WriteDivider["!=","==",2];
  ATmstr["! REFERENCES"];
  tls=PredetermineRefs[mdls];
  (*tls is of form { {WM,{WPH14,HPx-eos}}, {GRT,{WPH14}} }*)
  
  tkeys=Sort[DeleteDuplicates[Flatten[Map[#[[2]]&,tls]]]];
  Print[tkeys];
  refstr=GetTC2TDAXReferences[];
  Print[refstr];
  refsusedbyax={};
  Do[
    tls=Flatten[Select[refstr,#[[1]]==tkeys[[i]]&]];
    If[Length[tls]>1,AppendTo[refsusedbyax,tls];];
    ,{i,Length[tkeys]}
  ];
  Do[
    ATmstr["  "<>refsusedbyax[[i,1]]<>":"];
    Do[
      ATmstr[If[j==2,"      ","        "]<>refsusedbyax[[i,j]]];
      ,{j,2,Length[refsusedbyax[[i]]]}
    ];
    ,{i,Length[refsusedbyax]}
  ];
  WriteDivider["!=","==",2];
  ATmstr["\n\n"];
  
  
  (* 1. Write pures passed to func in pures argument. *)
  WriteDivider["!=","==",1];
  ATmstr["!  1. PURE SPECIES LISTED IN CONVERTED TC-AX FILE"];
  WriteDivider["!=","==",1];
  (*ATmstr[
    "! The following are pure phases to be considered in calculations. Only a subset is \n"<>
    "! considered here. The full HP 2011 dataset from file "<> tcds <>" is included at \n"<>
    "! the bottom of the file for you, but is currently deactivated. Copy/paste members to \n"<>
    "! this section if you want."
  ];*)
  TYPEKEY="RESET";
  WritePCs[DeleteDuplicates[Flatten[pures]],tcdsver];
  
  
  
  (*Write solution models and their pc's*)
  Block[{ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,nsitex,sitemult,
         esitemult,sitexls,aidlns,mdlmklns,dqflns,checklns,mdlmades,mdlmadefroms,
         verbatims,info,
         llen=96,ts,ts2,tlen,
         ansit,
         sitnames,
         sitmult,
         ells,
         nells,
         aidnorm,
         pcs,
         pcsitcomp,
         reqdpcs,
         tmstr,
         ntxttbls
         },
    (*1. Find all mdlmadefroms to put them first in td-ax file*)
    Do[
      
      (* sets 22 variables to mdls[[i]]; *)
      {ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,
       nsitex,sitemult,esitemult,sitexls,aidlns,mdlmklns,dqflns,checklns,mdlmades,mdlmadefroms,
       verbatims,info} = mdls[[i]];
      AppendTo[solnmades,{ph,mdlmades}]; 
      AppendTo[solnmadefroms,{ph,mdlmadefroms}];
      AppendTo[mstrmadefroms,mdlmadefroms];
      ,{i,Length[mdls]}
    ];
    (* info has form: {{"skip",False},{"name",ph},{"ref",""},{"note",""},{txttbl1,txttbl2}; *)
    (* all models that make it here from ParseCleanTCAxLines have skip=False. *)
    
    (* write madefroms minus phases *)
    mstrmadefroms=DeleteDuplicates[Flatten[mstrmadefroms]];
    ATmstr["\n\n\n"];
    WriteDivider["!=","==",1];
    ATmstr["!  2. MINUS SPECIES: Only used in COM statements; not considered for phase stability."];
    WriteDivider["!=","==",1];
    
    (*
    ATmstr["\n"];
    ATmstr["!  The following are 'minus' phases used in COM statements in solution models.\n"<>
           "!  They are only used for calculating G's of COM phases, and are not considered\n"<>
           "!  for phase stability.\n"
    ];
    *)
    
    ATmstr["*** MINERAL DATA ***"];
    Do[
      If[pcLoaded[mstrmadefroms[[jj]], tcdsver],
        tmstr = GetPCInTDFormat[mstrmadefroms[[jj]], tcdsver, "UseF3Comps"->usef3],
        tmstr = "!  MISSING " <> mstrmadefroms[[jj]]
      ];
      ATmstr[tmstr];
      ,{jj,Length[mstrmadefroms]}
    ];
  
  
    (*2. Loop through to write solutions.*)
    ATmstr["\n\n\n"];
    WriteDivider["!=","==",1];
    ATmstr["!  3. SOLUTION MODELS AND THEIR PHASE COMPONENTS"];
    WriteDivider["!=","==",1];
    ATmstr["!   In general, the phase components for the a-x model appear immediately \n"<>
           "!   before the mixing model.\n"
    ];
    Do[
    
      (* sets 22 variables to mdls[[i]]; *)
      {ph,npc,ver,pcnames,vars,guesses,ppnls,xcstype,wlist,alist,
       nsitex,sitemult,esitemult,sitexls,aidlns,mdlmklns,dqflns,checklns,mdlmades,mdlmadefroms,
       verbatims,info} = mdls[[i]];
       
      (*info has form: {{"skip",False},{"name",ph},{"ref",""},{"note",""},{txttbl1,txttbl2}};*)
      
      (* added extra flag {"ext",False}, which should now be at idx 5, so tbls should now be at idx 6. *)
      (*info has form: {{"skip",False},{"name",ph},{"ref",""},{"note",""},{"ext",False},{txttbl1,txttbl2}};*)
      
      (* ***** *)
      (* write ax model header *)
      ts="!";
      Do[ts=ts<>"-",{Range[llen-1]}]; ATmstr[ts];
      ts="!"; ts2=" AX MODEL " <> ToUpperCase[ph];
      tlen=llen - 1 - StringLength[ts2];
      Do[ts = ts <> "-", {Range[tlen]}];
      ts = ts <> ts2;
      ATmstr[ts];
      (*ATmstr[refline];*)
      
      (* MAKE MOD HERE TO LOOP IN CASE MORE THAN ONE REF *)
      ts="";
      Do[ts = ts<>info[[3,2,j]]<>"  ", {j,Length[info[[3,2]]]}];
      ATmstr["!   Ref:   "<>ts<>"\n"<>
             "!   note:  "<>info[[4,2]]<>"\n"<>
             "!   entry: Converted by tc2td.wl on " <> DateString[Date[]] <>"\n"<>
             "!   Converted from HPx-eos website files "<>FileNameTake[tcds]<>
             " & "<>FileNameTake[tcax]<>".\n"<>
             "!   Required dataset species for COM's: "<>
              TextString[DeleteDuplicates[mdlmadefroms]]<>"\n!"
      ];
      ATmstr[axsiteanal[ph][[-1]][[1]]];(*write ncomptbl, stored in axsiteanal dict*)
      (*Determine if any sites are effective sites; if so print tdtbl*)
      If[ AnyTrue[MapThread[#1[[2]]/#2[[2]]&,{esitemult,sitemult}], TrueQ[#!=1]&],
        ATmstr["!   ------"];
        ATmstr[axsiteanal[ph][[-1]][[2]]];
      ];
      
      
      WriteDivider["!-","--",1];
      ATmstr["*** MINERAL DATA ***"];
    
      (* write reqd pc's (normal ds species and made pc's; madefroms minus 
      phases already written above) *)
      reqdpcs=DeleteDuplicates[pcnames~Join~mdlmades];
      AppendTo[writtenpcs,reqdpcs];
      Do[
        If[pcLoaded[reqdpcs[[jj]], tcdsver],
          tmstr = GetPCInTDFormat[reqdpcs[[jj]], tcdsver, "UseF3Comps"->usef3],
          tmstr = "!  MISSING " <> reqdpcs[[jj]]
        ];
        ATmstr[tmstr];
        ,{jj,Length[reqdpcs]}
      ];
      WriteDivider["!."," .",1];
    
      (* ***** *)
      (*2022-04-04: need to send sm and esm vecs*)
      ansit = NewAnalSites[ph,sitexls, aidlns, sitemult, esitemult];
      {sitnames, sitmult, aidnorm, ells, nells, pcs, pcsitcomp, ntxttbls} = ansit;
      ATmstr["*** SOLUTION DATA ***"];
      
      (* write name and mixing sites line *)
      ts = ToUpperCase[ph] <> "    ";
      If[!TrueQ[info[[5,2]]],
        ts = ts <> "(-SITE,MARGULES)";
        ,
        ts = ts <> "(-EXT,MARGULES)";
      ];
      If[aidnorm != 1, ts = ts <> ToString[aidnorm, InputForm] <> "   "];
      If[!TrueQ[info[[5,2]]],
        ts = ts <> "   ";
        Do[ (* MARK - LOOP OVER Sites *)  
          If[j>1, ts = ts <> "- "];
          ts = ts <> sitnames[[j]] <> "(" <> ToString[sitmult[[j]]] <> "):";
          Do[
            ts = ts <> ells[[j, k]];
            If[k < Length[ells[[j]]], ts = ts <> ",", ts = ts <> " "];
            ,{k, 1, Length[ells[[j]]]}
          ];
          ,{j, 1, Length[sitnames]}
        ];
      ];
      ATmstr[ts];
      
      (* do ideal act lines + a terms for asf *)
      Do[(* j loop over pc's *)
        ts = "  " <> StringPadRight[ToString[pcs[[j]]], 14];
        
        
        
        If[!TrueQ[info[[5,2]]],
          Do[ (* k loop over sites sid *)
            Do[ (* ll loop over pcsitcomp[[j]] *)
              If[pcsitcomp[[j,ll,1]] == sitnames[[k]], 
                Do[ts = ts <> pcsitcomp[[j,ll,2]] <> ","
                  ,{lll,1,pcsitcomp[[j,ll,3]]}
                ]
              ]
              ,{ll,1,Length[pcsitcomp[[j]]]}
            ];
            ts = StringTake[ts, {1, Length[ts] - 2}];(* take of trail comma *)
            ts = ts <> " - ";
            ,{k,1,Length[sitnames]} 
          ];(* end loop over sites sid *)
          ts = StringTake[ts, {1, Length[ts] - 4}];(* take of trail " - " *)
          ,
          If[xcstype == "asf",
            ts = ts <> "    vterm    ";
          ];
        ];
        
        
        (* if asf, add a/v's *)
        If[xcstype == "asf",
          ts = ts <> "    ";
          ts = ts <> StringPadRight[ToString[alist[[j,2]]], 8];
          ts = ts <> StringPadRight[ToString[alist[[j,3]]], 8];
          ts = ts <> StringPadRight[ToString[alist[[j,4]]], 8];
          (* no models have cp in a terms yet, so don't add *)
          (*tstr = tstr <> "0";*)  (* commmented 2020-06-16 *)
        ];
        ATmstr[ts];
        ,{j,1,Length[pcs]}
      ];(* END MARK - loop over pc's *)
    
      (* add interaction parameters *) 
      If[xcstype == "asf" || xcstype == "sf",
        (*Print["itype = asf || itype = sf"];*)
        ts = "*** MARGULES PARAMETERS ***";
        ATmstr[ts];
        Do[
          ts = wlist[[j,1]] <> " - " <> wlist[[j,2]] <> "  ";
          ATmstr[ts];
          ts = StringPadRight["12", 9];
          (*Wh*)
          ts = ts <> StringPadLeft[N2S[ToEx[wlist[[j,3]]]*1000., 8], 9];
          (* Ws, -1 mult per TD convention *)
          ts = ts <> StringPadLeft[N2S[ToEx[wlist[[j,4]]]*(-1000.), 8], 9];
          (* Wv, don't mult by 1000, as kJ/kbar -> j/bar*)
          ts = ts <> StringPadLeft[N2S[ToEx[wlist[[j,5]]]*1.0, 8], 9];
          ATmstr[ts];
          ,{j,Length[wlist]}
        ];
      ];(* end if *)

      ATmstr["!END AX MODEL " <> ToUpperCase[ph]];
      WriteDivider["!-","--"];
      ATmstr["\n\n"];
    
      ,{i,Length[mdls]}
    ](* end loop through models *)
    
  ];(* end write solution models and their pc's block *)
  
  
  
  WriteDivider["!=","==",1];
  ATmstr["!  4. BUFFERS"];
  WriteDivider["!=","==",1];
  (*To Code: Get required pc's and make minus phases of them, then write.*)
  ATmstr[GetBufferStuff[]];
  ATmstr["\n\n"];
  
  
  
  WriteDivider["!=","==",1];
  ATmstr["!  5: REMAINING PURE SPECIES OF HP11 DATABASE NOT INCLUDED IN GROUPS 1 AND 3"];
  WriteDivider["!=","==",1];
  Block[{tmstr},
    Do[
      If[pcLoaded[remdbspcs[[jj]], tcdsver],
        tmstr = GetPCInTDFormat[remdbspcs[[jj]], tcdsver, "UseF3Comps"->usef3],
        tmstr = "!  MISSING " <> remdbspcs[[jj]]
      ];
      ATmstr[tmstr];
      ,{jj,Length[remdbspcs]}
    ];
  ];
  WriteDivider["!."," .",1];
  ATmstr["!  Nothing."];
  ATmstr["\n\n"];
  
  (*To Code: Remove all the pure unmodified species written so far*)
  
  
  
  
  WriteDivider["!=","==",2];
  ATmstr[
    "! All entries of HP database " <> FileNameTake[tcds] <> " to add above when needed.\n"<>
    "! I recommend adding pure unmodified species to section 5 unless they are components \n"<>
    "! of solution models or bufferes in sections 3 or 4. The aqueous species data are not \n"<>
    "! included here."
  ];
  WriteDivider["!-","--",2];
  ATmstr["*** MuNERAL DATA ***"];
  WritePCs[Flatten[dsabbls],tcdsver];
  
  timestr = DateString[];
  (* Write to file *)
  
  ostrm = OpenWrite[ofile];
  (*WriteString[ostrm,hdrstr];*)
  Map[WriteString[ostrm, # <> "\n"] &, mstr];
  Close[ostrm];
  
  If[Global`fctdbg==True, Print["..Leaving WriteTDAXFile"]];
  Return[mstr];
  
(*
!  ==============================================================================================
!  1. GROUP 1: PURE SPECIES LISTED IN CONVERTED TC-AX FILE
!  =============================================================================================

!  =============================================================================================
!  2. GROUP 2: MINUS SPECIES
!  ==============================================================================================
!  The following are 'minus' phases used in COM statements in solution models.
!  They are only used for calculating G's of COM phases, and are not considered
!  for phase stability. They are placed before the a-x models because the a-x
!  models use the COM feature to make pc's from these.
!  -----------------------------------------------------------------------------------------------

!  =============================================================================================
!  3. GROUP 3: SOLUTION MODELS AND THEIR PHASE COMPONENTS
!  ==============================================================================================
!  Additonal models not listed in the TC ax file are listed first
!  -----------------------------------------------------------------------------------------------

!  =============================================================================================
!  4. GROUP 4: BUFFERS
!  ==============================================================================================
!  
!  -----------------------------------------------------------------------------------------------

!  =============================================================================================
!  5. GROUP 5: REMAINING PURE SPECIES OF HP11 DATABASE NOT INCLUDED IN GROUPS 1 AND 3
!  ==============================================================================================
!  
!  -----------------------------------------------------------------------------------------------

!  =============================================================================================
!  6. GROUP 6: AQUEOUS SPECIES OF THE HP11 DATABASE (CANNOT USE AT THE MOMENT, AFAIK)
!  ==============================================================================================
!  
!  -----------------------------------------------------------------------------------------------
*)
];


(* ::Section:: *)
(*End*)


End[];

Protect[GetHP11TDHeader,GetTC2TDAXReferences,GetBufferStuff,MatchPatt,
GetEx,GenExListFromString,GenExListFromStringList,GetCleanTCAxLines,
ParseCleanTCAxLines,GenMinusSpeciesData,ProcessMakeVecTC350,
NewAnalSites,LoadTCAXModels,WriteTDAXFile];

SetAttributes[{GetHP11TDHeader,GetTC2TDAXReferences,GetBufferStuff,MatchPatt,
GetEx,GenExListFromString,GenExListFromStringList,GetCleanTCAxLines,
ParseCleanTCAxLines,GenMinusSpeciesData,ProcessMakeVecTC350,
NewAnalSites,LoadTCAXModels,WriteTDAXFile},
  ReadProtected];

Print["ngen`tc2td` loaded."];
EndPackage[];
