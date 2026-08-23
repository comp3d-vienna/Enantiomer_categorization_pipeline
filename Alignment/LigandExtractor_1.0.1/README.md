# LigandExtractor

LigandExtractor is a utility tool for extracting and enumerating all ligands of a protein-ligand
complex. When enumerating ligands, basic descriptors, such as the unique smiles, the
number of simple rotatable bonds, the heavy atom count, or the molecular weight are displayed.


## Table of Contents

- [LigandExtractor](#ligandextractor)
   * [Table of Contents](#table-of-contents)
   * [License](#license)
      + [Activation](#activation)
   * [Quick Start](#quick-start)
   * [Output](#output)
      + [Skip Reasons](#skip-reasons)
   * [Configuration](#configuration)
      + [Practical Tips](#practical-tips)
   * [Frequently Asked Questions](#frequently-asked-questions)
   * [Error Reporting](#error-reporting)


## License

LigandExtractor requires a license. Licenses are free for academic use. Non-academic users can apply
for a one-month evaluation period, after which a paid license is required. You can obtain a license
at:
https://software.zbh.uni-hamburg.de/


### Activation

After obtaining a license, you will need to activate LigandExtractor with that license. To do so,
open the license file, copy the content, and execute LigandExtractor as follows:
```
$ ./LigandExtractor --license <paste_your_license_here>
```


## Quick Start

LigandExtractor enumerates all ligands of a PDB or mmCIF file. Also, the ligands can be extracted.
To enumerate ligands use the '--enumeratePotentialLigands' option:
```
$ ./LigandExtractor --complex example_data/1r5g.pdb --enumeratePotentialLigands
```

To extract a ligand of your choice, use the '--ligandidentifier' option, define the ligand you
want to extract in the style of 'ResName_Chain_ResSeqNo', and define the output directory with
'--out':
```
$ ./LigandExtractor --complex example_data/1r5g.pdb --ligandidentifier AO1_A_501 --out results/
```
Please note that the output directory must always be defined.


## Usage

To run LigandExtractor, a protein file in the PDB or mmCIF format is mandatory and can be provided
with the option '--complex'. If you want to enumerate all ligands the option
'--enumeratePotentialLigands' is required. This will print out some information about the ligand
and the identifier (ResName_Chain_ResSeqNo) you can use to extract the ligand using the extraction
functionality of the tool. If you use the --ligandidentifier option, the output path is also
required ('--out').


## Output

LigandExtractor enumerates ligands with the following information:
pdb;name;skip_reason;hetcodes;nof_residues;usmiles;naomi_type;srotb;formula;heavy_atoms;mw

'pdb' is the PDB code of the given complex, 'name' is the ligand identifier to be used for
extraction, 'naomi_type' is an internal structure identification type of our software library.
Information on the 'skip_reason' can be found below.


### Skip Reasons

There are a number of so-called "Skip Reasons", which are originally used to skip molecules that
might not be handled correctly by our internal library 'NAOMI'.

A list of all skip reasons is given below:

| Skip Reason                           | Explanation                                                                        |
|---------------------------------------|------------------------------------------------------------------------------------|
| NotSkipped                            | ligand successfully parsed                                                         |
| LigandIsMetal                         | ligand is a metal ion                                                              |
| LigandIsUnknownMetal                  | ligand is a metal ion not listed in the HET record                                 |
| LigandContainsMetal                   | ligand contains a metal ion (e.g. HEM)                                             |
| LigandIsIncomplete                    | ligand has missing atoms (according to REMARK470/REMARK610)                        |
| LigandIsIncompleteNaomi               | ligand is split into multiple parts upon parsing                                   |
| LigandHasUnknownHetcode               | ligand has a FORMUL record but no formula string                                   |
| LigandHasUnexpectedFormula            | ligand does not correspond to the molecular formula according to the FORMUL record |
| LigandIsPolymerWithDifferentChains    | ligand is a polymer of residues from different chains                              |
| LigandIsPolymerWithChainBreak         | ligand is a polymer but seems to be part of a disconnected larger chain            |
| LigandIsPolymerWithMissingCoordinates | ligand is a polymer with missing coordinates                                       |
| LigandIsCovalentOrIncompletePolymer   | ligand has missing connections according to the LINK records                       |
| LigandIsCovalentOrHeteroResidue       | ligand is part of the protein                                                      |
| LigandWasSkippedByAltloc              | ligand atoms with different AltLoc skipped                                         |
| LigandCouldNotBeBuilt                 | ligand is in the HET record but coordinates are missing                            |
| LigandIsUnknownAtom                   | ligand is unknown (UNX, UNL)                                                       |
| LigandIsMissingInModel1               | ligand is not present in MODEL 1 that is used by default                           |
| InternalError                         | ligand could not be parsed                                                         |


## Configuration

For an overview of all settings, you can start LigandExtractor using the help command:
```
$ ./LigandExtractor --help
```


## Error Reporting

If you want to report a problem via software.zbh@uni-hamburg.de, please provide as much information
as possible. An error report should contain a short description of the problem, detailed
reproduction steps, and input data to reproduce the issue.
