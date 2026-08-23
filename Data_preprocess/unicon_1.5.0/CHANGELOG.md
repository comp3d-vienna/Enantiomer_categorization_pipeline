# UNICON

## Version 1.5.0

### New Features

- Reworked README
- In addition to the licence string, it is now also possible to activate the licence using a file
- Removed current user name from written SDF files
- When using the option 'split', the output file name is now set without a suffix in the middle of
  the output name, e.g., output.sdf_1.sdf is now output_1.sdf
- Fix typo in warning
- Output directories can now be provided using the '--output' option

### Bugfixes

- Added error handling when no output file suffix is given

## Version 1.4.2

### New Features

- Added flags to the -e/--stereoisomers option
- New 'All' flag allows to enumerate all stereocenters present in a molecule
- New 'Unspecified' flag allows only to enumerate unspecified stereocenters

## Version 1.4.1

### New Features

- If a molecule fails the bond length check (unreasonable long/short bonds), this molecule is not
  skipped anymore but a warning is issued and the molecule is further processed nevertheless
- The command line options for 'single' protomer was clarified