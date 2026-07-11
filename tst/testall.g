#
# CanonicalPcPres: Canonical Presentations of Finite Solvable Groups
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "CanonicalPcPres" );

ReadPackage( "CanonicalPcPres", "tst/utils.g" );

TestDirectory(DirectoriesPackageLibrary( "CanonicalPcPres", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
