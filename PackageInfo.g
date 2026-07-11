#
# CanonicalPcPres: Canonical Presentations of Finite Solvable Groups
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "CanonicalPcPres",
Subtitle := "Canonical Presentations of Finite Solvable Groups",
Version := "1.0",
Date := "11/07/2026", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Santiago",
    LastName := "Barrera Acevedo",
    WWWHome := "https://scholars.latrobe.edu.au/s2barreraace",
    Email := "S.BarreraAcevedo@latrobe.edu.au",
    IsAuthor := true,
    IsMaintainer := true,
    Place := "Melbourne, Australia",
    Institution := "La Trobe University",
  ),
  rec(
    FirstNames := "Heiko",
    LastName := "Dietrich",
    WWWHome := "http://users.monash.edu.au/~heikod/",
    Email := "heiko.dietrich@monash.edu",
    IsAuthor := true,
    IsMaintainer := true,
    #PostalAddress := Concatenation(
    #           "School of Mathematics\n",
    #           "Monash University\n",
    #           "VIC 3800\n",
    #           "Melbourne, Australia" ),
    Place := "Melbourne",
    Institution := "Monash University",
  ),
  rec(
    FirstNames := "Max",
    LastName := "Horn",
    WWWHome := "https://www.quendi.de/math",
    Email := "mhorn@rptu.de",
    IsAuthor := true,
    IsMaintainer := true,
    #PostalAddress := Concatenation(
    #           "Fachbereich Mathematik\n",
    #           "RPTU Kaiserslautern-Landau\n",
    #           "Gottlieb-Daimler-Straße 48\n",
    #           "67663 Kaiserslautern\n",
    #           "Germany" ),
    Place := "Kaiserslautern, Germany",
    Institution := "RPTU Kaiserslautern-Landau",
  ),
],

SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/gap-packages/CanonicalPcPres",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := "https://gap-packages.github.io/CanonicalPcPres/",
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

AbstractHTML := Concatenation(
  "The <span class=\"pkgname\">CanonicalPcPres</span> package computes ",
  "canonical polycyclic presentations of finite solvable groups, together ",
  "with explicit isomorphisms from the input groups to their canonical ",
  "presentations."
),

PackageDoc := rec(
  BookName  := "CanonicalPcPres",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Canonical Presentations of Finite Solvable Groups",
),

Dependencies := rec(
  GAP := ">= 4.13",
  NeededOtherPackages := [
     ["anupq",">= 3.3.0"],
     ["orb",">= 4.9"],
     ["images",">= 1.3"],
  ],
  SuggestedOtherPackages := [
  ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := [
  "canonical presentations",
  "solvable groups",
  "pc groups",
  "group isomorphisms",
],

));
