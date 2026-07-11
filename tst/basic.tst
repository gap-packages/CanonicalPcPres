gap> START_TEST( "basic.tst" );

# disable info messages
gap> SetInfoLevel(InfoCanonicalPcPres,0);
gap> SetInfoLevel(InfoCanonicalPcPresTimings,0);

#
gap> canform_testSingleGroup_solv(SmallGroup(200,43),20);
true

#
gap> canform_testSingleGroup_solv(SmallGroup(12100,156),4);
true

#
gap> G:=SymmetricGroup(4);
Sym( [ 1 .. 4 ] )
gap> iso := IsomorphismCanonicalPcGroup(G);;
gap> Source(iso) = G;
true
gap> IsBijective(iso);
true
gap> IsPcGroup(Range(iso));
true
gap> Range(iso) = CanonicalPcGroup(G);
true
gap> CodePcGroup(CanonicalPcGroup(G));
8281755524
gap> IdGroup(G);
[ 24, 12 ]
gap> CodePcGroup(CanonicalPcGroup(SmallGroup(24,12)));
8281755524

#
gap> canform_testSingleGroup_solv(SmallGroup(33708, 25),3);
true

#
gap> STOP_TEST( "basic.tst" );
