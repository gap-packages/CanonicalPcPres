#
# CanonicalPcPres: Canonical Presentations of Finite Solvable Groups
# by Santiago Barrera Acevedo, Heiko Dietrich, Max Horn
#

#! @Chapter Overview

#! @Section Introduction
#!
#! The &CanonicalPcPres; package computes canonical polycyclic
#! presentations of finite solvable groups, such that isomorphic input groups yield
#! equal canonical pc presentations.
#! For a finite solvable group
#! <M>G</M>, it constructs a pc group <M>can(G)</M> and an explicit
#! isomorphism from <M>G</M> to <M>can(G)</M>.
#!
#! The underlying algorithm, its mathematical justification, and references are 
#! all described in <Cite Key="BDH"/>. The package is an implementation intended for
#! experimentation and further development:
#! it works well for many groups, but some examples require expensive orbit
#! computations. Runtime depends strongly on the module ranks and orbit sizes,
#! and not only on the order of the input group.
#!
#! As explained in <Cite Key="BDH"/>, the canonical presentation of the group <M>G</M>
#! involves certain orbit computations and defining canonical elements in these orbits.
#! These orbit computations are currently the main bottleneck, and it is to be expected that
#! future versions of this package might implement different, more efficient methods for
#! these tasks. As a consequence, it can happen that the canonical presentation of <M>G</M>
#! in a future release will be different from the canonical presentation computed by an
#! earlier version of the package. When comparing groups via their canonical presentations,
#! it is therefore important that these presentations have been computed with the same
#! release version of the package.
#! 

#! @Section Installation and dependencies
#!
#! Install &CanonicalPcPres; as a GAP package, either from a release archive
#! or by cloning
#! <URL>https://github.com/gap-packages/CanonicalPcPres</URL> into a GAP
#! package directory. Then start GAP and load the package with
#! <C>LoadPackage("CanonicalPcPres");</C>.
#!
#! The package requires GAP 4.13 or newer and the GAP packages
#! <Package>ANUPQ</Package>, <Package>orb</Package>, and
#! <Package>images</Package>.

#! @Section Quick start
#!
#! The function <Ref Func="CanonicalPcGroup"/> returns the canonical pc group.
#! The attribute <C>IsomorphismCanonicalPcGroup</C> returns an explicit
#! isomorphism from the input finite solvable group to its canonical pc group.
#!
#! @BeginExample
#! gap> LoadPackage("CanonicalPcPres");;
#! gap> G := Group([ (9,10)(11,12)(13,14), (1,2), (3,4), (8,9,11,13,14,12,10) ]);;
#! gap> iso := IsomorphismCanonicalPcGroup(G);
#! [ (9,10)(11,12)(13,14), (1,2), (3,4), (8,9,11,13,14,12,10) ] -> 
#! [ f3, f2, f1, f4^5 ]
#! gap> Source(iso) = G;
#! true
#! gap> IsBijective(iso);
#! true
#! gap> Range(iso) = CanonicalPcGroup(G);
#! true
#! @EndExample

#! @Section Comparing canonical presentations
#!
#! Two canonical pc presentations can be compared with
#! <Ref Func="CodePcGroup" BookName="ref"/>. The pair consisting of the group
#! order and this integer code uniquely determines the pc presentation. Thus,
#! when the input groups have the same order, it is enough to compare their
#! codes. We note that while the canonical presentations of two isomorphic
#! groups are the same, as GAP objects they are different. Therefore one has
#! to verify equality of the underlying presentations by other means, for example
#! by <C>CodePcGroup</C> and a size comparison.
#!
#! @BeginExample
#! gap> G := Group([ (9,10)(11,12)(13,14), (1,2), (3,4), (8,9,11,13,14,12,10) ]);;
#! gap> H := SmallGroup(56, 12);;
#! gap> CanonicalPcGroup(G) = CanonicalPcGroup(H);
#! false
#! gap> Size(H) = Size(G);
#! true
#! gap> CodePcGroup(CanonicalPcGroup(G));
#! 3649541
#! gap> CodePcGroup(CanonicalPcGroup(H));
#! 3649541
#! gap> PrintPcPresentation(CanonicalPcGroup(G),true);
#! g1^2 = id
#! g2^2 = id
#! g3^2 = id
#! g4^7 = id
#! [g4,g3] = g4^5
#! true
#! gap> PrintPcPresentation(CanonicalPcGroup(H),true);
#! g1^2 = id
#! g2^2 = id
#! g3^2 = id
#! g4^7 = id
#! [g4,g3] = g4^5
#! true
#! gap> L := List([1..10], i->PcGroupCode(RandomSpecialPcgsCoded(H),Size(H)));;
#! gap> ForAll(L, U -> CodePcGroup(CanonicalPcGroup(U)) = 3649541);
#! true
#! @EndExample

#! @Chapter Reference

#! @Section Canonical presentations

#! @Arguments G
#! @Returns the canonical pc group of <A>G</A>
#! @Description
#! Let <A>G</A> be a finite solvable group. This function returns the
#! canonical pc presentation <M>can(G)</M>. Isomorphic input groups yield
#! equal pc presentations. If <A>G</A> is not already given as a pc group,
#! the package first converts it via
#! <Ref Func="IsomorphismSpecialPcGroup" BookName="ref"/>.
#! @BeginExample
#! gap> G := SmallGroup(36, 9);;
#! gap> H := CanonicalPcGroup(G);
#! <pc group of size 36 with 4 generators>
#! gap> Size(H);
#! 36
#! @EndExample
DeclareGlobalFunction("CanonicalPcGroup");

#! @Section Isomorphism attribute

#! @Arguments G
#! @Returns an isomorphism from <A>G</A> to its canonical pc presentation
#! @Description
#! Let <A>G</A> be a finite solvable group. This attribute returns an
#! isomorphism from <A>G</A> to <Ref Func="CanonicalPcGroup"/> applied to
#! <A>G</A>. If necessary, <A>G</A> is first converted to a pc group
#! using <Ref Func="IsomorphismSpecialPcGroup" BookName="ref"/>.
#! @BeginExample
#! gap> G := SmallGroup(24, 12);;
#! gap> iso := IsomorphismCanonicalPcGroup(G);;
#! gap> [ Source(iso) = G, Size(Range(iso)) ];
#! [ true, 24 ]
#! @EndExample
DeclareAttribute("IsomorphismCanonicalPcGroup", IsGroup);

#! @Section Information classes

#! @Description
#! This info class reports the stages of the canonical-presentation
#! computation. Level 0 suppresses these messages; higher levels provide
#! progressively more detail.
DeclareInfoClass("InfoCanonicalPcPres");

#! @Description
#! This info class reports timing measurements collected during the
#! computation. Set its level to 0 for quiet, reproducible examples.
DeclareInfoClass("InfoCanonicalPcPresTimings");

#! @Section Global variables

#! There are two global variables, <C>canform_USE_NEW_AUT</C> and <C>canform_USE_NEW_CP</C>,
#! that determine whether the required automorphism group and compatible pairs constructions
#! are performed using the standard GAP functions <C>AutomorphismGroup</C> and
#! <C>CompatiblePairs</C> or the package internal variants.
#! By default, both variables are set <C>true</C>, because
#! this often seems to improve the performance. However, there are groups (even of small
#! orders) where setting one / both variables to <C>false</C> would decrease the runtime.


#! @Chapter Limitations, credits, and citation

#! @Section Limitations
#!
#! The input must be a finite solvable group. Internally, the package works
#! with pc groups. If the input is not already given as one, it is converted
#! using <Ref Func="IsomorphismSpecialPcGroup" BookName="ref"/> before the
#! canonical-form computation starts. One can also use
#! <Ref Func="IsomorphismPcGroup" BookName="ref"/>, but a special pc
#! presentation can provide faster group arithmetic and therefore reduce the
#! overall runtime of the canonical-form computation.
#!
#! The implementation is intended as a framework for the algorithm in
#! <Cite Key="BDH"/>, not as a uniformly fast isomorphism test. The principal
#! bottlenecks are orbit computations for canonical module structures,
#! compatible pairs, and cohomology classes. Some small groups can therefore
#! be substantially harder than much larger groups.

#! @Section Credits and citation
#!
#! The implementation combines tools from several parts of computational
#! group theory. In particular, it uses <Package>ANUPQ</Package> to compute
#! standard presentations of finite <M>p</M>-groups when dealing with the
#! maximal nilpotent quotient of the input group. It uses
#! <Package>orb</Package> for orbit computations and hashing, and
#! <Package>images</Package> when computing smallest orbit representatives
#! under group actions. Please credit
#! those packages when their functionality is relevant to published
#! computations.
#!
#! The bibliography of this manual is intentionally brief. The algorithm and
#! its implementation build on many important results and computational tools;
#! for the full bibliography and appropriate credits, please consult the
#! accompanying paper <Cite Key="BDH"/>.
#!
#! If you use &CanonicalPcPres; in published work, please cite the
#! accompanying paper <Cite Key="BDH"/>.
