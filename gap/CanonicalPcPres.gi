#
# CanonicalPcPres: Canonical Presentations of Finite Solvable Groups
# by Santiago Barrera Acevedo, Heiko Dietrich, Max Horn
#
# This package provides the implementation of the algorithm developed in
# 
# [BADH] S. Barrera Acevedo, H. Dietrich, M. Horn.
#        Computing canonical labellings of finite solvable groups (2026)
#
#
# Some notes:
#
# 1) The namespace prefix is "canform_".
#
# 2) For efficiency, the code often uses NC variants; correctness checks are
#    concentrated in dedicated test code.
#
#
#
# Flags:
#
# The first flag is set true by default; in this case we compute Aut(H) directly using
# the information we have computed along the way. If set to false, then we use
# the GAP function AutomorphismGroup. Efficiency can vary depending on the flag.

canform_USE_NEW_AUT      := true;  # if we compute Aut(H) directly using our information
canform_USE_NEW_AUT_test := false; # tests if our Aut(H) is correct (just for debugging)

# The second flag is set true by default; in this case we compute the compatible pairs
# directly using the information we have computed along the way. If set to false, then we use
# the GAP function CompatiblePairs. Efficiency can vary depending on the flag.
canform_USE_NEW_CP      := true;  # if we compute Comp directly using our information
canform_USE_NEW_CP_test := false; # tests if our Comp is correct (just for debugging)


#############################################################################################
##
# main functions:
##
## IsomorphismCanonicalPcGroup(G)
##     input:  finite solvable group (given as pc-group, preferably with SpecialPcgs)
##     output: isomorphism G -> can(G) where can(G) is the canonical pc group isomorphic to G
##
##
## CanonicalPcGroup(G)
##     input:  finite solvable group (given as pc-group, preferably with SpecialPcgs)
##     output: the canonical pc group can(G) isomorphic to G
##
## There are some serious bottlenecks in this implementation
## that are mainly related to the computation of large orbits, see Sections 5+6 in [ABDH].




####
####
## PRELIMINARY FUNCTIONS
####
####


#############################################################################
##
##  input: dimension d, prime p
##  output: faithful permutation representation of GL(d,p)
##
canform_setup_perm_map := function(d,p)
  local g, f, o;
  g:=GL(d,p);
  f:=GF(p);
  # Acting on lines loses scalar information and changes the matrix
  # canonicalization problem. Use the nonzero vectors instead.
  o:=Filtered(AsList(f^d), v -> not IsZero(v));
  o:=Set(o, r -> ImmutableVector(f, r));
  return ActionHomomorphism(g,o,OnRight,"surjective");
end;

#############################################################################
##
##  input:  matrix m
##  output: true/false depending on m is a scalar matrix
##
canform_IsScalarMat := function(m)
  local x;
  Assert(1, NrRows(m) > 0 and NrRows(m) = NrCols(m));
  x := m[1,1];
  return ForAll([2..NrRows(m)], i -> m[i,i] = x) and IsDiagonalMat(m);
end;

#############################################################################
##
## input:  tuple T=(A_1,..,A_m) of  matrices in GL(d,p)
## output: B in GL(d,p) such that (B^-1 * A_1 * B, ..., B^-1 * A_m * B ) is
##         canonical element in GL(d,p) orbit on T (acting by conjugation)
##
canform_cantuples := function(perm_map,normalform_cache,centralizer_cache,T,p)
local lead, B, Bperm, new, C, opt, orb, orb2, bestidx, i, trace, gens, gensI, pt,
      j, active, trans, stab;

   # find the first non-central matrix in the list
   lead := PositionProperty(T, x -> not canform_IsScalarMat(x));
   if lead = fail then
      # all entries are central -> orbit has length 1 and we are done
      return rec(cantup:=T, trans:=T[1]^0);
   fi;

   # obtain the rational canonical form of T[lead]: either it is already
   # in the cache, or else we compute it and put it into the cache
   B := LookupDictionary(normalform_cache, T[lead]);
   if B = fail then
      B := RationalCanonicalFormTransform(T[lead]);
      ConvertToMatrixRep(B);
      AddDictionary(normalform_cache, T[lead], B);
   fi;

   if Length(T) = 1 then
      return rec(cantup:=[T[1]^B], trans:=B);
   fi;
   Info(InfoCanonicalPcPres, 3, "   canform_cantuples: tuple length ", Length(T));

   # map to permutation group
   Bperm := perm_map(B);

   # map to permutation group and apply "base change" B
   new   := List(T,x->perm_map(x)^Bperm);

   # compute the centralizer
   # note: element centralizers in GL(d,p) are known and could just be
   # "written down" if this becomes a bottleneck
   C := LookupDictionary(centralizer_cache, new[lead]);
   if C = fail then
      C := Centraliser(ImagesSource(perm_map), new[lead]);
      AddDictionary(centralizer_cache, new[lead], C);
   fi;

   # we want to compute the lexicographically minimal element under the
   # OnTuples action; we can do that by computing the orbit of each entry,
   # picking the minimal element in it, then moving on to the next. Since
   # we already stabilize entries 1 to lead, we start after that
   for active in [lead+1 .. Length(T)] do

      Info(InfoCanonicalPcPres, 4, "   canform_cantuples: centraliser size ", Size(C));
      if ForAll([active..Length(T)], i -> IsCentral(C, new[i])) then
         Info(InfoCanonicalPcPres, 4,
         "   canform_cantuples: orbit will be trivial, abort early");
         break;
      fi;

      if LargestMovedPoint(C) <= 256 and Size(C) > 25000 then

         # use images package to determine minimal element in orbit; but only
         # for "small" degree, and for groups that are "big enough", due to
         # some performance issues in version 1.3.3 and below of the images
         # package; see <https://github.com/gap-packages/images/issues/27>
         # for some background.
         pt := new[active];
         stab := Centralizer(C, pt);
         trans := MinimalImagePerm(C, pt, OnPoints, rec(stabilizer := stab));

         B := B * PreImagesRepresentative(perm_map, trans);
         new := OnTuples(new, trans);

         C := stab;

      else

         # brute-force orbit calculation on the truncated tuple
         opt := rec(permgens:=GeneratorsOfGroup(C));
         if active < Length(T) then
            # Unless we are in the final round, we should also compute a Schreier
            # vector so we can read off the stabilizer
            opt.schreier:=true;
         fi;

         Size(C);  # force computing size of C, helps `Orb`
         orb := Orb(C, new[active], OnPoints, opt);
         Enumerate(orb);

         Info(InfoCanonicalPcPres, 4, "   canform_cantuples:     orbit length ", Size(orb));

         # locate the minimal element in the orbit
         bestidx := PositionMinimum(orb);

         # compute a Schreier tree trace for the minimal element, and use that to
         # adjust new resp. B to map the original element to the minimal one
         trace := TraceSchreierTreeForward(orb, bestidx);
         trans := Product(trace, i -> C.(i), One(C));
         Assert(2, OnPoints(new[active], trans) = orb[bestidx]);
         Assert(2, ForAll([1..active-1], i -> new[i]^trans = new[i]));
         B := B * PreImagesRepresentative(perm_map, trans);
         new := OnTuples(new, trans);

         # compute the stabilizer, except in the last round
         if active = Length(T) then
            break;
         fi;
         C := Stabilizer(orb);
      fi;

      # abort early if C is trivial
      if IsTrivial(C) then
        break;
      fi;

      # now adjust C to be the stabilizer of the new element, by conjugating it
      C := C^trans;
      Assert(2, ForAll([1..active], i -> ForAll(GeneratorsOfGroup(C), x -> new[i]^x = new[i])));
   od;

   return rec(cantup:=OnTuples(T,B), trans:=B);

end;




##############################################################################
## input is H, Aut(H), tuple T corresponding to Pcgs(H) describing action on layer
## output is canonical tuple (describing action) and transversal element
##
canform_can_action := function(H, autH, T, p)
local cT, B, getmat, autact, orb, trans, stab, i, y, yi, yin, tt,
      perm_map, normalform_cache, centralizer_cache, posCache, Bstab,
      gensautH, gensautHexp, gpos, autact_brute, ps, m, M, CP, my, ttt;

   Info(InfoCanonicalPcPres, 2, "   canform_can_action: tuple length ", Length(T));

   # Compute a faithful permutation representation of GL(d,p) on nonzero vectors.
   perm_map := canform_setup_perm_map(NrRows(T[1]), p);

   # Caches for normal forms and centralizers.
   normalform_cache := NewDictionary(T[1], true);
   centralizer_cache := NewDictionary((1,2,3), true);

   cT := canform_cantuples(perm_map,normalform_cache,centralizer_cache,T,p);
   B  := cT.trans;
   cT := cT.cantup;

   getmat := function(T, exp)
   local new, i;
      new := T[1]^0;
      for i in [1..Size(exp)] do new := new * T[i]^exp[i]; od;
      return new;
   end;   

   ## avoid recomputing newgens for autH generators
   gensautH    := GeneratorsOfGroup(autH);
   gensautHexp := List(gensautH,aut->List(Pcgs(H),x->ExponentsOfPcElement(Pcgs(H),Image(aut^-1,x))));
   autact := function(T,gpos)
      return List(gensautHexp[gpos],exp->getmat(T,exp));
   end;   
 
   ## the version for general aut
   autact_brute := function(T,aut)
   local newgens;
      aut     := aut^-1;
      newgens := List(Pcgs(H),x->ExponentsOfPcElement(Pcgs(H),Image(aut,x)));
      return List(newgens,exp -> getmat(T,exp));
   end;

   orb  := [cT];
   trans:= [[One(autH),B]];
   stab := [];
   i      := 1;
   posCache := NewDictionary(orb[1], true);
   AddDictionary(posCache, Last(orb), Length(orb));
   while i <= Length(orb) do
      y := orb[i];
      Info(InfoCanonicalPcPres,2,"canform_can_action: iteration ", i, ", orb len = ",Size(orb));
      for gpos in [1..Size(gensautH)] do
         yi  := autact(y,gpos);

         yin := canform_cantuples(perm_map,normalform_cache,centralizer_cache,yi,p);
         B   := yin.trans;
         yin := yin.cantup;
         m := gensautH[gpos];
         ps := LookupDictionary(posCache, yin);
         if ps = fail then
            Add(orb, yin);
            Add(trans, [trans[i][1]*m,trans[i][2]*B]);
            AddDictionary(posCache, Last(orb), Length(orb));
         else
            Add(stab, [trans[i][1]*m*trans[ps][1]^-1,trans[i][2]*B*trans[ps][2]^-1]);
         fi;
      od;
      i := i+1;
   od;
   i  := PositionMinimum(orb);
   B  := trans[i];
   cT := orb[i];

   tt := List(autact_brute(T,B[1]),x-> B[2]^-1*x*B[2]);
   if not tt = cT then
      Error("canform_can_action: transported tuple does not match canonical tuple");
   fi;

   stab := Set(stab);  # Remove duplicate generators.

   ## we acted with Aut(H) on Aut(M)-canonical forms;
   ## need to act Stab_{Aut(M)}(T) gens to stab
   if canform_USE_NEW_CP then
      Bstab := GL(NrRows(T[1]),p);
      ttt := Runtime();
      for m in T do Bstab := Centraliser(Bstab,m); od;
      Append(stab, List(GeneratorsOfGroup(Bstab),x-> [One(autH),x]));
      Info(InfoCanonicalPcPresTimings,2,">>>> time comp pairs Aut(B) centraliser ",Runtime()-ttt);
   fi;

   ## in case we want to test correctness:
   if canform_USE_NEW_CP_test then
      ###### test stab original
      if not ForAll(stab, BB -> T = List(autact_brute(T,BB[1]),x-> BB[2]^-1*x*BB[2])) then
         Error("original stab elements wrong");
      fi;
   fi;
 

   ### Rewrite stabilizer generators so that they stabilize cT.
   stab := List(stab, x-> [B[1]^-1*x[1]*B[1],B[2]^-1*x[2]*B[2]]);

   ## in case we want to test correctness:
   if canform_USE_NEW_CP_test then
      ## test new stab gens
      if not ForAll(stab, BB -> cT = List(autact_brute(cT,BB[1]),x-> BB[2]^-1*x*BB[2])) then
         Error("new stab elements wrong");
      fi;

      ## test whether we get comp pairs
      M     := GModuleByMats(cT,GF(p));
      Info(InfoCanonicalPcPres, 3, "compute CP using GAPs function");
      CP    := Group(GeneratorsOfGroup(CompatiblePairs(H,M))); ##redefine CP, otherwise there was a bug!!!
      Info(InfoCanonicalPcPres, 3, "test size");
      my    := Group(Set(stab, DirectProductElement));
      if not ForAll(GeneratorsOfGroup(my), x-> x in CP) then
         Error("canform_can_action: computed stabilizer generators are not in CompatiblePairs(H, M)");
      fi;
      if not Size(CP) = Size(my) then
         Error("canform_can_action: compatible-pairs size check failed");
      fi;
      Info(InfoCanonicalPcPres, 3, "tested CP and all correct");
   fi;


   return rec(canact:=cT, stab:=stab, trans:=B);
end;


#############################################################################
##
## input:  m x d matrix T over GF(p)
## output: B in GL(d,p) such that T*B is a "canonical" elt in orbit T*GL(d,p)
##
canform_canmat := function(T,p)
local lrows, linds, d, mb, i, v, B;

## initialise linearly independent rows and rank
   lrows  := [];
   linds :=  [];
   d     := Length(T[1]);

## find max number of lin independent rows
   mb := MutableBasis( GF(p), [], Zero( T[1] ) );
   for i in [1..Length(T)] do
      if CloseMutableBasis( mb, T[i] ) then
         Add(linds,i);
         Add(lrows,T[i]);
      fi;
      if Length(lrows)=d then break; fi;
   od;

## complement B to basis: add std basis vectors until rank d
   if Length(lrows) < d then
      # TODO: we could do better, as we know the pivots in mb,
      # so know which rows of k to add
      for v in IdentityMat(d, GF(p)) do
         if CloseMutableBasis( mb, v ) then
            Add(lrows,v);
         fi;
         if Length(lrows)=d then break; fi;
      od;
   fi;

   B := lrows^-1;
   return rec(canmat:=T*B, trans:=B, inds:=linds);
end;





###################################################################################
######
###### this is essentially GAP's TwoCohomology (copied from GAP library), to see
###### what's going on. It's all deterministic, and since we assume
###### that H=can(G/G_i) and M is an H-module
###### this computation and choice of bases is canonical.
###### Also, the computation of preimages via cohom is deterministic
######
###### input:  pc-group H with corresponding GModuleByMats M
###### output: same as TwoCohomology, with added inverse map
######
canform_my_getH2 := function(H,M)
local C, d, z, Z2, co, B2, cb,  V, W, F, comp, gen, B, imb, canbas, zero,
      Bimgs, nathom, H2, mypreim, Wvectors, mb, pr, img, compl, cohominv;

   ## set up data structures; in particular, presentation of H etc
   C := CollectorSQ( H, M, false );;
   d := Length( C.orders );;
   d := d * (d + 1) / 2;;
   z := Flat( List( [ 1 .. d ], function ( x ) return C.mzero[1];  end ) );;

   ## now compute 2-coc and 2-cob tails
   Z2 := TwoCocyclesSQ( C, H, M );;
   co := VectorSpace( M.field, Z2, z );;
   B2 := TwoCoboundariesSQ( C, H, M );;
   cb := SubspaceNC( co, B2 );;
   pr := FpGroupPcGroupSQ( H );;

   ## now get map onto Z2/B2
   ## this is essentially NaturalHomomorphismBySubspaceOntoFullRowSpace
   V  := co;
   W  := cb;
   F  := M.field;
   if Dimension( V ) = Dimension( W ) then
        nathom   := ZeroMapping( V, FullRowModule( F, 0 ) );
        cohominv := ZeroMapping( W, V);
   else  
      Wvectors := BasisVectors( Basis( W ) );
      mb := MutableBasis( F, Wvectors, Zero( W ) );
      compl := [  ];
      for gen in BasisVectors( Basis( V ) ) do
         if CloseMutableBasis( mb, gen ) then
            Add( compl, gen );
         fi;
      od;
      B      := BasisNC( V, Concatenation( Wvectors, compl ) );
      img    := FullRowModule( F, Length( compl ) );
      canbas := CanonicalBasis( img );
      zero   := Zero( img );
      Bimgs  := Concatenation( List( Wvectors, function ( v ) return zero; end ), BasisVectors( canbas ) );
      nathom := LeftModuleHomomorphismByMatrix( B, Bimgs, canbas );
      SetIsSurjective( nathom, true );
      nathom!.basisimage := canbas;
      nathom!.preimagesbasisimage := Immutable( compl );
      SetKernelOfAdditiveGeneralMapping( nathom, W );
      UseFactorRelation( V, W, img );

      ## define what the inverse of nathom: Z2 --> Z2/B2 is doing
      ## this is essentially PreImagesRepresentatives
      cohominv := function(x)
      local res;
         res := Coefficients( nathom!.basisimage, x );
         return LinearCombination( nathom!.preimagesbasisimage, res );
      end;   

   fi;
   ##### return H2 data rec
   H2 := rec(
             group          := H,
             module         := M,
             collector      := C,
             isPcCohomology := true,
             cohom          := nathom,
             cohominv       := cohominv,
             presentation   := pr );;
   #####
return H2;
end;


#############################################################################
## this is mainly copied from GAP's MatrixOperationOfCPGroup to control what's
## happening
##
## input:  coh record cc, compatible pair g, 2-coc tailvector coc
## output: image tailvector coc^g
##
## NOTE: this is on tail vectors, not on classes in H^2
##       this function is not used for the orbit calculation
##       but for getting the 2-cob via coc^g - coc, which is 
##       necessary for writing down the isomorphism
##
canform_actWithCPonTail := function ( cc, g ,coc )
local pcgs, ords, imgs, n, d, fpgens, fprels, H,
      pcgsH, imgl, k, i, j, rel, tail, m, tails, field;

    pcgs := Pcgs( cc.group );
    ords := RelativeOrders( pcgs );
    imgs := List( [g], function ( x )
            return List( pcgs, function ( y )
                    return y ^ Inverse( x[1] );
                end );
        end );
    n      := Length( pcgs );
    d      := cc.module.dimension;
    field  := cc.module.field;
    fpgens := GeneratorsOfGroup( cc.presentation.group );
    fprels := cc.presentation.relators;

    H     := ExtensionSQ( cc.collector, cc.group, cc.module, coc );
    pcgsH := Pcgs( H );

    imgl := List( imgs[1], function ( x ) return MappedPcElement( x, pcgs, pcgsH ); end );
    if imgl <> pcgs then
       k := 0;
       tails := [  ];
       for i in [ 1 .. Length( pcgs ) ] do
          for j in [ 1 .. i ] do
             k := k + 1;
             rel := fprels[k];
             tail := MappedWord( rel, fpgens, imgl );
             if not IsBound( cc.module.isCentral ) or not cc.module.isCentral then
                if i = j then
                   m := imgl[i] ^ ords[i];
                else
                   m := imgl[i] ^ imgl[j];
                fi;
                tail := tail ^ m;
             fi;
             tail := ExponentsOfPcElement( pcgsH, tail, [ n + 1 .. n + d ] );
             tail := tail * g[2];
             ConvertToVectorRepNC( tail, field );
             if Length( tails ) = 0 then
                tails := tail;
             else
                Append( tails, tail );
             fi;
          od;
      od;
   else
      Error("canform_actWithCPonTail: image pc generators unexpectedly equal original pc generators");
   fi;
   return tails;
end;


######################################################################################
## input:  cohomology record cc, G, M, and and a 2-cob cob
## output: list v of exponents (wrt pcgs(M)), one for each Pcgs(G)=[g1,..,gn] element
## such that [g1 v[1], ..., gn v[n]] generates complement to M in Extension(G,M,cob)
## (here v[i] has to be replaced by Product([1..Size(Pcgs(M))],j-> Pcgs(M)[j]^v[i][j])
##
## cob can also be a list of cobs
##
## this is a modification of GAP's TwoCoboundariesSQ; this is used to write down the
## isomorphism: for a coboundary def by a map f we need the images f(g_i)
##
canform_getCOBfunctionImages := function( C, G, M, cob )
    local   n,  R,  MI,  j,  i,  x,  m,  e,  k,  r,  d, v, c, res;

    # start with zero matrix
    n := Length(Pcgs( G ));
    R := [];
    r := n*(n+1)/2;
    for i  in [ 1 .. n ]  do
        R[i] := [];
        for j in [ 1 .. r ] do
            R[i][j] := C.mzero;
        od;
    od;

    # compute inverse generators
    M  := M.generators;
    MI := List( M, x -> x^-1 );
    d  := Length(M[1]);

    # loop over all relators
    for j  in [ 1 .. n ]  do
        for i in  [ j .. n ]  do
            x := (i^2-i)/2 + j;

            # power relator
            if i = j  then
                m := C.mone;
                for e  in [ 1 .. C.orders[i] ]  do
                    R[i][x] := R[i][x] - m;  m := M[i] * m;
                od;

            # conjugate
            else
                R[i][x] := R[i][x] - M[j];
                R[j][x] := R[j][x] + MI[j]*M[i]*M[j] - C.mone;
            fi;

            # compute fox derivatives
            m := C.mone;
            r := C.relators[i][j];
            if r <> 0  then
                for k  in [ Length(r)-1, Length(r)-3 .. 1 ]  do
                    for e  in [ 1 .. r[k+1] ]  do
                        R[r[k]][x] := R[r[k]][x] + m;
                        m := M[r[k]] * m;
                    od;
                od;
            fi;
        od;
    od;

    # make one list
    m := [];
    r := n*(n+1)/2;
    for i  in [ 1 .. n ]  do
        for k  in [ 1 .. d ]  do
            e := [];
            for j  in [ 1 .. r ]  do
                Append( e, R[i][j][k] );
            od;
            Add( m, e );
        od;
    od;
   ## columns of m correspond to gens of G
   ## rows correspond to relations
   ## image of m is B^2 tail vectors
   ## so we want preiamge of cob \in image(m)
   if IsList(cob[1]) then
      res := [];
      for c in cob do
         v   := SolutionMat(m,c);
         if not Size(v) mod d = 0 then
            Error("canform_getCOBfunctionImages: solution vector length is not divisible by module dimension");
         fi;
         v   := List([1..Size(v)/d],ii-> List(v{[(ii-1)*d+1.. ii*d]},Int));
         Add(res,v);
      od;
      return res;
   else
      c := cob;
      v := SolutionMat(m,c);
      if not Size(v) mod d = 0 then
         Error("canform_getCOBfunctionImages: solution vector length is not divisible by module dimension");
      fi;
      v   := List([1..Size(v)/d],ii-> List(v{[(ii-1)*d+1.. ii*d]},Int));
      return v;
   fi;
end;





##############################################################################
## naive orbit-stab alg that keeps track of comp pairs
## input:  list of mats corresponding to cp action of comp pairs cps, and 
##         coc tailvector coc
## output: enumerated orbit orb, with transversal and stabiliser gens
##
canform_myorbit := function(mats,cps, coc)
local orb, i, trans, stab,  y, j, m, c, yi, p;

   orb   := [coc];
   trans := [[mats[1]^0,cps[1]^0]];
   stab  := [];
   i     := 1;
   while i<= Size(orb) do
      y := orb[i];
      for j in [1..Size(mats)] do
         m  := mats[j];
         c  := cps[j];
         yi := y*m;
         p  := Position(orb,yi);  # TODO: hashtable?
         if p = fail then
            Add(orb,yi);
            Add(trans,[trans[i][1]*m,trans[i][2]*c]);
         else
            Add(stab, [trans[i][1]*m*trans[p][1]^-1,trans[i][2]*c*trans[p][2]^-1]);
         fi;
      od;
      i:=i+1;
   od;

   stab := Set(stab);  # Remove duplicate generators.

   return rec(orb:=orb, trans:=trans, stab:=stab);
end;


##############################################################################
## here we assume that M is trivial H module, so Comp(H,M)=Aut(H) x Aut(M)
## naive orbit-stab alg that splits the action of Comp(H,M) = Aut(H) x Aut(M)
## such that Aut(H) acts on Aut(M)-canonical forms
## input:  Aut(H), cohom record H2, and tail vector
## output: Comp-orbit of corresponding cohom class (with transversal etc)
##
canform_myorbit_reducedAutM := function(autH, H2, tails)
local toH2, toTail, p, d, canM, mats, t, orb, trans, stab, i, j,
      y, yi, m, B, ps, matsCP, posCache;

  Info(InfoCanonicalPcPres, 1, "           canform_myorbit_reducedAutM: enter");

   toH2   := H2.cohom;
   toTail := H2.cohominv;
   p      := Size(H2.module.field);
   d      := Length(H2.collector.module[1][1]);

   # HACK: at this point, toH2 may be something like
   # <linear mapping by matrix, <vector space of dimension 20 over GF(2)> -> ( GF(2)^20 )>
   # and using that is expensive for no good reason. So replace it
   # by an identity map then
   if IsOne(toH2) then
      canM   := function(tt)
         tt := List([1..Length(tt)/d],i-> tt{[(i-1)*d+1 .. i*d]});
         tt := canform_canmat(tt,p);
         return rec(tail := Flat(tt.canmat), trans := tt.trans);
      end;

   else
      ## takes tail vector t in H2, preim in Z2, computed can_M, returns image in H2
      canM   := function(t)
         local tt;
         tt := toTail(t);
         tt := List([1..Length(tt)/d],i-> tt{[(i-1)*d+1 .. i*d]});
         tt := canform_canmat(tt,p);
         return rec(tail := toH2(Flat(tt.canmat)), trans := tt.trans);
      end;
   fi;

   Info(InfoCanonicalPcPres,1,"             ... get CP on H2 mats");

   mats     := GeneratorsOfGroup(autH);
   matsCP   := MatrixOperationOfCPGroup( H2, List(mats,x->[x,One(GL(d,p))]));;
   Info(InfoCanonicalPcPres,1,"             done; now orbit");

   t     := canM(toH2(tails));
   ConvertToVectorRep(t.tail, p);
   orb   := [t.tail];
   trans := [[mats[1]^0,t.trans]];
   stab  := [];
   posCache := NewDictionary(orb[1], true);
   AddDictionary(posCache, Last(orb), Length(orb));
   i     := 1;
   while i<= Length(orb) do
      Info(InfoCanonicalPcPres, 4, "      canform_myorbit_reducedAutM: orbit position ", i, " out of ", Length(orb));
      y := orb[i];
      for j in [1..Size(mats)] do
         m  := mats[j];
         yi := y*matsCP[j];
         yi := canM(yi);  # FIXME: bottleneck!
         B  := yi.trans; 
         yi := yi.tail;
         ps := LookupDictionary(posCache, yi);
         if ps = fail then
            Add(orb, yi);
            Add(trans, [trans[i][1]*m,trans[i][2]*B]);
            AddDictionary(posCache, Last(orb), Length(orb));
         else
            Add(stab, [trans[i][1]*m*trans[ps][1]^-1,trans[i][2]*B*trans[ps][2]^-1]);
         fi;
      od;
      i:=i+1;
   od;

   stab := Set(stab);  # Remove duplicate generators.

  # TODO: I don't think we need the full transversal. A Schreier vector  
  # is fine; in the end, the code locates the "canonical" element
  # in the orbit and uses its transversal representative, that's all
  
   return rec(orb:=orb, stab:= stab, trans:=List(trans,x->[1,x]));
end;


##############################################################################
##
## lift aut group from stab_comp(tails+B^2) and from Z1(H,M)
##
canform_getautgroup := function(H2, cps, Etails, tails, H, M)
local pct, n, cp, alphai, betai, d, hprei, mprei,  new, ii, el, j, iso, auts, cobs,
      vs, v, Z1;

   pct   := Pcgs(Etails);
   n     := Size(Pcgs(H));
   auts  := [];
   d     := M.dimension;

   ## first deal with auts induced by compatible pairs
   ##
   ## note: since we have several orbit/stab functions, there is a bit of
   ## inconsistency how we store the stabiliser in the comp pairs:
   ## sometimes these is just comp pairs (a,b), sometimes it is [m, (a,b)]
   ## where the matrix m describes (a,b) acting as a matrix. here we only
   ## need (a,b), so we need to filter this, if required
   if not cps in [ [], [[]] ] then
      if not IsGroupHomomorphism(cps[1][1]) then
         cps   := List(cps,x->x[2]);
      fi;

      cobs  := List(cps,c-> canform_actWithCPonTail(H2, c, tails)-tails);
      vs    := canform_getCOBfunctionImages(H2.collector, H, M, -cobs);

      for ii in [1..Size(cps)] do
         cp       := cps[ii];
         v        := vs[ii];
         alphai   := cp![1]^-1;
         betai    := List(cp![2]^-1,x->Concatenation(ListWithIdenticalEntries(n,0),List(x,Int)));
         hprei    := List( List(Pcgs(H),x-> Concatenation(ExponentsOfPcElement(Pcgs(H),Image(alphai,x)),
                                                     ListWithIdenticalEntries(M.dimension,0))),
                        ex -> PcElementByExponents(pct,ex));
         mprei    := List(betai,ex->PcElementByExponents(pct,ex));
         hprei    := Concatenation(hprei,mprei);

         pct  := Pcgs(Etails);
         new  := [];
         for ii in [1..Size(v)] do
            el := pct[ii];
            for j in [1..Size(v[ii])] do
               if v[ii][j]>0 then el := el * pct[n+j]^v[ii][j]; fi;
            od;
            Add(new,el);
         od;
         for ii in [n+1..Size(pct)] do Add(new,pct[ii]); od;
         iso := GroupHomomorphismByImagesNC(Etails,Etails,new,hprei);
         if iso = fail then
            Error("canform_getautgroup: failed to construct automorphism induced by a compatible pair");
         fi;
         Add(auts,iso);
      od;
   fi;

   ## now deal with Z1 automorphisms
   Z1 := BasisVectors(Basis(OneCocycles(Etails,Subgroup(Etails,pct{[n+1..n+d]})).oneCocycles));
   Z1 := List(Z1,v->List([1..Size(v)/d],ii-> List(v{[(ii-1)*d+1.. ii*d]},Int)));
   for v in Z1 do
       pct  := Pcgs(Etails);
       new  := [];
       for ii in [1..Size(v)] do
          el := pct[ii];
          for j in [1..Size(v[ii])] do
             if v[ii][j]>0 then el := el * pct[n+j]^v[ii][j]; fi;
          od;
          Add(new,el);
       od;
       for ii in [n+1..Size(pct)] do Add(new,pct[ii]); od;
       iso := GroupHomomorphismByImagesNC(Etails,Etails,new,pct);
       if iso = fail then
          Error("canform_getautgroup: failed to construct automorphism induced by a 1-cocycle");
       fi;
       Add(auts,iso);
   od;

   return auts;
end;


####################################################################################
## input: solvable group G; note that nilpotent groups are dealt with directly!
## output: LG series
##
## Can optionally be passed as second argument to canform_compute_canonical
##
canform_LGseries := function(G)
local s;
   s := SpecialPcgs(G);
   return List(LGFirst(s), x -> Subgroup(G, s{[x..Length(s)]}));
end;   


####################################################################################
## input:  pc group G; note that nilpotent groups are dealt with directly!
## output: char subgroup series (as discussed in paper)
##
## Can optionally be passed as second argument to canform_compute_canonical
canform_char_series := function(G)
local ser, i, j, primes, hom, char_series_ab, abs,new;



   if IsAbelian(G) then return [G,TrivialSubgroup(G)]; fi;

   char_series_ab := function(H)
   local ser, new, i, hom, abs, j,ss;
      ser  := List(SylowSystem(H),GeneratorsOfGroup);
      ser  := List([1..Size(ser)],i-> Subgroup(H,Concatenation(ser{[i..Size(ser)]})));
      Add(ser,Subgroup(H,[One(H)]));
      new  := [];
      for i in [1..Length(ser)-1] do
         hom := NaturalHomomorphismByNormalSubgroup(ser[i],ser[i+1]);
         abs := PCentralSeries(Image(hom));
         if i<Length(ser)-1 then abs := abs{[1..Size(abs)-1]}; fi;
         for j in abs do Add(new,PreImage(hom,j)); od;
      od;
      return new;
   end;

   ser := DerivedSeries(G);
   new := [ser[1]];
   for i in [2..Length(ser)-1] do
      hom := NaturalHomomorphismByNormalSubgroup(ser[i],ser[i+1]);
      abs := char_series_ab(Image(hom));
      if i<Length(ser)-1 then abs := abs{[1..Size(abs)-1]}; fi;
      for j in abs do Add(new,PreImage(hom,j)); od;
   od;

   new[Length(new)] := TrivialSubgroup(G);

   return new;
end;



####################################################################################
## input: solvable group G; note that nilpotent groups are dealt with directly!
## output: char subgroup series (as discussed in paper), first section nilpotent
##
## Can optionally be passed as second argument to canform_compute_canonical
## (and is currently the default if no second argument is passed)
canform_char_series_nilp := function(G)
local ser, i, j, primes, hom, char_series_ab, abs,new;



   if IsAbelian(G) then return [G,TrivialSubgroup(G)]; fi;

   char_series_ab := function(H)
   local ser, new, i, hom, abs, j,ss;
      ser  := List(SylowSystem(H),GeneratorsOfGroup);
      ser  := List([1..Size(ser)],i-> Subgroup(H,Concatenation(ser{[i..Size(ser)]})));
      Add(ser,Subgroup(H,[One(H)]));
      new  := [];
      for i in [1..Length(ser)-1] do
         hom := NaturalHomomorphismByNormalSubgroup(ser[i],ser[i+1]);
         abs := PCentralSeries(Image(hom));
         if i<Length(ser)-1 then abs := abs{[1..Size(abs)-1]}; fi;
         for j in abs do Add(new,PreImage(hom,j)); od;
      od;
      return new;
   end;

   ser := DerivedSeries(G);
   new := [ser[1]];
   for i in [2..Length(ser)-1] do
      hom := NaturalHomomorphismByNormalSubgroup(ser[i],ser[i+1]);
      abs := char_series_ab(Image(hom));
      if i<Length(ser)-1 then abs := abs{[1..Size(abs)-1]}; fi;
      for j in abs do Add(new,PreImage(hom,j)); od;
   od;

   new[Length(new)] := TrivialSubgroup(G);

   i := Last([1..Size(ser)],x->IsNilpotentGroup(G/ser[x]));
   ser := ser{Concatenation([1],[i..Size(ser)])};

   return new;
end;


####################################################################################
## input: finite nilpotent pc group G
## output: isom G -> can(G) by using stand pres for p-groups
##
canform_nilpotent := function(G)
local syls, cans, DP, gens, ims, iso, pcv;
   Assert(0, IsNilpotentGroup(G));
   syls := SylowSystem(G);
   cans := List(syls,x->EpimorphismStandardPresentation(x:Prime:=PrimePGroup(x)));
   pcv  := List(cans,x->PcGroupFpGroup(Range(x)));
   cans := List([1..Size(cans)],x->
           cans[x]*GroupHomomorphismByImagesNC(Range(cans[x]),pcv[x],GeneratorsOfGroup(Range(cans[x])),GeneratorsOfGroup(pcv[x])));
   DP   := DirectProduct(pcv);
   cans := List([1..Size(syls)],i-> cans[i]*Embedding(DP,i));
   gens := List(syls,x->GeneratorsOfGroup(x));
   ims  := Concatenation(List([1..Size(syls)],i->List(gens[i],x->Image(cans[i],x)  )));
   gens := Concatenation(gens);
   iso  := GroupHomomorphismByImagesNC(G,DP,gens,ims);

   return iso;
end;



#####################################################################################
# compatible pairs
#
# For now this is just a stripped down copy of GAP's CompatiblePairs,
# with lots of stuff we don't need removed, and an important performance
# bottleneck eliminated.
#
# This function is not used by default; we only use it if `canform_USE_NEW_CP`
# is set to `false`.
#
canform_compatible_pairs := function( G, M )
local Mgrp, oper, A, B, D, translate, gens, genimgs, triso, K, K1,
  K2, f, tmp, Ggens, pcgs, l, idx, u, tup,Dos,preimlist,pows,
  baspt,newimgs,i,j,basicact,neu,K1nontriv,epi,hf,pool,modulehom,test;

    Info( InfoCompPairs, 1, "    CompP: |G| = ", Size(G), ", dim(M) = ", M.dimension);
    Assert(0, IsPcGroup(G));
    A:=fail;
    Mgrp := GroupByGenerators( M.generators );
    Ggens:=Pcgs(G);
    oper:=fail;
    oper := GroupHomomorphismByImagesNC( G, Mgrp, Ggens, M.generators );

    # automorphism groups of G and M
    Info( InfoCompPairs, 1, "    CompP: compute aut group");
    A:=AutomorphismGroup(G);
    B := GL( M.dimension, Characteristic( M.field ) );
    D := DirectProduct( A, B );

    # the trivial case
    if IsBound( M.isCentral ) and M.isCentral then
        return D;
    fi;

    # do we translate D in a permutation group?
    translate:=EXPermutationActionPairs(D);
    if translate<>false then

      D:=translate.permgroup;
      gens:=translate.permgens;
      genimgs:=translate.pairgens;
      triso:=translate.isomorphism;
      translate:=true;
    else
      gens:=GeneratorsOfGroup(D);
      genimgs:=gens;
    fi;

    Dos:=Size(D);

    # get kernel of oper
    K := KernelOfMultiplicativeGeneralMapping( oper );
    Assert(0, IsPcGroup(K));

    # compute stabilizer of K in A
    if Size(K)>1 then

      # get its stabilizer
      K1:=CanonicalPcgsWrtFamilyPcgs(Centre(K));
      K1nontriv:=Length(K1)>0;
      K2:=CanonicalPcgsWrtFamilyPcgs(K);
      f := function( pt, a )
             return CanonicalPcgsWrtFamilyPcgs(Group(List(pt,i->Image( a[1], i ))));
           end;

      if K1nontriv and K1<>K2 then
        tmp := Stabilizer( D, K1,gens,genimgs, f );

        if Size(tmp)<Size(D) then
          Info( InfoMatOrb, 1, "    CompP: found orbit of centre of length ",
                Size(D)/Size( tmp ));
          D := tmp;
          if translate<>false then
            if HasIsSolvableGroup(D) and IsSolvableGroup(D) then
              gens:=Pcgs(D);
            else
              gens:=GeneratorsOfGroup(D);
            fi;
            genimgs:=List(gens,i->ImageElm(triso,i));
            translate:=rec(pairgens:=genimgs,
                           permgens:=gens,
                           isomorphism:=triso,
                           permgroup:=D);
            EXReducePermutationActionPairs(translate);
            gens:=translate.permgens;
            genimgs:=translate.pairgens;
            triso:=translate.isomorphism;
            D:=translate.permgroup;
          else
            gens:=GeneratorsOfGroup(D);
            genimgs:=gens;
          fi;
        fi;
        tmp:=false; # clear memory

      fi;

      tmp := Stabilizer( D, K2,gens,genimgs, f );

      if Size(tmp)<Size(D) then
        Info( InfoMatOrb, 1, "    CompP: found orbit of length ",
              Size(D)/Size(tmp));
        D := tmp;
        if translate<>false then
          if HasIsSolvableGroup(D) and IsSolvableGroup(D) then
            gens:=Pcgs(D);
          else
            gens:=GeneratorsOfGroup(D);
          fi;
          genimgs:=List(gens,i->ImageElm(triso,i));
          translate:=rec(pairgens:=genimgs,
                          permgens:=gens,
                          isomorphism:=triso,
                          permgroup:=D);
          EXReducePermutationActionPairs(translate);
          gens:=translate.permgens;
          genimgs:=translate.pairgens;
          triso:=translate.isomorphism;
          D:=translate.permgroup;
        else
          gens:=GeneratorsOfGroup(D);
          genimgs:=gens;
        fi;
      fi;
      tmp:=false; # clear memory

    fi;

    # compute stabilizer of M.generators in D

    basicact:=function( tup, elm )
    local gens;
      #gens := List( tup[1], x -> PreImagesRepresentative( elm[1], x ) );
      #gens := List( gens, x -> MappedPcElement( x, tup[1], tup[2] ) );
      gens := List( Ggens, x -> PreImagesRepresentative( elm[1], x ) );
      gens := List( gens, x -> MappedPcElement( x, Ggens, tup ) );
      gens := List( gens, x -> x ^ elm[2] );
      return gens;
      #return DirectProductElement( [tup[1], gens] );
    end;

    # build tails of the pcgs that are closed under automorphisms
    pcgs:=Pcgs(G);
    l:=Length(Pcgs(G))+1;
    repeat
      Unbind(tmp);
      repeat
        l:=l-1;
        idx:=[l..Length(pcgs)];
        u:=SubgroupNC(G,pcgs{idx});
      until ForAll(GeneratorsOfGroup(u),
        i->ForAll(GeneratorsOfGroup(A),j->Image(j,i) in u));

      Ggens:=InducedPcgsByPcSequence(pcgs,pcgs{idx});
      tup:=M.generators{idx};

      tmp := Stabilizer( D, tup,gens,genimgs, basicact );

      Info( InfoMatOrb, 1, "    CompP: ",l,"-tail found orbit of length ",
            Size(D)/Size(tmp));
      if Size(tmp)<Size(D) then
        D:=tmp;
        if IsPcgs(gens) then
          gens:=InducedPcgs(gens,tmp);
        else
          gens:=SmallGeneratingSet(tmp);
        fi;
        genimgs:=List(gens,i->ImageElm(triso,i));
        if translate<>false then
          translate:=rec(pairgens:=genimgs,
                          permgens:=gens,
                          isomorphism:=triso,
                          permgroup:=D);
          EXReducePermutationActionPairs(translate);
          gens:=translate.permgens;
          genimgs:=translate.pairgens;
          triso:=translate.isomorphism;
          D:=translate.permgroup;
        fi;
      fi;
    until l=1;

    if translate<>false then
      l:=Size(D);
      if Length(gens)>3 then
        # reduce generator number

        u:=SmallGeneratingSet(D);
        if IsSubset(gens,u) then
          Info( InfoMatOrb, 3, "Reduce generators subset");
          idx:=List(u,x->Position(gens,x));
          gens:=gens{idx};
          genimgs:=genimgs{idx};
        else
          Info( InfoMatOrb, 3, "Reduce generators new words");
          gens:=u;
          genimgs:=List(gens,i->ImageElm(triso,i));
        fi;
      fi;
      tmp:=SubgroupNC(Range(triso),genimgs);
      #SetIsGroupOfAutomorphismsFiniteGroup(tmp,true);
      SetSize(tmp,l);
      # cache the faithful permutation representation in case we need it
      # later
      tmp!.permrep:=rec(pairgens:=genimgs,
                      permgens:=gens,
                      permgroup:=D);
      D:=tmp;
    fi;
    Info( InfoMatOrb, 1, "Total index: ",Dos/Size(D));
    return D;
end;


####################################################################################
##
## compute a canonical pc presentation for a given pc group
##
## input: a pc group GG
##
## optional second argument is type of series, e.g. canform_char_series_nilp (default)
## or canform_char_series or canform_LGseries
## 
## output is an isomorphism GG -> can(GG)
##
BindGlobal("canform_compute_canonical", function(GG, useseries...)
local lcp, homs, G, P, hom, Q, H, N, isoHQ, isoPN, isoNP, M, H2, i, F, rel,
      pre, tails, C, mats, Mgrp, orb, can, p, pos, cob, v, ii, oldG, T,
      Etails, pct, n, cp, alphai, betai, d, hprei, mprei, new, el, j, isoCanToTail,
      isoHtoG, isoOldGtoQ, isoTailtoG, autgens, gap, q,
      ords, imgl, k, tail, jj, m, autH, isoQH, ttt, oldH;

   if Length(useseries)=0 then
      useseries := canform_char_series_nilp;
     #useseries :=  canform_char_series;              ## our series with first quotient G/G'
   else
      useseries := useseries[1]; ## e.g. canform_LGseries
   fi;   

   ## cover the trivial cases
   if not IsPcGroup(GG) then Error("input needs to be pc group"); fi;
   if Size(GG) = 1 then
      H := AbelianGroup([1]);
      return GroupHomomorphismByImagesNC(GG,H,[One(GG)],[One(H)]);
   fi;

 ## Abelian case: use GAP's generic isomorphism routine.
   if IsAbelian(GG) then
      n := AbelianInvariants(GG);
      H := AbelianGroup(n);
      return IsomorphismGroups(GG,H);
   fi;

   if IsNilpotent(GG) then
     return canform_nilpotent(GG);
   fi;


   ## get characteristic series for our iteration

   lcp := useseries(GG); 

   new := List([1..Length(lcp)-1],x->Collected(FactorsInt(Size(lcp[x])/Size(lcp[x+1])))[1]);
   Info(InfoCanonicalPcPres,1,"module sizes are:");
   for i in new do
      Info(InfoCanonicalPcPres,1,i[1],"^",i[2]);
   od;


   ## if group is not nilpotent then we need to iterate
   ## the first canonical factor is just GG/lcp[2] = C_p^r
   ##
   ## so G = GG/lcp[3] has Q=GG/lcp[2] as a quotient with kernel P=lcp[2]/lcp[3]
   ## we set up H = can(Q) and N\cong P \cong M and then consider 
   ## the extension of H by M that is isomorphic to G
   ##


   Info(InfoCanonicalPcPres,1,"get first round's quotients and isoms");
   ttt   := Runtime();


   ## depending on our series, first quotient is abelian or nilpotent
   if IsAbelian(GG/lcp[2]) then

      homs  := NaturalHomomorphismByNormalSubgroup(GG,lcp[3]);
      G     := Image(homs);;
      P     := Image(homs,lcp[2]);
      hom   := NaturalHomomorphismByNormalSubgroup(G,P);
      Q     := Image(hom);
      q     := AbelianInvariants(Q);
      H     := AbelianGroup(q);

      p     := PrimePGroup(P);
      N     := AbelianGroup(List([1..RankPGroup(P)],x->p));

      ## name indicates source and image
      isoHQ := IsomorphismGroups(H,Q);        ##brute force for abelian part
      isoQH := InverseGeneralMapping(isoHQ);

   elif IsNilpotentGroup(GG/lcp[2]) then

      homs  := NaturalHomomorphismByNormalSubgroup(GG,lcp[3]);
      G     := Image(homs);;
      P     := Image(homs,lcp[2]);
      hom   := NaturalHomomorphismByNormalSubgroup(G,P);
      Q     := Image(hom);

      isoQH := canform_nilpotent(Q);
      isoHQ := InverseGeneralMapping(isoQH);
      H     := Range(isoQH);

      p     := PrimePGroup(P);
      N     := AbelianGroup(List([1..RankPGroup(P)],x->p));
  
   fi;



   isoPN := GroupHomomorphismByImagesNC(P,N,Pcgs(P),Pcgs(N));
   isoNP := GroupHomomorphismByImagesNC(N,P,Pcgs(N),Pcgs(P));

   T     := List(Pcgs(H), x ->  PreImagesRepresentative(hom,Image(isoHQ,x)));
   T     := List(T, x-> List(Pcgs(P),a-> ExponentsOfPcElement(Pcgs(P),a^x))*One(GF(p)));
   List(T, ConvertToMatrixRep);
   autH  := AutomorphismGroup(H);

   Info(InfoCanonicalPcPresTimings,2,">>>> time for first quotient, isoms, autgrps: ",Runtime()-ttt);


   ####### get canonical module structure for the first iteration
   Info(InfoCanonicalPcPres,1,"get canonical module: group size ",Size(H), " autH size ",Size(autH)," module size ",Size(P));
 
   ttt   := Runtime();
   T     := canform_can_action(H, autH, T, p);

   Info(InfoCanonicalPcPresTimings,2,">>>> time canform_can_action: ",Runtime()-ttt);

   M     := GModuleByMats(T.canact,GF(p));
   ## recall that T.trans = [a,b] such that (a,b) * oldT = T.canact
   isoQH := isoQH*T.trans[1];
   isoHQ := InverseGeneralMapping(isoQH);
   isoPN := isoPN*GroupHomomorphismByImagesNC(N,N,Pcgs(N),List(T.trans[2],x->PcElementByExponents(Pcgs(N),List(x,Int))));
   isoNP := InverseGeneralMapping(isoPN);


   ##
   ## in the following loop we always have the following:
   ## recall GG is our original input group; below G is the current larger group, i.e.:
   ##
   ## G = GG/GG_{i+1} and Q=GG/GG_i and P = GG_i/GG_{i+1} with hom : G -> Q
   ## and H = can(GG/GG_i) with isoHQ : H -> Q and P is el-ab module
   ## and M is GModuleByMats
   ##

   Info(InfoCanonicalPcPres,1,"We need ", Length(lcp)-2," iterations\n");

   for i in [2..Size(lcp)-1] do

      Info(InfoCanonicalPcPres,1,"now start round ",i," group size ",Size(H)," module size ",Size(P));

      Info(InfoCanonicalPcPres,1,"         get 2-cohom");
      ttt:= Runtime();
      H2 := canform_my_getH2(H,M);
      Info(InfoCanonicalPcPres,1,"         done; get tail vector");

      ## now get tail of extension G/G_{i+1} of G/G_i \cong H by P\cong M \cong C_p^r
      ## as usual: lift gens from H=can(G/G_i) to G/G_{i+1} and evaluate relations
      ## make sure this is all moved between H=can(G/G_i) and Q=G/G_i etc

      F     := H2.presentation.group;;
      rel   := H2.presentation.relators;;
      pre   := List(Pcgs(H),x-> PreImagesRepresentative(hom,Image(isoHQ,x)));; #preims in G of gens of H


      ##### get tail vector; this rewriting is taken from MatrixOperationOfCPGroup
      ## evaluate preimages of gens in G/G_{i+1}
      ## need to adjust for non-central module due to structure of relations
      ## e.g., power relation a^4=b is encoded as a^4b^{-1}, and tail t will satisfy
      ##       a^4b^{-1}=t, so a^4=tb, but we need tails on the right, so a^4=bt^{b^-1}, etc
      ## 

      ords  := RelativeOrders(Pcgs(H));
      imgl := List(M.generators,x->x);
      k := 0;
      tails := [  ];
      for ii in [ 1 .. Length( pre ) ] do
         for jj in [ 1 .. ii ] do
            k := k + 1;
             tail := One(GF(p))*ExponentsOfPcElement(Pcgs(N),Image(isoPN,MappedWord( rel[k], GeneratorsOfGroup(F), pre )));
             if not IsBound( H2.module.isCentral ) or not H2.module.isCentral then
                if ii = jj then
                   m := imgl[ii] ^ ords[ii];
                else
                   m := imgl[ii] ^ imgl[jj];
                fi;
                tail := tail * m;
             fi;
             Add(tails,tail);
          od;
      od;
      tails := Flat(tails);
     ####


      ## get corresponding extension and set up isom between E(tails) and G
      Etails     := ExtensionSQ(H2.collector,H,M,tails);      
      isoTailtoG := GroupHomomorphismByImagesNC(Etails,G,Pcgs(Etails),
                                                Concatenation(pre,List(Pcgs(N),u->Image(isoNP,u ))));
      if isoTailtoG = fail then
         Error("round ",i," error with isoTailtoG");
      fi;

      Info(InfoCanonicalPcPresTimings,2,">>>> time H2 and tail: ",Runtime()-ttt);


      ## get canonical tail vector; catch the cases that H2 is trivial or comp-pairs
      ## are trivial, or tails is coboundary, or module is central,...
      if Size(Image(H2.cohom))=1 or Image(H2.cohom,tails)=Image(H2.cohom,0*tails) then
         Info(InfoCanonicalPcPres,1,"         done; have H^2=1");    
         #######
         Info(InfoCanonicalPcPres,1,"         now compute compatible pairs");
         ttt := Runtime();
         if canform_USE_NEW_CP then
             C := Group(Set(T.stab, DirectProductElement));
         else   
             C  := canform_compatible_pairs(H,M);;
         fi;
         Info(InfoCanonicalPcPresTimings,2,">>>> time CompPairs: ",Runtime()-ttt);
        ####### 

         orb := rec(orb:=[tails],
                    stab:=GeneratorsOfGroup(C),trans:=[[false,One(C)]]);
      else
        Info(InfoCanonicalPcPres,1,"         tail vector is not coboundary; H2 has dim ",Dimension(Image(H2.cohom)));   
        if ForAll(M.generators,IsOne)  then
           ttt := Runtime();
           Info(InfoCanonicalPcPres,1,"         central module, use AutM-orbit red");
           orb := canform_myorbit_reducedAutM(autH,H2, tails);
           Info(InfoCanonicalPcPresTimings,2,">>>> time can orbit reduced AutM: ",Runtime()-ttt);

        else
            Info(InfoCanonicalPcPres,1,"         now compute compatible pairs");
            ttt := Runtime();

         #######
            if canform_USE_NEW_CP then
               C := Group(Set(T.stab, DirectProductElement));
            else   
               C  := canform_compatible_pairs(H,M);;
            fi;
        ####### 
            Info(InfoCanonicalPcPresTimings,2,">>>> time CompPairs: ",Runtime()-ttt);
            Info(InfoCanonicalPcPres,1,"         done; now get matrix op and order");
           #this does exactly what canform_actWithCPonTail does, just maps it over with H2.cohom
            ttt  := Runtime();
            mats := MatrixOperationOfCPGroup( H2 , GeneratorsOfGroup( C ) );;
            Info(InfoCanonicalPcPresTimings,2,">>>> time CompPairs as mats: ",Runtime()-ttt);
            #Info(InfoCanonicalPcPres,1,"             order CP is ",Size(C));  # Size(C) can be slow here.
            Info(InfoCanonicalPcPres,1,"             now compute orbit");

            ## this used to be Size(C)=1... but Size(C) took ages for our comp pairs!
            if ForAll(GeneratorsOfGroup(C),x->x=One(C)) then  ## this can happen for p=2
               orb := rec(orb:=[tails],stab:=GeneratorsOfGroup(C),trans:=[[false,One(C)]]);
            else
               ttt  := Runtime();
               Mgrp := GroupByGenerators( mats );;
               orb  := canform_myorbit(mats,GeneratorsOfGroup(C),  Image( H2.cohom, tails ));
               Info(InfoCanonicalPcPresTimings,2,">>>> time InfoCanonicalPcPres_myorbit: ",Runtime()-ttt);
            fi;
         fi;

      fi; 



      Info(InfoCanonicalPcPres,1,"             done; orbit has size ",Size(orb.orb));
      ttt  := Runtime();
      pos  := PositionMinimum(orb.orb);
      can  := orb.orb[pos];
      cp   := orb.trans[pos][2];
      can  := H2.cohominv( can );;

      ## get 2-cob vector that is difference of starting coc tail and canonical tail can
      ## get isom from E(tails) to E(can)
      ## in contrast to paper, we map (mu^-1(g), \nu^-1(a)) -> (g,af(g))
      ##
      Info(InfoCanonicalPcPres,1,"             get coboundary and new isomorphisms");
      cob      := canform_actWithCPonTail(H2, cp, tails)-can;
      v        := canform_getCOBfunctionImages(H2.collector,H,M,-cob);

      pct      := Pcgs(Etails);
      n        := Size(Pcgs(H));
      alphai   := cp![1]^-1;
      betai    := List(cp![2]^-1,x->Concatenation(ListWithIdenticalEntries(n,0),List(x,Int)));
      d        := Length(betai[1]);
      hprei    := List( List(Pcgs(H),x-> Concatenation(ExponentsOfPcElement(Pcgs(H),Image(alphai,x)),
                                                     ListWithIdenticalEntries(M.dimension,0))),
                        ex -> PcElementByExponents(pct,ex));
      mprei    := List(betai,ex->PcElementByExponents(pct,ex));
      hprei    := Concatenation(hprei,mprei);

      ### now define it for the next round!
      ### now H = can(G/G_{i+1})
      oldH := H;
      H    := ExtensionSQ( H2.collector, H, M, can );;  
      ###
      ###

      #### continue with isom from E(tails) to H=E(can)
      pct  := Pcgs(H);
      new  := [];
      for ii in [1..Size(v)] do
         el := pct[ii];
         for j in [1..Size(v[ii])] do
            if v[ii][j]>0 then el := el * pct[n+j]^v[ii][j]; fi;
         od;
         Add(new,el);
      od;
      for ii in [n+1..Size(pct)] do Add(new,pct[ii]); od;
      isoCanToTail := GroupHomomorphismByImagesNC(H,Etails,new,hprei);
      if isoCanToTail = fail then
         Error("canform_compute_canonical: failed to construct isomorphism from canonical extension to tail extension");
      fi;




      isoHtoG := isoCanToTail*isoTailtoG; 
      oldG := G;

      Info(InfoCanonicalPcPresTimings,2,">>>> time get cob and isom: ",Runtime()-ttt);

      ## need another round:
      if i <= Size(lcp)-2 then     
         Info(InfoCanonicalPcPres,1,"             now get data for next step (isoms etc)");
         Info(InfoCanonicalPcPres,1,"             first, compute automorphism group");
         ttt  := Runtime();


         if not canform_USE_NEW_AUT then

             autH := AutomorphismGroup(H);

         else

            ## these are gens for Aut(Etails)
            autgens := canform_getautgroup(H2, orb.stab, Etails, tails,  oldH, M);
            ## now translate to gens for Aut(H)
            autgens := List(autgens, phi -> isoCanToTail*phi*InverseGeneralMapping(isoCanToTail));
            autH    := Group(autgens);
            AssignNiceMonomorphismAutomorphismGroup(autH, H);

            if canform_USE_NEW_AUT_test then
               Info(InfoCanonicalPcPres, 3, "aut done; we have ",Size(GeneratorsOfGroup(autH))," generators;  now start test");
               gap     := AutomorphismGroup(H);
               if not ForAll(autgens, x -> x in gap) then
                  Error("canform_compute_canonical: lifted generators are not all automorphisms");
               fi;
               if not Size(gap) = Size(autH) then
                  Error("canform_compute_canonical: lifted automorphism group has unexpected size");
               fi;
            fi;
         fi;

         Info(InfoCanonicalPcPres,1,"             done; number of Aut(H) generators: ",Size(GeneratorsOfGroup(autH)));
         Info(InfoCanonicalPcPresTimings,2,">>>> time Aut(H) ", Runtime()-ttt);
         homs  := NaturalHomomorphismByNormalSubgroup(GG,lcp[i+2]);
         G     := Image(homs);;
         P     := Image(homs,lcp[i+1]);
         hom   := NaturalHomomorphismByNormalSubgroup(G,P);
         Q     := Image(hom);

         isoOldGtoQ   := GroupHomomorphismByImagesNC(oldG,Q,Pcgs(oldG),Pcgs(Q));
         if isoOldGtoQ = fail then
            Error("canform_compute_canonical: failed to construct quotient map from previous round");
         fi;

         isoHQ := isoHtoG*isoOldGtoQ;
         isoQH := InverseGeneralMapping(isoHQ);
         p     := PrimePGroup(P);

         N     := AbelianGroup(List([1..RankPGroup(P)],x->p));
         isoPN := GroupHomomorphismByImagesNC(P,N,Pcgs(P),Pcgs(N));
         isoNP := GroupHomomorphismByImagesNC(N,P,Pcgs(N),Pcgs(P));

         Info(InfoCanonicalPcPres,1,"             second, get canonical module: group size ",Size(H), " module size ",Size(P));
      ##### canonical module
         ttt   := Runtime();
         T     := List(Pcgs(H), x ->  PreImagesRepresentative(hom,Image(isoHQ,x)));
         T     := List(T, x-> List(Pcgs(P),a-> ExponentsOfPcElement(Pcgs(P),a^x))*One(GF(p)));
         List(T, ConvertToMatrixRep);
         T     := canform_can_action(H, autH, T, p);
         Info(InfoCanonicalPcPresTimings,2,">>>> time next round canform_can_action: ",Runtime()-ttt);
         Info(InfoCanonicalPcPres,1,"             done;"); 
         M     := GModuleByMats(T.canact,GF(p));
       ##recall that T.trans = [a,b] such that (a,b) * oldT = T.canact
         isoQH := isoQH*T.trans[1];
         isoHQ := InverseGeneralMapping(isoQH);
         isoPN := isoPN*GroupHomomorphismByImagesNC(N,N,Pcgs(N),List(T.trans[2],x->PcElementByExponents(Pcgs(N),List(x,Int))));
         isoNP := InverseGeneralMapping(isoPN);


      else # no further iteration required

         isoHQ := isoHtoG;
         if not Source(isoHQ) = H or not Image(isoHQ) = GG or not Size(Kernel(isoHQ)) = 1 then
            Error("canform_compute_canonical: final isomorphism validation failed");
         fi;

      fi;

   od;


   return InverseGeneralMapping(isoHQ);
end);


####
#### now the main functions
####


InstallGlobalFunction(CanonicalPcGroup, function(G)
   local iso;
   iso := IsomorphismCanonicalPcGroup(G);
   Assert(0, IsBijective(iso));
   return Range(iso);
end);



InstallMethod(IsomorphismCanonicalPcGroup, [IsPcGroup],
function(G)
   return canform_compute_canonical(G);
end);

InstallMethod(IsomorphismCanonicalPcGroup, [IsGroup],
function(G)
   local isoToPc, isoPcToCan, res;

   if not IsFinite(G) then
      Error("<G> must be finite");
   fi;

   if not IsSolvableGroup(G) then
      Error("<G> must be solvable");
   fi;

   isoToPc := IsomorphismSpecialPcGroup(G);
   if isoToPc = fail then
     Error("IsomorphismCanonicalPcGroup: failed to convert input group to a special pc group");
   fi;

   isoPcToCan := IsomorphismCanonicalPcGroup(Image(isoToPc));
   res := GroupHomomorphismByImagesNC(G,Image(isoPcToCan),GeneratorsOfGroup(G),
              List(GeneratorsOfGroup(G),u->Image(isoPcToCan,Image(isoToPc,u))));
   return res;
end);
