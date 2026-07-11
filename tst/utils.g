# a bunch of examples and test functions


################################################################################
## input:  pc-group
## output: "random" isomorphic copy of the group
canform_random_copy := function(G)
  if not IsPcGroup(G) then return G; fi;
  return PcGroupCode(RandomSpecialPcgsCoded(G), Size(G));
end;

canform_random_copy_old := function(G)
local new, max, C, M, Gpc, g, hom, Q;

  if not IsPcGroup(G) then return G; fi;
  new := [];
  if Size(G) = 1 then return G; fi;
  max := MaximalNormalSubgroups(G);
  C   := G;
  M   := Random(max);
  while Size(M)>1 do
     hom    := NaturalHomomorphismByNormalSubgroup(C,M);
     Q      := Image(hom);
     repeat g := Random(Q); until not g = g^0;
     Add(new,PreImagesRepresentative(hom,g));
     C   := M;
     max := MaximalNormalSubgroups(M);
     M   := Random(max);
  od;
  repeat g:= Random(C); until not g in M;
  Add(new, g);
  new := PcgsByPcSequence(FamilyObj(One(G)),new);
  Assert(0, Product(RelativeOrders(new)) = Size(G));
  Gpc := PcGroupWithPcgs(new);
 #if not IdSmallGroup(Gpc)=IdSmallGroup(G) then Error("ups"); fi;
  return Gpc;
end;

################################################################################
## tests a single group G: makes 'random copies' and checks if pc code of can form
## is always the same (i.e. same can pc pres)
## input: group G, or group G and integer (nr) of tests 
##
canform_testSingleGroup_solv := function(arg)
local grps, H, pres, id, G, nr, r, i, code, hom, new, tt, lll,ll,stt;

   lll := InfoLevel(InfoCanonicalPcPres);
   ll  := InfoLevel(InfoCanonicalPcPresTimings);
   SetInfoLevel(InfoCanonicalPcPres,0);
   SetInfoLevel(InfoCanonicalPcPresTimings,0);

   Info(InfoCanonicalPcPres, 1, "get copies of group");
   if Length(arg)=1 then G:=arg[1]; nr:=10; fi;
   if Length(arg)=2 then G:=arg[1]; nr:=arg[2]; fi;
   grps := List([1..nr],x-> canform_random_copy_old(G));
   Info(InfoCanonicalPcPres, 1, "got copies, now get can form");
   code := false;
   if nr = 1 then
      Error("if nr=1, then there's no real test: we need at least 2 groups\n",
            "to see that they both have the same canonical form... \n",
            "please increase nr to be at least 2");
   fi;
   for i in [1..nr] do
      Info(InfoCanonicalPcPres, 1, "    start solvable group ",i," with code ",CodePcGroup(grps[i]));
      tt := Runtime();
      stt := State(GlobalMersenneTwister);
      hom := IsomorphismCanonicalPcGroup(grps[i]);
      if i = 1 then
         Info(InfoCanonicalPcPresTimings, 1, "    runtime for first cf: ",Runtime()-tt);
      fi;
      H   := Range(hom);
      hom := GroupHomomorphismByImages(grps[i],H,Pcgs(grps[i]),List(Pcgs(grps[i]),i->Image(hom,i)));
      if hom=fail or Size(Kernel(hom))>1 then Error("iso is not an iso!"); fi;
      new := CodePcGroup(H);
      Info(InfoCanonicalPcPres, 1, "    computed can form; got code ",new);
      if code = false then
         code:=new;
      else
         if not code=new then Error("CODE INCORRECT; rnd state stt"); fi;
      fi;
   od;
   SetInfoLevel(InfoCanonicalPcPres,lll);
   SetInfoLevel(InfoCanonicalPcPresTimings,ll);
   return true;
end;

################################################################################
## tests a single group G: makes one random copie and checks if pc code of can form
## is always the same (i.e. same can pc pres)
## input: group G, or group G and integer (nr) of tests 
##
## no prints, DISPLAY TIME of first test, do only 2 tests
##
canform_testSingleGroup_solv_time := function(G,max)
local grps, H, pres, id, nr, r, i, code, hom, new, tt, lll;
   lll := InfoLevel(InfoCanonicalPcPres);
   SetInfoLevel(InfoCanonicalPcPres,0);
   grps := List([1..2],x-> canform_random_copy(G));
   code := false;
   id   := IdSmallGroup(G);
   for i in [1..2] do
     #Print("    start solvable group ",i," with code ",CodePcGroup(grps[i]),"\n");
      tt := Runtime();
      hom := IsomorphismCanonicalPcGroup(grps[i]);
      if i = 1 then
         tt := Runtime()-tt;
         if tt>max[1] then max:=[tt,id]; fi;
         Info(InfoCanonicalPcPresTimings, 1, "runtime for ",id," is ", tt," and max is ",max,"\n");
      fi;
      H   := Range(hom);
      hom := GroupHomomorphismByImages(grps[i],H,Pcgs(grps[i]),List(grps[i],i->Image(hom,i)));
      if hom=fail or Size(Kernel(hom))>1 then Error("iso is not an iso!"); fi;
      new := CodePcGroup(H);
     #Print("    computed can form; got code ",new,"\n");
      if code = false then
         code:=new;
      else
         if not code=new then Error("CODE INCORRECT"); fi;
      fi;
   od;
   SetInfoLevel(InfoCanonicalPcPres,lll);
   return max;
end;


canform_test_time := function(start,finish)
local i, G, max, t,n;
   max := [0,0];
   for n in [start..finish] do
      if IsPrimePowerInt(n) then
         Info(InfoCanonicalPcPres, 1, "skip prime power");
      else
         for i in [1..NumberSmallGroups(n)] do
            G := SmallGroup(n,i);
            if IsSolvableGroup(G) then
               t := canform_testSingleGroup_solv_time(G,max);
               if t[1] > max[1] then max:=t; fi;
            fi;
         od;
       fi;
    od;
 end;   


################################################################################
## tests groups in the database
## input: order, or order and nr, or order and nr and start
##
canform_testAll_solv := function(arg)
local i, G, ord, nr, start;
   if Length(arg)=1 then
      ord:=arg[1]; nr:=3; start:=1;
   elif Length(arg)=2 then
      ord:=arg[1]; nr := arg[2]; start:=1;
   elif Length(arg)=3 then
      ord:=arg[1]; nr:=arg[2]; start:=arg[3];
   fi;
   for i in [start..NumberSmallGroups(ord)] do
      G := SmallGroup(ord,i);
      if IsSolvableGroup(G) then
         Info(InfoCanonicalPcPres, 1, "------------------------------------------- start: ",[ord,i]);
         canform_testSingleGroup_solv(G,nr);
      fi;
   od;
return true;
end;

canform_testAll_solv_skip_nilpotent:= function(ord,nr...)
local i, G;
   if IsPrimePowerInt(ord) then return; fi;
   if IsEmpty(nr) then nr := 2; else nr := nr[1]; fi;
   for i in [1..NumberSmallGroups(ord)] do
      G := SmallGroup(ord,i);
      if IsNilpotentGroup(G) then continue; fi;
      if IsSolvableGroup(G) then
         Info(InfoCanonicalPcPres, 1, "------------------------------------------- start: ",[ord,i]);
         canform_testSingleGroup_solv(G,nr);
      fi;
   od;
return true;
end;
