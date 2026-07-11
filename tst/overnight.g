#
# CanonicalPcPres overnight stress tests.
#
# This file is intentionally not part of tst/testall.g.  It is meant to be
# run manually, possibly many times in parallel with different seeds or order
# ranges.
#
# Example:
#   gap -q --quitonbreak --packagedirs . tst/overnight.g
#
# Useful environment variables:
#   CANFORM_OVERNIGHT_START=1
#   CANFORM_OVERNIGHT_FINISH=383
#   CANFORM_OVERNIGHT_CODE_REPEATS=100
#   CANFORM_OVERNIGHT_ISO_REPEATS=5
#   CANFORM_OVERNIGHT_SEED=123456789
#   CANFORM_OVERNIGHT_INCLUDE_NILPOTENT=false
#   CANFORM_OVERNIGHT_SKIP_MULTIPLES_OF_64=true
#   CANFORM_OVERNIGHT_SHARDS=1
#   CANFORM_OVERNIGHT_SHARD=1
#   CANFORM_OVERNIGHT_LIMIT=0
#   CANFORM_OVERNIGHT_PROGRESS_DIR=
#
# Set CANFORM_OVERNIGHT_CODE_REPEATS=0 or CANFORM_OVERNIGHT_ISO_REPEATS=0
# to disable one of the two passes.
#
# To use 32 independent shell jobs, run with CANFORM_OVERNIGHT_SHARDS=32
# and CANFORM_OVERNIGHT_SHARD set to 1, 2, ..., 32 in the individual jobs.
#

LoadPackage("CanonicalPcPres");
ReadPackage("CanonicalPcPres", "tst/utils.g");

canform_overnight_env := function(name, default)
  if IsBound(GAPInfo.SystemEnvironment.(name)) then
    return GAPInfo.SystemEnvironment.(name);
  fi;
  return default;
end;

canform_overnight_int_env := function(name, default)
local value;
  value := canform_overnight_env(name, fail);
  if value = fail or value = "" then
    return default;
  fi;
  value := Int(value);
  if value = fail then
    Error("environment variable ", name, " must be an integer");
  fi;
  return value;
end;

canform_overnight_bool_env := function(name, default)
local value;
  value := canform_overnight_env(name, fail);
  if value = fail or value = "" then
    return default;
  fi;
  value := LowercaseString(value);
  if value in ["1", "true", "yes", "y", "on"] then
    return true;
  fi;
  if value in ["0", "false", "no", "n", "off"] then
    return false;
  fi;
  Error("environment variable ", name, " must be boolean, not ", value);
end;

canform_overnight_config := rec(
  start := canform_overnight_int_env("CANFORM_OVERNIGHT_START", 1),
  finish := canform_overnight_int_env("CANFORM_OVERNIGHT_FINISH", 383),
  codeRepeats := canform_overnight_int_env("CANFORM_OVERNIGHT_CODE_REPEATS", 100),
  isoRepeats := canform_overnight_int_env("CANFORM_OVERNIGHT_ISO_REPEATS", 5),
  seed := canform_overnight_int_env("CANFORM_OVERNIGHT_SEED",
                                    NanosecondsSinceEpoch()),
  includeNilpotent := canform_overnight_bool_env(
                        "CANFORM_OVERNIGHT_INCLUDE_NILPOTENT", false),
  skipMultiplesOf64 := canform_overnight_bool_env(
                        "CANFORM_OVERNIGHT_SKIP_MULTIPLES_OF_64", true),
  shards := canform_overnight_int_env("CANFORM_OVERNIGHT_SHARDS", 1),
  shard := canform_overnight_int_env("CANFORM_OVERNIGHT_SHARD", 1),
  limit := canform_overnight_int_env("CANFORM_OVERNIGHT_LIMIT", 0),
  progressDir := canform_overnight_env("CANFORM_OVERNIGHT_PROGRESS_DIR", ""),
);

canform_overnight_copy_functions := [
  rec(name := "canform_random_copy", fun := canform_random_copy),
  rec(name := "canform_random_copy_old", fun := canform_random_copy_old),
];

canform_overnight_validate_config := function(config)
  if config.start < 1 then
    Error("CANFORM_OVERNIGHT_START must be positive");
  fi;
  if config.finish < config.start then
    Error("CANFORM_OVERNIGHT_FINISH must be at least START");
  fi;
  if config.codeRepeats <> 0 and config.codeRepeats < 2 then
    Error("CANFORM_OVERNIGHT_CODE_REPEATS must be 0 or at least 2");
  fi;
  if config.isoRepeats < 0 then
    Error("CANFORM_OVERNIGHT_ISO_REPEATS must be non-negative");
  fi;
  if config.codeRepeats = 0 and config.isoRepeats = 0 then
    Error("at least one overnight test pass must be enabled");
  fi;
  if config.shards < 1 then
    Error("CANFORM_OVERNIGHT_SHARDS must be positive");
  fi;
  if config.shard < 1 or config.shards < config.shard then
    Error("CANFORM_OVERNIGHT_SHARD must be in [1..SHARDS]");
  fi;
  if config.limit < 0 then
    Error("CANFORM_OVERNIGHT_LIMIT must be non-negative");
  fi;
end;

canform_overnight_reset_random_sources := function(seed)
  Reset(GlobalMersenneTwister, seed);
  Reset(GlobalRandomSource, seed);
end;

canform_overnight_fail := function(config, id, copyName, pass, iteration,
                                   state, message)
  Error("overnight test failed: ", message, "\n",
        "  pass: ", pass, "\n",
        "  group: ", id, "\n",
        "  copy function: ", copyName, "\n",
        "  iteration: ", iteration, "\n",
        "  seed: ", config.seed, "\n",
        "  GlobalMersenneTwister state before iteration: ", state, "\n");
end;

canform_overnight_checked_iso := function(config, id, copyName, iteration, iso)
local G, H, gens, checked, state;
  state := State(GlobalMersenneTwister);
  G := Source(iso);
  H := Range(iso);
  gens := GeneratorsOfGroup(G);
  checked := GroupHomomorphismByImages(G, H, gens,
               List(gens, g -> Image(iso, g)));

  if checked = fail then
    canform_overnight_fail(config, id, copyName, "iso", iteration, state,
                           "checked homomorphism construction returned fail");
  fi;
  if not IsBijective(checked) then
    canform_overnight_fail(config, id, copyName, "iso", iteration, state,
                           "checked homomorphism is not bijective");
  fi;
  if Size(Kernel(checked)) <> 1 then
    canform_overnight_fail(config, id, copyName, "iso", iteration, state,
                           "checked homomorphism has non-trivial kernel");
  fi;
  if Size(Image(checked)) <> Size(H) then
    canform_overnight_fail(config, id, copyName, "iso", iteration, state,
                           "checked homomorphism is not onto its range");
  fi;
  return true;
end;

canform_overnight_check_code := function(config, G, id, reference)
local copyData, i, copy, H, code, state;
  for copyData in canform_overnight_copy_functions do
    for i in [1..config.codeRepeats] do
      state := State(GlobalMersenneTwister);
      copy := copyData.fun(G);
      H := CanonicalPcGroup(copy);
      if Size(H) <> Size(G) then
        canform_overnight_fail(config, id, copyData.name, "code", i, state,
                               "canonical group has wrong size");
      fi;
      code := CodePcGroup(H);
      if code <> reference then
        canform_overnight_fail(config, id, copyData.name, "code", i, state,
                               Concatenation("expected code ",
                                             String(reference),
                                             ", got ", String(code)));
      fi;
    od;
  od;
  return true;
end;

canform_overnight_check_iso := function(config, G, id, reference)
local copyData, i, copy, iso, H, code, state;
  for copyData in canform_overnight_copy_functions do
    for i in [1..config.isoRepeats] do
      state := State(GlobalMersenneTwister);
      copy := copyData.fun(G);
      iso := IsomorphismCanonicalPcGroup(copy);
      H := Range(iso);
      if Size(H) <> Size(G) then
        canform_overnight_fail(config, id, copyData.name, "iso", i, state,
                               "canonical group has wrong size");
      fi;
      canform_overnight_checked_iso(config, id, copyData.name, i, iso);
      code := CodePcGroup(H);
      if code <> reference then
        canform_overnight_fail(config, id, copyData.name, "iso", i, state,
                               Concatenation("expected code ",
                                             String(reference),
                                             ", got ", String(code)));
      fi;
    od;
  od;
  return true;
end;

canform_overnight_group_is_selected := function(config, groupNr)
  return ((groupNr - config.shard) mod config.shards) = 0;
end;

canform_overnight_order_is_skipped := function(config, order)
  if not SmallGroupsAvailable(order) then
    return true;
  fi;
  if not config.includeNilpotent and IsPrimePowerInt(order) then
    return true;
  fi;
  if config.skipMultiplesOf64 and order mod 64 = 0 then
    return true;
  fi;
  return false;
end;

canform_overnight_worklist := function(config)
local worklist, order, number, i, G;
  worklist := [];
  for order in [config.start..config.finish] do
    if canform_overnight_order_is_skipped(config, order) then
      continue;
    fi;
    number := NumberSmallGroups(order);
    for i in [1..number] do
      G := SmallGroup(order, i);
      if not IsSolvableGroup(G) then
        continue;
      fi;
      if not config.includeNilpotent and IsNilpotentGroup(G) then
        continue;
      fi;
      Add(worklist, [order, i]);
    od;
  od;
  return worklist;
end;

canform_overnight_progress_path := function(config, subdir, name)
  return Concatenation(config.progressDir, "/", subdir, "/", name);
end;

canform_overnight_progress_tsv := function(config)
  return Concatenation(config.progressDir, "/shard-",
                      String(config.shard), ".tsv");
end;

canform_overnight_progress_append := function(config, event, id)
  if config.progressDir = "" then
    return;
  fi;
  AppendTo(canform_overnight_progress_tsv(config),
           Runtime(), "\t", event, "\t", id[1], "\t", id[2], "\n");
end;

canform_overnight_progress_current := function(config, id)
  if config.progressDir = "" then
    return;
  fi;
  PrintTo(canform_overnight_progress_path(config, "current",
            Concatenation("shard-", String(config.shard))),
          id[1], "\t", id[2], "\n");
  canform_overnight_progress_append(config, "start", id);
end;

canform_overnight_progress_done := function(config, id)
  if config.progressDir = "" then
    return;
  fi;
  PrintTo(canform_overnight_progress_path(config, "groups",
            Concatenation("grp-", String(id[1]), "-", String(id[2]), ".done")),
          "done\t", id[1], "\t", id[2], "\n");
  canform_overnight_progress_append(config, "done", id);
end;

canform_overnight_run := function(config)
local oldInfo, oldTimingInfo, groupNr, G, id, reference, tested, startRuntime,
      worklist, selected;
  canform_overnight_validate_config(config);
  canform_overnight_reset_random_sources(config.seed);

  oldInfo := InfoLevel(InfoCanonicalPcPres);
  oldTimingInfo := InfoLevel(InfoCanonicalPcPresTimings);
  SetInfoLevel(InfoCanonicalPcPres, 0);
  SetInfoLevel(InfoCanonicalPcPresTimings, 0);

  tested := 0;
  startRuntime := Runtime();

  Print("CanonicalPcPres overnight test configuration:\n");
  Print("  orders: ", config.start, "..", config.finish, "\n");
  Print("  code repeats per copy function: ", config.codeRepeats, "\n");
  Print("  iso repeats per copy function: ", config.isoRepeats, "\n");
  Print("  include nilpotent groups: ", config.includeNilpotent, "\n");
  Print("  skip prime-power orders: ", not config.includeNilpotent, "\n");
  Print("  skip multiples of 64: ", config.skipMultiplesOf64, "\n");
  if config.progressDir <> "" then
    Print("  progress directory: ", config.progressDir, "\n");
  fi;
  Print("  seed: ", config.seed, "\n");
  Print("  shard: ", config.shard, " of ", config.shards, "\n");
  if config.limit > 0 then
    Print("  limit: ", config.limit, " selected groups\n");
  fi;

  worklist := canform_overnight_worklist(config);
  selected := Filtered([1..Length(worklist)],
                       n -> canform_overnight_group_is_selected(config, n));
  Print("  eligible group ids: ", Length(worklist), "\n");
  Print("  group ids assigned to this shard: ", Length(selected), "\n");

  for groupNr in selected do
    id := worklist[groupNr];
    canform_overnight_progress_current(config, id);
    G := SmallGroup(id[1], id[2]);
    reference := CodePcGroup(CanonicalPcGroup(G));
    tested := tested + 1;
    Print("testing group ", id, " as selected group ", tested,
          " with reference code ", reference, "\n");
    if config.codeRepeats > 0 then
      canform_overnight_check_code(config, G, id, reference);
    fi;
    if config.isoRepeats > 0 then
      canform_overnight_check_iso(config, G, id, reference);
    fi;
    canform_overnight_progress_done(config, id);

    if config.limit > 0 and tested >= config.limit then
      break;
    fi;
  od;

  SetInfoLevel(InfoCanonicalPcPres, oldInfo);
  SetInfoLevel(InfoCanonicalPcPresTimings, oldTimingInfo);

  Print("CanonicalPcPres overnight tests completed: ", tested,
        " selected groups in ", Runtime() - startRuntime, " ms\n");
  return true;
end;

canform_overnight_run(canform_overnight_config);
