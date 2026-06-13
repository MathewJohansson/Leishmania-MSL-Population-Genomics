#!/usr/bin/env bash
# run_splitstree_all.sh
# Generates SplitsTree6 NeighbourNet .stree6 files for all NJ tree subsets.
# Step 1: Convert any missing tetraploid .txt distance matrices to .nex
# Step 2: Run workflow-run on each .nex to produce .stree6 output
#
# SplitsTree6 version: 6.7.5
# Workflow: neighbournet_workflow.stree6 (NeighbourNet from distance matrix)
# Output .stree6 files must be opened in SplitsTree6 GUI for image export.

set -e

# ─── PATHS ────────────────────────────────────────────────────────────────────

BASE="/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
DATA="${BASE}/Data/Phylogenetic_Data/NJ_Trees"
WORKFLOW="${BASE}/Scripts/Phylogenetic Scripts/NJ Trees/neighbournet_workflow.stree6"
PYTHON="python3"
NEI_SCRIPT="${BASE}/neis_txt_to_nex.py"
SPLITSTREE="/home/johansson/splitstreeapp/tools/workflow-run"

# ─── STEP 1: CONVERT .txt → .nex WHERE NEEDED ─────────────────────────────────

echo "=== Step 1: Converting distance matrices to NEXUS format ==="

# Whole-genome LO
$PYTHON "$NEI_SCRIPT" \
  "${DATA}/Whole-Genome/L.infantum-plus-outgroup/diploid_whole_genome_Linfantum_plus_outgroup_Neis_distance.txt" \
  "${DATA}/Whole-Genome/L.infantum-plus-outgroup/diploid_whole_genome_Linfantum_plus_outgroup_Neis_distance.nex"

# Whole-genome LI
$PYTHON "$NEI_SCRIPT" \
  "${DATA}/Whole-Genome/L.infantum-only/diploid_whole_genome_Linfantum_only_Neis_distance.txt" \
  "${DATA}/Whole-Genome/L.infantum-only/diploid_whole_genome_Linfantum_only_Neis_distance.nex"

# Whole-genome AM
$PYTHON "$NEI_SCRIPT" \
  "${DATA}/Whole-Genome/Americas-only/diploid_whole_genome_Americas_only_Neis_distance.txt" \
  "${DATA}/Whole-Genome/Americas-only/diploid_whole_genome_Americas_only_Neis_distance.nex"

# 100kb LI (tetraploid)
$PYTHON "$NEI_SCRIPT" \
  "${DATA}/100kb_windows/L.infantum-only/tetraploid_100kb_Linfantum_only_Neis_distance.txt" \
  "${DATA}/100kb_windows/L.infantum-only/tetraploid_100kb_Linfantum_only_Neis_distance.nex"

# 100kb AM (tetraploid)
$PYTHON "$NEI_SCRIPT" \
  "${DATA}/100kb_windows/Americas-only/tetraploid_100kb_Americas_only_Neis_distance.txt" \
  "${DATA}/100kb_windows/Americas-only/tetraploid_100kb_Americas_only_Neis_distance.nex"

# 20kb LI (tetraploid)
$PYTHON "$NEI_SCRIPT" \
  "${DATA}/20kb_windows/L.infantum-only/tetraploid_20kb_Linfantum_only_Neis_distance.txt" \
  "${DATA}/20kb_windows/L.infantum-only/tetraploid_20kb_Linfantum_only_Neis_distance.nex"

# 20kb AM (tetraploid)
$PYTHON "$NEI_SCRIPT" \
  "${DATA}/20kb_windows/Americas-only/tetraploid_20kb_Americas_only_Neis_distance.txt" \
  "${DATA}/20kb_windows/Americas-only/tetraploid_20kb_Americas_only_Neis_distance.nex"

echo "=== Step 1 complete ==="

# ─── STEP 2: RUN WORKFLOW-RUN ON EACH .nex ────────────────────────────────────

echo ""
echo "=== Step 2: Running SplitsTree6 NeighbourNet on all subsets ==="

run_splitstree() {
  local label="$1"
  local input_nex="$2"
  local output_stree6="$3"

  # Use symlinks to avoid spaces in paths
  local tmp_in="/tmp/st_input_$$.nex"
  local tmp_wf="/tmp/st_workflow_$$.stree6"

  ln -sf "$input_nex" "$tmp_in"
  ln -sf "$WORKFLOW" "$tmp_wf"

  echo "--- Running: $label ---"
  "$SPLITSTREE" -w "$tmp_wf" -i "$tmp_in" -o "$output_stree6"

  rm -f "$tmp_in" "$tmp_wf"
  echo "--- Done: $output_stree6 ---"
  echo ""
}

# Whole-genome LO
run_splitstree "Whole-genome L.infantum + outgroup" \
  "${DATA}/Whole-Genome/L.infantum-plus-outgroup/diploid_whole_genome_Linfantum_plus_outgroup_Neis_distance.nex" \
  "${DATA}/Whole-Genome/L.infantum-plus-outgroup/diploid_whole_genome_Linfantum_plus_outgroup_SplitsTree.stree6"

# Whole-genome LI
run_splitstree "Whole-genome L.infantum-only" \
  "${DATA}/Whole-Genome/L.infantum-only/diploid_whole_genome_Linfantum_only_Neis_distance.nex" \
  "${DATA}/Whole-Genome/L.infantum-only/diploid_whole_genome_Linfantum_only_SplitsTree.stree6"

# Whole-genome AM
run_splitstree "Whole-genome Americas-only" \
  "${DATA}/Whole-Genome/Americas-only/diploid_whole_genome_Americas_only_Neis_distance.nex" \
  "${DATA}/Whole-Genome/Americas-only/diploid_whole_genome_Americas_only_SplitsTree.stree6"

# Chr31 LI
run_splitstree "Chr31 L.infantum-only" \
  "${DATA}/Chr31/L.infantum-only/tetraploid_chr31_Linfantum_only_Neis_distance.nex" \
  "${DATA}/Chr31/L.infantum-only/tetraploid_chr31_Linfantum_only_SplitsTree.stree6"

# Chr31 AM
run_splitstree "Chr31 Americas-only" \
  "${DATA}/Chr31/Americas-only/tetraploid_chr31_Americas_Neis_distance.nex" \
  "${DATA}/Chr31/Americas-only/tetraploid_chr31_Americas_SplitsTree.stree6"

# 100kb LI
run_splitstree "100kb L.infantum-only" \
  "${DATA}/100kb_windows/L.infantum-only/tetraploid_100kb_Linfantum_only_Neis_distance.nex" \
  "${DATA}/100kb_windows/L.infantum-only/tetraploid_100kb_Linfantum_only_SplitsTree.stree6"

# 100kb AM
run_splitstree "100kb Americas-only" \
  "${DATA}/100kb_windows/Americas-only/tetraploid_100kb_Americas_only_Neis_distance.nex" \
  "${DATA}/100kb_windows/Americas-only/tetraploid_100kb_Americas_only_SplitsTree.stree6"

# 20kb LI
run_splitstree "20kb L.infantum-only" \
  "${DATA}/20kb_windows/L.infantum-only/tetraploid_20kb_Linfantum_only_Neis_distance.nex" \
  "${DATA}/20kb_windows/L.infantum-only/tetraploid_20kb_Linfantum_only_SplitsTree.stree6"

# 20kb AM
run_splitstree "20kb Americas-only" \
  "${DATA}/20kb_windows/Americas-only/tetraploid_20kb_Americas_only_Neis_distance.nex" \
  "${DATA}/20kb_windows/Americas-only/tetraploid_20kb_Americas_only_SplitsTree.stree6"

echo "=== All done! ==="
echo ""
echo "Open each .stree6 file in SplitsTree6 GUI and export as PDF/SVG/PNG."
