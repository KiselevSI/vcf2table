#!/usr/bin/env bash
# build_phage_db_conda.sh <taxid|assembly_accession> <snpEff_ID>
# One‑click builder for a custom SnpEff database in the current conda env.

set -euo pipefail
[[ $# -eq 2 ]] || { echo "Usage: $0 <taxid|assembly_accession> <snpEff_ID>"; exit 1; }
TAXID="$1"; DBID="$2"

# ── locate snpEff ───────────────────────────────────────────────────────────
SNPEFF_BIN=$(command -v snpEff || true)
[[ -x $SNPEFF_BIN ]] || { echo "[E] snpEff not in PATH" >&2; exit 1; }
CONDA_PREFIX=$(dirname $(dirname "$SNPEFF_BIN"))
shopt -s nullglob; snpeff_dirs=("$CONDA_PREFIX"/share/snpeff*); shopt -u nullglob
[[ ${#snpeff_dirs[@]} -gt 0 ]] || { echo "[E] snpEff jar dir not found under $CONDA_PREFIX/share" >&2; exit 2; }
SNPEFF_HOME="${snpeff_dirs[0]}"
CONFIG_FILE="$SNPEFF_HOME/snpEff.config"
[[ -f $CONFIG_FILE ]] || touch "$CONFIG_FILE"

# ── temp workspace ──────────────────────────────────────────────────────────
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT; cd "$TMP"

# ── download genome & annotation ────────────────────────────────────────────
mapfile -t ftp_paths < <(esearch -db assembly -query "$TAXID" | esummary | \
    xtract -pattern DocumentSummary -element FtpPath_RefSeq,FtpPath_GenBank | tr '\t' '\n' | grep ^ftp)
[[ ${#ftp_paths[@]} -gt 0 ]] || { echo "[E] No FTP path for '$TAXID'" >&2; exit 3; }
for ftp in "${ftp_paths[@]}"; do
  base=$(basename "$ftp")
  fna="${base}_genomic.fna.gz"; gff="${base}_genomic.gff.gz"; gbk="${base}_genomic.gbff.gz"
  wget -q "$ftp/$fna" || continue
  if wget -q "$ftp/$gff"; then ann=$gff; mode=gff3; break
  elif wget -q "$ftp/$gbk"; then ann=$gbk; mode=genbank; break; fi
done
[[ -f $fna && -f $ann ]] || { echo "[E] Download failed" >&2; exit 3; }
UNZIP="gzip -dc"; command -v pigz &>/dev/null && UNZIP="pigz -dc"
$UNZIP "$fna" > genome.fa
$UNZIP "$ann" > annotation.tmp
[[ $mode == gff3 ]] && head -n1 annotation.tmp | grep -q '^##gff' || mode=genbank

# ── install files ───────────────────────────────────────────────────────────
DB_DIR="$SNPEFF_HOME/data/$DBID"; mkdir -p "$DB_DIR"
cp genome.fa "$DB_DIR/sequences.fa"
[[ $mode == gff3 ]] && mv annotation.tmp "$DB_DIR/genes.gff" || mv annotation.tmp "$DB_DIR/genes.gbk"

# ── update config ───────────────────────────────────────────────────────────
sed -i "/^${DBID}\.codonTable/d" "$CONFIG_FILE"
if ! grep -q "^${DBID}\.genome" "$CONFIG_FILE"; then echo "${DBID}.genome : ${DBID}" >> "$CONFIG_FILE"; fi
if [[ $mode == gff3 ]]; then
  codon=$(grep -Pom1 'transl_table=\K[0-9]+' "$DB_DIR/genes.gff" || true)
  declare -A map=([11]=Bacterial_and_Plant_Plastid [1]=Standard)
  [[ -n ${map[$codon]:-} ]] && echo "${DBID}.codonTable : ${map[$codon]}" >> "$CONFIG_FILE"
fi

# ── build snpEff DB (skip heavy checks) ─────────────────────────────────────
cd "$SNPEFF_HOME"
snpEff build -noCheckCds -noCheckProtein -dataDir "$SNPEFF_HOME/data" ${mode:+-$mode} -v "$DBID"
[[ -f "$DB_DIR/snpEffectPredictor.bin" ]] || { echo "[E] snpEffectPredictor.bin missing" >&2; exit 4; }

echo "[OK] SnpEff database '$DBID' built at $DB_DIR"
