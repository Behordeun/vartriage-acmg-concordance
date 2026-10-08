#!/usr/bin/env bash
# Point-system eRepo re-classification to produce per-variant tags for the PM2 sensitivity.
# Writes a fresh classifications file WITHOUT clobbering the committed (old-era) one.
set -o pipefail

REPO=/Users/muhammad/Documents/DevProjects/personal_projects/Bioinformatics_Libraries
ENG=$REPO/vartriage
REF=$ENG/data/references
EREPO=$REPO/vartriage-streaming-acmg/data/erepo
OUT=$EREPO/erepo_classifications_pointsystem.json
RC=$REPO/vartriage-acmg-concordance/results/.pm2_rerun_rc

cd "$ENG" || { echo 999 > "$RC"; exit 999; }

uv run vartriage \
  --vcf "$EREPO/erepo_variants.vcf.gz" \
  --gene-annotation "$REF/gencode.v46.annotation.gtf" \
  --clinvar "$REF/clinvar.tsv" \
  --revel-scores "$EREPO/revel_erepo.tsv" \
  --gnomad-remote gnomad-exomes-v4-grch38 \
  --remote-cache-ttl -1 \
  --spliceai-db /Users/Shared/open-cravat/modules/annotators/spliceai/data/spliceai.sqlite \
  --output "$OUT" \
  --output-format json \
  --no-confirm \
  > "$REPO/vartriage-acmg-concordance/results/.pm2_rerun.log" 2>&1

echo $? > "$RC"
