# shellcheck shell=bash
# OpenCRAVAT annotation and matched-control read metrics.
# Sourced by CircuSNV.sh; uses its global settings.

# Annotate with OpenCRAVAT.
# Usage: annotate <input.txt (tab-delimited variants)> <run name> <output format>...
annotate() {
  local input=$1
  local name=$2
  shift 2
  echo -e "annotating $input with OpenCravat..."
  activate_env
  oc run "$input" \
    -l hg38 \
    -n "$name" --silent --mp "$THREADS" \
    -a cgc cadd_exome clinpred cosmic mutation_assessor thousandgenomes gnomad3 repeat \
    -t "$@"
  deactivate_env
}

# Function to add reads info from matched control to each mutation found in tumor sample 
add_matched_control_metrics() {
  local mutations="$1"
  local samplename="$2"
  local matched_control ctrl

  echo "[MATCHED CTRL] Adding matched control metrics for each mutation..."

  # find matched control (strip CR if present)
  matched_control=$(
    awk -F $'\t' -v s="$samplename" '$1==s{gsub(/\r/,"",$4); print $4; exit}' "$CrosslinkTable"
  )
  [[ -z "$matched_control" ]] && { echo "ERROR: no matched_control for $samplename" >&2; return 1; }

  # control BAM path
  ctrl="${Inputbam%/}/${matched_control}.bam"
  [[ ! -f "$ctrl" ]] && { echo "ERROR: control BAM not found: $ctrl" >&2; return 1; }

  # ensure BAM index exists (create if necessary)
  if [[ ! -f "${ctrl}.bai" && ! -f "${ctrl%.bam}.bai" ]]; then
    echo "[MATCHED CTRL] Index for $ctrl not found; creating..."
    samtools index "$ctrl" || { echo "ERROR: samtools index failed for $ctrl" >&2; return 1; }
  fi

  local out="${mutations%.txt}.control_reads.txt"
  local merged="${mutations%.txt}.with_control.txt"

  # make temporary files
  local bedfile
  local pileupfile
  bedfile=$(mktemp "${TMPDIR:-/tmp}/ctrl_pos.XXXXXX.bed")
  pileupfile=$(mktemp "${TMPDIR:-/tmp}/ctrl_pileup.XXXXXX.txt")

  # cleanup on exit
  trap 'rm -f "$bedfile" "$pileupfile"' RETURN INT TERM

  # create BED with 0-based start, 1-based end for each mutation position
  # (also remove potential CRs)
  awk -F $'\t' 'BEGIN{OFS="\t"} {gsub(/\r/,"",$0); if(NF>=2 && $1!="" && $2!="") print $1, $2-1, $2}' "$mutations" > "$bedfile"

  # run mpileup once for all positions (-l accepts a BED), write stdout to pileupfile
  # use -A -B to include anomalous reads & disable BAQ (same as before)
  echo "[MATCHED CTRL] Running samtools mpileup for all positions..."
  samtools mpileup -f "$Reference" -l "$bedfile" -d 10000000 -A -B "$ctrl" > "$pileupfile" 2>/dev/null || true

  # Now parse the pileup and the mutation file in one awk run:
  #  - first (NR==FNR) load pileup lines into an array keyed by chr:pos
  #  - then, for each mutation line, lookup the pileup entry and compute counts (ref/alt/depth)
  # The counting logic is similar to your previous awk, but runs only once per pileup line.
  awk -v OFS="\t" -v CONTROL="$matched_control" -v OUT="$out" -v MERGED="$merged" '
    function strip_start_end(s,   tmp) {
      # remove start-of-read markers ^X and read-end $ from the mpileup read_bases string
      gsub(/\^./, "", s)
      gsub(/\$/, "", s)
      return s
    }

    # count function: s is read_bases, REF and ALT provided
    function count_ref_alt(s, REF, ALT,   i, L, c, n_ref, n_alt, ref_len, alt_len, ins_seq, del_seq, sign, j, num, indlen, seq, base) {
      n_ref = 0; n_alt = 0
      ref_len = length(REF); alt_len = length(ALT)
      if (alt_len > ref_len) ins_seq = substr(ALT, ref_len+1)
      else ins_seq = ""
      if (ref_len > alt_len) del_seq = substr(REF, alt_len+1)
      else del_seq = ""

      L = length(s); i = 1
      while (i <= L) {
        c = substr(s, i, 1)
        if (c == "." || c == ",") {
          # check following token: could be indel marker +/-
          nextc = substr(s, i+1, 1)
          if (nextc == "+" || nextc == "-") {
            sign = nextc; j = i + 2; num = ""
            while (j <= L && substr(s, j, 1) ~ /[0-9]/) { num = num substr(s, j, 1); j++ }
            indlen = (num == "" ? 0 : num + 0)
            seq = substr(s, j, indlen)
            if (sign == "+") {
              if (ins_seq != "" && seq == ins_seq) n_alt++
            } else {
              if (del_seq != "" && seq == del_seq) n_alt++
            }
            i = j + indlen
          } else {
            n_ref++; i++
          }
        } else if (c == "+" || c == "-") {
          sign = c; j = i + 1; num = ""
          while (j <= L && substr(s, j, 1) ~ /[0-9]/) { num = num substr(s, j, 1); j++ }
          indlen = (num == "" ? 0 : num + 0)
          seq = substr(s, j, indlen)
          if (sign == "+") {
            if (ins_seq != "" && seq == ins_seq) n_alt++
          } else {
            if (del_seq != "" && seq == del_seq) n_alt++
          }
          i = j + indlen
        } else if (c == "*") {
          # placeholder for deletion in mpileup; ignore
          i++
        } else {
          base = toupper(c)
          if (ref_len == 1 && alt_len == 1) {
            if (base == toupper(ALT)) n_alt++
            else if (base == toupper(REF)) n_ref++
            else n_alt++    # conservative: other base = alt
          } else {
            # heuristic for indel/complex: match first base of ALT
            if (alt_len > 0 && base == toupper(substr(ALT,1,1))) n_alt++
            else n_ref++
          }
          i++
        }
      }
      return n_ref " " n_alt
    }

    # --- load pileup into array keyed by chr:pos ---
    FNR==NR {
      # pileup format: chr pos ref depth read_bases qual...
      chr = $1; pos = $2; depth = $4; read_bases = $5
      key = chr ":" pos
      read_bases = strip_start_end(read_bases)
      # store raw read_bases and depth for per-mutation counting later
      pile_rb[key] = read_bases
      pile_dp[key] = depth + 0
      next
    }

    # --- process mutations file lines ---
    {
      line = $0
      # preserve original fields reliably even if info contains tabs (we assume tabs separate the standard 6 columns)
      chrom = $1; pos = $2; ref = $3; alt = $4; sample = $5
      # info might be empty or could contain tabs; reconstruct it starting at $6 if present
      info = ""
      if (NF >= 6) {
        info = $6
        for (i = 7; i <= NF; i++) info = info OFS $i
      }

      key = chrom ":" pos
      if (key in pile_rb) {
        rb = pile_rb[key]
        dp = pile_dp[key] + 0
        # compute counts using actual REF/ALT
        cl = count_ref_alt(rb, ref, alt)
        split(cl, a, " ")
        ref_count = a[1] + 0
        alt_count = a[2] + 0
      } else {
        ref_count = 0; alt_count = 0; dp = 0
      }

      # write control reads output (control name only)
      print chrom, pos, ref, alt, sample, CONTROL, ref_count, alt_count, dp > OUT

      # build merged info: append CONTROL and counts to original info
      ctrl_append = "_CONTROL:" CONTROL "_CREF:" ref_count "_CALT:" alt_count "_CDP:" dp
      if (info == "") new_info = ctrl_append
      else new_info = info ctrl_append

      # write merged line: keep six-columns-like shape: chrom pos ref alt sample new_info
      print chrom, pos, ref, alt, sample, new_info > MERGED
    }
  ' "$pileupfile" "$mutations"

  echo "[MATCHED CTRL] Wrote control reads to: $out"
  echo "[MATCHED CTRL] Wrote merged input+control to: $merged"

  # cleanup (trap will remove tmp files)
  return 0
}
