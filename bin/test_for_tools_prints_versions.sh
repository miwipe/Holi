more test_tools_safe.sh 
#!/usr/bin/env bash
set -euo pipefail

TOOLS=(
  bamtools
  bowtie2
  fastANI
  fastp
  megahit
  mosdepth
  samtools
  htsfile
  seqkit
  seqtk
  sga
  spades.py
  vsearch
  filterBAM
  unicorn
)

echo "=== Tool availability + versions ==="
echo "Active conda env: ${CONDA_DEFAULT_ENV:-<none>}"
echo "CONDA_PREFIX: ${CONDA_PREFIX:-<none>}"
echo

# run command, capture first non-empty line from stdout/stderr
first_line() {
  "$@" 2>&1 | sed '/^[[:space:]]*$/d' | head -n 1 || true
}

# print where a command is coming from
where_cmd() {
  command -v "$1" 2>/dev/null || true
}

for t in "${TOOLS[@]}"; do
  printf "%-10s : " "$t"

  path="$(where_cmd "$t")"
  if [[ -z "$path" ]]; then
    echo "NOT FOUND on PATH"
    continue
  fi

  # show path, then version
  printf "%s | " "$path"

  case "$t" in
    bamtools)  first_line bamtools --version ;;
    bowtie2)   first_line bowtie2 --version ;;
    fastANI)   first_line fastANI --version || first_line fastANI -v ;;
    fastp)     first_line fastp --version ;;
    megahit)   first_line megahit --version || first_line megahit -v ;;
    mosdepth)  first_line mosdepth --version || first_line mosdepth -h ;;
    samtools)  first_line samtools --version ;;
    htsfile)   first_line htsfile --version ;;
    seqkit)    first_line seqkit version || first_line seqkit --version ;;
    seqtk)     first_line seqtk ;;  # prints help/version-ish header
    sga)       first_line sga --version || first_line sga ;;
    spades.py) first_line spades.py --version ;;
    vsearch)   first_line vsearch --version ;;
    filterBAM) first_line filterBAM --version || first_line filterBAM -h ;;
    *)         first_line "$t" --version || first_line "$t" -V || first_line "$t" -v || first_line "$t" -h ;;
  esac
done

echo
echo "=== Done ==="

