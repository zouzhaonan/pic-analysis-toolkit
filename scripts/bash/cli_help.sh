#!/usr/bin/env bash

show_main_help() {
  cat "${PIC_HELP_DIR}/help_main.txt"
}

print_registered_biomart_genomes() {
  local lookup_file
  local entries
  local biomart_list_awk

  lookup_file="$(pic_resolve_lib_dir)/register/biomart_lookup.tsv"
  biomart_list_awk="${PIC_ROOT}/scripts/awk/biomart_list_registered.awk"

  echo "Registered biomart genomes:"
  if [[ ! -f "$lookup_file" ]]; then
    echo "  (none)"
    return 0
  fi

  entries="$(awk -f "$biomart_list_awk" "$lookup_file" | sort -u)"
  if [[ -z "$entries" ]]; then
    echo "  (none)"
    return 0
  fi

  while IFS=$'\t' read -r genome dataset; do
    [[ -z "$genome" || -z "$dataset" ]] && continue
    echo "  - ${genome} (${dataset})"
  done <<< "$entries"
}

show_help_for_subcommand() {
  local subcommand="$1"
  local help_file="${PIC_HELP_DIR}/help_${subcommand}.txt"

  if [[ ! -f "$help_file" ]]; then
    show_main_help
    return 0
  fi

  cat "$help_file"

  case "$subcommand" in
    mapping|build-genome)
      echo
      print_registered_pic_genomes
      ;;
    deseq2|manage-biomart)
      echo
      print_registered_biomart_genomes
      ;;
  esac
}
