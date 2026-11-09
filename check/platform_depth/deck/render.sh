#!/usr/bin/env bash
# Render platform_depth_summary.pptx to PDF (LibreOffice, DengXian drawn as Liberation Sans) and to slide images.
# Usage: bash render.sh [img_dir]     (img_dir: where slide-NN.png go; default: no images)
set -euo pipefail
cd "$(dirname "$0")"
SOFFICE=/home/b/b37974/opt/lo/root/opt/libreoffice25.8/program/soffice
PROFILE=$(mktemp -d)
FONTCONFIG_FILE="$PWD/fonts.conf" "$SOFFICE" -env:UserInstallation=file://"$PROFILE" --headless \
  --convert-to pdf --outdir . platform_depth_summary.pptx >/dev/null 2>&1
rm -rf "$PROFILE"
if [ $# -ge 1 ]; then
  mkdir -p "$1"; rm -f "$1"/slide-*.png
  pdftoppm -png -r 80 platform_depth_summary.pdf "$1/slide"
fi
ls -la platform_depth_summary.pdf
