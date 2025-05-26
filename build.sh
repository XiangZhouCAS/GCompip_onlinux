#!/bin/bash
chmod +x gcompip
chmod +x ./bins/com
chmod +x ./bins/custom
chmod +x install.sh
comts_path=$(dirname "$(readlink -f "$0")")
echo 'export PATH="$PATH:'"$comts_path"'"' >> ~/.bashrc
bins_dir=$comts_path/$(echo bins)
sed -i "1s|.*|R_SCRIPTS_DIR=$bins_dir|" $comts_path/gcompip
sed -i "1s|.*|R_SCRIPTS_DIR=$bins_dir|" $comts_path/bins/com
sed -i "1s|.*|R_SCRIPTS_DIR=$bins_dir|" $comts_path/bins/custom
gunzip ter.dmnd.gz
R --no-save < ./bins/install.r
