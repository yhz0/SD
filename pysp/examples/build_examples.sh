python ../scripts/pysp2smps.py -m farmer/model.py -d farmer/sdinput/farmer_full_names_explicit -b farmer_full_names_explicit --symbolic-solver-labels
rm -rf farmer/sdinput/farmer_full_names_explicit/scenario_files

python ../scripts/pysp2smps.py -m farmer/model.py -d farmer/sdinput/farmer_short_names_explicit -b farmer_short_names_explicit
rm -rf farmer/sdinput/farmer_short_names_explicit/scenario_files

python ../scripts/pysp2smps.py -m baa99/baa99.py -d baa99/sdinput/baa99_full_names_explicit -b baa99_full_names_explicit --symbolic-solver-labels
rm -rf baa99/sdinput/baa99_full_names_explicit/scenario_files

python ../scripts/pysp2smps.py -m baa99/baa99.py -d baa99/sdinput/baa99_short_names_explicit -b baa99_short_names_explicit
rm -rf baa99/sdinput/baa99_short_names_explicit/scenario_files

python ../scripts/pysp2smps.py -m baa99/baa99_exotic.py -d baa99/sdinput/baa99_exotic_full_names_explicit -b baa99_exotic_full_names_explicit --symbolic-solver-labels
rm -rf baa99/sdinput/baa99_exotic_full_names_explicit/scenario_files

python ../scripts/pysp2smps.py -m baa99/baa99_exotic.py -d baa99/sdinput/baa99_exotic_short_names_explicit -b baa99_exotic_short_names_explicit
rm -rf baa99/sdinput/baa99_exotic_short_names_explicit/scenario_files
