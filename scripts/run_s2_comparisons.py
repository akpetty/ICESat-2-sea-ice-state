#!/bin/sh

#python S2_IS2_plots_APv4.py -b 'gt1r' -e 'Antarctic'
#python S2_IS2_plots_APv4.py -b 'gt2r' -e 'Antarctic'
#python S2_IS2_plots_APv4.py -b 'gt3r' -e 'Antarctic'
#python S2_IS2_plots_APv4.py -b 'gt1l' -e 'Antarctic'
#python S2_IS2_plots_APv4.py -b 'gt2l' -e 'Antarctic'
#python S2_IS2_plots_APv4.py -b 'gt3l' -e 'Antarctic'

python plot_s2_comparison.py -b 'gt1l' -e 'Arctic'
python plot_s2_comparison.py -b 'gt2l' -e 'Arctic'
python plot_s2_comparison.py -b 'gt3l' -e 'Arctic'
python plot_s2_comparison.py -b 'gt1r' -e 'Arctic'
python plot_s2_comparison.py -b 'gt2r' -e 'Arctic'
python plot_s2_comparison.py -b 'gt3r' -e 'Arctic'