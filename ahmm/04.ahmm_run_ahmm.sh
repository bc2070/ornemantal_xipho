input=APS_all_aim_input
ploidy=ahmm.ploidy
ahmm=/data/software/Ancestry_HMM/src/ancestry_hmm

export PATH=$PATH:/data/software/Ancestry_HMM/src/ancestry_hxmm

$ahmm -@ 12 -i $input -s $ploidy -a 3 0.8 0.18 0.02 -p 0 -20 -0.8 -p 1 -20 -0.18 -p 2 -20 -0.02 --tmax 50 -e 1e-3 >APS_out2.log 2>&1
