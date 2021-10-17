
#flip ancestral state

for i in {1..16} 
do 
/ohta1/apps/relate_v1.1.6_x86_64_static/bin/RelateFileFormats --mode FlipHapsUsingAncestor \
--haps Scaffold_${i}_phased.vcf.haps \
--sample Scaffold_${i}_phased.vcf.sample \
--ancestor AMATA_finishedtohypo_renamedfinal_palmacnestral_modheader.fasta \
-o palmaligned_${i}_flipped.haps

done

#FIRST ESTIMATION of trees
for i in {1..16}
do

/ohta1/apps/relate_v1.1.6_x86_64_static/scripts/RelateParallel/RelateParallel.sh \
      --mode All \
      -m 7e-9 \
      -N 500000 \
      --haps ../palmaligned_${i}_flipped.haps \
      --sample /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/Scaffold_${i}_shapeit.sample \
      --map ../../Scaffold_${i}_shapeit2.map \
      --seed 1 \
      -o palmref_chr${i} \
      --threads 20

done


#ESTIMATE POP SIZE
/ohta1/apps/relate_v1.1.6_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i palm_aligned_bys \
                --first_chr 1 \
                --last_chr 16 \
              -m 7.910025e-09 \
              --poplabels detailedpops.label \
              --seed 1 \
	      --years_per_gen 1 --bins 1,6,0.25 \
	      --threads 20 \
              -o /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/relate_runs/palm_aligned_bys_popsize


#ESTIMATE MUTATION RATE
/ohta1/apps/relate_v1.1.6_x86_64_static/bin/RelateMutationRate --mode Avg \
--first_chr 1 --last_chr 16 \
-i palm_aligned_bys \
-o /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/relate_runs/formatt_converged/palm_aligned_byregion_chr

#RESTIMATE BRANCH LENGTHS
for i in {1..16}
do

/ohta1/apps/relate_v1.1.6_x86_64_static/bin/RelateCoalescentRate --mode ReEstimateBranchLengths \
--input palm_aligned_bys \
--mrate /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/relate_runs/formatt_convergedpalm_aligned_byregion_chr_avg.rate \
--coal palm_aligned_bys_popsize.coal \
-o /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/relate_runs/formatt_converged/palm_aligned_bysubspecies_gen_updated \
--mutation_rate 7.910025e-09

done


###########################################
#for CLUES inference
###########################################
#subset haplotypes to include resistance haplotypes of focal origin + outgroup
cat ALS574_2 midwest_sus > subset_ALS574
cat ALS653_2 ontario_sus > subset_ALS653

comm -23 <(sed 's/..$//g' subset_ALS574 | uniq | sort ) <(cat <(sed 's/..$//g' ../vis_haps/ALS574_1 | uniq) <(sed 's/..$//g' ../vis_haps/ALS574_3 ) | sed '/^$/d' | sort | uniq) | uniq > uniqinds_ALS574
comm -23 <(sed 's/..$//g' subset_ALS653 | uniq | sort ) <(sed 's/..$//g' ../vis_haps/ALS653_1 | uniq | sed '/^$/d' | sort ) | uniq > uniqinds_ALS653


grep -vf uniqinds_ALS574 /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/Scaffold_11_shapeit.sample | tail -n+3 | cut -d" " -f1 > notALS574_origin2 #haploypes to exclude
grep -vf uniqinds_ALS653 /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/Scaffold_11_shapeit.sample | tail -n+3 | cut -d" " -f1 > notALS653_origin5 #haploypes to exclude

#subset datasets by mutational origin
for i in notALS574_origin2 notALS653_origin5
do
/ohta1/apps/relate_v1.1.6_x86_64_static/bin/RelateFileFormats \
                 --mode RemoveSamples \
                 --haps /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/bams/processed/shapeit4/remade.haps.haps \
                 --sample /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/bams/processed/shapeit4/remade.haps.mod.sample \
                 -i ${i} \
                 -o ${i} &
done

#rerun RELATE locally on subsetted samples
for i in notALS574_origin2 notALS653_origin5
do
/ohta1/apps/relate_v1.1.6_x86_64_static/scripts/RelateParallel/RelateParallel.sh \
      --mode All \
      -m 7.910025e-09 \
      --coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
      --haps ../${i}.haps \
      --sample ../${i}.sample \
      --map /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/Scaffold_11_shapeit2.map \
      --seed 1 \
      -N 500000 \
      -o ${i} \
      --threads 40 
done

#restimate branch lengths
for i in notALS574_origin2 notALS653_origin5
do
/ohta1/apps/relate_v1.1.6_x86_64_static/bin/RelateCoalescentRate --mode ReEstimateBranchLengths \
--input ${i} \
--mrate /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/relate_runs/formatt_converged/palm_aligned_byregion_chr_avg.rate \
--coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
-o ${i}_resamp \
--mutation_rate 7.910025e-09
done

for i in notALS653_origin5
do
/ohta1/apps/relate_v1.1.6_x86_64_static/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i ${i}_resamp \
        -o ${i}_resamp_branches \
        -m 7.910025e-09 --coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
        --format b \
        --num_samples 100 \
        --first_bp 24249814 --last_bp 24249814
done


for i in notALS574_origin2
do
/ohta1/apps/relate_v1.1.6_x86_64_static/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i ${i}_resamp \
        -o ${i}_resamp_branches \
        -m 7.910025e-09 --coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
        --format b \
        --num_samples 100 \
        --first_bp 24250051 --last_bp 24250051
done

for i in notALS574_origin2
do
python3 /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/CLUES/clues_new/inference.py  \
--times $[i}_resamp_branches \
--coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
--popFreq 0.31 \
--timeBins 10time.bins \
--out $[i}_clues_10gens

python3 /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/CLUES/clues_new/inference.py \
--times $[i}_resamp_branches \
--coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
--popFreq 0.31 \
--timeBins 30time.bins \
--out ${i}_clues_30gens
done

for i in notALS574_origin2
do
python3 /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/CLUES/clues_new/inference.py \
--times $[i}_resamp_branches --coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
--popFreq 0.29 \
--timeBins 10time.bins \
--out ${i}_clues_10gens

python3 /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/CLUES/clues_new/inference.py \
--times $[i}_resamp_branches \
--coal /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/processing/formatt/palm_aligned_bys_bys_popsize_025bins1yr.coal \
--popFreq 0.29 \
--timeBins 30time.bins \
--out ${i}_clues_30gens
done
