# Circular visualisation of ancient organelle DNA. 



![](figures/cp_example.png)

Create environment using .yaml file
```
conda env create -f org_visualization.yaml
```
Activate the environment

```
conda activate org_visualization
```

Make a directory where you want to install your programmes
```
mkdir $pwd/PROGRAM
cd $pwd/PROGRAM
```

Install bwa in directory
```
git clone https://github.com/lh3/bwa.git
cd bwa; make
cd ..
```

Install angsd
```
git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/angsd/angsd.git;
cd htslib
make
cd ../angsd
make HTSSRC=../htslib
cd ..
```

Install metaDMG
```
git clone https://github.com/metaDMG-dev/metaDMG-cpp.git
cd metaDMG-cpp
make HTTSRC=../htslib
cd ..
```

Use the script
```
./Alignscript.sh Reference.fasta Reads.fasta
```

The relevant outputs are pmd, bed, depth, and .mafs files. Move these to a directory where R can access them, then use the CircVis R markdown file to run the visualization. You can optionally download genes from your reference genome and add them to the visualization. You can use the Compare R markdown file to compare mutations in two sets of reads.

































