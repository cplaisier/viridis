# viridis
Tools for quantifying gene expression in diatoms.  
    
The natural analysis would follow as:   

#### cleaning reads
readsCleaner.py: script to call Trimmomatic and clean the reads.

#### mapping reads
readsMapper.py: script to call STAR and map the reads to the genome.

#### quantifying reads (counts)
readsCounter.py: script to call HTSeq to counts the reads per transcript.  
samplesCorrelationGrapher.py: script to compute and plot the correlation among samples.  

#### quantifying reads (FPKM)
cufflinksCaller.py: script to call cufflinks and quantify the mapped reads.  

#### finding differentially expressed transcripts
classifiersFinder.py: script that calls cuffdiff to find DETs between two conditions.  

#### visualizing samples
pcaGrapher.py: script to compute the PCA on expression values in FPKM.  
sampleMapper.py: script to map samples in diurnal, growth space trained in epochs 0 and 1.  
sampleMapper_epoch0.py: script to map samples in diurnal, growth space trained only in epoch 0.  

#### misc
logValuesMatrixCreator.py: script to transform absolute (FPKM) into relative expression (log2 fold change).

