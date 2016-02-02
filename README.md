# viridis
Tools for quantifying gene expression in diatoms.  
    
The natural analysis would follow as:   

## cleaning reads
readsCleaner.py: script to call Trimmomatic and clean the reads.

### mapping reads
readsMapper.py: script to call STAR and map the reads to the genome.

#### quantifying reads (counts)
readsCounter.py: script to call HTSeq to counts the reads per transcript.  
samplesCorrelationGrapher.py: script to compute and plot the correlation among samples.  

##### quantifying reads (FPKM)
cufflinksCaller.py: script to call cufflinks and quantify the mapped reads.  
pcaGrapher.py: script to compute the PCA on expression values in FPKM.  
sampleMapper.py: script to discriminate samples in the light-growth space.

logValuesMatrixCreator.py: script to transform absolute (FPKM) into relative expression (log2 fold change).

